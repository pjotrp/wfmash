/**
 * @file    computeAlignments.hpp
 * @brief   logic for generating alignments when given mashmap 
 *          mappings as input
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef COMPUTE_ALIGNMENTS_HPP 
#define COMPUTE_ALIGNMENTS_HPP

#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <zlib.h>
#include <cassert>
#include <thread>
#include <memory>
#include <htslib/faidx.h>

//Own includes
#include "align/include/align_types.hpp"
#include "align/include/align_parameters.hpp"
#include "map/include/base_types.hpp"
#include "map/include/commonFunc.hpp"

//External includes
#include "common/wflign/src/wflign.hpp"
#include "common/atomic_queue/atomic_queue.h"
#include "common/seqiter.hpp"
#include "common/progress.hpp"
#include "common/utils.hpp"

namespace align
{

  long double float2phred(long double prob) {
    if (prob == 1)
        return 255;  // guards against "-0"
    long double p = -10 * (long double) log10(prob);
    if (p < 0 || p > 255) // int overflow guard
        return 255;
    else
        return p;
  }

  struct seq_record_t {
      MappingBoundaryRow currentRecord;
      std::string mappingRecordLine;
      std::shared_ptr<std::string> qSequence;
      seq_record_t(const MappingBoundaryRow& c, const std::string& r, const shared_ptr<std::string>& q)
          : currentRecord(c)
          , mappingRecordLine(r)
          , qSequence(q)
          { }
  };
  // load into this
  typedef atomic_queue::AtomicQueue<seq_record_t*, 2 << 16> seq_atomic_queue_t;
  // results into this, write out
  typedef atomic_queue::AtomicQueue<std::string*, 2 << 16> paf_atomic_queue_t;

  /**
   * @class     align::Aligner
   * @brief     compute alignments and generate sam output
   *            from mashmap mappings
   */
  class Aligner
  {
    private:

      //algorithm parameters
      const align::Parameters &param;

      refSequenceMap_t refSequences;
      std::vector<faidx_t*> faidxs;

    public:
      /**
       * @brief                 destructor, cleans up faidx index
       */
      ~Aligner()
      {
          for (auto& faid : this->faidxs) {
              fai_destroy(faid);
          }
      }

      /**
       * @brief                 constructor, also reads reference sequences
       * @param[in] p           algorithm parameters
       */
      explicit Aligner(const align::Parameters &p) :
        param(p)
      {
          this->getRefSequences();
      }

      /**
       * @brief                 compute alignments
       */
      void compute()
      {
        this->computeAlignments();
      }

      /**
       * @brief       parse mashmap row sequence
       * @param[in]   mappingRecordLine
       * @param[out]  currentRecord
       */
      inline static void parseMashmapRow(const std::string &mappingRecordLine, MappingBoundaryRow &currentRecord)
      {
        std::stringstream ss(mappingRecordLine); // Insert the string into a stream
        std::string word; // Have a buffer string

        vector<std::string> tokens; // Create vector to hold our words

        while (ss >> word)
          tokens.push_back(word);

        //We expect and need at least these many values in a mashmap mapping
        assert(tokens.size() >= 9);

        // Extract the mashmap identity from the string
        const vector<string> mm_id_vec = skch::CommonFunc::split(tokens[12], ':');
        // if the estimated identity is missing, avoid assuming too low values
        const float mm_id = wfmash::is_a_number(mm_id_vec.back()) ? std::stof(mm_id_vec.back())/(float) 100.0 : skch::fixed::percentage_identity; // divide by 100 for consistency with block alignment

        //Save words into currentRecord
        {
            currentRecord.qId = tokens[0];
            currentRecord.qStartPos = std::stoi(tokens[2]);
            currentRecord.qEndPos = std::stoi(tokens[3]);
            currentRecord.strand = (tokens[4] == "+" ? skch::strnd::FWD : skch::strnd::REV);
            currentRecord.refId = tokens[5];
            currentRecord.rStartPos = std::stoi(tokens[7]);
            currentRecord.rEndPos = std::stoi(tokens[8]);
            currentRecord.mashmap_estimated_identity = mm_id;
        }
      }

  private:

      /**
       * @brief                 parse and save all the reference sequences
       */
      void getRefSequences() {
          auto& nthreads = param.threads;
          assert(param.refSequences.size() == 1);
          auto& filename = param.refSequences.front();
          for (int i = 0; i < nthreads; ++i) {
              auto faid = fai_load(filename.c_str());
              faidxs.push_back(faid);
          }
      }

      /**
       * @brief                 parse query sequences and mashmap mappings
       *                        to compute sequence alignments
       */
      void computeAlignments()
      {
          uint64_t total_seqs = 0;

          // Count the number of mapped bases to align
          uint64_t total_alignment_length = 0;
          {
              std::ifstream mappingListStream(param.mashmapPafFile);
              std::string mappingRecordLine;
              MappingBoundaryRow currentRecord;

              while(!mappingListStream.eof()) {
                  std::getline(mappingListStream, mappingRecordLine);
                  if (!mappingRecordLine.empty()) {
                      parseMashmapRow(mappingRecordLine, currentRecord);
                      total_alignment_length += currentRecord.qEndPos - currentRecord.qStartPos;
                  }
              }
          }

          progress_meter::ProgressMeter progress(total_alignment_length, "[wfmash::align::computeAlignments] aligned");

          // input atomic queue
          seq_atomic_queue_t seq_queue;
          // output atomic queues
          paf_atomic_queue_t paf_queue, tsv_queue;
          // flag when we're done reading
          std::atomic<bool> reader_done;
          reader_done.store(false);

          auto& nthreads = param.threads;
          //for (

          // atomics to record if we're working or not
          std::vector<std::atomic<bool>> working(nthreads);
          for (auto& w : working) {
              w.store(true);
          }

          // reader picks up candidate alignments from input
          auto reader_thread =
              [&]() {
                  //Parse query sequences
                  for(const auto &fileName : param.querySequences)
                  {
//#define DEBUG true
#ifdef DEBUG
                      std::cerr << "INFO, align::Aligner::computeAlignments, parsing query sequences in file " << fileName << std::endl;
#endif

                      //Open mashmap output file
                      std::ifstream mappingListStream(param.mashmapPafFile);
                      std::string mappingRecordLine;
                      MappingBoundaryRow currentRecord;

                      seqiter::for_each_seq_in_file(
                          fileName, {}, "", 
                          [&](const std::string& qSeqId,
                              const std::string& _seq) {
                              ++total_seqs;

                              uint64_t rank_mapping = 0;

                              // copy our input into a shared ptr
                              std::shared_ptr<std::string> seq(new std::string(_seq));
                              // todo: offset_t is an 32-bit integer, which could cause problems
                              skch::offset_t len = seq->length();
                              // upper-case our input and make sure it's canonical DNA (for WFA)
                              skch::CommonFunc::makeUpperCaseAndValidDNA((char*)seq->c_str(), len);
                              // todo maybe this should change to some kind of unique pointer?
                              // something where we can GC it when we're done aligning to it
                              //std::string qSequence = seq;
                              //std::cerr << seq << std::endl;

                              //Check if all mapping records are processed already
                              while(!mappingListStream.eof() && mappingRecordLine.empty()) {
                                  //Read first record from mashmap output file during first iteration
                                  std::getline(mappingListStream, mappingRecordLine);
                              }

                              if( !mappingRecordLine.empty() ) {
                                  parseMashmapRow(mappingRecordLine, currentRecord);

                                  //Check if mapping query id matches current query sequence id
                                  if(currentRecord.qId == qSeqId)
                                  {
                                      //Queue up this query record
                                      auto q = new seq_record_t(currentRecord, mappingRecordLine, seq);
                                      q->currentRecord.rankMapping = rank_mapping++;
                                      seq_queue.push(q);

                                      //Check if more mappings have same query sequence id
                                      while(std::getline(mappingListStream, mappingRecordLine))
                                      {
                                          parseMashmapRow(mappingRecordLine, currentRecord);

                                          if(currentRecord.qId != qSeqId)
                                          {
                                              //Break the inner loop to read query sequence
                                              break;
                                          }
                                          else
                                          {
                                              auto q = new seq_record_t(currentRecord, mappingRecordLine, seq);
                                              q->currentRecord.rankMapping = rank_mapping++;
                                              seq_queue.push(q);
                                          }
                                      }
                                  }
                              }
                          });

                      mappingListStream.close();

                  }
                  reader_done.store(true);
              };

          // helper to check if we're still aligning
          auto still_working =
              [&](const std::vector<std::atomic<bool>>& working) {
                  bool ongoing = false;
                  for (auto& w : working) {
                      ongoing = ongoing || w.load();
                  }
                  return ongoing;
              };

          // writer, picks output from queue and writes it to our output stream
          std::ofstream outstrm(param.pafOutputFile, ios::app);

          auto writer_thread =
              [&]() {
                  while (true) {
                      std::string* paf_lines = nullptr;
                      if (!paf_queue.try_pop(paf_lines)
                          && !still_working(working)) {
                          break;
                      } else if (paf_lines != nullptr) {
                          outstrm << *paf_lines;
                          delete paf_lines;
                      } else {
                          std::this_thread::sleep_for(100ns);
                      }
                  }
              };

#ifdef WFA_PNG_AND_TSV
          auto writer_thread_tsv =
                  [&]() {
              if (!param.tsvOutputPrefix.empty()) {
                  uint64_t num_alignments_completed = 0;

                  while (true) {
                      std::string* tsv_lines = nullptr;
                      if (!tsv_queue.try_pop(tsv_lines)
                      && !still_working(working)) {
                          break;
                      } else if (tsv_lines != nullptr) {
                          std::ofstream ofstream_tsv(param.tsvOutputPrefix + std::to_string(num_alignments_completed++) + ".tsv");
                          ofstream_tsv << *tsv_lines;
                          ofstream_tsv.close();

                          delete tsv_lines;
                      } else {
                          std::this_thread::sleep_for(100ns);
                      }
                  }
              }
          };
#endif

          // worker, takes candidate alignments and runs wfa alignment on them
          auto worker_thread = 
              [&](uint64_t tid,
                  std::atomic<bool>& is_working) {
                  is_working.store(true);
                  while (true) {
                      seq_record_t* rec = nullptr;
                      if (!seq_queue.try_pop(rec)
                          && reader_done.load()) {
                          break;
                      } else if (rec != nullptr) {
                          std::stringstream output;
#ifdef WFA_PNG_AND_TSV
                          std::stringstream output_tsv;
#endif
                          doAlignment(
                                  output,
#ifdef WFA_PNG_AND_TSV
                                  output_tsv,
#endif
                                  rec->currentRecord,
                                  rec->mappingRecordLine,
                                  rec->qSequence, tid);
                          progress.increment(rec->currentRecord.qEndPos - rec->currentRecord.qStartPos);

                          auto* paf_rec = new std::string(output.str());
                          if (!paf_rec->empty()) {
                              paf_queue.push(paf_rec);
                          } else {
                              delete paf_rec;
                          }

#ifdef WFA_PNG_AND_TSV
                          auto* tsv_rec = new std::string(output_tsv.str());
                          if (!tsv_rec->empty()) {
                              tsv_queue.push(tsv_rec);
                          } else {
                              delete tsv_rec;
                          }
#endif

                          delete rec;
                      } else {
                          std::this_thread::sleep_for(100ns);
                      }
                  }
                  is_working.store(false);
              };

          // launch reader
          std::thread reader(reader_thread);
          // launch PAF/SAM writer
          std::thread writer(writer_thread);
#ifdef WFA_PNG_AND_TSV
          // launch TSV writer
          std::thread writer_tsv(writer_thread_tsv);
#endif
          // launch workers
          std::vector<std::thread> workers; workers.reserve(nthreads);
          for (uint64_t t = 0; t < nthreads; ++t) {
              workers.emplace_back(worker_thread,
                                   t,
                                   std::ref(working[t]));
          }

          // wait for reader and workers to complete
          reader.join();
          for (auto& worker : workers) {
              worker.join();
          }
          // and finally the writer
          writer.join();
#ifdef WFA_PNG_AND_TSV
          writer_tsv.join();
#endif
          progress.finish();
          std::cerr << "[wfmash::align::computeAlignments] "
                    << "count of mapped reads = " << total_seqs
                    << ", total aligned bp = " << total_alignment_length << std::endl;
      }

        /**
       * @brief                           compute alignment using WFA
       * @param[in]   currentRecord       mashmap mapping parsed information
       * @param[in]   mappingRecordLine   mashmap mapping output raw string
       * @param[in]   qSequence           query sequence
       * @param[in]   outstrm             output stream
       */
      void doAlignment(
              std::stringstream& output,
#ifdef WFA_PNG_AND_TSV
              std::stringstream& output_tsv,
#endif
              MappingBoundaryRow &currentRecord,
              const std::string &mappingRecordLine,
              const std::shared_ptr<std::string> &qSequence,
              uint64_t tid) {

#ifdef DEBUG
        std::cerr << "INFO, align::Aligner::doAlignment, aligning mashmap record: " << mappingRecordLine << std::endl;
#endif

        //Obtain reference substring for this mapping
        // htslib caches are not threadsafe! so we use a thread-specific faidx_t
        faidx_t* faid = faidxs[tid];
        const int64_t ref_size = faidx_seq_len(faid, currentRecord.refId.c_str());

        // Take flanking sequences to support head/tail patching due to noisy (inaccurate) mapping boundaries
        const uint64_t head_padding = currentRecord.rStartPos >= param.wflign_max_len_minor ? param.wflign_max_len_minor : currentRecord.rStartPos;
        const uint64_t tail_padding = ref_size - currentRecord.rEndPos >= param.wflign_max_len_minor ? param.wflign_max_len_minor : ref_size - currentRecord.rEndPos;

        int64_t got_seq_len = 0;
        char * ref_seq = faidx_fetch_seq64(
                faid, currentRecord.refId.c_str(),
                currentRecord.rStartPos - head_padding,
                currentRecord.rEndPos + tail_padding,
                &got_seq_len
                );

        // hack to make it 0-terminated as expected by WFA
        ref_seq[got_seq_len] = '\0';

        // upper-case our input and make sure it's canonical DNA (for WFA)
        skch::CommonFunc::makeUpperCaseAndValidDNA(ref_seq, got_seq_len);

        // Shift the pointer to the currentRecord.rStartPos position
        ref_seq = ref_seq + head_padding;

        skch::offset_t refLen = currentRecord.rEndPos - currentRecord.rStartPos;

        //Define query substring for this mapping
        const char* queryRegion = qSequence->c_str();  //initially point to beginning
        const auto& querySize = qSequence->size();
        skch::offset_t queryLen = currentRecord.qEndPos - currentRecord.qStartPos;
        queryRegion += currentRecord.qStartPos;

        char* queryRegionStrand = new char[queryLen+1];

        if(currentRecord.strand == skch::strnd::FWD) {
          strncpy(queryRegionStrand, queryRegion, queryLen);    //Copy the same string
        } else {
          skch::CommonFunc::reverseComplement(queryRegion, queryRegionStrand, queryLen); //Reverse complement
        }

        assert(queryLen <= querySize);

        //Compute alignment
#ifdef DEBUG
        std::cerr << "INFO, align::Aligner::doAlignment, WFA execution starting, query region length = " << queryLen
          << ", reference region length= " << refLen << ", edit distance limit= " << editDistanceLimit << std::endl; 
#endif

        // To distinguish split alignment in SAM output format (currentRecord.rankMapping == 0 to avoid the suffix in there is just one alignment for the query)
        const std::string query_name_suffix = param.split && param.sam_format ? "_" + std::to_string(currentRecord.rankMapping) : "";

        wflign::wavefront::WFlign* wflign = new wflign::wavefront::WFlign(
                param.wflambda_segment_length,
                param.min_identity,
                param.wfa_mismatch_score,
                param.wfa_gap_opening_score,
                param.wfa_gap_extension_score,
                currentRecord.mashmap_estimated_identity,
                param.wflign_mismatch_score,
                param.wflign_gap_opening_score,
                param.wflign_gap_extension_score,
                param.wflign_max_mash_dist,
                param.wflign_min_wavefront_length,
                param.wflign_max_distance_threshold,
                param.wflign_max_len_major,
                param.wflign_max_len_minor,
                param.wflign_erode_k,
                param.chain_gap,
                param.wflign_max_patching_score);
        wflign->set_output(
                &output,
#ifdef WFA_PNG_AND_TSV
                !param.tsvOutputPrefix.empty(),
                &output_tsv,
                param.prefix_wavefront_plot_in_png,
                param.wfplot_max_size,
#endif
                true, // merge alignments
                param.emit_md_tag,
                !param.sam_format,
                param.no_seq_in_sam);
        wflign->wflign_affine_wavefront(
                currentRecord.qId + query_name_suffix,
                queryRegionStrand,
                querySize,
                currentRecord.qStartPos,
                queryLen,
                currentRecord.strand != skch::strnd::FWD,
                currentRecord.refId,
                ref_seq,
                ref_size,
                currentRecord.rStartPos,
                refLen);
        delete wflign;

        delete [] queryRegionStrand;

        // Re-shift the pointer to the malloc()-ed address
        ref_seq = ref_seq - head_padding;

        free(ref_seq);
      }
  };
}


#endif
