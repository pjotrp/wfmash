#ifndef WFLIGN_PATCH_HPP_
#define WFLIGN_PATCH_HPP_

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>
#include <sstream>
#include <functional>
#include <fstream>
#include <lodepng/lodepng.h>

#include "dna.hpp"
#include "rkmh.hpp"
#include "wflign.hpp"
#include "wflign_alignment.hpp"

/*
 * Configuration
 */
//#define WFLIGN_DEBUG true // for debugging messages

/*
 * Namespaces
 */
namespace wflign {
    namespace wavefront {

        bool do_wfa_segment_alignment(
                const std::string& query_name,
                const char* query,
                std::vector<rkmh::hash_t>*& query_sketch,
                const uint64_t& query_length,
                const int& j,
                const std::string& target_name,
                const char* target,
                std::vector<rkmh::hash_t>*& target_sketch,
                const uint64_t& target_length,
                const int& i,
                const uint16_t& segment_length_q,
                const uint16_t& segment_length_t,
                const uint16_t& step_size,
                wflign_extend_data_t* extend_data,
                alignment_t& aln);
        void do_wfa_patch_alignment(
                const char* query,
                const uint64_t& j,
                const uint64_t& query_length,
                const char* target,
                const uint64_t& i,
                const uint64_t& target_length,
                wfa::WFAlignerGapAffine& _wf_aligner,
                const wflign_penalties_t& affine_penalties,
                alignment_t& aln,
                const int64_t& chain_gap,
                const int& max_patching_score);
        void write_merged_alignment(
                std::ostream &out,
                const std::vector<alignment_t *> &trace,
                wfa::WFAlignerGapAffine& wf_aligner,
                const wflign_penalties_t& affine_penalties,
                const bool& emit_md_tag,
                const bool& paf_format_else_sam,
                const bool& no_seq_in_sam,
                const char* query,
                const std::string& query_name,
                const uint64_t& query_total_length,
                const uint64_t& query_offset,
                const uint64_t& query_length,
                const bool& query_is_rev,
                const char* target,
                const std::string& target_name,
                const uint64_t& target_total_length,
                const uint64_t& target_offset,
                const uint64_t& target_length,
                const float& min_identity,
                const long& elapsed_time_wflambda_ms,
                const uint64_t& num_alignments,
                const uint64_t& num_alignments_performed,
                const float& mashmap_estimated_identity,
                const uint64_t& wflign_max_len_major,
                const uint64_t& wflign_max_len_minor,
                const int& erode_k,
                const int64_t& chain_gap,
                const int& max_patching_score,
                const int& min_wf_length,
                const int& max_dist_threshold,
#ifdef WFA_PNG_AND_TSV
                const std::string* prefix_wavefront_plot_in_png,
                const uint64_t& wfplot_max_size,
                const bool& emit_patching_tsv,
                std::ostream* out_patching_tsv,
#endif
                const bool& with_endline = true);
        void write_alignment(
                std::ostream& out,
                const alignment_t& aln,
                const std::string& query_name,
                const uint64_t& query_total_length,
                const uint64_t& query_offset,
                const uint64_t& query_length,
                const bool& query_is_rev,
                const std::string& target_name,
                const uint64_t& target_total_length,
                const uint64_t& target_offset,
                const uint64_t& target_length, // unused
                const float& min_identity,
                const float& mashmap_estimated_identity,
                const bool& with_endline = true);
        double float2phred(const double& prob);
        void sort_indels(std::vector<char>& v);

    } /* namespace wavefront */

    void encodeOneStep(const char *filename, std::vector<unsigned char> &image, unsigned width, unsigned height);
} /* namespace wflign */

#endif /* WFLIGN_PATCH_HPP_ */