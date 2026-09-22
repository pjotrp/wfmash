// Vendored from wfmash v0.14.1 (branch v0.14.1, commit 9b2a7388) for the
// low-divergence engine (src/engine/low-divergence/).  Files are renamed
// with an lde_ prefix and the mashmap/yeet/align namespaces are prefixed
// lde_ so the 0.14-lineage engine code cannot collide with the mainline
// 0.24 engine.  Provenance: waveygang/wfmash.
/**
 * @file    align.cpp
 * @ingroup src
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include <sstream>
#include <fstream>
#include <iostream>
#include <ctime>
#include <chrono>
#include <functional>
#include <cstdio>

#include "engine/low-divergence/lde_map_parameters.hpp"
#include "engine/low-divergence/lde_base_types.hpp"
#include "engine/low-divergence/lde_winSketch.hpp"
#include "engine/low-divergence/lde_computeMap.hpp"
#include "engine/low-divergence/lde_parseCmdArgs.hpp"

#include "engine/low-divergence/lde_parse_args.hpp"

#include "engine/low-divergence/lde_align_parameters.hpp"
#include "engine/low-divergence/lde_parseCmdArgs.hpp"



//External includes
#include "common/args.hxx"
#include "engine/low-divergence/lde_ALeS.hpp"
#include "engine/low-divergence/lde_agc_index.hpp"

int lde_main(int argc, char** argv) {
    /*
     * Make sure env variable MALLOC_ARENA_MAX is unset 
     * for efficient multi-thread execution
     */
    unsetenv((char *)"MALLOC_ARENA_MAX");

    // get our parameters from the command line
    lde_skch::Parameters map_parameters;
    lde_align::Parameters align_parameters;
    lde_yeet::Parameters yeet_parameters;
    lde_yeet::parse_args(argc, argv, map_parameters, align_parameters, yeet_parameters);

    //parameters.refSequences.push_back(ref);

    //lde_skch::parseandSave(argc, argv, cmd, parameters);
    if (!yeet_parameters.remapping) {
        lde_skch::printCmdOptions(map_parameters);

        auto t0 = lde_skch::Time::now();

        if (map_parameters.use_spaced_seeds) {
          std::cerr << "[wfmash::map] Generating spaced seeds" << std::endl;
          uint32_t seed_weight = map_parameters.spaced_seed_params.weight;
          uint32_t seed_count = map_parameters.spaced_seed_params.seed_count;
          float similarity = map_parameters.spaced_seed_params.similarity;
          uint32_t region_length = map_parameters.spaced_seed_params.region_length;

          lde_ales::spaced_seeds sps = lde_ales::generate_spaced_seeds(seed_weight, seed_count, similarity, region_length);
          std::chrono::duration<double> time_spaced_seeds = lde_skch::Time::now() - t0;
          std::cerr << "[wfmash::map] Time spent generating spaced seeds " << time_spaced_seeds.count()  << " seconds" << std::endl;
          map_parameters.spaced_seed_sensitivity = sps.sensitivity;
          map_parameters.spaced_seeds =  sps.seeds;
          lde_ales::printSpacedSeeds(map_parameters.spaced_seeds);
          std::cerr << "[wfmash::map] Spaced seed sensitivity " << sps.sensitivity << std::endl;
        }

        //Build the sketch for reference
        lde_skch::Sketch referSketch(map_parameters);

        std::chrono::duration<double> timeRefSketch = lde_skch::Time::now() - t0;
        std::cerr << "[wfmash::map] time spent computing the reference index: " << timeRefSketch.count() << " sec" << std::endl;

        if (referSketch.minmerIndex.size() == 0)
        {
            std::cerr << "[wfmash::map] ERROR, reference sketch is empty. Reference sequences shorter than the segment length are not indexed" << std::endl;
            return 1;
        }

        //Map the sequences in query file
        t0 = lde_skch::Time::now();

        lde_skch::Map mapper = lde_skch::Map(map_parameters, referSketch);

        std::chrono::duration<double> timeMapQuery = lde_skch::Time::now() - t0;
        std::cerr << "[wfmash::map] time spent mapping the query: " << timeMapQuery.count() << " sec" << std::endl;
        std::cerr << "[wfmash::map] mapping results saved in: " << map_parameters.outFileName << std::endl;

        if (yeet_parameters.approx_mapping) {
            return 0;
        }
    }

    // [low-divergence engine] The alignment stage of the 0.14 engine
    // (the full wflign aligner) is not vendored into this build yet;
    // this engine is mapping-only for now.
    std::cerr << "[wfmash::lde] ERROR: the low-divergence engine currently"
                 " supports mapping only; pass -m/--approx-map" << std::endl;
    return 1;
}
