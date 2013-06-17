// ==========================================================================
//                         Mason - A Read Simulator
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Given a genome, simulate fragments.
// ==========================================================================

// TODO(holtgrew): Make methylation levels optional.

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <seqan/arg_parse.h>

#include "fragment_generation.h"

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    MasonFragmentsOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cerr << "MASON FRAGMENT SIMULATOR\n"
              << "========================\n\n";
    
    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cerr << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY     \t" << options.verbosity << '\n'
                  << "\n"
                  << "SEED          \t" << options.seed << '\n'
                  << "FRAGMENT COUNT\t" << options.numFragments << '\n'
                  << "\n"
                  << "OUTPUT FILE   \t" << options.outputFilename << "\n"
                  << "INPUT FILE    \t" << options.inputFilename << "\n"
                  << "\n"
                  << "DISTRIBUTION  \t" << distributionStr(options.distribution) << "\n"
                  << "MIN SIZE      \t" << options.minFragmentSize << "\n"
                  << "MAX SIZE      \t" << options.maxFragmentSize << "\n"
                  << "MEAN SIZE     \t" << options.meanFragmentSize << "\n"
                  << "SIZE STD DEV  \t" << options.stdDevFragmentSize << "\n"
                  << "\n"
                  << "BS-SEQ SIMULATION\n"
                  << "  ENABLED     \t" << yesNo(options.bsSimEnabled) << "\n"
                  << "  PROTOCOL    \t" << bsProtocolStr(options.bsProtocol) << "\n"
                  << "  METH LEVELS \t" << options.methLevelsInFasta << "\n"
                  << "  CONVERSION  \t" << options.bsConversionRate << "\n"
                  << "\n";
    }

    std::cerr << "__PREPARATION________________________________________________________________\n"
              << "\n";

    // The Genome allows loading of all information for a contig, one contig at a time.
    Genome genome(options.bsSimEnabled);

    // Loading sequence FAI Index.
    std::cerr << "Loading Sequence        " << options.inputFilename << " ...";
    if (read(genome.seqFaiIndex, toCString(options.inputFilename)) != 0)
    {
        std::cerr << " FAILED (not fatal, we can just build it)\n";
        std::cerr << "Building FAI          " << options.inputFilename << ".fai ...";
        if (build(genome.seqFaiIndex, toCString(options.inputFilename)) != 0)
        {
            std::cerr << "Could not build FAI index.\n";
            return 1;
        }
        std::cerr << " OK\n";
        seqan::CharString faiPath = options.inputFilename;
        append(faiPath, ".fai");
        std::cerr << "Reference Index       " << faiPath << " ...";
        if (write(genome.seqFaiIndex, toCString(faiPath)) != 0)
        {
            std::cerr << "Could not write FAI index we just built.\n";
            return 1;
        }
        std::cerr << " OK (" << numSeqs(genome.seqFaiIndex) << " seqs)\n";
    }
    else
    {
        std::cerr << " OK (" << numSeqs(genome.seqFaiIndex) << " seqs)\n";
    }

    // Loading methylation levels.
    if (options.bsSimEnabled)
    {
        std::cerr << "Loading Meth. Level FAI " << options.methLevelsInFasta << " ...";
        if (read(genome.lvlFaiIndex, toCString(options.methLevelsInFasta)) != 0)
        {
            std::cerr << " FAILED (not fatal, we can just build it)\n";
            std::cerr << "Building FAI          " << options.methLevelsInFasta << ".fai ...";
            if (build(genome.lvlFaiIndex, toCString(options.methLevelsInFasta)) != 0)
            {
                std::cerr << "Could not build FAI index.\n";
                return 1;
            }
            std::cerr << " OK\n";
            seqan::CharString faiPath = options.methLevelsInFasta;
            append(faiPath, ".fai");
            std::cerr << "Reference Index       " << faiPath << " ...";
            if (write(genome.lvlFaiIndex, toCString(faiPath)) != 0)
            {
                std::cerr << "Could not write FAI index we just built.\n";
                return 1;
            }
            std::cerr << " OK (" << length(genome.lvlFaiIndex) << " seqs)\n";
        }
        else
        {
            std::cerr << " OK (" << length(genome.lvlFaiIndex) << " seqs)\n";
        }
    }

    // Open output file.
    std::cerr << "Output File             " << options.outputFilename << " ...";
    seqan::SequenceStream outStream;
    open(outStream, toCString(options.outputFilename), seqan::SequenceStream::WRITE, seqan::SequenceStream::FASTA);
    if (!isGood(outStream))
    {
        std::cerr << "\nERROR: Could not open output file " << options.outputFilename << "\n";
        return 1;
    }
    std::cerr << " OK\n";

    // Perform genome simulation.
    std::cerr << "\n__SIMULATING FRAGMENTS_______________________________________________________\n"
              << "\n";

    TRng rng(options.seed);

    FragmentOptions fragOptions;
    fragOptions.model = (options.distribution == MasonFragmentsOptions::NORMAL) ?
            FragmentOptions::NORMAL : FragmentOptions::UNIFORM;
    fragOptions.numFragments = options.numFragments;
    fragOptions.minFragmentSize = options.minFragmentSize;
    fragOptions.maxFragmentSize = options.maxFragmentSize;
    fragOptions.meanFragmentSize = options.meanFragmentSize;
    fragOptions.stdDevFragmentSize = options.stdDevFragmentSize;
    fragOptions.embedSamplingInfo = options.embedSamplingInfo;

    fragOptions.bsSimEnabled = options.bsSimEnabled;
    fragOptions.bsConversionRate = options.bsConversionRate;
    fragOptions.bsProtocol = (options.bsProtocol == MasonFragmentsOptions::DIRECTIONAL) ?
            FragmentOptions::DIRECTIONAL : FragmentOptions::UNDIRECTIONAL;

    FragmentSimulator fragSim(rng, outStream, genome, fragOptions);

    std::cerr << "Simulating Fragments ...";
    fragSim.simulate();
    std::cerr << " OK\n";

    std::cerr << "\nDone.\n";    
    return 0;
}
