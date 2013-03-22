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

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <seqan/arg_parse.h>

#include "fragment_generation.h"

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class MasonFragmentsOption
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct MasonFragmentsOptions
{
    // Enum for selecting size distribution.
    enum Distribution
    {
        NORMAL,
        UNIFORM
    };

    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The input file name.
    seqan::CharString inputFilename;
    // The output file name.
    seqan::CharString outputFilename;

    // The seed to use for the RNG.
    int seed;
    // The number of fragments.
    int numFragments;
    // Whether or not to embed sampling information in FASTA header.
    int embedSamplingInfo;

    // The distribution to use for fragment size simulation.
    Distribution distribution;

    // The minimal fragment size, used when distribution is UNIFORM.
    int minFragmentSize;
    // The maximal fragment size, used when distribution is UNIFORM.
    int maxFragmentSize;
    // The mean fragment size, used when distribution is NORMAL.
    int meanFragmentSize;
    // The fragment size standard deviation, used when distribution is NORMAL.
    int stdDevFragmentSize;

    MasonFragmentsOptions() :
            verbosity(1), seed(0), numFragments(0), embedSamplingInfo(0), distribution(NORMAL), minFragmentSize(0),
            maxFragmentSize(0), meanFragmentSize(0), stdDevFragmentSize(0)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(MasonFragmentsOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("mason_fragments");
    // Set short description, version, and date.
    setShortDescription(parser, "Random Fragment Simulation");
    setVersion(parser, "2.1");
    setDate(parser, "March 2013");
    setCategory(parser, "Simulators");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-n\\fP \\fINUM\\fP \\fB-i\\fP \\fIIN.fa\\fP \\fB-o\\fP \\fIOUT.fa\\fP");
    addDescription(parser, "Generate \\fINUM\\fP fragments from the \\fIIN.fa\\fP file to the \\fIOUT.fa\\fP file.");

    // We require one argument.
    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    addSection(parser, "Fragment Simulation Options");

    addOption(parser, seqan::ArgParseOption("s", "seed", "The seed to use for the random number generator.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "seed", "0");

    addOption(parser, seqan::ArgParseOption("n", "num-fragments", "The number of fragments to simulate.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setMinValue(parser, "num-fragments", "1");
    setRequired(parser, "num-fragments");

    addOption(parser, seqan::ArgParseOption("", "embed-sampling-info", "Whether or not to embed sampling information "
                                            "in FASTA id."));

    addOption(parser, seqan::ArgParseOption("", "distribution", "The size distribution.",
                                            seqan::ArgParseOption::STRING, "DIST"));
    setValidValues(parser, "distribution", "uniform normal");
    setDefaultValue(parser, "distribution", "normal");

    addOption(parser, seqan::ArgParseOption("", "min-size", "The minimal fragment size, used if "
                                            "\\fB--distribution\\fP is \\fIuniform\\fP.",
                                            seqan::ArgParseOption::INTEGER, "LENGTH"));
    setMinValue(parser, "min-size", "1");
    setDefaultValue(parser, "min-size", "150");

    addOption(parser, seqan::ArgParseOption("", "max-size", "The maximal fragment size, used if "
                                            "\\fB--distribution\\fP is \\fIuniform\\fP.",
                                            seqan::ArgParseOption::INTEGER, "LENGTH"));
    setMinValue(parser, "max-size", "1");
    setDefaultValue(parser, "max-size", "450");

    addOption(parser, seqan::ArgParseOption("", "mean-size", "The mean fragment size, used if "
                                            "\\fB--distribution\\fP is \\fInormal\\fP.",
                                            seqan::ArgParseOption::INTEGER, "LENGTH"));
    setMinValue(parser, "mean-size", "1");
    setDefaultValue(parser, "mean-size", "300");

    addOption(parser, seqan::ArgParseOption("", "size-stddev", "The fragment size standard deviation, used if "
                                            "\\fB--distribution\\fP is \\fInormal\\fP.",
                                            seqan::ArgParseOption::INTEGER, "LENGTH"));
    setMinValue(parser, "size-stddev", "1");
    setDefaultValue(parser, "size-stddev", "40");

    addSection(parser, "Input / Output Options");

    addOption(parser, seqan::ArgParseOption("i", "in-file", "Input file.",
                                            seqan::ArgParseOption::INPUTFILE, "FILE"));
    setValidValues(parser, "in-file", "fa fasta");
    setRequired(parser, "in-file");

    addOption(parser, seqan::ArgParseOption("o", "out-file", "Output file.",
                                            seqan::ArgParseOption::OUTPUTFILE, "FILE"));
    setValidValues(parser, "out-file", "fa fasta");
    setRequired(parser, "out-file");

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBmason_genome\\fP \\fB-l\\fP 1000 \\fB-l\\fP 4000 \\fB-o\\fP \\fIgenome.fa\\fP",
                "Simulate a genome with two contigs of lengths 1000 and 4000 and write it to genome.fa.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValue(options.inputFilename, parser, "in-file");
    getOptionValue(options.outputFilename, parser, "out-file");

    getOptionValue(options.seed, parser, "seed");
    getOptionValue(options.numFragments, parser, "num-fragments");
    options.embedSamplingInfo = isSet(parser, "embed-sampling-info");

    seqan::CharString tmp;
    getOptionValue(tmp, parser, "distribution");
    if (tmp == "normal")
        options.distribution = MasonFragmentsOptions::NORMAL;
    else
        options.distribution = MasonFragmentsOptions::UNIFORM;
    
    getOptionValue(options.minFragmentSize, parser, "min-size");
    getOptionValue(options.maxFragmentSize, parser, "max-size");
    getOptionValue(options.meanFragmentSize, parser, "mean-size");
    getOptionValue(options.stdDevFragmentSize, parser, "size-stddev");

    if (options.minFragmentSize > options.maxFragmentSize)
    {
        std::cerr << "ERROR: --min-size must not be greater than --max-size.";
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

char const * distributionStr(MasonFragmentsOptions::Distribution d)
{
    return (d == MasonFragmentsOptions::NORMAL) ? "NORMAL" : "UNIFORM";
}

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
                  << "SIZE STD DEV  \t" << options.stdDevFragmentSize << "\n\n";
    }

    std::cerr << "__PREPARATION________________________________________________________________\n"
              << "\n";

    // Open reference FAI index.
    std::cerr << "Loading reference     " << options.inputFilename << " ...";
    seqan::SequenceStream inGenome(toCString(options.inputFilename));
    if (!isGood(inGenome))
    {
        std::cerr << " ERROR\n"
                  << "Could not open " << options.inputFilename << "\n";
        return 1;
    }

    Genome genome;
    if (readAll(genome.ids, genome.seqs, inGenome))
    {
        std::cerr << " ERROR\n"
                  << "Could not load genome.\n";
        return 1;
    }
    genome.trimIds();
    std::cerr << " OK\n";

    // Open output file.
    std::cerr << "\nOutput File           " << options.outputFilename << " ...";
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

    FragmentSimulator fragSim(rng, outStream, genome, fragOptions);

    std::cerr << "Simulating Fragments ...";
    fragSim.simulate();
    std::cerr << " OK\n";

    std::cerr << "\nDone.\n";    
    return 0;
}
