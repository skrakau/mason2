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

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class MasonFragmentsOption
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
struct MasonFragmentsOptions
{
    // Enum for selecting size distribution.
    enum Distribution
    {
        NORMAL,
        UNIFORM
    };

    // Enum for selecting the BS-seq protocol (directional vs. undirectional).
    enum BSProtocol
    {
        DIRECTIONAL,
        UNDIRECTIONAL
    };

    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // -----------------------------------------------------------------------
    // Input / Output Options
    // -----------------------------------------------------------------------

    // The input file name.
    seqan::CharString inputFilename;
    // The output file name.
    seqan::CharString outputFilename;

    // -----------------------------------------------------------------------
    // Fragment Simulation Options
    // -----------------------------------------------------------------------

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

    // -----------------------------------------------------------------------
    // BS-Seq Simulation Options
    // -----------------------------------------------------------------------

    // Whether or not to BS-seq methylation simulation.  Set to !empty(methLevelsInFasta).
    bool bsSimEnabled;
    // Path to FASTA file with methylation levels.  The levels of 0-100% are encoded in steps of 1.25% in ASCII
    // characters starting with '!'.  '>' is skipped and may not occur in the input.
    seqan::CharString methLevelsInFasta;
    // The rate for unmethylated Cs to become Ts.
    double bsConversionRate;
    // The protocol to use for the simulation.
    BSProtocol bsProtocol;

    MasonFragmentsOptions() :
            verbosity(1), seed(0), numFragments(0), embedSamplingInfo(0), distribution(NORMAL), minFragmentSize(0),
            maxFragmentSize(0), meanFragmentSize(0), stdDevFragmentSize(0), bsSimEnabled(false), bsConversionRate(0),
            bsProtocol(DIRECTIONAL)
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

    // BS-Seq Options
    addSection(parser, "BS-Seq Options");

    addOption(parser, seqan::ArgParseOption("", "meth-levels-in", "Methylation levels FASTA input.  "
                                            "If given, BS-Seq simulation is enabled.",
                                            seqan::ArgParseOption::INPUTFILE, "FILE"));
    setValidValues(parser, "meth-levels-in", "fa fasta");
    setRequired(parser, "meth-levels-in");

    addOption(parser, seqan::ArgParseOption("", "bs-conversion-rate", "Methylation conversion rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "bs-conversion-rate", "0.0");
    setMaxValue(parser, "bs-conversion-rate", "1.0");
    setDefaultValue(parser, "bs-conversion-rate", "0.98");

    addOption(parser, seqan::ArgParseOption("", "bs-protocol", "Methylation protocol.",
                                            seqan::ArgParseOption::STRING, "STR"));
    setValidValues(parser, "bs-protocol", "directional undirectional");
    setDefaultValue(parser, "bs-protocol", "directional");

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBmason_fragments\\fP \\fB-n\\fP 1000 \\fB-i\\fP \\fIgenome.fa\\fP \\fB-o\\fP \\fIfragments.fa\\fP",
                "Simulate 1000 fragments of file \\fIgenome.fa\\fP and write them to \\fIfragments.fa\\fP.  The "
                "fragments will be simulated using the default configuration.");

    // Add BS-Seq Section.
    addTextSection(parser, "BS-Seq Simulation");
    addText(parser,
            "If \\fB--bs-levels-in\\fP is given, methylation levels are loaded for each element in "
            "\\fB--in-file\\fP.  The methylation levels for the TOP strand are expected to have the "
            "name of the contig with the suffix \"/TOP\", the name of the bottom strand is expected "
            "to have the suffix \"/BOT\".");
    addText(parser,
            "The sequence characters encode the values 0..80 with an offset of 33 (\"!\" = 0).  "
            "The character \">\" is skipped.");
    addText(parser,
            "When simulating, we simulate for each C whether it is methylated.  Methylated Cs are "
            "kept intact.  Unmethylated Cs are converted into Ts with a probability given by "
            "\\fB--bs-conversion-rate\\fP.");

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

    getOptionValue(options.methLevelsInFasta, parser, "meth-levels-in");
    options.bsSimEnabled = !empty(options.methLevelsInFasta);
    getOptionValue(options.bsConversionRate, parser, "bs-conversion-rate");

    getOptionValue(tmp, parser, "bs-protocol");
    options.bsProtocol = (tmp == "directional") ?
            MasonFragmentsOptions::DIRECTIONAL : MasonFragmentsOptions::UNDIRECTIONAL;

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

char const * yesNo(bool b)
{
    return b ? "YES" : "NO";
}

char const * bsProtocolStr(MasonFragmentsOptions::BSProtocol p)
{
    switch (p)
    {
        case MasonFragmentsOptions::DIRECTIONAL:
            return "DIRECTIONAL";
        default:
            return "UNDIRECTIONAL";
    }
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
