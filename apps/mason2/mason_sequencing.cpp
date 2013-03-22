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
// Simulate sequencing process on fragments.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <seqan/arg_parse.h>

#include "sequencing.h"

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class MasonSequencingOptions
// --------------------------------------------------------------------------

struct MasonSequencingOptions
{
    // The sequencing technology to simulate.
    enum Technology
    {
        ILLUMINA,
        SANGER,
        ROCHE_454
    };

    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The input file name.
    seqan::CharString inputFilename;
    // The output file name.
    seqan::CharString outputFilename;
    // The output file name for the right reads.
    seqan::CharString outputFilenameRight;

    // The seed to use for the RNG.
    int seed;

    // The selected technology to simulate.
    Technology technology;

    // The SequencingOptions subclass to use for the configuration.
    std::auto_ptr<SequencingOptions> optionsImpl;

    MasonSequencingOptions() : verbosity(1), seed(0), technology(ILLUMINA)
    {}

    // Set technology by creating optionsImpl.
    void initTechnology(seqan::CharString const & str)
    {
        if (str == "illumina")
        {
            technology = ILLUMINA;
            optionsImpl.reset(new IlluminaSequencingOptions());
        }
        else if (str == "454")
        {
            technology = ROCHE_454;
            optionsImpl.reset(new Roche454SequencingOptions());
        }
        else
        {
            technology = SANGER;
            optionsImpl.reset(new SangerSequencingOptions());
        }
    }

    void print(std::ostream & out)
    {
        out << "  INPUT FILENAME                \t" << inputFilename << "\n"
            << "  OUTPUT FILENAME               \t" << outputFilename << "\n"
            << "  OUTPUT FILENAME (RIGHT MATES) \t" << outputFilenameRight << "\n"
            << "  SEED                          \t" << seed << "\n"
            << "\n";
        if (optionsImpl.get())
            optionsImpl->print(out);
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(MasonSequencingOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("mason_seqencing");
    // Set short description, version, and date.
    setShortDescription(parser, "Sequencing Simulation");
    setVersion(parser, "2.1");
    setDate(parser, "March 2013");
    setCategory(parser, "Simulators");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-i\\fP \\fIIN.fa\\fP \\fB-o\\fP \\fIOUT.{fa,fq}\\fP "
                 "[\\fB-r\\fP \\fIOUT2.{fa,fq}\\fP]");
    addDescription(parser, "Given a FASTA file with fragments, simulate sequencing thereof.");

    // -----------------------------------------------------------------------
    // Global Options Options
    // -----------------------------------------------------------------------

    // We require one argument.
    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // -----------------------------------------------------------------------
    // Input / Output Options
    // -----------------------------------------------------------------------

    addSection(parser, "Input / Output Options");

    addOption(parser, seqan::ArgParseOption("i", "in-file", "Input file.",
                                            seqan::ArgParseOption::INPUTFILE, "FILE"));
    setValidValues(parser, "in-file", "fa fasta");
    setRequired(parser, "in-file");

    addOption(parser, seqan::ArgParseOption("o", "out-file", "Output file, when doing paired read simulation, this is "
                                            "the output file for the left reads.  Paired read simulation is activated "
                                            "by giving \\fB-r\\fP/\\fB--out-file-right\\fP",
                                            seqan::ArgParseOption::OUTPUTFILE, "FILE"));
    setValidValues(parser, "out-file", "fa fasta fq fastq");
    setRequired(parser, "out-file");

    addOption(parser, seqan::ArgParseOption("r", "out-file-right", "Output file for right reads when doing mate-pair "
                                            "simulation, this is the output file for the right reads.  Paired read "
                                            "simulation is activated by giving \\fB-r\\fP/\\fB--out-file-right\\fP.",
                                            seqan::ArgParseOption::OUTPUTFILE, "FILE"));
    setValidValues(parser, "out-file-right", "fa fasta fq fastq");

    // -----------------------------------------------------------------------
    // Global / Generic Sequencing Options
    // -----------------------------------------------------------------------

    addSection(parser, "Global Simulation Options");

    addOption(parser, seqan::ArgParseOption("", "seed", "Seed to use for random number generator.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "seed", "0");

    addOption(parser, seqan::ArgParseOption("t", "technology", "Set sequencing technology to simulate.",
                                            seqan::ArgParseOption::STRING, "TECHNOLOGY"));
    setValidValues(parser, "technology", "illumina 454 sanger");
    setDefaultValue(parser, "technology", "illumina");

    addOption(parser, seqan::ArgParseOption("", "mate-orientation", "Orientation for paired reads.  See section Read "
                                            "Orientation below."));
    setValidValues(parser, "mate-orientation", "FR RF FF FF2");
    setDefaultValue(parser, "mate-orientation", "FR");

    addOption(parser, seqan::ArgParseOption("", "strands", "Strands to simulate from, only applicable to paired "
                                            "sequencing simulation.", seqan::ArgParseOption::STRING, "STRAND"));
    setValidValues(parser, "strands", "forward reverse both");
    setDefaultValue(parser, "strands", "both");

    // -----------------------------------------------------------------------
    // Illumina Sequencing Options
    // -----------------------------------------------------------------------

    addSection(parser, "Illumina Options");

    addOption(parser, seqan::ArgParseOption("", "illumina-read-length", "Read length for Illumina simulation.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "illumina-read-length", "1");
    setDefaultValue(parser, "illumina-read-length", "100");

    addOption(parser, seqan::ArgParseOption("", "illumina-prob-insert",
                                            "Insert per-base probability for insertion in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "illumina-prob-insert", "0");
    setMaxValue(parser, "illumina-prob-insert", "1");
    setDefaultValue(parser, "illumina-prob-insert", "0.001");

    addOption(parser, seqan::ArgParseOption("", "illumina-prob-deletion",
                                            "Insert per-base probability for deletion in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "illumina-prob-insert", "0");
    setMaxValue(parser, "illumina-prob-insert", "1");
    setDefaultValue(parser, "illumina-prob-insert", "0.001");

    addOption(parser, seqan::ArgParseOption("", "illumina-prob-mismatch-scale",
                                            "Scaling factor for Illumina mismatch probability.",
                                            seqan::ArgParseOption::DOUBLE, "FACTOR"));
    setMinValue(parser, "illumina-prob-mismatch-scale", "0");
    setDefaultValue(parser, "illumina-prob-mismatch-scale", "1.0");

    addOption(parser, seqan::ArgParseOption("", "illumina-prob-mismatch",
                                            "Average per-base mismatch probability in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "illumina-prob-mismatch", "0.0");
    setMaxValue(parser, "illumina-prob-mismatch", "1.0");
    setDefaultValue(parser, "illumina-prob-mismatch", "0.004");

    addOption(parser, seqan::ArgParseOption("", "illumina-prob-mismatch-begin",
                                            "Per-base mismatch probability of first base in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "illumina-prob-mismatch-begin", "0.0");
    setMaxValue(parser, "illumina-prob-mismatch-begin", "1.0");
    setDefaultValue(parser, "illumina-prob-mismatch-begin", "0.002");

    addOption(parser, seqan::ArgParseOption("", "illumina-prob-mismatch-end",
                                            "Per-base mismatch probability of last base in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "illumina-prob-mismatch-end", "0.0");
    setMaxValue(parser, "illumina-prob-mismatch-end", "1.0");
    setDefaultValue(parser, "illumina-prob-mismatch-end", "0.012");

    addOption(parser, seqan::ArgParseOption("", "illumina-position-raise",
                                            "Point where the error curve raises in relation to read length.",
                                            seqan::ArgParseOption::DOUBLE, "FLOAT"));
    setMinValue(parser, "illumina-position-raise", "0.0");
    setMaxValue(parser, "illumina-position-raise", "1.0");
    setDefaultValue(parser, "illumina-position-raise", "0.66");

    addOption(parser, seqan::ArgParseOption("", "illumina-quality-mean-begin",
                                            "Mean PHRED quality for non-mismatch bases of first base in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-quality-mean-begin", "40.0");

    addOption(parser, seqan::ArgParseOption("", "illumina-quality-mean-end",
                                            "Mean PHRED quality for non-mismatch bases of last base in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-quality-mean-end", "39.5");

    addOption(parser, seqan::ArgParseOption("", "illumina-quality-stddev-begin",
                                            "Standard deviation of PHRED quality for non-mismatch bases of first base "
                                            "in Illumina sequencing.", seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-quality-stddev-begin", "0.05");

    addOption(parser, seqan::ArgParseOption("", "illumina-quality-stddev-end",
                                            "Standard deviation of PHRED quality for non-mismatch bases of last base "
                                            "in Illumina sequencing.", seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-quality-stddev-end", "10.0");

    addOption(parser, seqan::ArgParseOption("", "illumina-mismatch-quality-mean-begin",
                                            "Mean PHRED quality for mismatch bases of first base in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-mismatch-quality-mean-begin", "40.0");

    addOption(parser, seqan::ArgParseOption("", "illumina-mismatch-quality-mean-end",
                                            "Mean PHRED quality for mismatch bases of last base in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-mismatch-quality-mean-end", "39.5");

    addOption(parser, seqan::ArgParseOption("", "illumina-mismatch-quality-stddev-begin",
                                            "Standard deviation of PHRED quality for mismatch bases of first base "
                                            "in Illumina sequencing.", seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-mismatch-quality-stddev-begin", "0.05");

    addOption(parser, seqan::ArgParseOption("", "illumina-mismatch-quality-stddev-end",
                                            "Standard deviation of PHRED quality for mismatch bases of last base "
                                            "in Illumina sequencing.", seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-mismatch-quality-stddev-end", "10.0");

    // -----------------------------------------------------------------------
    // 454 Sequencing Options
    // -----------------------------------------------------------------------

    addSection(parser, "454 Sequencing Options");

    addOption(parser, seqan::ArgParseOption("", "454-read-length-model", "The model to use for sampling the 454 read length.",
                                            seqan::ArgParseOption::STRING, "MODEL"));
    setValidValues(parser, "454-read-length-model", "normal uniform");
    setDefaultValue(parser, "454-read-length-model", "normal");

    addOption(parser, seqan::ArgParseOption("", "454-read-length-min", "The minimal read length when the read length is "
                                            "sampled uniformly.", seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "454-read-length-min", "0");
    setDefaultValue(parser, "454-read-length-min", "10");

    addOption(parser, seqan::ArgParseOption("", "454-read-length-max", "The maximal read length when the read length is "
                                            "sampled uniformly.", seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "454-read-length-max", "0");
    setDefaultValue(parser, "454-read-length-max", "600");

    addOption(parser, seqan::ArgParseOption("", "454-read-length-mean", "The mean read length when the read length is "
                                            "sampled with normal distribution.", seqan::ArgParseOption::DOUBLE, "LEN"));
    setMinValue(parser, "454-read-length-mean", "0");
    setDefaultValue(parser, "454-read-length-mean", "400");

    addOption(parser, seqan::ArgParseOption("", "454-read-length-stddev", "The read length standard deviation when the "
                                            "read length is sampled uniformly.", seqan::ArgParseOption::DOUBLE, "LEN"));
    setMinValue(parser, "454-read-length-stddev", "0");
    setDefaultValue(parser, "454-read-length-stddev", "40");

    addOption(parser, seqan::ArgParseOption("", "454-no-sqrt-in-std-dev", "For error model, if set then "
                                            "(sigma = k * r)) is used, otherwise (sigma = k * sqrt(r))."));

    addOption(parser, seqan::ArgParseOption("", "454-k", "Proportionality factor for calculating the standard "
                                            "deviation proportional to the read length.",
                                            seqan::ArgParseOption::DOUBLE, "FLOAT"));
    setMinValue(parser, "454-k", "0");
    setDefaultValue(parser, "454-k", "0.15");

    addOption(parser, seqan::ArgParseOption("", "454-background-noise-mean", "Mean of lognormal distribution to use for "
                                            "the noise.", seqan::ArgParseOption::DOUBLE, "MEAN"));
    setMinValue(parser, "454-background-noise-mean", "0");
    setDefaultValue(parser, "454-background-noise-mean", "0.23");

    addOption(parser, seqan::ArgParseOption("", "454-background-noise-stddev", "Standard deviation of lognormal "
                                            "distribution to use for the noise.", seqan::ArgParseOption::DOUBLE,
                                            "SIGMA"));
    setMinValue(parser, "454-background-noise-stddev", "0");
    setDefaultValue(parser, "454-background-noise-stddev", "0.15");

    // -----------------------------------------------------------------------
    // Sanger Sequencing Options
    // -----------------------------------------------------------------------
    
    addSection(parser, "Sanger Sequencing Options");

    addOption(parser, seqan::ArgParseOption("", "sanger-read-length-model", "The model to use for sampling the Sanger "
                                            "read length.", seqan::ArgParseOption::STRING, "MODEL"));
    setValidValues(parser, "sanger-read-length-model", "normal uniform");
    setDefaultValue(parser, "sanger-read-length-model", "normal");

    addOption(parser, seqan::ArgParseOption("", "sanger-read-length-min", "The minimal read length when the read length is "
                                            "sampled uniformly.", seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "sanger-read-length-min", "0");
    setDefaultValue(parser, "sanger-read-length-min", "400");

    addOption(parser, seqan::ArgParseOption("", "sanger-read-length-max", "The maximal read length when the read length is "
                                            "sampled uniformly.", seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "sanger-read-length-max", "0");
    setDefaultValue(parser, "sanger-read-length-max", "600");

    addOption(parser, seqan::ArgParseOption("", "sanger-read-length-mean", "The mean read length when the read length is "
                                            "sampled with normal distribution.", seqan::ArgParseOption::DOUBLE, "LEN"));
    setMinValue(parser, "sanger-read-length-mean", "0");
    setDefaultValue(parser, "sanger-read-length-mean", "400");

    addOption(parser, seqan::ArgParseOption("", "sanger-read-length-stddev", "The read length standard deviation when the "
                                            "read length is sampled uniformly.", seqan::ArgParseOption::DOUBLE, "LEN"));
    setMinValue(parser, "sanger-read-length-stddev", "0");
    setDefaultValue(parser, "sanger-read-length-stddev", "40");

    addOption(parser, seqan::ArgParseOption("", "sanger-prob-mismatch-scale",
                                            "Scaling factor for Sanger mismatch probability.",
                                            seqan::ArgParseOption::DOUBLE, "FACTOR"));
    setMinValue(parser, "sanger-prob-mismatch-scale", "0");
    setDefaultValue(parser, "sanger-prob-mismatch-scale", "1.0");

    addOption(parser, seqan::ArgParseOption("", "sanger-prob-mismatch-begin",
                                            "Per-base mismatch probability of first base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "sanger-prob-mismatch-begin", "0.0");
    setMaxValue(parser, "sanger-prob-mismatch-begin", "1.0");
    setDefaultValue(parser, "sanger-prob-mismatch-begin", "0.005");

    addOption(parser, seqan::ArgParseOption("", "sanger-prob-mismatch-end",
                                            "Per-base mismatch probability of last base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "sanger-prob-mismatch-end", "0.0");
    setMaxValue(parser, "sanger-prob-mismatch-end", "1.0");
    setDefaultValue(parser, "sanger-prob-mismatch-end", "0.001");

    addOption(parser, seqan::ArgParseOption("", "sanger-prob-insertion-begin",
                                            "Per-base insertion probability of first base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "sanger-prob-insertion-begin", "0.0");
    setMaxValue(parser, "sanger-prob-insertion-begin", "1.0");
    setDefaultValue(parser, "sanger-prob-insertion-begin", "0.0025");

    addOption(parser, seqan::ArgParseOption("", "sanger-prob-insertion-end",
                                            "Per-base insertion probability of last base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "sanger-prob-insertion-end", "0.0");
    setMaxValue(parser, "sanger-prob-insertion-end", "1.0");
    setDefaultValue(parser, "sanger-prob-insertion-end", "0.005");

    addOption(parser, seqan::ArgParseOption("", "sanger-prob-deletion-begin",
                                            "Per-base deletion probability of first base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "sanger-prob-deletion-begin", "0.0");
    setMaxValue(parser, "sanger-prob-deletion-begin", "1.0");
    setDefaultValue(parser, "sanger-prob-deletion-begin", "0.0025");

    addOption(parser, seqan::ArgParseOption("", "sanger-prob-deletion-end",
                                            "Per-base deletion probability of last base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "sanger-prob-deletion-end", "0.0");
    setMaxValue(parser, "sanger-prob-deletion-end", "1.0");
    setDefaultValue(parser, "sanger-prob-deletion-end", "0.005");

    // -----------------------------------------------------------------------
    // Read Orientation Information
    // -----------------------------------------------------------------------

    // Add Examples Section.
    addTextSection(parser, "Read Orientation");
    addText(parser,
            "You can use the \\fB--mate-orientation\\fP to set the relative orientation when doing paired-end "
            "sequencing.  The valid values are given in the following.");
    addListItem(parser, "FR",
                "Reads are inward-facing, the same as Illumina paired-end reads: R1 --> <-- R2.");
    addListItem(parser, "RF",
                "Reads are outward-facing, the same as Illumina mate-pair reads: R1 <-- --> R2.");
    addListItem(parser, "FF",
                "Reads are on the same strand: R1 --> --> R2.");
    addListItem(parser, "FF2",
                "Reads are on the same strand but the \"right\" reads are sequenced to the left of the \"left\" reads, "
                "same as 454 paired: R2 --> --> R1.");

    // -----------------------------------------------------------------------
    // Examples
    // -----------------------------------------------------------------------

    // Add Examples Section.
    addTextSection(parser, "Examples");

    addListItem(parser, "\\fBmason_sequencing\\fP \\fB-i\\fP \\fIfragments.fa\\fP \\fB-o\\fP \\fIreads.fa\\fP",
                "Given the fragments from \\fIfragments.fa\\fP, perform single-end sequencing and write the "
                "results to the FASTA file \\fIreads.fa\\fP.  Illumina sequencing is simulated using its default "
                "settings.");

    addListItem(parser, "\\fBmason_sequencing\\fP \\fB-i\\fP \\fIfragments.fa\\fP \\fB-o\\fP \\fIreads_1.fq\\fP "
                "\\fB-r\\fP \\fIreads_2.fq\\fP",
                "Given the fragments from \\fIfragments.fa\\fP, perform paired-end sequencing and write the "
                "results to the FASTQ files \\fIreads_1.fq\\fP and \\fIreads_2.fq\\fP.  Illumina sequencing is "
                "simulated using its default settings.");

    addListItem(parser, "\\fBmason_sequencing\\fP \\fB-i\\fP \\fIfragments.fa\\fP \\fB-o\\fP \\fIreads_1.fq\\fP "
                "\\fB-r\\fP \\fIreads_2.fq\\fP \\fB-t\\fP \\fI454\\fP",
                "Given the fragments from \\fIfragments.fa\\fP, perform paired-end sequencing and write the "
                "results to the FASTQ files \\fIreads_1.fq\\fP and \\fIreads_2.fq\\fP.  454 sequencing is used "
                "with its default settings.");

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
    getOptionValue(options.outputFilenameRight, parser, "out-file-right");
    getOptionValue(options.seed, parser, "seed");

    seqan::CharString tmpTech;
    getOptionValue(tmpTech, parser, "technology");
    options.initTechnology(tmpTech);

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    MasonSequencingOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cerr << "MASON SEQUENCING SIMULATOR\n"
              << "==========================\n\n";
    
    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cerr << "__OPTIONS____________________________________________________________________\n"
                  << '\n';
        options.print(std::cerr);
    }

    std::cerr << "\n__PREPARATION________________________________________________________________\n"
              << "\n";

    // // Open reference FAI index.
    // std::cerr << "Loading reference     " << options.inputFilename << " ...";
    // seqan::SequenceStream inGenome(toCString(options.inputFilename));
    // if (!isGood(inGenome))
    // {
    //     std::cerr << " ERROR\n"
    //               << "Could not open " << options.inputFilename << "\n";
    //     return 1;
    // }

    // Genome genome;
    // if (readAll(genome.ids, genome.seqs, inGenome))
    // {
    //     std::cerr << " ERROR\n"
    //               << "Could not load genome.\n";
    //     return 1;
    // }
    // genome.trimIds();
    // std::cerr << " OK\n";

    // // Open output file.
    // std::cerr << "\nOutput File           " << options.outputFilename << " ...";
    // seqan::SequenceStream outStream;
    // open(outStream, toCString(options.outputFilename), seqan::SequenceStream::WRITE, seqan::SequenceStream::FASTA);
    // if (!isGood(outStream))
    // {
    //     std::cerr << "\nERROR: Could not open output file " << options.outputFilename << "\n";
    //     return 1;
    // }
    // std::cerr << " OK\n";

    // // Perform genome simulation.
    // std::cerr << "\n__SIMULATING FRAGMENTS_______________________________________________________\n"
    //           << "\n";

    // TRng rng(options.seed);

    // FragmentOptions fragOptions;
    // fragOptions.model = (options.distribution == MasonFragmentsOptions::NORMAL) ?
    //         FragmentOptions::NORMAL : FragmentOptions::UNIFORM;
    // fragOptions.numFragments = options.numFragments;
    // fragOptions.minFragmentSize = options.minFragmentSize;
    // fragOptions.maxFragmentSize = options.maxFragmentSize;
    // fragOptions.meanFragmentSize = options.meanFragmentSize;
    // fragOptions.stdDevFragmentSize = options.stdDevFragmentSize;
    // fragOptions.embedSamplingInfo = options.embedSamplingInfo;

    // FragmentSimulator fragSim(rng, outStream, genome, fragOptions);

    // std::cerr << "Simulating Fragments ...";
    // fragSim.simulate();
    // std::cerr << " OK\n";

    // std::cerr << "\nDone.\n";    
    return 0;
}
