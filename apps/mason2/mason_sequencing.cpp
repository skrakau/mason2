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
public:
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The read name prefix.
    seqan::CharString readNamePrefix;

    // The input file name.
    seqan::CharString inputFilename;
    // The output file name.
    seqan::CharString outputFilename;
    // The output file name for the right reads.
    seqan::CharString outputFilenameRight;

    // The seed to use for the RNG.
    int seed;

    // The selected technology to simulate.
    SequencingSimulatorFactory::Technology technology;

    // The SequencingOptions subclass to use for the configuration.
    std::auto_ptr<SequencingOptions> optionsImpl;

    MasonSequencingOptions() : verbosity(1), seed(0), technology(SequencingSimulatorFactory::ILLUMINA)
    {}

    // Set technology by creating optionsImpl.
    void initTechnology(seqan::CharString const & str)
    {
        if (str == "illumina")
        {
            technology = SequencingSimulatorFactory::ILLUMINA;
            optionsImpl.reset(new IlluminaSequencingOptions());
        }
        else if (str == "454")
        {
            technology = SequencingSimulatorFactory::ROCHE_454;
            optionsImpl.reset(new Roche454SequencingOptions());
        }
        else
        {
            technology = SequencingSimulatorFactory::SANGER;
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
                                            "Orientation below.", seqan::ArgParseOption::STRING, "ORIENTATION"));
    setValidValues(parser, "mate-orientation", "FR RF FF FF2");
    setDefaultValue(parser, "mate-orientation", "FR");

    addOption(parser, seqan::ArgParseOption("", "strands", "Strands to simulate from, only applicable to paired "
                                            "sequencing simulation.", seqan::ArgParseOption::STRING, "STRAND"));
    setValidValues(parser, "strands", "forward reverse both");
    setDefaultValue(parser, "strands", "both");

    addOption(parser, seqan::ArgParseOption("", "embed-read-info", "Whether or not to embed information abou the "
                                            "sequencing process in the read meta data."));

    addOption(parser, seqan::ArgParseOption("", "read-name-prefix", "Prefix for the read name.",
                                            seqan::ArgParseOption::STRING, "NAME"));
    setDefaultValue(parser, "read-name-prefix", "reads.");

    // -----------------------------------------------------------------------
    // Illumina Sequencing Options
    // -----------------------------------------------------------------------

    addSection(parser, "Illumina Options");

    addOption(parser, seqan::ArgParseOption("", "illumina-read-length", "Read length for Illumina simulation.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "illumina-read-length", "1");
    setDefaultValue(parser, "illumina-read-length", "100");

    addOption(parser, seqan::ArgParseOption("", "illumina-error-profile-file",
                                            "Path to file with Illumina error profile.  The file must be a text file "
                                            "with floating point numbers separated by space, each giving a positional "
                                            "error rate.",
                                            seqan::ArgParseOption::INPUTFILE, "FILE"));
    setValidValues(parser, "illumina-error-profile-file", "txt");

    addOption(parser, seqan::ArgParseOption("", "illumina-prob-insert",
                                            "Insert per-base probability for insertion in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "illumina-prob-insert", "0");
    setMaxValue(parser, "illumina-prob-insert", "1");
    setDefaultValue(parser, "illumina-prob-insert", "0.001");

    addOption(parser, seqan::ArgParseOption("", "illumina-prob-deletion",
                                            "Insert per-base probability for deletion in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "illumina-prob-deletion", "0");
    setMaxValue(parser, "illumina-prob-deletion", "1");
    setDefaultValue(parser, "illumina-prob-deletion", "0.001");

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
                                            "read length is sampled with normal distribution.",
                                            seqan::ArgParseOption::DOUBLE, "LEN"));
    setMinValue(parser, "454-read-length-stddev", "0");
    setDefaultValue(parser, "454-read-length-stddev", "40");

    addOption(parser, seqan::ArgParseOption("", "454-no-sqrt-in-std-dev", "For error model, if set then "
                                            "(sigma = k * r)) is used, otherwise (sigma = k * sqrt(r))."));

    addOption(parser, seqan::ArgParseOption("", "454-proportionality-factor", "Proportionality factor for calculating the standard "
                                            "deviation proportional to the read length.",
                                            seqan::ArgParseOption::DOUBLE, "FLOAT"));
    setMinValue(parser, "454-proportionality-factor", "0");
    setDefaultValue(parser, "454-proportionality-factor", "0.15");

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

    addOption(parser, seqan::ArgParseOption("", "sanger-read-length-error", "The read length standard deviation when the "
                                            "read length is sampled uniformly.", seqan::ArgParseOption::DOUBLE, "LEN"));
    setMinValue(parser, "sanger-read-length-error", "0");
    setDefaultValue(parser, "sanger-read-length-error", "40");

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

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-match-start-mean",
                                            "Mean PHRED quality for non-mismatch bases of first base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-match-start-mean", "40.0");

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-match-end-mean",
                                            "Mean PHRED quality for non-mismatch bases of last base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-match-end-mean", "39.5");

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-match-start-stddev",
                                            "Mean PHRED quality for non-mismatch bases of first base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-match-start-stddev", "0.1");

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-match-end-stddev",
                                            "Mean PHRED quality for non-mismatch bases of last base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-match-end-stddev", "2");

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-error-start-mean",
                                            "Mean PHRED quality for errorneous bases of first base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-error-start-mean", "30");

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-error-end-mean",
                                            "Mean PHRED quality for errorneous bases of last base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-error-end-mean", "20");

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-error-start-stddev",
                                            "Mean PHRED quality for errorneous bases of first base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-error-start-stddev", "2");

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-error-end-stddev",
                                            "Mean PHRED quality for errorneous bases of last base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-error-end-stddev", "5");

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
    getOptionValue(options.readNamePrefix, parser, "read-name-prefix");

    seqan::CharString tmpTech;
    getOptionValue(tmpTech, parser, "technology");
    options.initTechnology(tmpTech);

    seqan::CharString tmpOrientation;
    getOptionValue(tmpOrientation, parser, "mate-orientation");
    if (tmpOrientation == "FR")
        options.optionsImpl->mateOrientation = SequencingOptions::FORWARD_REVERSE;
    else if (tmpOrientation == "RF")
        options.optionsImpl->mateOrientation = SequencingOptions::REVERSE_FORWARD;
    else if (tmpOrientation == "FF")
        options.optionsImpl->mateOrientation = SequencingOptions::FORWARD_FORWARD;
    else if (tmpOrientation == "FF2")
        options.optionsImpl->mateOrientation = SequencingOptions::FORWARD_FORWARD2;

    options.optionsImpl->verbosity = options.verbosity;
    options.optionsImpl->embedReadInfo = isSet(parser, "embed-read-info");

    if (tmpTech == "illumina")
    {
        IlluminaSequencingOptions * ptr = static_cast<IlluminaSequencingOptions *>(options.optionsImpl.get());

        getOptionValue(ptr->readLength, parser, "illumina-read-length");

        getOptionValue(ptr->probabilityMismatchFile, parser, "illumina-error-profile-file");

        getOptionValue(ptr->probabilityInsert, parser, "illumina-prob-insert");
        getOptionValue(ptr->probabilityDelete, parser, "illumina-prob-deletion");

        getOptionValue(ptr->probabilityMismatchScale, parser, "illumina-prob-mismatch-scale");

        getOptionValue(ptr->probabilityMismatch, parser, "illumina-prob-mismatch");
        getOptionValue(ptr->probabilityMismatchBegin, parser, "illumina-prob-mismatch-begin");
        getOptionValue(ptr->probabilityMismatchEnd, parser, "illumina-prob-mismatch-end");
        getOptionValue(ptr->positionRaise, parser, "illumina-position-raise"); 

        getOptionValue(ptr->meanQualityBegin, parser, "illumina-quality-mean-begin");
        getOptionValue(ptr->meanQualityEnd, parser, "illumina-quality-mean-end");
        getOptionValue(ptr->stdDevQualityBegin, parser, "illumina-quality-stddev-begin");
        getOptionValue(ptr->stdDevQualityEnd, parser, "illumina-quality-stddev-end");

        getOptionValue(ptr->meanMismatchQualityBegin, parser, "illumina-mismatch-quality-mean-begin");
        getOptionValue(ptr->meanMismatchQualityEnd, parser, "illumina-mismatch-quality-mean-end");
        getOptionValue(ptr->stdDevMismatchQualityBegin, parser, "illumina-mismatch-quality-stddev-begin");
        getOptionValue(ptr->stdDevMismatchQualityEnd, parser, "illumina-mismatch-quality-stddev-end");
    }
    else if (tmpTech == "454")
    {
        Roche454SequencingOptions * ptr = static_cast<Roche454SequencingOptions *>(options.optionsImpl.get());

        seqan::CharString tmp;
        getOptionValue(tmp, parser, "454-read-length-model");
        ptr->lengthModel = (tmp == "uniform") ? Roche454SequencingOptions::UNIFORM : Roche454SequencingOptions::NORMAL;

        getOptionValue(ptr->minReadLength, parser, "454-read-length-min");
        getOptionValue(ptr->maxReadLength, parser, "454-read-length-max");
        getOptionValue(ptr->meanReadLength, parser, "454-read-length-mean");
        getOptionValue(ptr->stdDevReadLength, parser, "454-read-length-stddev");
        ptr->sqrtInStdDev = !isSet(parser, "454-no-sqrt-in-std-dev");
        getOptionValue(ptr->k, parser, "454-proportionality-factor");
        getOptionValue(ptr->backgroundNoiseMean, parser, "454-background-noise-mean");
        getOptionValue(ptr->backgroundNoiseStdDev, parser, "454-background-noise-stddev");
    }
    else  // tmpTech == "sanger"
    {
        SangerSequencingOptions * ptr = static_cast<SangerSequencingOptions *>(options.optionsImpl.get());

        seqan::CharString tmp;
        getOptionValue(tmp, parser, "sanger-read-length-model");
        ptr->readLengthIsUniform = (tmp == "uniform");

        getOptionValue(ptr->readLengthMean, parser, "sanger-read-length-mean");
        getOptionValue(ptr->readLengthError, parser, "sanger-read-length-error");
        getOptionValue(ptr->readLengthMin, parser, "sanger-read-length-min");
        getOptionValue(ptr->readLengthMax, parser, "sanger-read-length-max");
        
        getOptionValue(ptr->probabilityMismatchBegin, parser, "sanger-prob-mismatch-begin");
        getOptionValue(ptr->probabilityMismatchEnd, parser, "sanger-prob-mismatch-end");
        getOptionValue(ptr->probabilityInsertBegin, parser, "sanger-prob-insertion-begin");
        getOptionValue(ptr->probabilityInsertEnd, parser, "sanger-prob-insertion-end");
        getOptionValue(ptr->probabilityDeleteBegin, parser, "sanger-prob-deletion-begin");
        getOptionValue(ptr->probabilityDeleteEnd, parser, "sanger-prob-deletion-end");

        getOptionValue(ptr->qualityMatchStartMean, parser, "sanger-quality-match-start-mean");
        getOptionValue(ptr->qualityMatchEndMean, parser, "sanger-quality-match-end-mean");
        getOptionValue(ptr->qualityMatchStartStdDev, parser, "sanger-quality-match-start-stddev");
        getOptionValue(ptr->qualityMatchEndStdDev, parser, "sanger-quality-match-end-stddev");
        getOptionValue(ptr->qualityErrorStartMean, parser, "sanger-quality-error-start-mean");
        getOptionValue(ptr->qualityErrorEndMean, parser, "sanger-quality-error-end-mean");
        getOptionValue(ptr->qualityErrorStartStdDev, parser, "sanger-quality-error-start-stddev");
        getOptionValue(ptr->qualityErrorEndStdDev, parser, "sanger-quality-error-end-stddev");
    }

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

    // Open fragments FASTA file.
    std::cerr << "Opening fragments       " << options.inputFilename << " ...";
    seqan::SequenceStream inFragments(toCString(options.inputFilename));
    if (!isGood(inFragments))
    {
        std::cerr << " ERROR\n"
                  << "Could not open " << options.inputFilename << "\n";
        return 1;
    }
    std::cerr << " OK\n";

    // Open reads output file.
    std::cerr << "Opening output file (L) " << options.outputFilename << " ...";
    seqan::SequenceStream outReads(toCString(options.outputFilename), seqan::SequenceStream::WRITE);
    if (!isGood(outReads))
    {
        std::cerr << " ERROR\n"
                  << "Could not open " << options.outputFilename << "\n";
        return 1;
    }
    std::cerr << " OK\n";

    // Open output file for the right reads.
    seqan::SequenceStream outReadsRight;
    if (!empty(options.outputFilenameRight))
    {
        std::cerr << "Opening output file (R) " << options.outputFilenameRight << " ...";
        open(outReadsRight, toCString(options.outputFilenameRight), seqan::SequenceStream::WRITE);
        if (!isGood(outReadsRight))
        {
            std::cerr << " ERROR\n"
                      << "Could not open " << options.outputFilenameRight << "\n";
            return 1;
        }
        std::cerr << " OK\n";
    }

    // Configure output streams to write out each sequence in a single line.
    outReads.outputOptions.lineLength = 0;
    outReadsRight.outputOptions.lineLength = 0;

    // Perform genome simulation.
    std::cerr << "\n__SIMULATING READS___________________________________________________________\n"
              << "\n"
              << "Simulating reads ...";

    TRng rng(options.seed);

    SequencingSimulatorFactory factory(rng, &outReads, &outReadsRight, &inFragments, *options.optionsImpl, options.technology);
    std::auto_ptr<SequencingSimulator> sim = factory.make();

    // Buffers for reading in fragments.
    seqan::CharString fragId;
    seqan::Dna5String fragSeq;
    // Buffers for simulated reads.
    seqan::Dna5String seqL, seqR;
    seqan::CharString qualsL, qualsR;
    // The information for storing the simulation info.
    SequencingSimulationInfo simInfoL, simInfoR;

    // We will use these string streams to generate the read identifier strings.
    std::stringstream ssL, ssR;

    for (unsigned readId = 1; !atEnd(inFragments); ++readId)
    {

        // Reset the string streams.
        ssL.str("");
        ssL.clear();
        ssR.str("");
        ssR.clear();

        // Read fragment to simulate from.
        if (readRecord(fragId, fragSeq, inFragments) != 0)
            return 1;

        if (empty(options.outputFilenameRight))  // Single-end sequencing.
        {
            sim->simulateSingleEnd(seqL, qualsL, simInfoL, infix(fragSeq, 0, length(fragSeq)));
            ssL << options.readNamePrefix << readId;
            if (options.optionsImpl->embedReadInfo)
            {
                ssL << ' ';
                simInfoL.serialize(ssL);
            }
            if (writeRecord(outReads, ssL.str(), seqL, qualsL) != 0)
            {
                std::cerr << "ERROR writing to " << options.outputFilename << "\n";
                return 1;
            }
        }
        else  // Paired sequencing.
        {
            sim->simulatePairedEnd(seqL, qualsL, simInfoL, seqR, qualsR, simInfoR, infix(fragSeq, 0, length(fragSeq)));
            ssL << options.readNamePrefix << readId;
            ssR << options.readNamePrefix << readId;
            if (options.optionsImpl->embedReadInfo)
            {
                ssL << ' ';
                simInfoL.serialize(ssL);
                ssR << ' ';
                simInfoR.serialize(ssR);
            }

            // std::cerr << seqL << "\t" << qualsL << "\n"
            //           << seqR << "\t" << qualsR << "\n\n";

            if (writeRecord(outReads, ssL.str(), seqL, qualsL) != 0)
            {
                std::cerr << "ERROR writing to " << options.outputFilename << "\n";
                return 1;
            }
            if (writeRecord(outReadsRight, ssR.str(), seqR, qualsR) != 0)
            {
                std::cerr << "ERROR writing to " << options.outputFilenameRight << "\n";
                return 1;
            }
        }
    }
    
    std::cerr << " OK\n";

    std::cerr << "\nDONE.\n";
    return 0;
}
