// ==========================================================================
//                         Mason - A Read Simulator
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// We define the options of all Mason components in a central place.
//
//
// There is one *Options struct for each simulation core component.  Each
// program has its own Mason*Options struct for storing its configuration.
//
// Each struct has a function addOptions() and getOptions() function for
// adding options to an ArgumentParser and getting the values back from the
// parser after parsing.
//
// Also, some helper functions are provided for converting booleans and enums
// to human readable strings.
// ==========================================================================

#ifndef SANDBOX_MASON2_APPS_MASON2_MASON_OPTIONS_H_
#define SANDBOX_MASON2_APPS_MASON2_MASON_OPTIONS_H_

#include <seqan/arg_parse.h>
#include <seqan/sequence.h>

#include "mason_types.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class MaterializerOptions
// ----------------------------------------------------------------------------

// Configuration for the contig- and haplotype-wise materialization of variants from a VCF file to a reference sequence
// FASTA file.

struct MaterializerOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Path to reference file.  Required.
    seqan::CharString fastaFileName;
    // Path to VCF file.  No variation is applied if empty.
    seqan::CharString vcfFileName;

    // Path to output sequence files for left (and single end) and right reads.
    seqan::CharString outFileNameLeft, outFileNameRight;

    MaterializerOptions() : verbosity(0)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser.
    void addTextSection(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);
};

// ----------------------------------------------------------------------------
// Class FragmentSamplerOptions
// ----------------------------------------------------------------------------

// Configuration for the sampling of uniformly or normally distributed fragments from a sequence.

struct FragmentSamplerOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Type for selecting the fragment size distribution: uniformly or normally distributed.
    enum FragmentSizeModel
    {
        UNIFORM,
        NORMAL
    };

    // Smallest fragment size if uniformly distributed.
    int minFragmentSize;
    // Maximal fragment size if uniformly distributed.
    int maxFragmentSize;

    // Mean fragment size if normally distributed.
    int meanFragmentSize;
    // Standard deviation of fragment size if normally distributed.
    int stdDevFragmentSize;

    // The model to use for the fragment size.
    FragmentSizeModel model;

    FragmentOptions() :
            verbosity(0), minFragmentSize(0), maxFragmentSize(0), meanFragmentSize(0), stdDevFragmentSize(0),
            model(UNIFORM)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser.
    void addTextSection(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);
};

// ----------------------------------------------------------------------------
// Class SequencingOptions
// ----------------------------------------------------------------------------

// Configuration for generic read simulation.
//
// Generic configuration for read simulation, such as mate orientation, enabling or disabling quality and paired
// simulation, and selecting specific strands.

struct SequencingOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Enum for selecting the relative read orientation.
    enum MateOrientation
    {
        FORWARD_REVERSE = 0,  // R1 --> <-- R2
        REVERSE_FORWARD = 1,  // R1 <-- --> R2
        FORWARD_FORWARD = 2,  // R1 --> --> R2
        FORWARD_FORWARD2 = 3  // R2 --> --> R1
    };

    // Enum for selecting forward and/or reverse strand.
    enum SourceStrands
    {
        BOTH,
        FORWARD,
        REVERSE
    };

    // Enum for selecting the sequencing technology to simulate.
    enum SequencingTechnology
    {
        ILLUMINA,
        ROCHE_454,
        SANGER
    };

    // Verbosity level.
    int verbosity;

    // Whether or not to simulate qualities.
    bool simulateQualities;
    // Whether or not to simulate mate pairs.
    bool simulateMatePairs;
    // Mate orientation.
    MateOrientation mateOrientation;
    // Whether to simulate from forward/reverse strand or both.
    SourceStrands strands;
    // The sequencing technology.
    SequencingTechnology sequencingTechnology;

    SequencingOptions() :
            verbosity(0), simulateQualities(false), simulateMatePairs(false),
            mateOrientation(FORWARD_REVERSE), strands(BOTH), sequencingTechnology(ILLUMINA)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);
};

// ----------------------------------------------------------------------------
// Class IlluminaSequencingOptions
// ----------------------------------------------------------------------------

// Configuration for Illumina read simulation.
//
// Configuration specific to the Illumina sequencing model.

// TODO(holtgrew): Allow for giving a FASTQ file as the input for qualities and N-patterns.

struct IlluminaSequencingOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Length of the reads to simulate.
    unsigned readLength;

    // Path to file with positional error probabilities.
    seqan::CharString probabilityMismatchFile;

    // The default orientation for Illumina paired-end reads.
    SequencingOptions::MateOrientation defaultOrientation;

    // -----------------------------------------------------------------------
    // Base Calling Error Model Parameters.
    // -----------------------------------------------------------------------

    // Probability of an insertion.
    double probabilityInsert;
    // Probability of a deletion.
    double probabilityDelete;

    // TODO(holtgrew): Load probabilities from file or use FASTQ file as quality/N template.

    // Scale factor to apply to mismatch probabilities,
    double probabilityMismatchScale;

    // Probability of a mismatch (single-base polymorphism).
    double probabilityMismatch;
    // Probability for a mismatch in the first base.
    double probabilityMismatchBegin;
    // Probability for a mismatch in the last base.
    double probabilityMismatchEnd;
    // Relative position in the read between 0 and 1 where the steeper curve begins.
    double positionRaise;

    // // If set then no Ns will be introduced into the read.
    // bool illuminaNoN;

    // -----------------------------------------------------------------------
    // Base Calling Quality Model Parameters.
    // -----------------------------------------------------------------------

    // Mean quality for non-mismatches at the first base.
    double meanQualityBegin;
    // Mean quality for non-mismatches at last base.
    double meanQualityEnd;
    // Standard deviation quality for non-mismatches at the first base.
    double stdDevQualityBegin;
    // Standard deviation quality for non-mismatches at the last base.
    double stdDevQualityEnd;

    // Mean quality for mismatches at the first base.
    double meanMismatchQualityBegin;
    // Mean quality for mismatches at last base.
    double meanMismatchQualityEnd;
    // Standard deviation quality for mismatches at the first base.
    double stdDevMismatchQualityBegin;
    // Standard deviation quality for mismatches at the last base.
    double stdDevMismatchQualityEnd;

    IlluminaSequencingOptions() :
            verbosity(0),
            readLength(0),
            defaultOrientation(SequencingOptions::FORWARD_REVERSE),
            // Base Calling Error Model Parameters
            probabilityInsert(0.001),
            probabilityDelete(0.001),
            probabilityMismatchScale(1.0),
            probabilityMismatch(0.004),
            probabilityMismatchBegin(0.002),
            probabilityMismatchEnd(0.012),
            positionRaise(0.66),
            // illuminaNoN(false),
            // Base Calling Quality Model Parameters
            meanQualityBegin(40),
            meanQualityEnd(39.5),
            stdDevQualityBegin(0.05),
            stdDevQualityEnd(10),
            meanMismatchQualityBegin(39.5),
            meanMismatchQualityEnd(30),
            stdDevMismatchQualityBegin(3),
            stdDevMismatchQualityEnd(15)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);
};

// ----------------------------------------------------------------------------
// Class SangerSequencingOptions
// ----------------------------------------------------------------------------

// Configuration for Sanger read simulation.
//
// Configuration specific to the Sanger sequencing model.

struct SangerSequencingOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The default orientation for Illumina paired-end reads.
    SequencingOptions::MateOrientation defaultOrientation;

    // Read Length Parameters.

    // Iff true, read lengths follow a uniform distribution, otherwise a
    // standard distribution will be used.
    bool readLengthIsUniform;

    // Average read length for normal distribution.
    double readLengthMean;

    // For standard distributed read lengths, this is the standard deviation.
    double readLengthError;

    // Minimal and maximal read lenght in case of uniform distribution.
    double readLengthMin;
    double readLengthMax;

    // Base Calling Error Model Parameters.

    // Mismatch probability ramp.
    double probabilityMismatchBegin;
    double probabilityMismatchEnd;

    // Insert probability ramp.
    double probabilityInsertBegin;
    double probabilityInsertEnd;

    // Delete probability ramp.
    double probabilityDeleteBegin;
    double probabilityDeleteEnd;

    // Quality mean and standard deviation ramp specification for matches.
    double qualityMatchStartMean;
    double qualityMatchEndMean;
    double qualityMatchStartStdDev;
    double qualityMatchEndStdDev;

    // Quality mean and standard deviation ramp specification for errors.
    double qualityErrorStartMean;
    double qualityErrorEndMean;
    double qualityErrorStartStdDev;
    double qualityErrorEndStdDev;

    SangerSequencingOptions() :
            verbosity(0),
            defaultOrientation(SequencingOptions::FORWARD_REVERSE),
            readLengthIsUniform(false),
            readLengthMean(400),
            readLengthError(40),
            readLengthMin(100),
            readLengthMax(200),
            probabilityMismatchBegin(0.005),
            probabilityMismatchEnd(0.01),
            probabilityInsertBegin(0.0025),
            probabilityInsertEnd(0.005),
            probabilityDeleteBegin(0.0025),
            probabilityDeleteEnd(0.005),
            qualityMatchStartMean(40),
            qualityMatchEndMean(39),
            qualityMatchStartStdDev(0.1),
            qualityMatchEndStdDev(2),
            qualityErrorStartMean(30),
            qualityErrorEndMean(20),
            qualityErrorStartStdDev(2),
            qualityErrorEndStdDev(5)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);
};

// ----------------------------------------------------------------------------
// Class Roche454SequencingOptions
// ----------------------------------------------------------------------------

// Configuration for Roche 454 read simulation.
//
// Configuration specific to the Roche 454 sequencing model.

struct Roche454SequencingOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The default orientation for Illumina paired-end reads.
    SequencingOptions::MateOrientation defaultOrientation;

    // Enum for selecting the read length model.
    enum ReadLengthModel
    {
        UNIFORM,
        NORMAL
    };

    // Read Length Parameters

    // The model for the read length.
    ReadLengthModel lengthModel;
    // The minimal and maximal read length if uniformly distributed.
    int minReadLength, maxReadLength;
    // The mean read length if normally distributed.
    double meanReadLength;
    // The read length standard deviation if normally distributed.
    double stdDevReadLength;

    // Base-Calling Error Model Parameters

    // If true then $\sigma = k * \sqrt(r)$, otherwise $\sigma = k * r$ is used.
    bool sqrtInStdDev;

    // Proportionality factor for calculating standard deviation proportional
    // to sqrt(homopolymer length).
    double k;

    // Noise parameters.  We take the default values 0.23 and 0.15 from Metasim.
    
    // The mean of the lognormal distribution for the noise.
    double backgroundNoiseMean;

    // The standard deviation of the lognormal distribution for the noise.
    double backgroundNoiseStdDev;

    Roche454SequencingOptions() :
            verbosity(0),
            defaultOrientation(SequencingOptions::FORWARD_FORWARD2),
            lengthModel(UNIFORM),
            minReadLength(0),
            maxReadLength(0),
            meanReadLength(0),
            stdDevReadLength(0),
            sqrtInStdDev(true),
            k(0),
            backgroundNoiseMean(0),
            backgroundNoiseStdDev(0)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);
};

// ----------------------------------------------------------------------------
// Class MasonSimulatorOptions
// ----------------------------------------------------------------------------

// Configuration for the program mason_simulator.
//
// mason_simulator reads in a VCF file and a genome file and then materializes the haplotype using the core components
// of mason_materializer.  Then, fragments are sampled from the materialized file and reads are sequencing is simulated
// to yield the reads and alignmetns.  For this, the core components of mason_fragments and mason_seq_fragments are
// used.
//
// Thus, we store only a few more settings besides re-using the configuration.

struct MasonSimulatorOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;
    // The seed for the random number generator.
    int seed;
    // The spacing of the see when using multi-threading.  Thread i (beginning with 0) will get (seed + i) as its
    // initial seed.  While not crytographically safe, this should be OK for read simulation.
    int seedSpacing;
    // The number of threads to use for the simulation.
    int numThreads;

    // Configuration for the reading of the reference and application of the variants from the VCF file.
    MaterializerOptions matOptions;

    // Configuration for sampling references.
    FragmentSamplerOptions fragSamplerOptions;

    // Generic sequencing configuration.
    SequencingOptions seqOptions;
    // Configuration of the Illumina read simulation.
    IlluminaSequencingOptions illuminaOptions;
    // Configuration of the Sanger read simulation.
    SangerSequencingOptions sangerOptions;
    // Configuration of the Roche 454 read simulation.
    Roche454SequencingOptions rocheOptions;

    MasonSimulatorOptions() : verbosity(0), seed(0), seedSpacing(2048), numThreads(1)
    {}

    // Add options to the argument parser.  Calls addOptions() on the nested *Options objects.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser, calls addTextSections() on *Options objects.
    void addTextSection(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.  Calls getOptionValues() on the nested *Option objects.
    void getOptionValues(seqan::ArgumentParser const & parser);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function MaterializerOptions::addOptions()
// ----------------------------------------------------------------------------

void MaterializerOptions::addOptions(seqan::ArgumentParser & parser) const
{
    addSection(parser, "Apply VCF Variants to Reference");

    addOption(parser, seqan::ArgParseOption("ir", "input-reference", "Path to FASTA file to read the reference from.",
                                            seqan::ArgParseOption::INPUTFILE, "IN.fa"));
    setValidValues(parser, "input-reference", "fa fasta");
    setRequired(parser, "input-reference");

    addOption(parser, seqan::ArgParseOption("iv", "input-vcf", "Path to the VCF file with variants to apply.",
                                            seqan::ArgParseOption::STRING, "IN.vcf"));
    setValidValues(parser, "input-vcf");
}

// ----------------------------------------------------------------------------
// Function MaterializerOptions::addTextSections()
// ----------------------------------------------------------------------------

void MaterializerOptions::addTextSections(seqan::ArgumentParser & parser) const
{
    addTextSection(parser, "VCF Variant Notes");

    addText(parser,
            "If the option \\fB--input-vcf\\fP/\\fB-iv\\fP is given then the given VCF file is read and the variants "
            "are applied to the input reference file.  If it is not given then the input reference file is taken "
            "verbatimly for simulating reads.");
    addText(parser,
            "There are some restrictions on the VCF file and the application of the variants to the reference will "
            "fail if the VCF file is non-conforming.  VCF files from the \\fImason_variator\\fP program are "
            "guaranteed to be read.");
     addText(parser,
             "Only the haplotypes of the first individual will be generated.");
}

// ----------------------------------------------------------------------------
// Function MaterializerOptions::getOptionValues()
// ----------------------------------------------------------------------------

void MaterializerOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    getOptionValue(fastaFileName, parser, "input-reference");
    getOptionValue(vcfFileName, parser, "input-vcf");
}

// ----------------------------------------------------------------------------
// Function FragmentSamplerOptions::addOptions()
// ----------------------------------------------------------------------------

void FragmentSamplerOptions::addOptions(seqan::ArgumentParser & parser) const
{
    addSection(parser, "Fragment Size (Insert Size) Options");

    addOption(parser, seqan::ArgParseOption("", "fragment-size-model", "The model to use for the fragment size simulation.",
                                            seqan::ArgParseOption::STRING, "MODEL"));
    setDefaultValue(parser, "fragment-size-model", "normal");
    setValidValues(parser, "fragment-size-model", "normal uniform");

    addOption(parser, seqan::ArgParseOption("", "fragment-min-size", "Smallest fragment size to use when using "
                                            "uniform fragment size simulation.",
                                            seqan::ArgParseOption::INTEGER, "SIZE"));
    setDefaultValue(parser, "fragment-min-size", "100");
    setMinValue(parser, "fragment-min-size", "1");

    addOption(parser, seqan::ArgParseOption("", "fragment-max-size", "Largest fragment size to use when using "
                                            "uniform fragment size simulation.",
                                            seqan::ArgParseOption::INTEGER, "SIZE"));
    setDefaultValue(parser, "fragment-max-size", "400");
    setMaxValue(parser, "fragment-max-size", "1");

    addOption(parser, seqan::ArgParseOption("", "fragment-mean-size", "Mean fragment size for normally distributed "
                                            "fragment size simulation.",
                                            seqan::ArgParseOption::INTEGER, "SIZE"));
    setDefaultValue(parser, "fragment-mean-size", "300");
    setMaxValue(parser, "fragment-mean-size", "1");

    addOption(parser, seqan::ArgParseOption("", "fragment-size-std-dev", "Fragment size standard deviation when using "
                                            "normally distributed fragment size simulation.",
                                            seqan::ArgParseOption::INTEGER, "SIZE"));
    setDefaultValue(parser, "fragment-size-std-dev", "30");
    setMaxValue(parser, "fragment-size-std-dev", "1");
}

// ----------------------------------------------------------------------------
// Function FragmentSamplerOptions::addTextSections()
// ----------------------------------------------------------------------------

void FragmentSamplerOptions::addTextSections(seqan::ArgumentParser & parser) const
{
    addTextSection(parser, "Fragment Size (Insert Size) Simulation");

    addText(parser,
            "You can choose between a normal and a uniform distribution of fragment lengths.  When sequencing "
            "these fragments from both sides in a paired protocol, the fragment size will become the insert "
            "size.");
}

// ----------------------------------------------------------------------------
// Function FragmentSamplerOptions::getOptionValues()
// ----------------------------------------------------------------------------

void FragmentSamplerOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    seqan::CharString tmp;
    getOptionValue(tmp, parser, "fragment-size-model");
    if (tmp == "normal")
        fragSamplerOptions.model = FragmentSamplerOptions::NORMAL;
    else
        fragSamplerOptions.model = FragmentSamplerOptions::UNIFORM;

    getOptionValue(fragSamplerOptions.minFragmentSize, parser, "fragment-min-size");
    getOptionValue(fragSamplerOptions.maxFragmentSize, parser, "fragment-max-size");
    getOptionValue(fragSamplerOptions.meanFragmentSize, parser, "fragment-mean-size");
    getOptionValue(fragSamplerOptions.stdDevFragmentSize, parser, "fragment-size-std-dev");
}

// ----------------------------------------------------------------------------
// Function SequencingOptions::addOptions()
// ----------------------------------------------------------------------------

void SequencingOptions::addOptions(seqan::ArgumentParser & parser) const
{
    addSection(parser, "Global Read Simulation Options");

    addOption(parser, seqan::ArgParseOption("", "seq-technology", "Set sequencing technology to simulate.",
                                            seqan::ArgParseOption::STRING, "TECH"));
    setValidValues(parser, "seq-technology", "illumina 454 sanger");
    setDefaultValue(parser, "seq-technology", "illumina");

    addOption(parser, seqan::ArgParseOption("", "seq-mate-orientation", "Orientation for paired reads.  See section Read "
                                            "Orientation below.", seqan::ArgParseOption::STRING, "ORIENTATION"));
    setValidValues(parser, "seq-mate-orientation", "FR RF FF FF2");
    setDefaultValue(parser, "seq-mate-orientation", "FR");

    addOption(parser, seqan::ArgParseOption("", "seq-strands", "Strands to simulate from, only applicable to paired "
                                            "sequencing simulation.", seqan::ArgParseOption::STRING, "STRAND"));
    setValidValues(parser, "seq-strands", "forward reverse both");
    setDefaultValue(parser, "seq-strands", "both");
}

// ----------------------------------------------------------------------------
// Function SequencingOptions::addTextSections()
// ----------------------------------------------------------------------------

void SequencingOptions::addTextSections(seqan::ArgumentParser & parser) const
{
    addTextSection(parser, "Sequencing Simulation");

    addText(parser,
            "Simulation of base qualities is disabled when writing out FASTA files.  Simulation of paired-end "
            "sequencing is enabled when specifying two output files.");

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
}

// ----------------------------------------------------------------------------
// Function SequencingOptions::getOptionValues()
// ----------------------------------------------------------------------------

void SequencingOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    seqan::CharString tmp;

    getOptionValue(tmp, parser, "seq-mate-orientation");
    if (tmp == "FR")
        mateOrientation = FORWARD_REVERSE;
    else if (tmp == "RF")
        mateOrientation = REVERSE_FORWARD;
    else if (tmp == "FF")
        mateOrientation = FORWARD_FORWARD;
    else if (tmp == "FF2")
        mateOrientaiton = FORWARD_FORWARD2;

    getOptionValue(tmp, parser, "seq-strands");
    if (tmp == "forward")
        strands = FORWARD;
    else if (tmp == "reverse")
        strands = REVERSE;
    else
        strands = BOTH;

    getOptionValue(tmp, parser, "seq-technology");
    if (tmp == "illumina")
        sequencingTechnology = SequencingOptions::ILLUMINA;
    else if (tmp == "454")
        sequencingTechnology = SequencingOptions::ROCHE_454;
    else
        sequencingTechnology = SequencingOptions::SANGER;
}

// ----------------------------------------------------------------------------
// Function IlluminaSequencingOptions::addOptions()
// ----------------------------------------------------------------------------

void IlluminaSequencingOptions::addOptions(seqan::ArgumentParser & parser) const
{
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
}

// ----------------------------------------------------------------------------
// Function IlluminaSequencingOptions::addTextSections()
// ----------------------------------------------------------------------------

void IlluminaSequencingOptions::addTextSections(seqan::ArgumentParser & parser) const
{}

// ----------------------------------------------------------------------------
// Function IlluminaSequencingOptions::getOptionValues()
// ----------------------------------------------------------------------------

void IlluminaSequencingOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    getOptionValue(readLength, parser, "illumina-read-length");

    getOptionValue(probabilityMismatchFile, parser, "illumina-error-profile-file");

    getOptionValue(probabilityInsert, parser, "illumina-prob-insert");
    getOptionValue(probabilityDelete, parser, "illumina-prob-deletion");

    getOptionValue(probabilityMismatchScale, parser, "illumina-prob-mismatch-scale");

    getOptionValue(probabilityMismatch, parser, "illumina-prob-mismatch");
    getOptionValue(probabilityMismatchBegin, parser, "illumina-prob-mismatch-begin");
    getOptionValue(probabilityMismatchEnd, parser, "illumina-prob-mismatch-end");
    getOptionValue(positionRaise, parser, "illumina-position-raise"); 

    getOptionValue(meanQualityBegin, parser, "illumina-quality-mean-begin");
    getOptionValue(meanQualityEnd, parser, "illumina-quality-mean-end");
    getOptionValue(stdDevQualityBegin, parser, "illumina-quality-stddev-begin");
    getOptionValue(stdDevQualityEnd, parser, "illumina-quality-stddev-end");

    getOptionValue(meanMismatchQualityBegin, parser, "illumina-mismatch-quality-mean-begin");
    getOptionValue(meanMismatchQualityEnd, parser, "illumina-mismatch-quality-mean-end");
    getOptionValue(stdDevMismatchQualityBegin, parser, "illumina-mismatch-quality-stddev-begin");
    getOptionValue(stdDevMismatchQualityEnd, parser, "illumina-mismatch-quality-stddev-end");
}

// ----------------------------------------------------------------------------
// Function SangerSequencingOptions::addOptions()
// ----------------------------------------------------------------------------

void SangerSequencingOptions::addOptions(seqan::ArgumentParser & parser) const
{
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
}

// ----------------------------------------------------------------------------
// Function SangerSequencingOptions::addTextSections()
// ----------------------------------------------------------------------------

void SangerSequencingOptions::addTextSections(seqan::ArgumentParser & parser) const
{}

// ----------------------------------------------------------------------------
// Function SangerSequencingOptions::getOptionValues()
// ----------------------------------------------------------------------------

void SangerSequencingOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    seqan::CharString tmp;
    getOptionValue(tmp, parser, "sanger-read-length-model");
    readLengthIsUniform = (tmp == "uniform");

    getOptionValue(readLengthMean, parser, "sanger-read-length-mean");
    getOptionValue(readLengthError, parser, "sanger-read-length-error");
    getOptionValue(readLengthMin, parser, "sanger-read-length-min");
    getOptionValue(readLengthMax, parser, "sanger-read-length-max");
        
    getOptionValue(probabilityMismatchBegin, parser, "sanger-prob-mismatch-begin");
    getOptionValue(probabilityMismatchEnd, parser, "sanger-prob-mismatch-end");
    getOptionValue(probabilityInsertBegin, parser, "sanger-prob-insertion-begin");
    getOptionValue(probabilityInsertEnd, parser, "sanger-prob-insertion-end");
    getOptionValue(probabilityDeleteBegin, parser, "sanger-prob-deletion-begin");
    getOptionValue(probabilityDeleteEnd, parser, "sanger-prob-deletion-end");

    getOptionValue(qualityMatchStartMean, parser, "sanger-quality-match-start-mean");
    getOptionValue(qualityMatchEndMean, parser, "sanger-quality-match-end-mean");
    getOptionValue(qualityMatchStartStdDev, parser, "sanger-quality-match-start-stddev");
    getOptionValue(qualityMatchEndStdDev, parser, "sanger-quality-match-end-stddev");
    getOptionValue(qualityErrorStartMean, parser, "sanger-quality-error-start-mean");
    getOptionValue(qualityErrorEndMean, parser, "sanger-quality-error-end-mean");
    getOptionValue(qualityErrorStartStdDev, parser, "sanger-quality-error-start-stddev");
    getOptionValue(qualityErrorEndStdDev, parser, "sanger-quality-error-end-stddev");
}

// ----------------------------------------------------------------------------
// Function Roche454SequencingOptions::addOptions()
// ----------------------------------------------------------------------------

void Roche454SequencingOptions::addOptions(seqan::ArgumentParser & parser) const
{
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
}

// ----------------------------------------------------------------------------
// Function Roche454SequencingOptions::addTextSections()
// ----------------------------------------------------------------------------

void Roche454SequencingOptions::addTextSections(seqan::ArgumentParser & parser) const
{}

// ----------------------------------------------------------------------------
// Function Roche454SequencingOptions::getOptionValues()
// ----------------------------------------------------------------------------

void Roche454SequencingOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    seqan::CharString tmp;
    getOptionValue(tmp, parser, "454-read-length-model");
    lengthModel = (tmp == "uniform") ? Roche454SequencingOptions::UNIFORM : Roche454SequencingOptions::NORMAL;

    getOptionValue(minReadLength, parser, "454-read-length-min");
    getOptionValue(maxReadLength, parser, "454-read-length-max");
    getOptionValue(meanReadLength, parser, "454-read-length-mean");
    getOptionValue(stdDevReadLength, parser, "454-read-length-stddev");
    sqrtInStdDev = !isSet(parser, "454-no-sqrt-in-std-dev");
    getOptionValue(k, parser, "454-proportionality-factor");
    getOptionValue(backgroundNoiseMean, parser, "454-background-noise-mean");
    getOptionValue(backgroundNoiseStdDev, parser, "454-background-noise-stddev");
}

// ----------------------------------------------------------------------------
// Function MasonSimulatorOptions::addOptions()
// ----------------------------------------------------------------------------

void MasonSimulatorOptions::addOptions(seqan::ArgumentParser & parser) const
{
    // Add top-level options.

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Low verbosity."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Higher verbosity."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Highest verbosity."));

    addOption(parser, seqan::ArgParseOption("", "seed", "Seed to use for random number generator.",
                                            seqan::ArgParseOptions::INTEGER, "NUM"));
    setDefaultValue(parser, "seed", "0");

    addOption(parser, seqan::ArgParseOption("", "seed-spacing", "Offset for seeds to use when multi-threading.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setDefaultValue(parser, "seed-spacing", "2048");

    addOption(parser, seqan::ArgParseOption("", "num-threads", "Number of threads to use.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setMinValue(parser, "num-threads", "1");
    setDefaultValue(parser, "num-threads", "1");

    addOption(parser, seqan::ArgParseOption("o", "out", "Output of single-end/left end reads.",
                                            seqan::ArgParseOption::OUTFILE, "OUT"));
    setRequired(parser, "out");
    setValidValues(parser, "out", "fa fasta fq fastq");

    addOption(parser, seqan::ArgParseOption("or", "out-right", "Output of right reads.  Giving this options enables "
                                            "paired-end simulation.", seqan::ArgParseOption::OUTFILE, "OUT2"));
    setRequired(parser, "out-right");
    setValidValues(parser, "out-right", "fa fasta fq fastq");

    // Add options of the component options.
    matOptions.addOptions(parser);
    fragSamplerOptions.addOptions(parser);
    seqOptions.addOptions(parser);
    illuminaOptions.addOptions(parser);
    sangerOptions.addOptions(parser);
    rocheOptions.addOptions(parser);
}

// ----------------------------------------------------------------------------
// Function MasonSimulatorOptions::addTextSections()
// ----------------------------------------------------------------------------

void MasonSimulatorOptions::addTextSections(seqan::ArgumentParser & parser) const
{
    // Add top-level text sections.
    addTextSection(parser, "Simulation Overview");

    addText(parser,
            "The first step is the application of VCF variants to the input reference file.");

    addText(parser,
            "After the generation of the haplotypes, fragments are sampled from the sequence.  These fragments "
            "correspond to the fragments in the real preparation step for the sequencing.  They are later sequenced "
            "from one or both sides depending on whether a single-end or a paired protocol is used.");

    addTextSection(parser, "Important Parameters");

    addText(parser, "For most users, the following options are most important.");

    addListItem(parser, "Paired-End Simulation",
                "Use --fragment-length-model to switch between normally and uniformly distributed insert sizes. "
                "Use the --fragment-* options for configuring the insert size simulation.");

    addTextSection(parser, "Multi-Threading");

    addText(parser,
            "When using multi-threading, each thread gets its own random number generator (RNG).  The RNG of thread "
            "i is initialized with the value of \\fB--seed\\fP plus i."):

    // Add text sections of the component options.
    matOptions.addTextSections(parser);
    fragSamplerOptions.addTextSections(parser);
    seqOptions.addTextSections(parser);
    illuminaOptions.addTextSections(parser);
    sangerOptions.addTextSections(parser);
    rocheOptions.addTextSections(parser);
}

// ----------------------------------------------------------------------------
// Function MasonSimulatorOptions::getOptionValues()
// ----------------------------------------------------------------------------

void MasonSimulatorOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    // Get top-level options.
    getOptionValues(verbosity, parser, "verbosity");
    getOptionValues(seed, parser, "seed");
    getOptionValues(seedSpacing, parser, "seed-spacing");
    getOptionValues(numThreads, parser, "num-threads");

    // Get options for the other components that we use.
    matOptions.getOptionValues(parser);
    fragSamplerOptions.getOptionValues(parser);
    seqOptions.getOptionValues(parser);
    illuminaOptions.getOptionValues(parser);
    sangerOptions.getOptionValues(parser);
    rocheOptions.getOptionValues(parser);

    // Copy in the verbosity flag into the component options.
    matOptions.verbosity = verbosity;
    fragSamplerOptions.verbosity = verbosity;
    seqOptions.verbosity = verbosity;
    illuminaOptions.verbosity = verbosity;
    sangerOptions.verbosity = verbosity;
    rocheOptions.verbosity = verbosity;

    // Configure simulation of pairs and mates depending on output files.
    seqOptions.simulateQualities = (endsWith(outFileNameLeft, ".fastq") || endsWith(outFileNameLeft, ".fq"));
    seqOptions.simulateMatePairs = !empty(outFileNameRight);
}

#endif  // #ifndef SANDBOX_MASON2_APPS_MASON2_MASON_OPTIONS_H_
