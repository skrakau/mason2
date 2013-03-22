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
// Simulator for the sequencing process.
// ==========================================================================

#ifndef SANDBOX_MASON2_APPS_MASON2_SEQUENCING_H_
#define SANDBOX_MASON2_APPS_MASON2_SEQUENCING_H_

#include <seqan/sequence.h>
#include <seqan/random.h>

// ============================================================================
// Forwards
// ============================================================================

typedef seqan::Dna5String TRead;
typedef seqan::CharString TQualities;
typedef seqan::Rng<seqan::MersenneTwister> TRng;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Sequencing options that are relevant for all sequencing technologies.

struct SequencingOptions
{
    // Enum for selecting the relative read orientation.
    enum MateOrientation
    {
        FORWARD_REVERSE,  // R1 --> <-- R2
        REVERSE_FORWARD,  // R1 <-- --> R2
        FORWARD_FORWARD,  // R1 --> --> R2
        FORWARD_FORWARD2  // R2 --> --> R1
    };

    // Enum for selecting forward and/or reverse strand.
    enum SourceStrands
    {
        BOTH,
        FORWARD,
        REVERSE
    };

    // Number of reads.
    int numReads;
    // Whether or not to simulate qualities.
    bool simulateQualities;
    // Whether or not to simulate mate pairs.
    bool simulateMatePairs;
    // Mate orientation.
    MateOrientation mateOrientation;
    // Whether to simulate from forward/reverse strand or both.
    SourceStrands strands;

    SequencingOptions() :
            numReads(0), simulateQualities(false), simulateMatePairs(false), mateOrientation(FORWARD_REVERSE),
            strands(BOTH)
    {}

    static const char * yesNo(bool b)
    {
        return b ? "YES" : "NO";
    }

    static const char * mateOrientationStr(MateOrientation o)
    {
        switch (o)
        {
            case FORWARD_REVERSE:
                return "FR (R1 --> <-- R2)";
            case REVERSE_FORWARD:
                return "RF (R1 <-- --> R2)";
            case FORWARD_FORWARD:
                return "FF (R1 --> --> R2)";
            case FORWARD_FORWARD2:
                return "FF2 (R2 --> --> R1)";
            default:
                return "INVALID";
        }
    }            

    static const char * strandsStr(SourceStrands s)
    {
        switch (s)
        {
            case BOTH:
                return "BOTH";
            case FORWARD:
                return "FORWARD";
            case REVERSE:
                return "REVERSE";
            default:
                return "INVALID";
        }
    }

    virtual void print(std::ostream & out)
    {
        out << "SIMULATION OPTIONS\n"
            << "\n";
        if (numReads > 0)  // is 0 if doing simulatin from fragment (induces numReads)
            out << "  NUMBER OF READS               \t" << numReads << "\n";
        out << "  SIMULATE QUALITIES            \t" << yesNo(simulateQualities) << "\n"
            << "  SIMULATE MATE PAIRS           \t" << yesNo(simulateMatePairs) << "\n"
            << "  MATE ORIENTATION              \t" << mateOrientationStr(mateOrientation) << "\n"
            << "  STRANDS                       \t" << strandsStr(strands) << "\n";
    }
};

// Sequencing options that are relevant for Illumina sequencing.

struct IlluminaSequencingOptions : SequencingOptions
{
    // Length of the reads to simulate.
    unsigned readLength;

    // Base Calling Error Model Parameters.

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

    // If set then no Ns will be introduced into the read.
    bool illuminaNoN;

    // Base Calling Quality Model Parameters.

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
            SequencingOptions(),
            readLength(0),
            // Base Calling Error Model Parameters
            probabilityInsert(0.001),
            probabilityDelete(0.001),
            probabilityMismatchScale(1.0),
            probabilityMismatch(0.004),
            probabilityMismatchBegin(0.002),
            probabilityMismatchEnd(0.012),
            positionRaise(0.66),
            illuminaNoN(false),
            // Base Calling Quality Model Parameters
            meanQualityBegin(40),
            meanQualityEnd(39.5),
            stdDevQualityBegin(0.05),
            stdDevQualityEnd(10),
            meanMismatchQualityBegin(39.5),
            meanMismatchQualityEnd(30),
            stdDevMismatchQualityBegin(3),
            stdDevMismatchQualityEnd(15)
    {
        this->mateOrientation = FORWARD_REVERSE;
    }

    virtual void print(std::ostream & out)
    {
        SequencingOptions::print(out);

        out << "\n"
            << "  SIMULATED TECHNOLOGY          \tillumina\n"
            << "\n"
            << "ILLUMINA OPTIONS\n"
            << "\n"
            << "  PROB INSERTION                \t" << probabilityInsert << "\n"
            << "  PROB DELETION                 \t" << probabilityDelete << "\n"
            << "  PROB MISMATCH SCALE           \t" << probabilityMismatchScale << "\n"
            << "  PROB MISMATCH                 \t" << probabilityMismatch << "\n"
            << "  PROB MISMATCH BEGIN           \t" << probabilityMismatchBegin << "\n"
            << "  PROB MISMATCH END             \t" << probabilityMismatchEnd << "\n"
            << "  POSITION RAISE                \t" << positionRaise << "\n"
            << "  QUALITY MEAN BEGIN            \t" << meanQualityBegin << "\n"
            << "  QUALITY MEAN END              \t" << meanQualityEnd << "\n"
            << "  QUALITY STD DEV BEGIN         \t" << stdDevQualityBegin << "\n"
            << "  QUALITY STD DEV END           \t" << stdDevQualityEnd << "\n"
            << "  MISMATCH QUALITY MEAN BEGIN   \t" << meanMismatchQualityBegin << "\n"
            << "  MISMATCH QUALITY MEAN END     \t" << meanMismatchQualityEnd << "\n"
            << "  MISMATCH QUALITY STD DEV BEGIN\t" << stdDevMismatchQualityBegin << "\n"
            << "  MISMATCH QUALITY STD DEV END  \t" << stdDevMismatchQualityEnd << "\n";
    }
};

// Sequencing options that are relevant for sanger sequencing.

struct SangerSequencingOptions : SequencingOptions
{
    // Read Length Parameters.

    // Iff true, read lengths follow a uniform distribution, otherwise a
    // standard distribution will be used.
    bool readLengthIsUniform;

    // Average read length.
    double readLengthMean;

    // For standard distributed read lengths, this is the standard deviation,
    // for uniform read length the interval around the average to use for
    // picking the read lengths.
    double readLengthError;

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

    SangerSequencingOptions()
            : readLengthIsUniform(false),
              readLengthMean(400),
              readLengthError(40),
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
    {
        this->mateOrientation = FORWARD_REVERSE;
    }

    virtual void print(std::ostream & out)
    {
        SequencingOptions::print(out);

        out << "\n"
            << "  SIMULATED TECHNOLOGY          \tsanger\n"
            << "\n"
            << "SANGER OPTIONS\n"
            << "\n"
            << "  UNIFORM READ LENGTH           \t" << yesNo(readLengthIsUniform) << "\n"
            << "  MEAN READ LENGTH              \t" << readLengthMean << "\n"
            << "  READ LENGTH ERROR             \t" << readLengthError << "\n"
            << "\n"
            << "  PROB MISMATCH BEGIN           \t" << probabilityMismatchBegin << "\n"
            << "  PROB MISMATCH END             \t" << probabilityMismatchEnd << "\n"
            << "  PROB INSERTION BEGIN          \t" << probabilityInsertBegin << "\n"
            << "  PROB INSERTION END            \t" << probabilityInsertEnd << "\n"
            << "  PROB DELETION BEGIN           \t" << probabilityDeleteBegin << "\n"
            << "  PROB DELETION END             \t" << probabilityDeleteEnd << "\n"
            << "\n"
            << "  QUALITY MEAN MATCH BEGIN      \t" << qualityMatchStartMean << "\n"
            << "  QUALITY MEAN MATCH END        \t" << qualityMatchEndMean << "\n"
            << "  QUALITY STD DEV MATCH BEGIN   \t" << qualityMatchStartStdDev << "\n"
            << "  QUALITY STD DEV MATCH END     \t" << qualityMatchEndStdDev << "\n"
            << "\n"
            << "  QUALITY MEAN MISMATCH BEGIN   \t" << qualityErrorStartMean << "\n"
            << "  QUALITY MEAN MISMATCH END     \t" << qualityErrorEndMean << "\n"
            << "  QUALITY STD DEV MISMATCH BEGIN\t" << qualityErrorStartStdDev << "\n"
            << "  QUALITY STD DEV MISMATCH END  \t" << qualityErrorEndStdDev << "\n";
    }
};

// Sequencing options that are relevant for 454 sequencing.

struct Roche454SequencingOptions : SequencingOptions
{
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
    int meanReadLength;
    // The read length standard deviation if normally distributed.
    int stdDevReadLength;

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
            SequencingOptions(), 
            lengthModel(UNIFORM), minReadLength(0), maxReadLength(0),
            meanReadLength(0), stdDevReadLength(0), sqrtInStdDev(true), k(0),
            backgroundNoiseMean(0), backgroundNoiseStdDev(0)
    {
        this->mateOrientation = FORWARD_FORWARD2;
    }

    static const char * lengthModelStr(ReadLengthModel m)
    {
        switch (m)
        {
            case UNIFORM:
                return "UNIFORM";
            case NORMAL:
                return "NORMAL";
            default:
                return "INVALID";
        }
    }

    virtual void print(std::ostream & out)
    {
        SequencingOptions::print(out);

        out << "\n"
            << "  SIMULATED TECHNOLOGY          \t454\n"
            << "\n"
            << "454 OPTIONS\n"
            << "\n"
            << "  LENGTH MODEL                  \t" << lengthModelStr(lengthModel) << "\n"
            << "  MIN READ LENGTH               \t" << minReadLength << "\n"
            << "  MAX LENGTH ERROR              \t" << maxReadLength << "\n"
            << "  MEAN READ LENGTH              \t" << meanReadLength << "\n"
            << "  READ LENGTH ERROR             \t" << stdDevReadLength << "\n"
            << "\n"
            << "  SQRT IN STD DEV               \t" << yesNo(sqrtInStdDev) << "\n"
            << "  K                             \t" << k << "\n"
            << "  BACKGROUND NOISE MEAN         \t" << backgroundNoiseMean << "\n"
            << "  BACKGROUND NOISE STD DEV      \t" << backgroundNoiseStdDev << "\n";
    }
};

// Responsible for the simulation of sequencing.
//
// We pick the read length independently of the error since the sequencing simulator might have context dependent
// errors.

class SequencingSimulator
{
public:
    // Sequencing direction from left or right side.
    enum Direction
    {
        LEFT,
        RIGHT
    };

    // Strand for sequencing.
    enum Strand
    {
        FORWARD,
        REVERSE
    };
    
    // The random number generator.
    TRng & rng;
    // Overall sequencing options.
    SequencingOptions options;

    SequencingSimulator(TRng & rng, SequencingOptions const & options) : rng(rng), options(options)
    {}

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength() = 0;

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, Direction dir, Strand strand) = 0;
};

// Illumina read simulation.

class IlluminaSequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Illumina sequencing.
    IlluminaSequencingOptions illuminaOptions;

    IlluminaSequencingSimulator(TRng & rng,
                                SequencingOptions const & options,
                                IlluminaSequencingOptions const & illuminaOptions) :
            SequencingSimulator(rng, options), illuminaOptions(illuminaOptions)
    {}

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength();

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, Direction dir, Strand strand);
};

// 454 read simulation.

class Roche454SequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Roche454 sequencing.
    Roche454SequencingOptions roche454Options;

    Roche454SequencingSimulator(TRng & rng,
                                SequencingOptions const & options,
                                Roche454SequencingOptions const & roche454Options) :
            SequencingSimulator(rng, options), roche454Options(roche454Options)
    {}

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength();

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, Direction dir, Strand strand);
};

// Sanger read simulation.

class SangerSequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Sanger sequencing.
    SangerSequencingOptions sangerOptions;

    SangerSequencingSimulator(TRng & rng,
                              SequencingOptions const & options,
                              SangerSequencingOptions const & sangerOptions) :
            SequencingSimulator(rng, options), sangerOptions(sangerOptions)
    {}

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength();

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, Direction dir, Strand strand);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_MASON2_APPS_MASON2_SEQUENCING_H_
