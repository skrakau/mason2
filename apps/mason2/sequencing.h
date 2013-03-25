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

#include <seqan/bam_io.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/random.h>

// ============================================================================
// Forwards
// ============================================================================

typedef seqan::Dna5String TRead;
typedef seqan::CharString TQualities;
typedef seqan::Rng<seqan::MersenneTwister> TRng;
typedef seqan::Infix<seqan::Dna5String const>::Type TFragment;
typedef seqan::String<seqan::CigarElement<> > TCigarString;

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

    // Verbosity level.
    int verbosity;

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
    // Whether or not to embed meta data into the read identifier.
    bool embedReadInfo;

    SequencingOptions() :
            verbosity(0), numReads(0), simulateQualities(false), simulateMatePairs(false),
            mateOrientation(FORWARD_REVERSE), strands(BOTH), embedReadInfo(false)
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

    // // If set then no Ns will be introduced into the read.
    // bool illuminaNoN;

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
            << "  QUALITY MEAN ERROR BEGIN      \t" << qualityErrorStartMean << "\n"
            << "  QUALITY MEAN ERROR END        \t" << qualityErrorEndMean << "\n"
            << "  QUALITY STD DEV ERROR BEGIN   \t" << qualityErrorStartStdDev << "\n"
            << "  QUALITY STD DEV ERROR END     \t" << qualityErrorEndStdDev << "\n";
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

// Composition of verbose information for sequencing simulation.
//
// Objects of this type store information such as the CIGAR string and the original sampled sequence for debug and
// evaluation purposes.

struct SequencingSimulationInfo
{
    // Originally sampled sequence, that together with the errors introduced by sequencing gives the read sequence.
    TRead sampleSequence;

    // The CIGAR string.
    TCigarString cigar;

    // Whether or not this comes from the forward strand.
    bool isForward;

    SequencingSimulationInfo() : isForward(false)
    {}

    template <typename TStream>
    void serialize(TStream & stream) const
    {
        stream << "SAMPLE_SEQUENCE=" << sampleSequence << " CIGAR=";
        for (unsigned i = 0; i < length(cigar); ++i)
            stream << cigar[i].count << cigar[i].operation;
        stream << " STRAND=" << (isForward ? 'F' : 'R');
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

    static Direction direction(Direction s)
    {
        return (s == LEFT) ? RIGHT : LEFT;
    }

    // Strand for sequencing.
    enum Strand
    {
        FORWARD,
        REVERSE
    };

    static Strand toggle(Strand s)
    {
        return (s == FORWARD) ? REVERSE : FORWARD;
    }
    
    // The random number generator.
    TRng & rng;
    // Overall sequencing options.
    SequencingOptions options;

    SequencingSimulator(TRng & rng, SequencingOptions const & options) : rng(rng), options(options)
    {}

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength() = 0;

    // Simulate paired-end sequencing from a fragment.
    void simulatePairedEnd(TRead & seqL, TQualities & qualsL, SequencingSimulationInfo & infoL,
                           TRead & seqR, TQualities & qualsR, SequencingSimulationInfo & infoR,
                           TFragment const & frag)
    {
        bool isForward = (pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, 1)) == 1);
        Strand strand = isForward ? FORWARD : REVERSE;
        this->simulateRead(seqL, qualsL, infoL, frag, LEFT, strand);
        this->simulateRead(seqR, qualsR, infoR, frag, RIGHT, toggle(strand));
    }

    // Simulate single-end sequencing from a fragment.
    void simulateSingleEnd(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                           TFragment const & frag)
    {
        bool isForward = (pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, 1)) == 1);
        Strand strand = isForward ? FORWARD : REVERSE;
        this->simulateRead(seq, quals, info, frag, LEFT, strand);
    }

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    //
    // seq -- target sequence of the read to simulate
    // quals -- target qualities of the read to simulate
    // frag -- source fragment
    // strand -- the strand of the fragment, coordinates are relative to forward and will be affected by this flag
    // dir -- whether this is the left or right read
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand) = 0;
};

// The internal Illumina model representation.

class IlluminaModel
{
public:
    // Probabilities for a mismatch at a given position.
    seqan::String<double> mismatchProbabilities;

    // Standard deviations for the normal distributions of base
    // qualities for the mismatch case.
    seqan::String<double> mismatchQualityMeans;
    // Standard deviations for the normal distributions of base
    // qualities for the mismatch case.
    seqan::String<double> mismatchQualityStdDevs;

    // Standard deviations for the normal distributions of base
    // qualities for the non-mismatch case.
    seqan::String<double> qualityMeans;
    // Standard deviations for the normal distributions of base
    // qualities for the non-mismatch case.
    seqan::String<double> qualityStdDevs;

    IlluminaModel()
    {}
};

// Illumina read simulation.

class IlluminaSequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Illumina sequencing.
    IlluminaSequencingOptions illuminaOptions;

    // Storage for the Illumina simulation.
    IlluminaModel model;

    IlluminaSequencingSimulator(TRng & rng,
                                IlluminaSequencingOptions const & illuminaOptions) :
            SequencingSimulator(rng, options), illuminaOptions(illuminaOptions)
    {
        this->_initModel();
    }

    void _initModel()
    {
        // Compute mismatch probabilities, piecewise linear function.
        resize(model.mismatchProbabilities, illuminaOptions.readLength);
        // Compute probability at raise point.
        double y_r = 2 * illuminaOptions.probabilityMismatch - illuminaOptions.positionRaise * illuminaOptions.probabilityMismatchBegin - illuminaOptions.probabilityMismatchEnd + illuminaOptions.probabilityMismatchEnd * illuminaOptions.positionRaise;
        if (illuminaOptions.verbosity >= 2)
        {
            std::cerr << "Illumina error curve:\n"
                      << "  (0, " << illuminaOptions.probabilityMismatchBegin << ") -- (" << illuminaOptions.positionRaise << ", " << y_r << ") -- (1, " << illuminaOptions.probabilityMismatchEnd << ")\n";
        }
        // std::cout << "y_r = " << y_r << std::endl;
        // Compute mismatch probability at each base.
        /*if (illuminaOptions.probabilityMismatchFromFile)
        {
            // Open file.
            std::fstream file;
            file.open(toCString(illuminaOptions.probabilityMismatchFile), std::ios_base::in);
            if (!file.is_open())
            {
                std::cerr << "Failed to load mismatch probabilities from " << illuminaOptions.probabilityMismatchFile << std::endl;
                return 1;
            }
            // Load probabilities.
            double x;
            file >> x;
            unsigned i;
            for (i = 0; i < illuminaOptions.readLength && !file.eof(); ++i) {
                model.mismatchProbabilities[i] = x;
                file >> x;
            }
            if (i != illuminaOptions.readLength)
            {
                std::cerr << "Not enough mismatch probabilites in " << illuminaOptions.probabilityMismatchFile << " (" << i << " < " << illuminaOptions.readLength << ")!" << std::endl;
                return 1;
            }
        } else */{
            // Use piecewise linear function for mismatch probability simulation.
            for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
                double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
                if (x < illuminaOptions.positionRaise) {
                    double b = illuminaOptions.probabilityMismatchBegin;
                    double m = (y_r - illuminaOptions.probabilityMismatchBegin) / illuminaOptions.positionRaise;
                    model.mismatchProbabilities[i] = m * x + b;
                    // std::cout << "model.mismatchProbabilities[" << i << "] = " << model.mismatchProbabilities[i] << std::endl;
                } else {
                    double b = y_r;
                    double m = (illuminaOptions.probabilityMismatchEnd - y_r) / (1 - illuminaOptions.positionRaise);
                    x -= illuminaOptions.positionRaise;
                    model.mismatchProbabilities[i] = m * x + b;
                    // std::cout << "model.mismatchProbabilities[" << i << "] = " << model.mismatchProbabilities[i] << std::endl;
                }
            }
        }
        if (illuminaOptions.probabilityMismatchScale != 1.0) {
            for (unsigned i = 0; i < illuminaOptions.readLength; ++i)
                model.mismatchProbabilities[i] *= illuminaOptions.probabilityMismatchScale;
        }

        // Compute match/mismatch means and standard deviations.
        resize(model.mismatchQualityMeans, illuminaOptions.readLength);
        for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
            double b = illuminaOptions.meanMismatchQualityBegin;
            double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
            double m = (illuminaOptions.meanMismatchQualityEnd - illuminaOptions.meanMismatchQualityBegin);
            model.mismatchQualityMeans[i] = m * x + b;
            // std::cout << "model.mismatchQualityMeans[" << i << "] = " << model.mismatchQualityMeans[i] << std::endl;
        }
        resize(model.mismatchQualityStdDevs, illuminaOptions.readLength);
        for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
            double b = illuminaOptions.stdDevMismatchQualityBegin;
            double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
            double m = (illuminaOptions.stdDevMismatchQualityEnd - illuminaOptions.stdDevMismatchQualityBegin);
            model.mismatchQualityStdDevs[i] = m * x + b;
            // std::cout << "model.mismatchQualityStdDevs[" << i << "] = " << model.mismatchQualityStdDevs[i] << std::endl;
        }
        resize(model.qualityMeans, illuminaOptions.readLength);
        for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
            double b = illuminaOptions.meanQualityBegin;
            double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
            double m = (illuminaOptions.meanQualityEnd - illuminaOptions.meanQualityBegin);
            model.qualityMeans[i] = m * x + b;
            // std::cout << "model.qualityMeans[" << i << "] = " << model.qualityMeans[i] << std::endl;
        }
        resize(model.qualityStdDevs, illuminaOptions.readLength);
        for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
            double b = illuminaOptions.stdDevQualityBegin;
            double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
            double m = (illuminaOptions.stdDevQualityEnd - illuminaOptions.stdDevQualityBegin);
            model.qualityStdDevs[i] = m * x + b;
            // std::cout << "model.qualityStdDevs[" << i << "] = " << model.qualityStdDevs[i] << std::endl;
        }
    }

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength()
    {
        return illuminaOptions.readLength;
    }

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand)
    {
        // Simulate sequencing operations.
        TCigarString cigar;
        _simulateCigar(cigar);

        // TODO(holtgrew): Check that the CIGAR string does not need more sequence than we have in frag.

        // Simulate sequence (materialize mismatches and insertions).
        typedef seqan::ModifiedString<seqan::ModifiedString<TFragment, seqan::ModView<seqan::FunctorComplement<seqan::Dna5> > >, seqan::ModReverse> TRevCompFrag;
        if ((dir == LEFT) == (strand == FORWARD))
            _simulateSequence(seq, frag, cigar);
        else
            _simulateSequence(seq, TRevCompFrag(frag), cigar);

        // Simulate qualities.
        _simulateQualities(quals, cigar);
        SEQAN_ASSERT_EQ(length(seq), length(quals));

        // Reverse-complement sequence and reverse qualities if necessary.
        if (strand == REVERSE)
            reverseComplement(seq);

        // Write out sequencing information info if configured to do so.
        if (options.embedReadInfo)
        {
            info.cigar = cigar;
            unsigned len = 0;
            _getLengthInRef(cigar, len);
            if ((dir == LEFT) == (strand == FORWARD))
                info.sampleSequence = prefix(frag, len);
            else
                info.sampleSequence = suffix(frag, length(frag) - len);
            info.isForward = (strand == FORWARD);
            if (strand == REVERSE)
                reverseComplement(info.sampleSequence);
        }
    }

    // Simulate PHRED qualities from the CIGAR string.
    void _simulateQualities(TQualities & quals, TCigarString const & cigar)
    {
        clear(quals);

        unsigned pos = 0;
        for (unsigned i = 0; i < length(cigar); ++i)
        {
            for (unsigned j = 0; j < cigar[i].count; ++j, ++pos)
            {
                int q = 0;
                if (cigar[i].operation == 'M')
                {
                    seqan::Pdf<seqan::Normal> pdf(model.qualityMeans[pos], model.qualityStdDevs[pos]);
                    q = static_cast<int>(pickRandomNumber(rng, pdf));
                }
                else if (cigar[i].operation == 'I' || cigar[i].operation == 'X')
                {
                    seqan::Pdf<seqan::Normal> pdf(model.mismatchQualityMeans[pos], model.mismatchQualityStdDevs[pos]);
                    q = static_cast<int>(pickRandomNumber(rng, pdf));
                }
                else
                {
                    // Deletion/padding, no quality required.
                    continue;
                }
                q = std::max(0, std::min(q, 40));  // limit quality to 0..40
                appendValue(quals, (char)('!' + q));
            }
        }
    }

    // Simulate the characters that polymorphisms turn into and inserted characters.
    //
    // Through the usage of ModifiedString, we will always go from the left to the right end.
    template <typename TFrag>
    void _simulateSequence(TRead & read, TFrag const & frag,
                           TCigarString const & cigar)
    {
        clear(read);

        typedef typename seqan::Iterator<TFrag>::Type TFragIter;
        TFragIter it = begin(frag, seqan::Standard());

        for (unsigned i = 0; i < length(cigar); ++i)
        {
            unsigned numSimulate = 0;
            if (cigar[i].operation == 'M')
            {
                for (unsigned j = 0; j < cigar[i].count; ++j, ++it)
                    appendValue(read, *it);
                continue;
            }
            else if (cigar[i].operation == 'D')
            {
                it += cigar[i].count;
                continue;
            }

            // Otherwise, we have insertions or mismatches.
            for (unsigned j = 0; j < cigar[i].count; ++j)
            {
                // Pick a value between 0 and 1.
                double x = 1.0;
                while (x == 1.0)
                    x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double> >(0, 1));
                int num = x / 0.25;

                // NOTE: We can only insert CGAT, but we can have a polymorphism to N.

                if (cigar[i].operation == 'I')
                    appendValue(read, seqan::Dna5(num));
                else
                    appendValue(read, seqan::Dna5(num + (num == ordValue(*it))));
            }
        }
    }

    // Simulate CIGAR string.  We can do this with position specific parameters only and thus independent of any
    // context.
    void _simulateCigar(TCigarString & cigar)
    {
        clear(cigar);
        unsigned len = this->readLength();

        for (unsigned i = 0; i < len;)
        {
            double x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double> >(0, 1));
            double pMismatch = model.mismatchProbabilities[i];
            double pInsert   = illuminaOptions.probabilityInsert;
            double pDelete   = illuminaOptions.probabilityDelete;
            double pMatch    = 1.0 - pMismatch - pInsert - pDelete;

            // Simulate mutation/insertion/deletion events.  If possible we reuse the last CIGAR entry.  Adjacent
            // insertion/deletion pairs cancel each other out.

            if (x < pMatch)  // match
            {
                ++i;
                if (empty(cigar) || back(cigar).operation != 'M')
                    appendValue(cigar, seqan::CigarElement<>('M', 1));
                else
                    back(cigar).count += 1;
            }
            else if (x < pMatch + pMismatch)  // point polymorphism
            {
                ++i;
                if (empty(cigar) || back(cigar).operation != 'X')
                    appendValue(cigar, seqan::CigarElement<>('X', 1));
                else
                    back(cigar).count += 1;
            }
            else if (x < pMatch + pMismatch + pInsert) // insertion
            {
                if (!empty(cigar) && back(cigar).operation == 'D')
                {
                    // This is an insertion following a deletion, they swallow each other up.
                    if (back(cigar).count > 1u)
                        back(cigar).count -= 1;
                    else
                        eraseBack(cigar);
                }
                else if (empty(cigar) || back(cigar).operation != 'I')
                {
                    ++i;
                    appendValue(cigar, seqan::CigarElement<>('I', 1));
                }
                else
                {
                    ++i;
                    back(cigar).count += 1;
                }
            }
            else  // deletion
            {
                if (!empty(cigar) && back(cigar).operation == 'I')
                {
                    // This is an deletion following an insertion, they swallow each other up.
                    if (back(cigar).count > 1u)
                        back(cigar).count -= 1;
                    else
                        eraseBack(cigar);
                }
                else if (empty(cigar) || back(cigar).operation != 'D')
                {
                    ++i;
                    appendValue(cigar, seqan::CigarElement<>('D', 1));
                }
                else
                {
                    ++i;
                    back(cigar).count += 1;
                }
            }
        }
    }
};

// 454 read simulation.

class Roche454SequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Roche454 sequencing.
    Roche454SequencingOptions roche454Options;

    Roche454SequencingSimulator(TRng & rng,
                                Roche454SequencingOptions const & roche454Options) :
            SequencingSimulator(rng, options), roche454Options(roche454Options)
    {}

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength()
    { return 0; }

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand)
    {}
};

// Sanger read simulation.

class SangerSequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Sanger sequencing.
    SangerSequencingOptions sangerOptions;

    SangerSequencingSimulator(TRng & rng,
                              SangerSequencingOptions const & sangerOptions) :
            SequencingSimulator(rng, options), sangerOptions(sangerOptions)
    {}

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength()
    { return 0; }

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand)
    {}
};

// Factory for SequencingSimulator objects.

class SequencingSimulatorFactory
{
public:

    enum Technology
    {
        ILLUMINA,
        SANGER,
        ROCHE_454
    };

    TRng & rng;

    seqan::SequenceStream * in;
    seqan::SequenceStream * out;
    seqan::SequenceStream * outR;

    SequencingOptions const * options;

    Technology tech;

    SequencingSimulatorFactory(TRng & rng,
                               seqan::SequenceStream * out,
                               seqan::SequenceStream * outR,
                               seqan::SequenceStream * in,
                               SequencingOptions const & options,
                               Technology tech) :
            rng(rng), in(in), out(out), outR(outR), tech(tech), options(&options)
    {}

    std::auto_ptr<SequencingSimulator> make()
    {
        std::auto_ptr<SequencingSimulator> res;
        
        switch (tech)
        {
            case ILLUMINA:
                res.reset(new IlluminaSequencingSimulator(rng, *dynamic_cast<IlluminaSequencingOptions const *>(options)));
                break;
            case SANGER:
                res.reset(new SangerSequencingSimulator(rng, *dynamic_cast<SangerSequencingOptions const *>(options)));
                break;
            case ROCHE_454:
                res.reset(new Roche454SequencingSimulator(rng, *dynamic_cast<Roche454SequencingOptions const *>(options)));
                break;
        }

        return res;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_MASON2_APPS_MASON2_SEQUENCING_H_
