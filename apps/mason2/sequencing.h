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

#include <stdexcept>

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

// Utility.

// Returns (a, b) where a is the difference in resulting read length and b is the difference in used up input sequence
// length.

inline std::pair<int, int> appendOperation(TCigarString & cigar, char op)
{
    // Canceling out of events.  This can only happen if the CIGAR string is not empty.
    if (!empty(cigar) && ((back(cigar).operation == 'I' && op == 'D') ||
                          (back(cigar).operation == 'D' && op == 'I')))
    {
        if (back(cigar).count > 1)
            back(cigar).count -= 1;
        else
            eraseBack(cigar);
        return std::make_pair(-(op == 'I'), -(op == 'D'));
    }

    // No canceling out of events.  The read length increases by one if the operation is no deletion and one base of
    // input sequence is used up if the operation is not an insertion.
    if (!empty(cigar) && back(cigar).operation == op)
        back(cigar).count += 1;
    else
        appendValue(cigar, seqan::CigarElement<>(op, 1));
    return std::make_pair((op != 'D'), (op != 'I'));
}

// Sequencing options that are relevant for Illumina sequencing.

struct IlluminaSequencingOptions : SequencingOptions
{
    // Length of the reads to simulate.
    unsigned readLength;

    // Path to file with positional error probabilities.
    seqan::CharString probabilityMismatchFile;

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

    SangerSequencingOptions()
            : readLengthIsUniform(false),
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
    SequencingOptions const * options;

    SequencingSimulator(TRng & rng, SequencingOptions const & _options) : rng(rng), options(&_options)
    {}

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength() = 0;

    // Simulate paired-end sequencing from a fragment.
    void simulatePairedEnd(TRead & seqL, TQualities & qualsL, SequencingSimulationInfo & infoL,
                           TRead & seqR, TQualities & qualsR, SequencingSimulationInfo & infoR,
                           TFragment const & frag)
    {
        bool isForward = (pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, 1)) == 1);
        // TODO(holtgrew): Use a table for this to simplify things?
        Strand leftStrand = isForward ? FORWARD : REVERSE;
        switch (options->mateOrientation)
        {
            case SequencingOptions::FORWARD_REVERSE:
                if (leftStrand == FORWARD)
                {
                    this->simulateRead(seqL, qualsL, infoL, frag, LEFT, FORWARD);
                    this->simulateRead(seqR, qualsR, infoR, frag, RIGHT, REVERSE);
                }
                else
                {
                    this->simulateRead(seqL, qualsL, infoL, frag, RIGHT, REVERSE);
                    this->simulateRead(seqR, qualsR, infoR, frag, LEFT, FORWARD);
                }
                break;
            case SequencingOptions::REVERSE_FORWARD:
                if (leftStrand == FORWARD)
                {
                    this->simulateRead(seqL, qualsL, infoL, frag, LEFT, REVERSE);
                    this->simulateRead(seqR, qualsR, infoR, frag, RIGHT, FORWARD);
                }
                else
                {
                    this->simulateRead(seqL, qualsL, infoL, frag, RIGHT, FORWARD);
                    this->simulateRead(seqR, qualsR, infoR, frag, LEFT, REVERSE);
                }
                break;
            case SequencingOptions::FORWARD_FORWARD:
                if (leftStrand == FORWARD)
                {
                    this->simulateRead(seqL, qualsL, infoL, frag, LEFT, FORWARD);
                    this->simulateRead(seqR, qualsR, infoR, frag, RIGHT, FORWARD);
                }
                else
                {
                    this->simulateRead(seqL, qualsL, infoL, frag, RIGHT, REVERSE);
                    this->simulateRead(seqR, qualsR, infoR, frag, LEFT, REVERSE);
                }
                break;
            case SequencingOptions::FORWARD_FORWARD2:
                if (leftStrand == FORWARD)
                {
                    this->simulateRead(seqL, qualsL, infoL, frag, RIGHT, FORWARD);
                    this->simulateRead(seqR, qualsR, infoR, frag, LEFT, FORWARD);
                }
                else
                {
                    this->simulateRead(seqL, qualsL, infoL, frag, LEFT, REVERSE);
                    this->simulateRead(seqR, qualsR, infoR, frag, RIGHT, REVERSE);
                }
                break;
        }
        // std::cerr << "\n";
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
            SequencingSimulator(rng, illuminaOptions), illuminaOptions(illuminaOptions)
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
        if (!empty(illuminaOptions.probabilityMismatchFile))
        {
            // Open file.
            std::fstream file;
            file.open(toCString(illuminaOptions.probabilityMismatchFile), std::ios_base::in);
            if (!file.is_open())
            {
                std::cerr << "Failed to load mismatch probabilities from " << illuminaOptions.probabilityMismatchFile << std::endl;
                // return 1;
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
                // return 1;
            }
        } else {
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
        // std::cerr << "simulateRead(" << (char const *)(dir == LEFT ? "L" : "R") << ", " << (char const *)(strand == FORWARD ? "-->" : "<--") << ")\n";
        // Simulate sequencing operations.
        TCigarString cigar;
        _simulateCigar(cigar);
        unsigned lenInRef = 0;
        _getLengthInRef(cigar, lenInRef);

        // TODO(holtgrew): Check that the CIGAR string does not need more sequence than we have in frag.

        // Simulate sequence (materialize mismatches and insertions).
        typedef seqan::ModifiedString<seqan::ModifiedString<TFragment, seqan::ModView<seqan::FunctorComplement<seqan::Dna5> > >, seqan::ModReverse> TRevCompFrag;
        if ((dir == LEFT) && (strand == FORWARD))
            _simulateSequence(seq, prefix(frag, lenInRef), cigar);
        else if ((dir == LEFT) && (strand == REVERSE))
            _simulateSequence(seq, TRevCompFrag(prefix(frag, lenInRef)), cigar);
        else if ((dir == RIGHT) && (strand == FORWARD))
            _simulateSequence(seq, suffix(frag, length(frag) - lenInRef), cigar);
        else  // ((dir == RIGHT) && (strand == REVERSE))
            _simulateSequence(seq, TRevCompFrag(suffix(frag, length(frag) - lenInRef)), cigar);

        // Simulate qualities.
        _simulateQualities(quals, cigar);
        SEQAN_ASSERT_EQ(length(seq), length(quals));

        // Reverse qualities if necessary.
        if (strand == REVERSE)
            reverse(quals);

        // Write out sequencing information info if configured to do so.
        if (illuminaOptions.embedReadInfo)
        {
            info.cigar = cigar;
            unsigned len = 0;
            _getLengthInRef(cigar, len);
            if (dir == LEFT)
                info.sampleSequence = prefix(frag, len);
            else
                info.sampleSequence = suffix(frag, length(frag) - len);
            info.isForward = (strand == FORWARD);
            if (strand == REVERSE)
                reverseComplement(info.sampleSequence);
        }
        // std::cerr << "  frag=" << frag << "\n";
        // std::cerr << "  " << seq << "\t" << quals << "\n";
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

        for (int i = 0; i < (int)len;)
        {
            double x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double> >(0, 1));
            double pMismatch = model.mismatchProbabilities[i];
            double pInsert   = illuminaOptions.probabilityInsert;
            double pDelete   = illuminaOptions.probabilityDelete;
            double pMatch    = 1.0 - pMismatch - pInsert - pDelete;

            // Simulate mutation/insertion/deletion events.  If possible we reuse the last CIGAR entry.  Adjacent
            // insertion/deletion pairs cancel each other out.

            // TODO(holtgrew): No indel at beginning or ending! Same for other simulators!

            if (x < pMatch)  // match
                i += appendOperation(cigar, 'M').first;
            else if (x < pMatch + pMismatch)  // point polymorphism
                i += appendOperation(cigar, 'X').first;
            else if (x < pMatch + pMismatch + pInsert) // insertion
                i += appendOperation(cigar, 'I').first;
            else  // deletion
                i += appendOperation(cigar, 'D').first;
        }
    }
};

// Maximal homopolymer length we will observe.

const unsigned MAX_HOMOPOLYMER_LEN = 40;

// 454 Model

class ThresholdMatrix
{
public:
    // The scaling parameter k.
    double _k;
    // Whether or not to use the sqrt for the std deviation computation.
    bool _useSqrt;
    // Mean of the log normally distributed noise.
    double _noiseMu;
    // Standard deviation of the log normally distributed noise.
    double _noiseSigma;
    // The edge length of the matrix.
    mutable unsigned _size;
    // The data of the matrix.
    mutable seqan::String<double> _data;

    ThresholdMatrix()
            : _k(0), _useSqrt(false), _noiseMu(0), _noiseSigma(0), _size(0)
    {}
    
    ThresholdMatrix(double k, bool useSqrt, double noiseMu, double noiseSigma)
            : _k(k), _useSqrt(useSqrt), _noiseMu(noiseMu), _noiseSigma(noiseSigma), _size(0)
    {}

    ThresholdMatrix(double k, bool useSqrt, double noiseMean, double noiseStdDev, seqan::MeanStdDev const & /*tag*/)
            : _k(k), _useSqrt(useSqrt),
              _noiseMu(::std::log(noiseMean) - 0.5 * ::std::log(1.0 + noiseStdDev * noiseStdDev / noiseMean / noiseMean)),
              _noiseSigma(::std::sqrt(::std::log(1.0 + noiseStdDev * noiseStdDev / noiseMean / noiseMean))),
              _size(0)
    {}
};

inline double
normalDensityF(double x, double mu, double sigma)
{
    const double PI = 3.14159265;
    double sigma2 = sigma * sigma;
    return exp(- (x - mu) * (x - mu) / (2 * sigma2)) / sqrt(2 * PI * sigma2);
}

inline double
lognormalDensityF(double x, double mu, double sigma)
{
    if (x <= 0)
        return 0;
    const double PI = 3.14159265;
    double sigma2 = sigma * sigma;
    double log_mu2 = (log(x) - mu) * (log(x) - mu);
    return exp(-log_mu2 / (2 * sigma2)) / (x * sigma * sqrt(2 * PI));
}

inline double
dispatchDensityFunction(ThresholdMatrix const & matrix, unsigned r, double x)
{
    if (r == 0) {
        return lognormalDensityF(x, matrix._noiseMu, matrix._noiseSigma);
    } else {
        double rd = static_cast<double>(r);
        return normalDensityF(x, rd, (matrix._useSqrt ? sqrt(rd) : rd));
    }
}

inline double
computeThreshold(ThresholdMatrix const & matrix, unsigned r1, unsigned r2)
{
    if (r1 > r2)
        return computeThreshold(matrix, r2, r1);
    // The epsilon we use for convergence detection.
    const double EPSILON = 0.00001;

    // In i, we will count the number of iterations so we can limit the maximal
    // number of iterations.
    unsigned i = 0;

    // f1 is the density function for r1 and f2 the density function for r2.

    // Pick left such that f1(left) > f2(left).
    double left = r1;
    if (left == 0) left = 0.23;
    while (dispatchDensityFunction(matrix, r1, left) <= dispatchDensityFunction(matrix, r2, left))
        left /= 2.0;
    // And pick right such that f1(right) < f2(right).
    double right = r2;
    if (right == 0) right = 0.5;
    while (dispatchDensityFunction(matrix, r1, right) >= dispatchDensityFunction(matrix, r2, right))
        right *= 2.;

    // Now, search for the intersection point.
    while (true)
    {
        SEQAN_ASSERT_LT_MSG(i, 1000u, "Too many iterations (%u)! r1 = %u, r2 = %u.", i, r1, r2);
        i += 1;

        double center = (left + right) / 2;
        double fCenter1 = dispatchDensityFunction(matrix, r1, center);
        double fCenter2 = dispatchDensityFunction(matrix, r2, center);
        double delta = fabs(fCenter1 - fCenter2);
        if (delta < EPSILON)
            return center;

        if (fCenter1 < fCenter2)
            right = center;
        else
            left = center;
    }
}

inline void
extendThresholds(ThresholdMatrix const & matrix, unsigned dim)
{
    // Allocate new data array for matrix.  Then compute values or copy
    // over existing ones.
    seqan::String<double> newData;
    resize(newData, dim * dim);
    for (unsigned i = 0; i < dim; ++i) {
        for (unsigned j = 0; j < dim; ++j) {
            if (i == j)
                continue;
            if (i < matrix._size && j < matrix._size)
                newData[i * dim + j] = matrix._data[i * matrix._size + j];
            else
                newData[i * dim + j] = computeThreshold(matrix, i, j);
        }
    }
    // Update matrix.
    assign(matrix._data, newData);
    matrix._size = dim;
}

inline double
getThreshold(ThresholdMatrix const & matrix, unsigned r1, unsigned r2)
{
    if (matrix._size <= r1 || matrix._size <= r2)
        extendThresholds(matrix, std::max(r1, r2) + 1);
    return matrix._data[r1 * matrix._size + r2];
}

inline void
setK(ThresholdMatrix & matrix, double k)
{
    matrix._k = k;
}

inline void
setUseSqrt(ThresholdMatrix & matrix, bool useSqrt)
{
    matrix._useSqrt = useSqrt;
}

inline void
setNoiseMu(ThresholdMatrix & matrix, double mu)
{
    matrix._noiseMu = mu;
}

inline void
setNoiseSigma(ThresholdMatrix & matrix, double sigma)
{
    matrix._noiseSigma = sigma;
}

inline void
setNoiseMeanStdDev(ThresholdMatrix & matrix, double mean, double stdDev)
{
    matrix._noiseMu = ::std::log(mean) - 0.5 * ::std::log(1.0 + stdDev * stdDev / mean / mean);
    matrix._noiseSigma = ::std::sqrt(::std::log(1.0 + stdDev * stdDev / mean / mean));
}

// Stores the threshold matrix.

struct Roche454Model
{
    ThresholdMatrix thresholdMatrix;
};

// 454 read simulation.

class Roche454SequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Roche454 sequencing.
    Roche454SequencingOptions roche454Options;

    // Precomputed model data for 454 Sequencing.
    Roche454Model model;

    Roche454SequencingSimulator(TRng & rng,
                                Roche454SequencingOptions const & roche454Options) :
            SequencingSimulator(rng, roche454Options), roche454Options(roche454Options)
    {
        _initModel();
    }

    // Initialize the threshold matrix.
    void _initModel()
    {
        setK(model.thresholdMatrix, roche454Options.k);
        setUseSqrt(model.thresholdMatrix, roche454Options.sqrtInStdDev);
        setNoiseMeanStdDev(model.thresholdMatrix, roche454Options.backgroundNoiseMean, roche454Options.backgroundNoiseStdDev);
    }

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength()
    {
        if (roche454Options.lengthModel == Roche454SequencingOptions::UNIFORM)
        {
            // Pick uniformly.
            double minLen = roche454Options.minReadLength;
            double maxLen = roche454Options.maxReadLength;
            double len = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double> >(minLen, maxLen));
            return static_cast<unsigned>(round(len));
        }
        else
        {
            // Pick normally distributed.
            double len = pickRandomNumber(rng, seqan::Pdf<seqan::Normal>(roche454Options.meanReadLength,
                                                                         roche454Options.stdDevReadLength));
            return static_cast<unsigned>(round(len));
        }
    }

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand)
    {
        clear(seq);
        clear(quals);

        // Compute read length and check whether it fits in fragment.
        unsigned readLength = this->readLength();
        if (readLength > length(frag))
        {
            throw std::runtime_error("454 read is too long, increase fragment length");
        }

        // Get a copy of the to be sequenced base stretch.
        TRead haplotypeInfix;
        if (dir == LEFT)
            haplotypeInfix = prefix(frag, readLength);
        else
            haplotypeInfix = suffix(frag, length(frag) - readLength);
        if (strand == REVERSE)
            reverseComplement(haplotypeInfix);

        // In the flow cell simulation, we will simulate light intensities which will be stored in observedIntensities.
        seqan::String<double> observedIntensities;
        reserve(observedIntensities, 4 * readLength);
        seqan::Dna5String observedBases;
        // We also store the real homopolymer length.
        seqan::String<unsigned> realBaseCount;

        // Probability density function to use for the background noise.
        seqan::Pdf<seqan::LogNormal> noisePdf(roche454Options.backgroundNoiseMean,
                                              roche454Options.backgroundNoiseStdDev, seqan::MeanStdDev());

        // Initialize information about the current homopolymer length.
        unsigned homopolymerLength = 0;
        seqan::Dna homopolymerType = haplotypeInfix[0];
        while (homopolymerLength < length(haplotypeInfix) && haplotypeInfix[homopolymerLength] == homopolymerType)
            ++homopolymerLength;

        // Simulate flowcell.
        for (unsigned i = 0, j = 0; i < readLength; ++j, j = j % 4)  // i indicates first pos of current homopolymer, j indicates flow phase
        {
            if (ordValue(homopolymerType) == j)
            {
                // Simulate positive flow observation.
                double l = homopolymerLength;
                double sigma = roche454Options.k * (roche454Options.sqrtInStdDev ? sqrt(l) : l);
                double intensity = pickRandomNumber(rng, seqan::Pdf<seqan::Normal>(homopolymerLength, sigma));
                intensity += pickRandomNumber(rng, noisePdf);  // Add noise.
                appendValue(observedIntensities, intensity);
                appendValue(realBaseCount, homopolymerLength);
                // Get begin pos and length of next homopolymer.
                i += homopolymerLength;
                if (i < length(haplotypeInfix))
                {
                    homopolymerType = haplotypeInfix[i];
                    homopolymerLength = 0;
                    while (((i + homopolymerLength) < length(haplotypeInfix)) && haplotypeInfix[i + homopolymerLength] == homopolymerType)
                        ++homopolymerLength;
                }
            }
            else
            {
                // Simulate negative flow observation.
                //
                // Constants taken from MetaSim paper which have it from the
                // original 454 publication.
                double intensity = std::max(0.0, pickRandomNumber(rng, noisePdf));
                appendValue(observedIntensities, intensity);
                appendValue(realBaseCount, 0);
            }
        }

        seqan::String<seqan::CigarElement<> > cigar;

        // Call bases, from this build the edit string and maybe qualities.  We only support the "inter" base calling
        // method which was published by the MetaSim authors in the PLOS paper.
        typedef seqan::Iterator<seqan::String<double>, seqan::Standard>::Type IntensitiesIterator;
        int i = 0;  // Flow round, Dna(i % 4) gives base.
        for (IntensitiesIterator it = begin(observedIntensities); it != end(observedIntensities); ++it, ++i)
        {
            double threshold = getThreshold(model.thresholdMatrix, static_cast<unsigned>(floor(*it)), static_cast<unsigned>(ceil(*it)));
            unsigned calledBaseCount = static_cast<unsigned>(*it < threshold ? floor(*it) : ceil(*it));
            // Add any matches.
            unsigned j = 0;
            for (; j < std::min(calledBaseCount, realBaseCount[i]); ++j)
            {
                appendOperation(cigar, 'M');
                appendValue(seq, seqan::Dna(i % 4));
            }
            // Add insertions, if any.
            for (; j < calledBaseCount; ++j)
            {
                appendOperation(cigar, 'I');
                appendValue(seq, seqan::Dna(i % 4));
            }
            // Add deletions, if any.
            for (; j < realBaseCount[i]; ++j)
                appendOperation(cigar, 'D');
            // Simulate qualities if configured to do so.
            if (roche454Options.simulateQualities)
            {
                // Compute likelihood for calling the bases, given this intensity and the Phred score from this.
                double densitySum = 0;
                for (unsigned j = 0; j <= std::max(4u, 2 * MAX_HOMOPOLYMER_LEN); ++j)  // Anecdotally through plot in maple: Enough to sum up to 4 or 2 times the maximal homopolymer length.
                    densitySum += dispatchDensityFunction(model.thresholdMatrix, j, *it);
                double x = 0;  // Probability of seeing < (j+1) bases.
                for (unsigned j = 0; j < calledBaseCount; ++j) {
                    x += dispatchDensityFunction(model.thresholdMatrix, j, *it);
                    int q = -static_cast<int>(10 * ::std::log10(x / densitySum));
                    q = std::max(0, std::min(40, q));
                    appendValue(quals, (char)('!' + q));
                }
            }
        }

        // Write out sequencing information info if configured to do so.
        if (roche454Options.embedReadInfo)
        {
            info.cigar = cigar;
            unsigned len = 0;
            _getLengthInRef(cigar, len);
            if (dir == LEFT)
                info.sampleSequence = prefix(frag, len);
            else
                info.sampleSequence = suffix(frag, length(frag) - len);
            info.isForward = (strand == FORWARD);
            if (strand == REVERSE)
                reverseComplement(info.sampleSequence);
        }
    }
};

// Sanger read simulation.

class SangerSequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Sanger sequencing.
    SangerSequencingOptions sangerOptions;

    SangerSequencingSimulator(TRng & rng,
                              SangerSequencingOptions const & sangerOptions) :
            SequencingSimulator(rng, sangerOptions), sangerOptions(sangerOptions)
    {}

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength()
    {
        if (sangerOptions.readLengthIsUniform)
        {
            // Pick uniformly.
            double minLen = sangerOptions.readLengthMin;
            double maxLen = sangerOptions.readLengthMax;
            double len = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double> >(minLen, maxLen));
            return static_cast<unsigned>(round(len));
        }
        else
        {
            // Pick normally distributed.
            double len = pickRandomNumber(rng, seqan::Pdf<seqan::Normal>(sangerOptions.readLengthMean, sangerOptions.readLengthError));
            return static_cast<unsigned>(round(len));
        }
    }

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand)
    {
        // TODO(holtgrew): Pick better names for read length here.
        // Pick sampled length.
        unsigned readLength = this->readLength();

        if (readLength > length(frag))
        {
            throw std::runtime_error("Sanger read is too long, increase fragment length");
        }

        // Simulate CIGAR string.
        TCigarString cigar;
        this->_simulateCigar(cigar, readLength);

        // Simulate sequence (materialize mismatches and insertions).
        typedef seqan::ModifiedString<seqan::ModifiedString<TFragment, seqan::ModView<seqan::FunctorComplement<seqan::Dna5> > >, seqan::ModReverse> TRevCompFrag;
        if ((dir == LEFT) && (strand == FORWARD))
            _simulateSequence(seq, prefix(frag, readLength), cigar);
        else if ((dir == LEFT) && (strand == REVERSE))
            _simulateSequence(seq, TRevCompFrag(prefix(frag, readLength)), cigar);
        else if ((dir == RIGHT) && (strand == FORWARD))
            _simulateSequence(seq, suffix(frag, length(frag) - readLength), cigar);
        else  // ((dir == RIGHT) && (strand == REVERSE))
            _simulateSequence(seq, TRevCompFrag(suffix(frag, length(frag) - readLength)), cigar);

        // Simulate Qualities.
        this->_simulateQualities(quals, cigar, readLength);

        // Reverse qualities if necessary.
        if (strand == REVERSE)
            reverse(quals);

        // Write out sequencing information info if configured to do so.
        if (sangerOptions.embedReadInfo)
        {
            info.cigar = cigar;
            unsigned len = 0;
            _getLengthInRef(cigar, len);
            if (dir == LEFT)
                info.sampleSequence = prefix(frag, len);
            else
                info.sampleSequence = suffix(frag, length(frag) - len);
            info.isForward = (strand == FORWARD);
            if (strand == REVERSE)
                reverseComplement(info.sampleSequence);
        }
    }

    // Simulate CIGAR string.  We can do this with position specific parameters only and thus independent of any
    // context.
    void _simulateCigar(TCigarString & cigar, unsigned readLength)
    {
        clear(cigar);

        for (unsigned i = 0; i < readLength;)
        {
            double x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double> >(0, 1));
            double pos = 1.0 * i / (readLength - 1);
            double pMismatch = sangerOptions.probabilityMismatchBegin + pos * (sangerOptions.probabilityMismatchEnd - sangerOptions.probabilityMismatchBegin);
            double pInsert   = sangerOptions.probabilityInsertBegin + pos * (sangerOptions.probabilityInsertEnd - sangerOptions.probabilityInsertBegin);
            double pDelete   = sangerOptions.probabilityDeleteBegin + pos * (sangerOptions.probabilityDeleteEnd - sangerOptions.probabilityDeleteBegin);
            double pMatch    = 1.0 - pMismatch - pInsert - pDelete;

            // Simulate mutation/insertion/deletion events.  If possible we reuse the last CIGAR entry.  Adjacent
            // insertion/deletion pairs cancel each other out.  We count i up to the input read length, thus using the
            // "second" member of appendOperation()'s result.

            // TODO(holtgrew): No indels at beginning or end of read.

            if (x < pMatch)  // match
                i += appendOperation(cigar, 'M').second;
            else if (x < pMatch + pMismatch)  // point polymorphism
                i += appendOperation(cigar, 'X').second;
            else if (x < pMatch + pMismatch + pInsert) // insertion
                i += appendOperation(cigar, 'I').second;
            else  // deletion
                i += appendOperation(cigar, 'D').second;
        }
    }

    // TODO(holtgrew): Copy-and-paste from IlluminaSequencingSimulator.
    //
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

    void _simulateQualities(TQualities & quals, TCigarString const & cigar, unsigned sampleLength)
    {
        clear(quals);

        unsigned pos = 0;   // Position in result.
        unsigned rPos = 0;  // Position in fragment.
        for (unsigned i = 0; i < length(cigar); ++i)
        {
            for (unsigned j = 0; j < cigar[i].count; ++j, ++pos)
            {
                double mean = 0.0, stdDev = 0.0;
                double relPos = 1.0 * rPos / (1.0 * sampleLength);

                rPos += (cigar[i].operation != 'D');

                if (cigar[i].operation == 'D')
                {
                    continue;  // No quality to give out.
                }
                else if (cigar[i].operation == 'I' || cigar[i].operation == 'X')
                {
                    mean = sangerOptions.qualityMatchStartMean + pos * (sangerOptions.qualityMatchEndMean - sangerOptions.qualityMatchStartMean);
                    stdDev = sangerOptions.qualityMatchStartStdDev + pos * (sangerOptions.qualityMatchEndStdDev - sangerOptions.qualityMatchStartStdDev);
                }
                else  // cigar[i].operation == 'M'
                {
                    mean = sangerOptions.qualityErrorStartMean + pos * (sangerOptions.qualityErrorEndMean - sangerOptions.qualityErrorStartMean);
                    stdDev = sangerOptions.qualityErrorStartStdDev + pos * (sangerOptions.qualityErrorEndStdDev - sangerOptions.qualityErrorStartStdDev);
                }

                seqan::Pdf<seqan::Normal> pdf(mean, stdDev);
                int q = static_cast<int>(pickRandomNumber(rng, pdf));
                q = std::max(0, std::min(40, q));
                appendValue(quals, (char)('!' + q));
            }
        }
    }
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
