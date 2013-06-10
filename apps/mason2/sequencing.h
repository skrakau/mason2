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
#include <seqan/seq_io.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/random.h>

// ============================================================================
// Forwards
// ============================================================================

class IlluminaModel;
class Roche454Model;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

typedef seqan::Dna5String TRead;
typedef seqan::CharString TQualities;
typedef seqan::Rng<seqan::MersenneTwister> TRng;
typedef seqan::Infix<seqan::Dna5String const>::Type TFragment;
typedef seqan::String<seqan::CigarElement<> > TCigarString;

// ----------------------------------------------------------------------------
// Class SequencingOptions
// ----------------------------------------------------------------------------

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

    // Default constructor.
    SequencingOptions() :
            verbosity(0), numReads(0), simulateQualities(false), simulateMatePairs(false),
            mateOrientation(FORWARD_REVERSE), strands(BOTH), embedReadInfo(false)
    {}

    // Print the options to the given stream.
    virtual void print(std::ostream & out);

    // Convert the boolean to a string for output in print().
    static const char * yesNo(bool b);

    // Convert the MateOrientation to a string for output in print().
    static const char * mateOrientationStr(MateOrientation o);

    // Convert the SourceStrand to a string for output in print().
    static const char * strandsStr(SourceStrands s);
};

// ----------------------------------------------------------------------------
// Class IlluminaSequencingOptions
// ----------------------------------------------------------------------------

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

    virtual void print(std::ostream & out);
};

// ----------------------------------------------------------------------------
// Class SangerSequencingOptions
// ----------------------------------------------------------------------------

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

    virtual void print(std::ostream & out);
};

// ----------------------------------------------------------------------------
// Class Roche454SequencingOptions
// ----------------------------------------------------------------------------

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

    static const char * lengthModelStr(ReadLengthModel m);

    virtual void print(std::ostream & out);
};

// ----------------------------------------------------------------------------
// Class SequencingSimulationInfo
// ----------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------
// Class SequencingSimulator
// ----------------------------------------------------------------------------

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
                           TFragment const & frag);

    // Simulate single-end sequencing from a fragment.
    void simulateSingleEnd(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                           TFragment const & frag);

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

// ----------------------------------------------------------------------------
// Class IlluminaSequencingSimulator
// ----------------------------------------------------------------------------

// Illumina read simulation.

class IlluminaSequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Illumina sequencing.
    IlluminaSequencingOptions illuminaOptions;

    // Storage for the Illumina simulation.
    std::auto_ptr<IlluminaModel> model;

    IlluminaSequencingSimulator(TRng & rng, IlluminaSequencingOptions const & illuminaOptions);

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength()
    {
        return illuminaOptions.readLength;
    }

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand);

private:
    // Initialize the model.
    void _initModel();

    // Simulate PHRED qualities from the CIGAR string.
    void _simulateQualities(TQualities & quals, TCigarString const & cigar);

    // Simulate CIGAR string.  We can do this with position specific parameters only and thus independent of any
    // context.
    void _simulateCigar(TCigarString & cigar);
};

// ----------------------------------------------------------------------------
// Class Roche454SequencingSimulator
// ----------------------------------------------------------------------------

// 454 read simulation.

class Roche454SequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Roche454 sequencing.
    Roche454SequencingOptions roche454Options;

    // Precomputed model data for 454 Sequencing.
    std::auto_ptr<Roche454Model> model;

    Roche454SequencingSimulator(TRng & rng, Roche454SequencingOptions const & roche454Options);

    // Initialize the threshold matrix.
    void _initModel();

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength();

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand);
};

// ----------------------------------------------------------------------------
// Class SangerSequencingSimulator
// ----------------------------------------------------------------------------

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
    virtual unsigned readLength();

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand);

    // Simulate CIGAR string.  We can do this with position specific parameters only and thus independent of any
    // context.
    void _simulateCigar(TCigarString & cigar, unsigned readLength);

    void _simulateQualities(TQualities & quals, TCigarString const & cigar, unsigned sampleLength);
};

// ----------------------------------------------------------------------------
// Class SequencingSimulatorFactory
// ----------------------------------------------------------------------------

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
            rng(rng), in(in), out(out), outR(outR), options(&options), tech(tech)
    {}

    std::auto_ptr<SequencingSimulator> make();
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function appendOrientation()
// ----------------------------------------------------------------------------

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

#endif  // #ifndef SANDBOX_MASON2_APPS_MASON2_SEQUENCING_H_
