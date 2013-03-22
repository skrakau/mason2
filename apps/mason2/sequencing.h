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
    enum ReadOrientation
    {
        FORWARD_FORWARD,  // --> -->
        FORWARD_REVERSE,  // --> <--
        REVERSE_FORWARD   // <-- -->
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
    // Whether to simulate from forward/reverse strand or both.
    SourceStrands strands;
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
            probabilityInsert(0),
            probabilityDelete(0),
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
    {}
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
            SequencingOptions(), lengthModel(UNIFORM), minReadLength(0), maxReadLength(0),
            meanReadLength(0), stdDevReadLength(0), sqrtInStdDev(true), k(0),
            backgroundNoiseMean(0), backgroundNoiseStdDev(0)
    {}
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
    SequenceingOptions options;

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
