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
// Simulation of the fragmentation step.
//
// DNA is shattered into fragments after PCR.  These fragments are then
// sequenced from one or both sides to form single/paired reads.  At the
// moment, we assume that no errors are introduced in PCR and also ignore any
// biases for PCR or fragment selection.
// ==========================================================================

// TODO(holtgrew): We should not load the whole genome into memory, especially bad if we have multiple haplotypes and methylation levels.

#ifndef SANDBOX_MASON2_APPS_MASON2_FRAGMENT_GENERATION_H_
#define SANDBOX_MASON2_APPS_MASON2_FRAGMENT_GENERATION_H_

#include <sstream>
#include <fstream>
#include <vector>
#include <memory>

#include <seqan/random.h>
#include <seqan/sequence.h>

#include "fstream_temp.h"

// ============================================================================
// Forwards
// ============================================================================

typedef seqan::Rng<seqan::MersenneTwister> TRng;

inline void trimAfterSpace(seqan::CharString & s);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class FragmentOptions
// ----------------------------------------------------------------------------

// Configuration for the simulation of fragments.

struct FragmentOptions
{
    // Type for selecting the fragment size distribution: uniformly or normally distributed.
    enum FragmentSizeModel
    {
        UNIFORM,
        NORMAL
    };

    // The BS simulation type to use.
    enum BSProtocol
    {
        DIRECTIONAL,
        UNDIRECTIONAL
    };

    // Number of fragments to generate.
    int numFragments;
    // Whether or not to embed sampling information in id.
    bool embedSamplingInfo;
    
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

    // Wether or not to simulate BS-seq.
    bool bsSimEnabled;
    // The rate that unmethylated Cs to become Ts.
    double bsConversionRate;
    // The protocol to use for the simulation.
    BSProtocol bsProtocol;

    FragmentOptions() :
            numFragments(0), embedSamplingInfo(false), minFragmentSize(0), maxFragmentSize(0), meanFragmentSize(0),
            stdDevFragmentSize(0), model(UNIFORM), bsSimEnabled(false), bsConversionRate(0), bsProtocol(DIRECTIONAL)
    {}
};

// --------------------------------------------------------------------------
// Class Genome
// --------------------------------------------------------------------------

// Container for genomic information.

class Genome
{
public:
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;

    // Forward and reverse methylation levels.
    seqan::StringSet<seqan::CharString> fwdMethLevels, revMethLevels;

    void trimIds()
    {
        for (unsigned i = 0; i < length(ids); ++i)
            trimAfterSpace(ids[i]);
    }

    // Returns methylation level in [0.0, 1.0] for given reference at given position on forward/reverse strand.
    float levelAt(int rId, int pos, bool reverse) const
    {
        char c = '\0';
        if (reverse)
            c = revMethLevels[rId][pos];
        else
            c = fwdMethLevels[rId][pos];
        SEQAN_ASSERT_NEQ(c, '>');
        if (c > '>')
            --c;
        return (c - '!') / 80.0 * 1.25;
    }

    int loadMethLevels(seqan::FaiIndex const & faiIndex)
    {
        resize(fwdMethLevels, length(ids));
        resize(revMethLevels, length(ids));
        
        for (unsigned i = 0; i < length(ids); ++i)
        {
            seqan::CharString topId = ids[i];
            append(topId, "/TOP");
            unsigned idx = 0;
            if (!getIdByName(faiIndex, topId, idx))
            {
                std::cerr << "ERROR: Methylation levels \"" << topId << "\" not found for " << ids[i] << "\n";
                return 1;
            }
            if (readSequence(fwdMethLevels[i], faiIndex, idx) != 0)
            {
                std::cerr << "ERROR: Problem reading \"" << topId << "\"\n";
                return 1;
            }

            seqan::CharString bottomId = ids[i];
            append(bottomId, "/BOT");
            if (!getIdByName(faiIndex, bottomId, idx))
            {
                std::cerr << "ERROR: Methylation levels \"" << bottomId << "\" not found for " << ids[i] << "\n";
                return 1;
            }
            if (readSequence(revMethLevels[i], faiIndex, idx) != 0)
            {
                std::cerr << "ERROR: Problem reading \"" << bottomId << "\"\n";
                return 1;
            }
        }
    }
};

// ----------------------------------------------------------------------------
// Class Fragment
// ----------------------------------------------------------------------------

class Fragment
{
public:
    // Index of contig to simulate on.
    int rId;
    // Begin and end position of the fragment.
    int beginPos, endPos;

    Fragment() : rId(-1), beginPos(0), endPos(0)
    {}
};

// ----------------------------------------------------------------------------
// Class FragmentGeneratorImpl
// ----------------------------------------------------------------------------

// Abstract base class for the fragment generation classes.

class FragmentGeneratorImpl
{
public:
    virtual void generate(Fragment & frag, int rId, unsigned contigLength) = 0;
};

// ----------------------------------------------------------------------------
// Class UniformFragmentGeneratorImpl
// ----------------------------------------------------------------------------

// Generation of fragments with uniform length.

class UniformFragmentGeneratorImpl : public FragmentGeneratorImpl
{
public:
    // Minimal and maximal fragment length.
    int minLength;
    int maxLength;

    // The random number generator to use.
    TRng & rng;
    // The probability density function for the simulation.
    seqan::Pdf<seqan::Uniform<int> > pdf;

    UniformFragmentGeneratorImpl(TRng & rng, int minLength, int maxLength) :
            rng(rng), minLength(minLength), maxLength(maxLength), pdf(minLength, maxLength)
    {}

    virtual void generate(Fragment & frag, int rId, unsigned contigLength);
};

// ----------------------------------------------------------------------------
// Class NormalFragmentGeneratorImpl
// ----------------------------------------------------------------------------

// Generation of fragments with normally distributed length.

class NormalFragmentGeneratorImpl : public FragmentGeneratorImpl
{
public:
    // Mean length and length standard deviation.
    int meanLength;
    int stdDevLength;

    // The random number generator to use.
    TRng & rng;
    // The probability density function for the simulation.
    seqan::Pdf<seqan::Normal> pdf;

    NormalFragmentGeneratorImpl(TRng & rng, int meanLength, int stdDevLength) :
            rng(rng), meanLength(meanLength), stdDevLength(stdDevLength), pdf(meanLength, stdDevLength)
    {}

    virtual void generate(Fragment & frag, int rId, unsigned contigLength);
};

// ----------------------------------------------------------------------------
// Class FragmentGenerator
// ----------------------------------------------------------------------------

// Generator for sequence fragments.

class FragmentGenerator
{
public:
    // Configuration for the generator.
    FragmentOptions options;

    // The actual generator implementation to use.
    std::auto_ptr<FragmentGeneratorImpl> impl;

    FragmentGenerator(TRng & rng, FragmentOptions const & options) : options(options)
    {
        if (options.model == FragmentOptions::UNIFORM)
            impl.reset(new UniformFragmentGeneratorImpl(rng, options.minFragmentSize,
                                                        options.maxFragmentSize));
        else
            impl.reset(new NormalFragmentGeneratorImpl(rng, options.meanFragmentSize,
                                                       options.stdDevFragmentSize));
    }

    // TODO(holtgrew): Ultimately, the generator should know about the contig and simulate biases etc.
    
    // Generate a fragment given a contig id and the length of the contig.
    void generate(Fragment & frag, int rId, unsigned contigLength)
    {
        impl->generate(frag, rId, contigLength);
    }
};

class FragmentSimulator
{
public:
    // The type to use for storing parts of Genome::seqs.
    typedef seqan::Value<seqan::StringSet<seqan::CharString> >::Type TGenomeSeq;

    // The random number generator.
    TRng & rng;
    // The SequenceStream to write the resulting fragment sequences to.
    seqan::SequenceStream & outputStream;
    // The genome to sequence from;
    Genome const & genome;
    // Simulation options for FragmentGenerator.
    FragmentOptions fragOptions;

    // Partial sums of the contig lengths.  Used for picking the contig to simulate from.
    seqan::String<int> lengthSums;

    FragmentGenerator fragGenerator;

    FragmentSimulator(TRng & rng,
                      seqan::SequenceStream & outStream,
                      Genome const & genome,
                      FragmentOptions const & fragOptions) :
            rng(rng), outputStream(outStream), genome(genome), fragOptions(fragOptions),
            fragGenerator(rng, fragOptions)
    {
        // Build the partial sums.
        resize(lengthSums, length(genome.seqs), 0);
        for (unsigned i = 0; i < length(genome.seqs); ++i)
        {
            appendValue(lengthSums, length(genome.seqs[i]));
            if (i > 0u)
                lengthSums[i] += lengthSums[i - 1];
        }
    }

    // Simulate the fragments.
    void simulate()
    {
        seqan::CharString id;
        seqan::Dna5String contig, seq;

        std::stringstream ss;

        for (int i = 0; i < fragOptions.numFragments; ++i)
        {
            ss.str("");
            ss.clear();
            ss << (i + 1);

            Fragment frag;
            int rId = _pickContig();
            fragGenerator.generate(frag, rId, length(genome.seqs[rId]));

            if (fragOptions.embedSamplingInfo)
                ss << " REF=" << genome.ids[rId] << " BEGIN=" << frag.beginPos << " END=" << frag.endPos;

            seqan::Pdf<seqan::Uniform<int> > pdfCoin(0, 1);
            bool reverse = pickRandomNumber(rng, pdfCoin);
            if (reverse)
                ss << " STRAND=FWD";
            else
                ss << " STRAND=REV";

            // Get fragment from forward/reverse strand.
            TGenomeSeq fragSeq = infix(genome.seqs[frag.rId], frag.beginPos, frag.endPos);

            // In case of BS-Seq, compute result of BS treatment now.
            if (fragOptions.bsSimEnabled)
                _simulateBSTreatment(fragSeq, frag, reverse);

            // Write out resulting sequence after reverse-complementation.
            if (reverse)
                reverseComplement(fragSeq);

            seqan::SequenceOutputOptions outOptions(0);
            // TODO(holtgrew): Check return value?
            writeRecord(outputStream, ss.str(), fragSeq, outOptions);
        }

        return;
    }

    // Simulate bisulphite treatment of the given sequence with the location information given by the fragment.
    void _simulateBSTreatment(TGenomeSeq & fragSeq, Fragment const & frag, bool reverse)
    {
        for (int pos = 0; pos != frag.endPos - frag.beginPos; ++pos)
        {
            if ((!reverse && fragSeq[pos] != 'C') || (reverse && fragSeq[pos] != 'G'))  // Skip all non-cyteline chars
            {
                SEQAN_ASSERT_EQ_MSG(genome.levelAt(frag.rId, pos + frag.beginPos, reverse), 0.0,
                                    "Methylation for non-C should be 0 (frag.rId=%d, pos+frag.beginPos=%d, reverse=%d",
                                    frag.rId, pos + frag.beginPos, reverse);
                continue;
            }

            // Decide whether fragSeq[pos] is methylated.  If this is the case then we leave it untouched.
            seqan::Pdf<seqan::Uniform<double> > pdf(0, 1);
            if (pickRandomNumber(rng, pdf) < genome.levelAt(frag.rId, pos + frag.beginPos, reverse))
                continue;

            // Otherwise, pick whether we will convert.
            if (pickRandomNumber(rng, pdf) < fragOptions.bsConversionRate)
                fragSeq[pos] = reverse ? 'A' : 'T';
        }
    }

    // Pick the contig to simulate the fragment from.
    //
    // The contig is picked in relation to their length.
    int _pickContig()
    {
        if (length(lengthSums) == 1u)
            return 0;
        seqan::Pdf<seqan::Uniform<int> > pdf(0, length(lengthSums) - 1);
        int result = 0;
        int x = pickRandomNumber(rng, pdf);
        for (unsigned i = 0; i < length(lengthSums); ++i)
        {
            if (x < lengthSums[i])
                result = i;
            if (x >= lengthSums[i])
                break;
        }
        return result;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function trimAfterSpace()
// --------------------------------------------------------------------------

// TODO(holtgrew): Put into header or SeqAn library.

// Trim after the first whitespace.
inline
void trimAfterSpace(seqan::CharString & s)
{
    unsigned i = 0;
    for (; i < length(s); ++i)
        if (isspace(s[i]))
            break;
    resize(s, i);
}

// TODO(holtgrew): The virtual functions below are the same!

void UniformFragmentGeneratorImpl::generate(Fragment & frag, int rId, unsigned contigLength)
{
    // TODO(holtgrew): Test whether there are too many tries.
    int fragLength = 0;
    while (true)
    {
        fragLength = pickRandomNumber(rng, pdf);
        if (fragLength <= 0 || fragLength > (int)contigLength)
            continue;  // Try again

        seqan::Pdf<seqan::Uniform<int> > posPdf(0, contigLength - fragLength);
        int beginPos = pickRandomNumber(rng, posPdf);

        frag.rId = rId;
        frag.beginPos = beginPos;
        frag.endPos = beginPos + fragLength;
        SEQAN_ASSERT_LEQ(frag.endPos, (int)contigLength);
        break;
    }
}

void NormalFragmentGeneratorImpl::generate(Fragment & frag, int rId, unsigned contigLength)
{
    // TODO(holtgrew): Test whether there are too many tries.
    int fragLength = 0;
    while (true)
    {
        fragLength = pickRandomNumber(rng, pdf);
        if (fragLength <= 0 || fragLength > (int)contigLength)
            continue;  // Try again

        seqan::Pdf<seqan::Uniform<int> > posPdf(0, contigLength - fragLength);
        int beginPos = pickRandomNumber(rng, posPdf);

        frag.rId = rId;
        frag.beginPos = beginPos;
        frag.endPos = beginPos + fragLength;
        SEQAN_ASSERT_LEQ(frag.endPos, (int)contigLength);
        break;
    }
}

#endif  // #ifndef SANDBOX_MASON2_APPS_MASON2_FRAGMENT_GENERATION_H_
