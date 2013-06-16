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
#include <seqan/seq_io.h>

#include "fstream_temp.h"
#include "external_split_merge.h"
#include "mason_options.h"

// ============================================================================
// Forwards
// ============================================================================

typedef seqan::Rng<seqan::MersenneTwister> TRng;

inline void trimAfterSpace(seqan::CharString & s);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

#if 0
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
//
// Allows the loading of all information for one contig at a time.

class Genome
{
public:
    // FAI index for loading the sequences.
    seqan::FaiIndex seqFaiIndex;
    // FAI index for loading the methylation levels.
    seqan::FaiIndex lvlFaiIndex;

    // ID of currently loaded contig.
    unsigned rID;
    // ID and sequence of currently loaded contig.
    seqan::CharString id;
    seqan::Dna5String seq;

    // Whether or not BS-Seq simulation is enabled or not.
    bool bsSimEnabled;
    // Forward and reverse strand methylation levels.
    seqan::CharString fwdMethLevels, revMethLevels;

    // Constructor.
    Genome(bool bsSimEnabled) : rID(seqan::maxValue<unsigned>()), bsSimEnabled(bsSimEnabled)
    {}

    // Loads the given contig from the genome.
    int loadContig(unsigned idx)
    {
        id = sequenceName(seqFaiIndex, idx);
        trimAfterSpace(id);

        // Load sequence.
        if (readSequence(seq, seqFaiIndex, idx) != 0)
            return 1;

        // Load methylation levels if BS-Seq is enabled.
        if (bsSimEnabled)
        {
            if (_readMethLevels(fwdMethLevels, "/TOP") != 0)
                return 1;
            if (_readMethLevels(revMethLevels, "/BOT") != 0)
                return 1;
        }
        return 0;
    }

    // Load methylation levels (/TOP or /BOT) into lvls.
    int _readMethLevels(seqan::CharString & lvls, char const * suffix)
    {
        // Load methylation levels.
        seqan::CharString fullId = id;
        append(fullId, suffix);
        unsigned idx = 0;
        if (!getIdByName(lvlFaiIndex, fullId, idx))
        {
            std::cerr << "ERROR: Methylation levels \"" << fullId << "\" not found for " << id << "\n";
            return 1;
        }
        if (readSequence(lvls, lvlFaiIndex, idx) != 0)
        {
            std::cerr << "ERROR: Problem reading \"" << fullId << "\"\n";
            return 1;
        }
        return 0;
    }

    // Returns methylation level in [0.0, 1.0] for current reference at given position on forward/reverse strand.
    float levelAt(int pos, bool reverse) const
    {
        char c = '\0';
        if (reverse)
            c = revMethLevels[pos];
        else
            c = fwdMethLevels[pos];
        SEQAN_ASSERT_NEQ(c, '>');
        if (c > '>')
            --c;
        return (c - '!') / 80.0 * 1.25;
    }
};
#endif  // #if 0

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
// Class FragmentSamplerImpl
// ----------------------------------------------------------------------------

// Abstract base class for the fragment generation classes.

class FragmentSamplerImpl
{
public:
    virtual void generate(Fragment & frag, int rId, unsigned contigLength) = 0;
    virtual void generateMany(std::vector<Fragment> & frags, int rId, unsigned contigLength,
                              unsigned count) = 0;
};


// ----------------------------------------------------------------------------
// Class UniformFragmentSamplerImpl
// ----------------------------------------------------------------------------

// Generation of fragments with uniform length.

class UniformFragmentSamplerImpl : public FragmentSamplerImpl
{
public:
    // Minimal and maximal fragment length.
    int minLength;
    int maxLength;

    // The random number generator to use.
    TRng & rng;
    // The probability density function for the simulation.
    seqan::Pdf<seqan::Uniform<int> > pdf;

    UniformFragmentSamplerImpl(TRng & rng, int minLength, int maxLength) :
            minLength(minLength), maxLength(maxLength), rng(rng), pdf(minLength, maxLength)
    {}

    virtual void generate(Fragment & frag, int rId, unsigned contigLength);

    virtual void generateMany(std::vector<Fragment> & frags, int rId, unsigned contigLength, unsigned count);

    void _generate(Fragment & frag, int rId, unsigned contigLength);
};

// ----------------------------------------------------------------------------
// Class NormalFragmentSamplerImpl
// ----------------------------------------------------------------------------

// Generation of fragments with normally distributed length.

class NormalFragmentSamplerImpl : public FragmentSamplerImpl
{
public:
    // Mean length and length standard deviation.
    int meanLength;
    int stdDevLength;

    // The random number generator to use.
    TRng & rng;
    // The probability density function for the simulation.
    seqan::Pdf<seqan::Normal> pdf;

    NormalFragmentSamplerImpl(TRng & rng, int meanLength, int stdDevLength) :
            meanLength(meanLength), stdDevLength(stdDevLength), rng(rng), pdf(meanLength, stdDevLength)
    {}

    virtual void generate(Fragment & frag, int rId, unsigned contigLength);

    virtual void generateMany(std::vector<Fragment> & frags, int rId, unsigned contigLength, unsigned count);

    void _generate(Fragment & frag, int rId, unsigned contigLength);
};

// ----------------------------------------------------------------------------
// Class FragmentSampler
// ----------------------------------------------------------------------------

// Generator for sequence fragments.

class FragmentSampler
{
public:
    // Configuration for the generator.
    FragmentSamplerOptions options;

    // The actual generator implementation to use.
    std::auto_ptr<FragmentSamplerImpl> impl;

    FragmentSampler(TRng & rng, FragmentSamplerOptions const & options) : options(options)
    {
        if (options.model == FragmentSamplerOptions::UNIFORM)
            impl.reset(new UniformFragmentSamplerImpl(rng, options.minFragmentSize,
                                                        options.maxFragmentSize));
        else
            impl.reset(new NormalFragmentSamplerImpl(rng, options.meanFragmentSize,
                                                       options.stdDevFragmentSize));
    }

    // TODO(holtgrew): Ultimately, the generator should know about the contig and simulate biases etc.
    
    // Generate a fragment given a contig id and the length of the contig.
    void generate(Fragment & frag, int rId, unsigned contigLength)
    {
        impl->generate(frag, rId, contigLength);
    }

    // Generate multiple fragments.
    void generateMany(std::vector<Fragment> & frags, int rId, unsigned contigLength, unsigned count)
    {
        impl->generateMany(frags, rId, contigLength, count);
    }
};

#if 0
class FragmentSimulator
{
public:
    // The type to use for storing parts of Genome::seqs.
    typedef seqan::Value<seqan::StringSet<seqan::CharString> >::Type TGenomeSeq;

    // The random number generator.
    TRng & rng;
    // The SequenceStream to write the resulting fragment sequences to.
    seqan::SequenceStream & outputStream;
    // The genome to access sequence with.
    Genome & genome;
    // Simulation options for FragmentSampler.
    FragmentOptions fragOptions;

    // Partial sums of the contig lengths.  Used for picking the contig to simulate from.
    seqan::String<int> lengthSums;

    FragmentSampler fragGenerator;

    FragmentSimulator(TRng & rng,
                      seqan::SequenceStream & outStream,
                      Genome & genome,
                      FragmentOptions const & fragOptions) :
            rng(rng), outputStream(outStream), genome(genome), fragOptions(fragOptions),
            fragGenerator(rng, fragOptions)
    {
        // Build the partial sums.
        clear(lengthSums);
        for (unsigned i = 0; i < numSeqs(genome.seqFaiIndex); ++i)
        {
            appendValue(lengthSums, sequenceLength(genome.seqFaiIndex, i));
            if (i > 0u)
                lengthSums[i] += lengthSums[i - 1];
        }
    }

    // Simulate the fragments.
    void simulate();

    // Simulate bisulphite treatment of the given sequence with the location information given by the fragment.
    void _simulateBSTreatment(TGenomeSeq & fragSeq, Fragment const & frag, bool reverse)
    {
        for (int pos = 0; pos != frag.endPos - frag.beginPos; ++pos)
        {
            if ((!reverse && fragSeq[pos] != 'C') || (reverse && fragSeq[pos] != 'G'))  // Skip all non-cyteline chars
            {
                SEQAN_ASSERT_EQ_MSG(genome.levelAt(pos + frag.beginPos, reverse), 0.0,
                                    "Methylation for non-C should be 0 (frag.rId=%d, pos+frag.beginPos=%d, reverse=%d",
                                    frag.rId, pos + frag.beginPos, reverse);
                continue;
            }

            // Decide whether fragSeq[pos] is methylated.  If this is the case then we leave it untouched.
            seqan::Pdf<seqan::Uniform<double> > pdf(0, 1);
            if (pickRandomNumber(rng, pdf) < genome.levelAt(pos + frag.beginPos, reverse))
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
        seqan::Pdf<seqan::Uniform<int> > pdf(0, back(lengthSums) - 1);
        int result = 0;
        int x = pickRandomNumber(rng, pdf);
        for (unsigned i = 0; i < length(lengthSums); ++i)
        {
            if (x >= lengthSums[i])
                result = i + 1;
            if (x < lengthSums[i])
                break;
        }
        return result;
    }
};
#endif  // #if 0

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#if 0

void FragmentSimulator::simulate()
{
    // ------------------------------------------------------------------------
    // (1) Distribute Fragments to Contigs (ids only).
    // ------------------------------------------------------------------------

    // Distributing of ids to per-contig files.
    IdSplitter idSplitter(numSeqs(genome.seqFaiIndex));
    idSplitter.open();

    for (int fragID = 0; fragID < fragOptions.numFragments; ++fragID)
    {
        int contigID = _pickContig();
        if (fwrite(&fragID, 1, sizeof(int), idSplitter.files[contigID]) != sizeof(int))
        {
            std::cerr << "ERROR: Problem writing fragment id in first step for splitting.\n";
            exit(1);
        }
    }

    // ------------------------------------------------------------------------
    // (2) Perform Simulation for Each Contig.
    // ------------------------------------------------------------------------

    idSplitter.reset();

    // Collect the fragments in a per-contig file.
    IdSplitter fragSplitter(numSeqs(genome.seqFaiIndex));
    fragSplitter.open();

    // String stream for building fragment ids.
    std::stringstream ss;

    // We read chunk-wise from the id file.
    unsigned const CHUNK_SIZE = 4096;
    std::vector<int> fragIDs;
    for (unsigned contigID = 0; contigID < numSeqs(genome.seqFaiIndex); ++contigID)
    {
        if (genome.loadContig(contigID) != 0)
        {
            std::cerr << "ERROR: Could not switch to contig " << contigID << "\n";
            exit(1);
        }
        
        fragIDs.resize(CHUNK_SIZE);
        while (!feof(idSplitter.files[contigID]))
        {
            unsigned numRead = fread(&fragIDs[0], sizeof(int), CHUNK_SIZE, idSplitter.files[contigID]);
            fragIDs.resize(numRead);
            for (std::vector<int>::const_iterator it = fragIDs.begin(); it != fragIDs.end(); ++it)
            {
                // Reset string stream and start with id.
                ss.str("");
                ss.clear();
                ss << (*it + 1);

                Fragment frag;
                fragGenerator.generate(frag, contigID, length(genome.seq));

                if (fragOptions.embedSamplingInfo)
                    ss << " REF=" << genome.id << " BEGIN=" << frag.beginPos
                       << " END=" << frag.endPos;

                seqan::Pdf<seqan::Uniform<int> > pdfCoin(0, 1);
                bool reverse = pickRandomNumber(rng, pdfCoin);
                if (fragOptions.embedSamplingInfo)
                {
                    if (reverse)
                        ss << " STRAND=FWD";
                    else
                        ss << " STRAND=REV";
                }
                
                // Get fragment from forward/reverse strand.
                TGenomeSeq fragSeq = infix(genome.seq, frag.beginPos, frag.endPos);
                
                // In case of BS-Seq, compute result of BS treatment now.
                if (fragOptions.bsSimEnabled)
                    _simulateBSTreatment(fragSeq, frag, reverse);
                
                // Write out resulting sequence after reverse-complementation.
                if (reverse)
                    reverseComplement(fragSeq);
                
                // TODO(holtgrew): Check return value?
                seqan::SequenceOutputOptions outOptions(0);
                writeRecord(fragSplitter.files[contigID], ss.str(), fragSeq, seqan::Fasta(), outOptions);
            }
        }
    }

    // ------------------------------------------------------------------------
    // (2) Join Fragments on Names (up to first whitespace).
    // ------------------------------------------------------------------------

    fragSplitter.reset();

    // Buffer for ID and sequence to load data into.
    seqan::CharString id, seq;

    // Join the split FASTA data by id.
    FastaJoiner fragJoiner(fragSplitter);
    while (!fragJoiner.atEnd())
    {
        if (fragJoiner.get(id, seq) != 0)
        {
            std::cerr << "ERROR: Problem when joining.\n";
            exit(1);
        }

        seqan::SequenceOutputOptions outOptions(0);
        if (writeRecord(outputStream, id, seq, outOptions) != 0)
        {
            std::cerr << "ERROR: Problem when writing.\n";
            exit(1);
        }
    }
}

#endif  // #if 0

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

// TODO(holtgrew): Too much redundancy here.

void UniformFragmentSamplerImpl::generate(Fragment & frag, int rId, unsigned contigLength)
{
    _generate(frag, rId, contigLength);
}

void UniformFragmentSamplerImpl::generateMany(std::vector<Fragment> & frags, int rId,
                                                unsigned contigLength, unsigned count)
{
    frags.resize(count);
    for (unsigned i = 0; i < count; ++i)
        _generate(frags[0], rId, contigLength);
}

void UniformFragmentSamplerImpl::_generate(Fragment & frag, int rId, unsigned contigLength)
{
    int fragLength = 0;
    unsigned const MAX_TRIES = 1000;
    for (unsigned tryNo = 0; tryNo < MAX_TRIES; ++tryNo)
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

void NormalFragmentSamplerImpl::generate(Fragment & frag, int rId, unsigned contigLength)
{
    _generate(frag, rId, contigLength);
}

void NormalFragmentSamplerImpl::generateMany(std::vector<Fragment> & frags, int rId,
                                               unsigned contigLength, unsigned count)
{
    frags.resize(count);
    for (unsigned i = 0; i < count; ++i)
        _generate(frags[i], rId, contigLength);
}

void NormalFragmentSamplerImpl::_generate(Fragment & frag, int rId, unsigned contigLength)
{
    int fragLength = 0;
    unsigned const MAX_TRIES = 1000;
    for (unsigned tryNo = 0; tryNo < MAX_TRIES; ++tryNo)
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
