// ==========================================================================
//                            methylation_levels.h
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

#ifndef SANDBOX_MASON2_APPS_MASON2_METHYLATION_LEVELS_H_
#define SANDBOX_MASON2_APPS_MASON2_METHYLATION_LEVELS_H_

#include <seqan/sequence.h>

// ============================================================================
// Forwards
// ============================================================================

typedef seqan::Rng<> TRng;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Class MethylationLevels
// --------------------------------------------------------------------------

// Stores methylation levels separately for forward and reverse strand.

struct MethylationLevels
{
    // Forward and reverse levels, encoded as round(level / 0.125) + 33.
    seqan::CharString forward, reverse;

    void resize(unsigned len)
    {
        seqan::resize(forward, len, '!');
        seqan::resize(reverse, len, '!');
    }

    void clear()
    {
        seqan::clear(forward);
        seqan::clear(reverse);
    }

    // Translate character in forward/reverse to level (0..80).
    inline char charToLevel(char c)
    {
        if (c < '>')  // '>' cannot be used as value
            return c - 33;
        else
            return c - 34;
    }

    // Translate level (0..80) to character in forward/reverse.
    inline char levelToChar(char c)
    {
        if (c + '!' < '>')
            return c + 33;
        else
            return c + 34;
    }

    // Returns methylation level for forward strand at position i.
    inline float levelF(unsigned i)
    {
        return (charToLevel(forward[i]) / 80.0) / 0.0125;
    }

    // Sets methylation level for forward strand at position i.
    inline void setLevelF(unsigned i, float level)
    {
        SEQAN_ASSERT_GEQ(level, 0.0);
        SEQAN_ASSERT_LEQ(level, 1.0);
        // std::cerr << "forward[i] = " << levelToChar(round(level / 0.0125)) << " ~ " << (level / 0.0125) << " ~ " << level << "\n";
        char c = levelToChar(round(level / 0.0125));
        SEQAN_ASSERT_NEQ((int)c, (int)'>');
        forward[i] = std::max(reverse[i], c);
    }

    // Returns methylation level for reverse strand at position i.
    inline float levelR(unsigned i)
    {
        return (charToLevel(reverse[i]) / 80.0) / 0.0125;
    }

    // Sets methylation level for reverse strand at position i if level is larger than current.
    inline void setLevelR(unsigned i, float level)
    {
        SEQAN_ASSERT_GEQ(level, 0.0);
        SEQAN_ASSERT_LEQ(level, 1.0);
        // std::cerr << "reverse[" << i << "] = " << levelToChar(round(level / 0.0125)) << " ~ " << (level / 0.0125) << " ~ " << level << "\n";
        char c = levelToChar(round(level / 0.0125));
        SEQAN_ASSERT_NEQ((int)c, (int)'>');
        reverse[i] = std::max(reverse[i], c);
    }
};

// --------------------------------------------------------------------------
// Class MethylationLevelSimulatorOptions
// --------------------------------------------------------------------------

struct MethylationLevelSimulatorOptions
{
    // Enable simulation of methylation levels.
    bool simulateMethylationLevels;
    // Path to write methylation rates to.  There will be one sequence entry for the reference and one entry for each
    // haplotype.  The probability will be encoded in ASCII characters 33-114, at a resolution of 1.25%, ">" is ignored.
    seqan::CharString methFastaOutFile;
    // Median and standard deviation for picking methylation level for all Cs.
    double methMuC, methSigmaC;
    // Median and standard deviation for picking methylation level for CpGs.
    double methMuCG, methSigmaCG;
    // Median and standard deviation for picking methylation level for CHGs.
    double methMuCHG, methSigmaCHG;
    // Median and standard deviation for picking methylation level for CHHs.
    double methMuCHH, methSigmaCHH;

    MethylationLevelSimulatorOptions() :
            simulateMethylationLevels(false), methMuC(0), methSigmaC(0), methMuCG(0), methSigmaCG(0),
            methMuCHG(0), methSigmaCHG(0), methMuCHH(0), methSigmaCHH(0)
    {}
};

// --------------------------------------------------------------------------
// Class MethylationLevelSimulator
// --------------------------------------------------------------------------

// Simulate methylation levels for a Dna sequence/contig on forward and reverse strand.

class MethylationLevelSimulator
{
public:
    // Options for the mu/sigma values.
    MethylationLevelSimulatorOptions const & options;

    // Random number generator to use.
    TRng & rng;

    // Beta probability density functions for level generation.
    seqan::Pdf<seqan::Beta> pdfC, pdfCG, pdfCHG, pdfCHH;

    MethylationLevelSimulator(TRng & rng, MethylationLevelSimulatorOptions const & options) :
            options(options), rng(rng),
            pdfC(options.methMuC, options.methSigmaC, seqan::MeanStdDev()),
            pdfCG(options.methMuCG, options.methSigmaCG, seqan::MeanStdDev()),
            pdfCHG(options.methMuCHG, options.methSigmaCHG, seqan::MeanStdDev()),
            pdfCHH(options.methMuCHH, options.methSigmaCHH, seqan::MeanStdDev())
    {}

    // Simulate methylation levels for the sequence in contig.  The results are stored in levels.
    void run(MethylationLevels & levels, seqan::Dna5String const & contig)
    {
        levels.resize(length(contig));

        typedef seqan::Iterator<seqan::Dna5String const>::Type TContigIter;
        TContigIter it = begin(contig, seqan::Standard());
        TContigIter itEnd = end(contig, seqan::Standard()) - 3;

        // We will go over the contig with hashes to search for patterns efficiently.
        seqan::Shape<seqan::Dna5> shape2, shape3;
        handleOneMer(levels, 0, ordValue(contig[0]));
        if (levels.forward[0] != '!') SEQAN_ASSERT_EQ(contig[0], 'C');
        if (levels.reverse[0] != '!') SEQAN_ASSERT_EQ(contig[0], 'G');
        if (length(contig) >= 2u)
        {
            resize(shape2, 2);
            hash(shape2, it);
            handleTwoMer(levels, 0, value(shape2));
            if (levels.forward[1] != '!') SEQAN_ASSERT_EQ(contig[1], 'C');
            if (levels.reverse[1] != '!') SEQAN_ASSERT_EQ(contig[1], 'G');
        }
        if (length(contig) >= 3u)
        {
            resize(shape3, 3);
            hash(shape3, it);
            handleThreeMer(levels, 0, value(shape3));
            if (levels.forward[2] != '!') SEQAN_ASSERT_EQ(contig[2], 'C');
            if (levels.reverse[2] != '!') SEQAN_ASSERT_EQ(contig[2], 'G');
        }
        ++it;
        unsigned pos = 1;
        for (; (pos + 3 < length(contig)) && (it != itEnd); ++it, ++pos)
        {
            hashNext(shape2, it);
            hashNext(shape3, it);
            handleOneMer(levels, pos, *it);
            handleTwoMer(levels, pos, value(shape2));
            handleThreeMer(levels, pos, value(shape3));
            if (levels.forward[pos] != '!') SEQAN_ASSERT_EQ(contig[pos], 'C');
            if (levels.reverse[pos] != '!') SEQAN_ASSERT_EQ(contig[pos], 'G');
        }
        if (pos < length(contig))
            handleOneMer(levels, pos, ordValue(*it));
        if (pos + 1 < length(contig))
        {
            hashNext(shape2, it++);
            handleTwoMer(levels, pos++, value(shape2));
        }
        if (pos < length(contig))
            handleOneMer(levels, pos, ordValue(*it));
    }

    // Handle 3mer, forward case.
    void handleThreeMer(MethylationLevels & levels, unsigned pos, unsigned hashValue)
    {
        // seqan::Dna5String dbg;
        // unhash(dbg, hashValue, 3);
        switch (hashValue)
        {
            case 27:  // CHG, symmetric
            case 32:
            case 42:
                // std::cerr << "CHG       \t" << dbg << "\n";
                levels.setLevelF(pos, pickRandomNumber(rng, pdfCHG));
                levels.setLevelR(pos + 2, pickRandomNumber(rng, pdfCHG));
                break;
            case 25:  // CHH
            case 26:
            case 28:
            case 30:
            case 31:
            case 33:
            case 40:
            case 41:
            case 43:
                // std::cerr << "CHH       \t" << dbg << "\n";
                levels.setLevelF(pos, pickRandomNumber(rng, pdfCHH));
                break;
            case 2:  // rc(CHH)
            case 12:
            case 17:
            case 52:
            case 62:
            case 67:
            case 77:
            case 87:
            case 92:
                // std::cerr << "rc(CHH)   \t" << dbg << "\n";
                levels.setLevelR(pos + 2, pickRandomNumber(rng, pdfCHH));
                break;
            default:
                // nop
                break;
        }
    }

    // Handle 2mer.
    void handleTwoMer(MethylationLevels & levels, unsigned pos, unsigned hashValue)
    {
        if (hashValue == 7)  // CpG forward (symmetric, also reverse)
        {
            // seqan::Dna5String dbg;
            // unhash(dbg, hashValue, 2);
            // std::cerr << "CpG     \t" << dbg << "\n";
            levels.setLevelF(pos, pickRandomNumber(rng, pdfCG));
            levels.setLevelR(pos + 1, pickRandomNumber(rng, pdfCG));
        }
    }

    // Handle 1mer.
    void handleOneMer(MethylationLevels & levels, unsigned pos, unsigned val)
    {
        // if (val == 1 || val == 2)
        //     std::cerr << "C/G     \t" << seqan::Dna5(val) << "\n";

        if (val == 1)   // C forward
            levels.setLevelF(pos, pickRandomNumber(rng, pdfC));
        else if (val == 2)  // C reverse (G)
            levels.setLevelR(pos, pickRandomNumber(rng, pdfC));
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function VariantMaterializer::_fixVariationLevels()
// ----------------------------------------------------------------------------

// Fix variation levels on contig given the points (second == true -> SNP, second == false -> breakpoint).

void fixVariationLevels(MethylationLevels & levels,
                        TRng & rng,
                        seqan::Dna5String const & contig,
                        seqan::String<std::pair<int, bool> > const & varPoints,
                        MethylationLevelSimulatorOptions const & options)
{
    MethylationLevelSimulator methSim(rng, options);
    seqan::Shape<seqan::Dna5> shape2, shape3;
    resize(shape2, 2);
    resize(shape3, 3);

    for (unsigned i = 0; i < length(varPoints); ++i)
    {
        int pos = varPoints[i].first;
        if (varPoints[i].second)  // is SNP
        {
            if (pos > 2)
            {
                levels.forward[pos - 2] = levels.reverse[pos - 2] = '!';
                methSim.handleOneMer(levels, pos - 2, ordValue(contig[pos - 2]));
                methSim.handleTwoMer(levels, pos - 2, hash(shape2, iter(contig, pos - 2, seqan::Standard())));
                methSim.handleThreeMer(levels, pos - 2, hash(shape3, iter(contig, pos - 2, seqan::Standard())));
            }
            if (pos > 1)
            {
                levels.forward[pos - 1] = levels.reverse[pos - 1] = '!';
                methSim.handleOneMer(levels, pos - 1, ordValue(contig[pos - 1]));
                methSim.handleTwoMer(levels, pos - 1, hash(shape2, iter(contig, pos - 1, seqan::Standard())));
            }
            levels.forward[pos] = levels.reverse[pos] = '!';
            methSim.handleOneMer(levels, pos, ordValue(contig[pos]));
            if (pos + 1 < (int)length(contig))
            {
                levels.forward[pos + 1] = levels.reverse[pos + 1] = '!';
                methSim.handleOneMer(levels, pos + 1, ordValue(contig[pos + 1]));
                methSim.handleTwoMer(levels, pos, hash(shape2, iter(contig, pos, seqan::Standard())));
            }
            if (pos + 2 < (int)length(contig))
            {
                levels.forward[pos + 2] = levels.reverse[pos + 2] = '!';
                methSim.handleOneMer(levels, pos + 2, ordValue(contig[pos + 2]));
                methSim.handleTwoMer(levels, pos + 1, hash(shape2, iter(contig, pos + 1, seqan::Standard())));
                methSim.handleThreeMer(levels, pos, hash(shape3, iter(contig, pos, seqan::Standard())));
            }
        }
        else  // is no SNP but breakpoint
        {
            // TODO(holtgrew): Double-check for correctness, might recompute too much around breakpoints.
            if (pos > 2)
            {
                levels.forward[pos - 2] = levels.reverse[pos - 2] = '!';
                methSim.handleOneMer(levels, pos - 2, ordValue(contig[pos - 2]));
                methSim.handleTwoMer(levels, pos - 2, hash(shape2, iter(contig, pos - 2, seqan::Standard())));
                methSim.handleThreeMer(levels, pos - 2, hash(shape3, iter(contig, pos - 2, seqan::Standard())));
            }
            if (pos > 1)
            {
                levels.forward[pos - 1] = levels.reverse[pos - 1] = '!';
                methSim.handleOneMer(levels, pos - 1, ordValue(contig[pos - 1]));
                methSim.handleTwoMer(levels, pos - 1, hash(shape2, iter(contig, pos - 1, seqan::Standard())));
            }
            levels.forward[pos] = levels.reverse[pos] = '!';
            methSim.handleOneMer(levels, pos, ordValue(contig[pos]));
            if (pos + 1 < (int)length(contig))
            {
                methSim.handleTwoMer(levels, pos, hash(shape2, iter(contig, pos, seqan::Standard())));
                levels.forward[pos + 1] = levels.reverse[pos + 1] = '!';
                methSim.handleOneMer(levels, pos + 1, ordValue(contig[pos + 2]));
            }
            if (pos + 2 < (int)length(contig))
            {
                methSim.handleTwoMer(levels, pos + 1, hash(shape2, iter(contig, pos + 1, seqan::Standard())));
                methSim.handleThreeMer(levels, pos, hash(shape3, iter(contig, pos, seqan::Standard())));
            }
        }
    }
}

#endif  // #ifndef SANDBOX_MASON2_APPS_MASON2_METHYLATION_LEVELS_H_
