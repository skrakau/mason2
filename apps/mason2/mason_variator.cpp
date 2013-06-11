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
// Given a genome, create variations thereof and export them VCF and
// materialized into a FASTA file.  Variations can also be imported from a VCF
// file and materialized into a FASTA file.
// ==========================================================================

// TODO(holtgrew): Currently, inserted sequence is picked at random, we could also give an insertion database.
// TODO(holtgrew): Currently, there only is support for left-to-right translocations.
// TODO(holtgrew): Allow inversion in translocation.
// TODO(holtgrew): Add support for parsing VCF.
// TODO(holtgrew): Simulate different SNPs/small variations for duplications, input for repeat separation.

#include <seqan/arg_parse.h>
#include <seqan/random.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/vcf_io.h>
#include <seqan/sequence_journaled.h>
#include <seqan/index.h>  // for Shape<>

#include "variation_size_tsv.h"
#include "genomic_variants.h"

// ==========================================================================
// Forwards
// ==========================================================================

// ==========================================================================
// Classes
// ==========================================================================

typedef seqan::Rng<> TRng;
typedef seqan::JournalEntries<seqan::JournalEntry<unsigned, int>, seqan::SortedArray> TJournalEntries;

// --------------------------------------------------------------------------
// Class MasonVariatorOptions
// --------------------------------------------------------------------------

struct MasonVariatorOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Seed for RNG.
    int seed;

    // ----------------------------------------------------------------------
    // Input / Output Options
    // ----------------------------------------------------------------------

    // VCF file to import.
    seqan::CharString vcfInFile;
    // FASTA file to import.
    seqan::CharString fastaInFile;
    // VCF file to write out.
    seqan::CharString vcfOutFile;
    // FASTA file to write out with variations.
    seqan::CharString fastaOutFile;

    // Path to a TSV file where the first two columns giving the type of the SV to simulate and the size of the SV.
    // This overrides the simulation of SV from the sv*Rate parameters.
    seqan::CharString inputSVSizeFile;
    // Path to TSV file to write the resulting breakpoints in variant genomes to.
    seqan::CharString outputBreakpointFile;

    // ----------------------------------------------------------------------
    // Haplotype / Allele Configuration
    // ----------------------------------------------------------------------

    // The number of haplotypes to simulate.
    int numHaplotypes;

    // The string to use to separate the haplotype identifier from the chromosome name in the output FASTA ifle.
    seqan::CharString haplotypeSep;

    // ----------------------------------------------------------------------
    // Variation Simulation
    // ----------------------------------------------------------------------

    // Per-base probability for SNPs and small-scale indels.
    double snpRate;
    double smallIndelRate;

    // Minimal and maximal size for small indels.  Indels will be simulated uniformly in this range.  The range is
    // stored internally as [min, max) but given as [min, max] from the command line.
    int minSmallIndelSize;
    int maxSmallIndelSize;
    
    // Per-base probability for having a structural variation.
    double svIndelRate;
    double svInversionRate;
    double svTranslocationRate;
    double svDuplicationRate;

    // Minimal and maximal size for structural variations.  SVs will be simulated uniformly in this range.  The range is
    // stored internally as [min, max) but given as [min, max] from the command line.
    int minSVSize;
    int maxSVSize;

    // ----------------------------------------------------------------------
    // Methylation Simulation
    // ----------------------------------------------------------------------

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

    MasonVariatorOptions() :
            verbosity(1), seed(0),
            snpRate(0), smallIndelRate(0), minSmallIndelSize(0), maxSmallIndelSize(0), svIndelRate(0),
            svInversionRate(0), svTranslocationRate(0), svDuplicationRate(0), minSVSize(0), maxSVSize(0),
            simulateMethylationLevels(false), methMuC(0), methSigmaC(0), methMuCG(0), methSigmaCG(0),
            methMuCHG(0), methSigmaCHG(0), methMuCHH(0), methSigmaCHH(0)
    {}
};

void print(std::ostream & out, MasonVariatorOptions const & options)
{
    out << "__OPTIONS_____________________________________________________________________\n"
        << "\n"
        // << "VCF IN               \t" << options.vcfInFile << "\n"
        << "FASTA IN             \t" << options.fastaInFile << "\n"
        << "SV SIZE TSV IN       \t" << options.inputSVSizeFile << "\n"
        << "VCF OUT              \t" << options.vcfOutFile << "\n"
        << "FASTA OUT            \t" << options.fastaOutFile << "\n"
        << "BREAKPOINT TSV OUT   \t" << options.outputBreakpointFile << "\n"
        << "\n"
        << "NUM HAPLOTYPES       \t" << options.numHaplotypes << "\n"
        << "HAPLOTYPE SEP        \t\"" << options.haplotypeSep << "\"\n"
        << "\n"
        << "SNP RATE             \t" << options.snpRate << "\n"
        << "SMALL INDEL RATE     \t" << options.smallIndelRate << "\n"
        << "\n"
        << "MIN SMALL INDEL SIZE \t" << options.minSmallIndelSize << "\n"
        << "MAX SMALL INDEL SIZE \t" << options.maxSmallIndelSize << "\n"
        << "\n"
        << "SV INDEL RATE        \t" << options.svIndelRate << "\n"
        << "SV INVERSION RATE    \t" << options.svInversionRate << "\n"
        << "SV TRANSLOCATION RATE\t" << options.svTranslocationRate << "\n"
        << "SV DUPLICATION RATE  \t" << options.svDuplicationRate << "\n"
        << "\n"
        << "MIN SV SIZE          \t" << options.minSVSize << "\n"
        << "MAX SV SIZE          \t" << options.maxSVSize << "\n"
        << "\n"
        << "SIM. METHYL. LEVELS  \t" << options.simulateMethylationLevels << "\n"
        << "METHYLATION OUT FILE \t" << options.methFastaOutFile << "\n"
        << "METHYLATION LEVELS\n"
        << "  C   MU             \t" << options.methMuC << "\n"
        << "  C   SIGMA          \t" << options.methSigmaC << "\n"
        << "  CG  MU             \t" << options.methMuCG << "\n"
        << "  CG  SIGMA          \t" << options.methSigmaCG << "\n"
        << "  CHG MU             \t" << options.methMuCHG << "\n"
        << "  CHG SIGMA          \t" << options.methSigmaCHG << "\n"
        << "  CHH MU             \t" << options.methMuCHH << "\n"
        << "  CHH SIGMA          \t" << options.methSigmaCHH << "\n"
        << "\n";
}

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
// Class MethylationLevelSimulator
// --------------------------------------------------------------------------

// Simulate methylation levels for a Dna sequence/contig on forward and reverse strand.

class MethylationLevelSimulator
{
public:
    // Options for the mu/sigma values.
    MasonVariatorOptions const & options;

    // Random number generator to use.
    TRng & rng;

    // Beta probability density functions for level generation.
    seqan::Pdf<seqan::Beta> pdfC, pdfCG, pdfCHG, pdfCHH;

    MethylationLevelSimulator(TRng & rng, MasonVariatorOptions const & options) :
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

// --------------------------------------------------------------------------
// Class StructuralVariantSimulator
// --------------------------------------------------------------------------

// Simulation of structural variants given error rates.

class StructuralVariantSimulator
{
public:
    // Random number generator
    TRng & rng;

    // FAI Index for loading sequence contig-wise.
    seqan::FaiIndex const & faiIndex;

    // The variator options.
    MasonVariatorOptions options;

    // Structural variation records.
    seqan::String<VariationSizeRecord> const & variationSizeRecords;
    seqan::String<int> variationToContig;

    StructuralVariantSimulator(TRng & rng, seqan::FaiIndex const & faiIndex,
                               seqan::String<VariationSizeRecord> const & variationSizeRecords,
                               MasonVariatorOptions const & options) :
            rng(rng), faiIndex(faiIndex), options(options), variationSizeRecords(variationSizeRecords)
    {
        _distributeVariations();
    }

    // Distribute variations to contigs.
    void _distributeVariations()
    {
        // Build prefix sume for distributing variations to contig proportional to the length.
        seqan::String<__int64> limits;
        appendValue(limits, 0);
        for (unsigned i = 0; i < numSeqs(faiIndex); ++i)
            appendValue(limits, back(limits) + sequenceLength(faiIndex, i));
        __int64 lengthSum = back(limits);

        if (options.verbosity >= 3)
        {
            for (unsigned i = 0; i < length(limits); ++i)
                std::cerr << "limit\t" << i << "\t" << limits[i] << "\n";
            std::cerr << "length sum\t" << lengthSum << "\n";
        }

        for (unsigned i = 0; i < length(variationSizeRecords); ++i)
        {
            __int64 x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<__int64> >(0, lengthSum - 1));
            if (options.verbosity >= 3)
                std::cerr << "  x == " << x << "\n";
            for (unsigned j = 0; j + 1 < length(limits); ++j)
                if (x >= limits[j] && x < limits[j + 1])
                {
                    if (options.verbosity >= 3)
                        std::cerr << "==> distributing " << i << " to " << j << "\n";
                    appendValue(variationToContig, j);
                    break;
                }
        }
    }

    void simulateContig(Variants & variants, unsigned rId, int haploCount)
    {
        seqan::CharString seq;
        if (readSequence(seq, faiIndex, rId) != 0)
        {
            std::cerr << "ERROR: Could not read sequence " << rId << " from FASTA file!\n";
            return;
        }

        if (!empty(options.inputSVSizeFile))
            _simulateFromSizes(variants, rId, haploCount, seq);
        else
            _simulateFromRates(variants, rId, haploCount, seq);
    }

    // Simulate the variants given variation types and kinds.
    void _simulateFromSizes(Variants & variants, unsigned rId, int haploCount, seqan::CharString const & seq)
    {
        // Picking SVs in a non-overlapping manner without any biases is too complicated to implement for a simulator.
        // Instead, we simulate the position uniformly at random and rerun picking the positions if we overlap with one
        // of the existing SVs within one base pair.  We impose a maximal number of retries of 1000 and stop for this
        // chromosome with a warning.  This leads to a quadratic running time but should be OK since we only need to
        // read hundreds of SV records from TSV at most.
        
        for (unsigned i = 0; i < length(variationSizeRecords); ++i)
        {
            VariationSizeRecord const & record = variationSizeRecords[i];
            int const MAX_TRIES = 1000;
            int tries = 0;
            for (; tries < MAX_TRIES; ++tries)
            {
                int pos = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, length(seq) - 1));

                switch (record.kind)
                {
                    case VariationSizeRecord::INDEL:
                        if (!simulateSVIndel(variants, haploCount, rId, pos, record.size))
                            continue;
                        break;
                    case VariationSizeRecord::INVERSION:
                        if (!simulateInversion(variants, haploCount, rId, pos, record.size))
                            continue;
                        break;
                    case VariationSizeRecord::TRANSLOCATION:
                        if (!simulateTranslocation(variants, haploCount, rId, pos, record.size))
                            continue;
                        break;
                    case VariationSizeRecord::DUPLICATION:
                        if (!simulateDuplication(variants, haploCount, rId, pos, record.size))
                            continue;
                        break;
                    default:
                        SEQAN_FAIL("Invalid SV type from TSV file!");
                }

                // Check whether the variant fits.
                bool variantFits = true;
                for (unsigned j = 0; j + 1 < length(variants.svRecords); ++j)
                    if (back(variants.svRecords).overlapsWith(variants.svRecords[j]))
                    {
                        variantFits = false;
                        break;
                    }
                if (!variantFits)
                    eraseBack(variants.svRecords);
                else
                    break;
            }
            if (tries == MAX_TRIES)
            {
                std::cerr << "WARNING: Could not place variant " << i << " for contig " << rId << " giving up for this contig.\n";
                return;
            }
        }

        // Sort simulated variants.
        std::sort(begin(variants.svRecords, seqan::Standard()), end(variants.svRecords, seqan::Standard()));
    }

    // Simulate the variants given per-position error rates.
    void _simulateFromRates(Variants & variants, unsigned rId, int haploCount, seqan::CharString const & seq)
    {
        // For each base, compute the whether to simulate a SNP and/or small indel.
        for (unsigned pos = 0; pos < length(seq); ++pos)
        {
            seqan::Pdf<seqan::Uniform<double> > pdf(0, 1);

            // Pick variant type if any.
            bool isIndel = (pickRandomNumber(rng, pdf) < options.svIndelRate);
            bool isInversion = (pickRandomNumber(rng, pdf) < options.svInversionRate);
            bool isTranslocation = (pickRandomNumber(rng, pdf) < options.svTranslocationRate);
            bool isDuplication = (pickRandomNumber(rng, pdf) < options.svDuplicationRate);
            while (isIndel + isInversion + isTranslocation + isDuplication > 1)
            {
                isIndel = (pickRandomNumber(rng, pdf) < options.svIndelRate);
                isInversion = (pickRandomNumber(rng, pdf) < options.svInversionRate);
                isTranslocation = (pickRandomNumber(rng, pdf) < options.svTranslocationRate);
                isDuplication = (pickRandomNumber(rng, pdf) < options.svDuplicationRate);
            }
            if (!isIndel && !isInversion && !isTranslocation && !isDuplication)
                continue;  // no variant picked

            int size = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(options.minSVSize,
                                                                              options.maxSVSize));
            if (isIndel)
            {
                if (!simulateSVIndel(variants, haploCount, rId, pos, size))
                    continue;
                if (back(variants.svRecords).size < 0)
                    pos += -back(variants.svRecords).size + 1;
            }
            else if (isInversion)
            {
                if (!simulateInversion(variants, haploCount, rId, pos, size))
                    continue;
                pos += back(variants.svRecords).size + 1;
            }
            else if (isTranslocation)
            {
                if (!simulateTranslocation(variants, haploCount, rId, pos, size))
                    continue;
                pos = back(variants.svRecords).targetPos + 1;
            }
            else if (isDuplication)
            {
                if (!simulateDuplication(variants, haploCount, rId, pos, size))
                    continue;
                pos = back(variants.svRecords).targetPos + 1;
            }

            SEQAN_ASSERT_LT(back(variants.svRecords).pos, (int)length(seq));
            if (back(variants.svRecords).targetPos != -1)
                SEQAN_ASSERT_LT(back(variants.svRecords).targetPos, (int)length(seq));

            if (options.verbosity >= 3)
                std::cerr << back(variants.svRecords) << "\n";
        }
    }

    bool simulateSVIndel(Variants & variants, int haploCount, int rId, unsigned pos, int size)
    {
        // Indels are simulated for one haplotype only.
        if (options.verbosity >= 2)
            std::cerr << "Simulating SV INDEL size = " << size << '\n';
        int hId = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, haploCount - 1));
        seqan::CharString indelSeq;
        reserve(indelSeq, options.maxSVSize);
        bool deletion = (size < 0);
        if (deletion && (pos + size) > sequenceLength(faiIndex, rId))
            return false;  // not enough space at the end
        seqan::Pdf<seqan::Uniform<int> > pdf(0, 3);
        for (int i = 0; i < size; ++i)  // not executed in case of deleted sequence
            appendValue(indelSeq, seqan::Dna5(pickRandomNumber(rng, pdf)));
        appendValue(variants.svRecords, StructuralVariantRecord(
                StructuralVariantRecord::INDEL, hId, rId, pos, size));
        back(variants.svRecords).seq = indelSeq;
        return true;
    }

    bool simulateInversion(Variants & variants, int haploCount, int rId, unsigned pos, int size)
    {
        if (pos + size >= sequenceLength(faiIndex, rId))
            return false;
        int hId = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, haploCount - 1));
        appendValue(variants.svRecords, StructuralVariantRecord(
                StructuralVariantRecord::INVERSION, hId, rId, pos, size));
        return true;
    }

    bool simulateTranslocation(Variants & variants, int haploCount, int rId, unsigned pos, int size)
    {
        int hId = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, haploCount - 1));
        int tPos = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(pos + size + options.minSVSize,
                                                                          pos + size + options.maxSVSize));
        if (tPos >= (int)sequenceLength(faiIndex, rId))
            return false;
        appendValue(variants.svRecords, StructuralVariantRecord(
                StructuralVariantRecord::TRANSLOCATION, hId, rId, pos, size, rId, tPos));
        return true;
    }

    bool simulateDuplication(Variants & variants, int haploCount, int rId, unsigned pos, int size)
    {
        if (!simulateTranslocation(variants, haploCount, rId, pos, size))
            return false;
        back(variants.svRecords).kind = StructuralVariantRecord::DUPLICATION;
        return true;
    }
};

// --------------------------------------------------------------------------
// Class SmallVariantSimulator
// --------------------------------------------------------------------------

// Simulation of small variants.

class SmallVariantSimulator
{
public:
    // Random number generator.
    TRng & rng;

    // FAI Index for loading sequence contig-wise.
    seqan::FaiIndex const & faiIndex;

    // The variator options.
    MasonVariatorOptions options;

    SmallVariantSimulator(TRng & rng, seqan::FaiIndex const & faiIndex, MasonVariatorOptions const & options) :
            rng(rng), faiIndex(faiIndex), options(options)
    {}

    // Perform simulation for one contig.
    //
    // variants should already contain structural variations for the current contig.  No small variants will be simulate
    // within 1 base to each side of the SV breakends.
    //
    // The variants in variants.svRecords must be non-overlapping.
    void simulateContig(Variants & variants, unsigned rId, int haploCount)
    {
        seqan::CharString seq;
        if (readSequence(seq, faiIndex, rId) != 0)
        {
            std::cerr << "ERROR: Could not read sequence " << rId << " from FASTA file!\n";
            return;
        }

        // Index into variants.svRecords.
        unsigned svIdx = 0;
        StructuralVariantRecord svRecord;
        if (!empty(variants.svRecords))
            svRecord = variants.svRecords[0];

        // For each base, compute the whether to simulate a SNP and/or small indel.
        for (unsigned pos = 0; pos < length(seq); ++pos)
        {
            // Seek next possible SV record that could have pos close to its breakends.
            bool skip = false;  // marker in case we switch SV records
            while ((int)pos > svRecord.endPosition())
            {
                // Skip if near breakend.
                skip = svRecord.nearBreakend(pos);

                if (options.verbosity >= 3)
                    std::cerr << " FROM " << svRecord;
                
                svIdx += 1;
                if (svIdx < length(variants.svRecords))
                    svRecord = variants.svRecords[svIdx];
                else
                    svRecord.pos = -1;  // mark as sentinel

                if (options.verbosity >= 3)
                    std::cerr << " TO " << svRecord << "\n";
            }
            // Skip if pos is near the SV breakend.
            if (skip || svRecord.nearBreakend(pos))
            {
                if (options.verbosity >= 3)
                    std::cerr << "pos " << pos << " is near " << svRecord << " (or previous)\n";
                continue;
            }
            
            seqan::Pdf<seqan::Uniform<double> > pdf(0, 1);

            // Perform experiment for SNP and small indel.
            bool isSnp = (pickRandomNumber(rng, pdf) < options.snpRate);
            double isIndel = (pickRandomNumber(rng, pdf) < options.smallIndelRate);
            int const MAX_TRIES = 1000;
            int tryNo = 0;
            for (; isSnp && isIndel && (tryNo < MAX_TRIES); ++tryNo)
            {
                isSnp = (pickRandomNumber(rng, pdf) < options.snpRate);
                isIndel = (pickRandomNumber(rng, pdf) < options.smallIndelRate);
                if (pos == 0)
                    isIndel = false;  // No indel at beginning, complex VCF case.
            }
            if (tryNo == MAX_TRIES)  // picked SNP and indel for MAX_TRIES time, pick none
                isSnp = (isIndel == false);

            // Simulate either SNP or indel.  In the case of a deletion, advance position such that there
            // is no variation in the deleted sequence.
            if (isSnp)
            {
                if (options.verbosity >= 3)
                    std::cerr << "Simulating SNP at (" << rId << ", " << pos << ")\n";
                simulateSnp(variants, seq, haploCount, rId, pos);
            }
            else if (isIndel)
            {
                if (!simulateSmallIndel(variants, haploCount, rId, pos))
                    continue;
                if (back(variants.smallIndels).size < 0)
                    pos += -back(variants.smallIndels).size + 1;
            }
        }
    }

    // Return true if SNP could be simulated.
    bool simulateSnp(Variants & variants, seqan::CharString & seq, int haploCount, int rId, unsigned pos)
    {
        // We simulate an alternative base for each haplotype.

        seqan::Dna5 from = seq[pos];
        for (int hId = 0; hId < haploCount; ++hId)
        {
            int toInt = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, 2));
            if (ordValue(from) <= toInt)
                toInt += 1;
            // std::cerr << hId << "\t" << rId << "\t" << pos << "\t" << from << "\t" << seqan::Dna5(toInt) << "\n";
            SEQAN_ASSERT_NEQ((int)ordValue(from), toInt);
            seqan::Dna5 to(toInt);
            appendValue(variants.snps, SnpRecord(hId, rId, pos, to));
        }

        return true;
    }

    // Return true if SNP could be simulated.
    bool simulateSmallIndel(Variants & variants, int haploCount, int rId, unsigned pos)
    {
        // Indels are simulated for one haplotype only.
        int hId = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, haploCount - 1));
        seqan::CharString indelSeq;
        reserve(indelSeq, options.maxSmallIndelSize);
        int indelSize = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(options.minSmallIndelSize,
                                                                               options.maxSmallIndelSize));
        bool deletion = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, 1));
        if (deletion && (pos + indelSize) > sequenceLength(faiIndex, rId))
            return false;  // not enough space at the end
        indelSize = deletion ? -indelSize : indelSize;
        seqan::Pdf<seqan::Uniform<int> > pdf(0, 3);
        for (int i = 0; i < indelSize; ++i)  // not executed in case of deleted sequence
            appendValue(indelSeq, seqan::Dna5(pickRandomNumber(rng, pdf)));
        appendValue(variants.smallIndels, SmallIndelRecord(hId, rId, pos, indelSize, indelSeq));
        return true;
    }
};

// --------------------------------------------------------------------------
// Class MasonVariatorApp
// --------------------------------------------------------------------------

class MasonVariatorApp
{
public:
    TRng & rng;

    MasonVariatorOptions options;

    seqan::VcfStream vcfStream;
    seqan::SequenceStream outSeqStream;

    seqan::SequenceStream outMethLevelStream;

    // FAI Index for loading sequence contig-wise.
    seqan::FaiIndex const & faiIndex;

    // Variation size record.
    seqan::String<VariationSizeRecord> variationSizeRecords;

    // File to write breakpoints to.
    std::fstream breakpointsOut;

    MasonVariatorApp(TRng & rng, seqan::FaiIndex const & faiIndex,
                     MasonVariatorOptions const & options) :
            rng(rng), options(options), faiIndex(faiIndex)
    {}

    int run()
    {
        // Open output breakpoints TSV file.
        if (!empty(options.outputBreakpointFile))
        {            
            breakpointsOut.open(toCString(options.outputBreakpointFile), std::ios::binary | std::ios::out);
            if (!breakpointsOut.good())
            {
                std::cerr << "ERROR: Could not open " << options.outputBreakpointFile << " for writing.\n";
                return 1;
            }
            breakpointsOut << "#ref\tpos\n";
        }

        // Open VCF stream to write to.
        open(vcfStream, toCString(options.vcfOutFile), seqan::VcfStream::WRITE);
        if (!isGood(vcfStream))
        {
            std::cerr << "Could not open " << options.vcfOutFile << " for writing.\n";
            return 1;
        }
        // Create header.
        appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord("fileformat", "VCFv4.1"));
        appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord("source", "mason_variator"));
        appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord("reference", options.fastaInFile));
        appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
                "INFO", "<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"));
        appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
                "INFO", "<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"));
        appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
                "INFO", "<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"));
        appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
                "INFO", "<ID=TARGETPOS,Number=1,Type=String,Description=\"Target position for duplications.\">"));
        appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
                "ALT", "<ID=INV,Description=\"Inversion\">"));
        appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
                "ALT", "<ID=DUP,Description=\"Duplication\">"));
        // We don't need DEL and INS here since we report exact one with the sequence.
        // appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
        //         "ALT", "<ID=DEL,Description=\"Deletion\">"));
        // appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
        //         "ALT", "<ID=INS,Description=\"Insertion of novel sequence\">"));

        // Copy over sequence names.
        for (unsigned i = 0; i < numSeqs(faiIndex); ++i)
        {
            seqan::CharString contigStr = "<ID=";
            append(contigStr, sequenceName(faiIndex, i));
            append(contigStr, ",length=");
            std::stringstream ss;
            ss << sequenceLength(faiIndex, i);
            append(contigStr, ss.str());
            append(contigStr, ">");
            appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord("contig", contigStr));
            appendName(*vcfStream._context.sequenceNames,
                       sequenceName(faiIndex, i),
                       vcfStream._context.sequenceNamesCache);
        }
        // Copy over sample names.
        appendName(*vcfStream._context.sampleNames, "simulated", vcfStream._context.sampleNamesCache);

        // Open output FASTA file if necessary.
        if (!empty(options.fastaOutFile))
        {
            open(outSeqStream, toCString(options.fastaOutFile), seqan::SequenceStream::WRITE,
                 seqan::SequenceStream::FASTA);
            if (!isGood(outSeqStream))
            {
                std::cerr << "ERROR: Could not open " << options.fastaOutFile << " for writing!\n";
                return 1;
            }
        }

        // Open methylation level output file if necessary.
        if (options.simulateMethylationLevels && !empty(options.methFastaOutFile))
        {
            open(outMethLevelStream, toCString(options.methFastaOutFile), seqan::SequenceStream::WRITE,
                 seqan::SequenceStream::FASTA);
            if (!isGood(outMethLevelStream))
            {
                std::cerr << "ERROR: Could not open " << options.methFastaOutFile << " for writing!\n";
                return 1;
            }
        }

        // Read in variant size TSV if path is given.
        if (_readVariationSizes() != 0)
            return 1;

        // Actually perform the variant simulation.
        if (options.verbosity >= 1)
            std::cerr << "\nSimulation...\n";
        StructuralVariantSimulator svSim(rng, faiIndex, variationSizeRecords, options);
        SmallVariantSimulator smallSim(rng, faiIndex, options);
        for (int rId = 0; rId < (int)numSeqs(faiIndex); ++rId)  // ref seqs
        {
            // Simulate methylation levels if configured to do so and write out for reference.  We always pass the
            // levels down into _simulateContigs() but it can be empty and ignored if methylation levels are not of
            // interest.
            MethylationLevels methLevels;
            if (options.simulateMethylationLevels)
                _simulateMethLevels(methLevels, rId);
            // Simulate contigs.
            _simulateContig(svSim, smallSim, methLevels, options, rId);
        }
        if (options.verbosity >= 1)
            std::cerr << "OK.\n\n";

        return 0;
    }

    // Simulate methylation levels.
    int _simulateMethLevels(MethylationLevels & levels, int rId)
    {
        MethylationLevelSimulator methSim(rng, options);
        seqan::Dna5String contig;
        if (!readSequence(contig, faiIndex, rId))
            methSim.run(levels, contig);

        return 0;
    }

    // Write out methylation levels to output file.
    //
    // levels -- levels
    // hId -- haplotype id, -1 for original
    // rId -- reference id
    int _writeMethylationLevels(MethylationLevels const & levels, int hId, int rId)
    {
        std::stringstream idTop;
        idTop << sequenceName(faiIndex, rId);
        if (hId != -1)
            idTop << options.haplotypeSep << hId;
        idTop << options.haplotypeSep << "TOP";
        if (writeRecord(outMethLevelStream, idTop.str(), levels.forward) != 0)
        {
            std::cerr << "ERROR: Problem writing to " << options.methFastaOutFile << "\n";
            return 1;
        }

        std::stringstream idBottom;
        idBottom << sequenceName(faiIndex, rId);
        if (hId != -1)
            idBottom << options.haplotypeSep << hId;
        idBottom << options.haplotypeSep << "BOT";
        if (writeRecord(outMethLevelStream, idBottom.str(), levels.reverse) != 0)
        {
            std::cerr << "ERROR: Problem writing to " << options.methFastaOutFile << "\n";
            return 1;
        }

        return 0;
    }

    // Read variation size TSV file if any is given.
    int _readVariationSizes()
    {
        if (empty(options.inputSVSizeFile))
            return 0;  // Nothing to do

        if (options.verbosity >= 1)
            std::cerr << "Variation Sizes " << options.inputSVSizeFile << " ...";

        std::fstream inF(toCString(options.inputSVSizeFile), std::ios::in | std::ios::binary);
        if (!inF.good())
        {
            std::cerr << "ERROR: Could not open " << options.inputSVSizeFile << "\n";
            return 1;
        }

        seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(inF);
        VariationSizeRecord record;
        while (!atEnd(reader))
        {
            if (value(reader) == '#')
            {
                skipLine(reader);
                continue;  // Skip comment.
            }

            if (readRecord(record, reader, VariationSizeTsv()) != 0)
            {
                std::cerr << "ERROR: Problem reading from " << options.inputSVSizeFile << "\n";
                return 1;
            }
            appendValue(variationSizeRecords, record);
        }

        if (options.verbosity >= 1)
            std::cerr << " OK\n";
        
        return 0;
    }

    // Perform simulation of one contig.
    //
    // If options.simulateMethylationLevels then methLevels will be used, otherwise it can be (and is) empty.
    //
    // svSim -- simulator for structural variants
    // smallSim -- simulator for small variants
    // methLevels -- methylation level information for reference
    // options -- configuration
    // rId -- ID of reference sequence that we are using now
    int _simulateContig(StructuralVariantSimulator & svSim,
                        SmallVariantSimulator & smallSim,
                        MethylationLevels const & methLevels,
                        MasonVariatorOptions const & options,
                        int rId)
    {
        if (options.verbosity >= 1)
            std::cerr << "  " << sequenceName(faiIndex, rId) << "\n";

        // Simulate variants.
        Variants variants;
        svSim.simulateContig(variants, rId, options.numHaplotypes);
        std::sort(begin(variants.svRecords, seqan::Standard()),
                  end(variants.svRecords, seqan::Standard()));
        smallSim.simulateContig(variants, rId, options.numHaplotypes);
        if (options.verbosity >= 1)
            std::cerr << "  snps:                " << length(variants.snps) << "\n"
                      << "  small indels:        " << length(variants.smallIndels) << "\n"
                      << "  structural variants: " << length(variants.svRecords) << "\n";

        // Load contig seq.
        // TODO(holtgrew): Pass from outside, could reuse sequence from methylation levels.
        seqan::Dna5String contig;
        if (readSequence(contig, faiIndex, rId) != 0)
        {
            std::cerr << "Could not read contig seq " << rId << "\n";
            return 1;
        }

        // Write out variants for contig to VCF file.
        if (_writeVcf(contig, variants, rId) != 0)
            return 1;

        // Apply variants to contigs and write out.
        if (!empty(options.fastaOutFile))
        {
            // Write out methylation levels for reference.
            if (options.simulateMethylationLevels && !empty(options.methFastaOutFile))
                if (_writeMethylationLevels(methLevels, -1, rId) != 0)
                    return 1;
            // Apply variations to contigs and write out.
            for (int hId = 0; hId < options.numHaplotypes; ++hId)
            {
                if (_writeContigs(contig, variants, methLevels, rId, hId) != 0)
                    return 1;
            }
        }

        return 0;
    }

    int _writeContigs(seqan::Dna5String const & contig, Variants const & variants, MethylationLevels const & levels, int rId, int hId)
    {
        // Build sequence id.
        seqan::CharString id = sequenceName(faiIndex, rId);
        append(id, options.haplotypeSep);
        char buffer[20];
        snprintf(buffer, 19, "%d", hId);
        append(id, buffer);

        // Apply small variants.  We get a sequence with the small variants and a journal of the difference to contig.
        seqan::Dna5String seqSmallVariants;
        TJournalEntries journal;
        MethylationLevels levelsSmallVariants;
        if (_materializeSmallVariants(seqSmallVariants, journal, levelsSmallVariants, contig, variants, levels, rId, hId) != 0)
            return 1;

        // Apply structural variants.
        seqan::Dna5String seqLargeVariants;
        MethylationLevels levelsLargeVariants;
        if (_materializeLargeVariants(seqLargeVariants, levelsLargeVariants, journal, seqSmallVariants, variants, levelsSmallVariants, rId, hId, id) != 0)
            return 1;

        // Write out methylation levels if necessary.
        if (options.simulateMethylationLevels && !empty(options.methFastaOutFile))
            if (_writeMethylationLevels(levelsLargeVariants, hId, rId) != 0)
                return 1;

        return writeRecord(outSeqStream, id, seqLargeVariants);
    }

    // Apply large structural variants from variants into seq.  The input is given as the sequence including small
    // variants and a journal for translating coordinates in contig into coordinates of the underlying sequence which is
    // used as the coordinate system in variants.
    int _materializeLargeVariants(seqan::Dna5String & seq,
                                  MethylationLevels & levelsLargeVariants,
                                  TJournalEntries const & journal,
                                  seqan::Dna5String const & contig,
                                  Variants const & variants,
                                  MethylationLevels const & levels,
                                  int /*rId*/,
                                  int hId,
                                  seqan::CharString const & ref)
    {
        if (options.verbosity >= 2)
            std::cerr << "\nMATERIALIZING LARGE VARIANTS FOR HAPLOTYPE " << hId << "\n\n";
        // Clear output methylation levels.
        levelsLargeVariants.clear();
        // Store variation points.  We reuse the _fixVariationLevels() function from small indel/snp simulation and thus
        // have to store a bool that is always set to false.
        seqan::String<std::pair<int, bool> > varPoints;
        
        // Track last position from contig appended to seq so far.
        int lastPos = 0;
        if (options.verbosity >= 3)
            std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";

        // Number of bytes written out so far/current position in variant.
        int currentPos = 0;

        for (unsigned i = 0; i < length(variants.svRecords); ++i)
        {
            if (variants.svRecords[i].haplotype != hId)  // Ignore all but the current contig.
                continue;
            // We obtain a copy of the current SV record since we translate its positions below.
            StructuralVariantRecord svRecord = variants.svRecords[i];

            // Translate positions and lengths of SV record.
            if (options.verbosity >= 2)
                std::cerr << "  Translating SvRecord\n  " << svRecord << '\n';
            svRecord.pos = hostToVirtualPosition(journal, svRecord.pos);
            SEQAN_ASSERT_LT(svRecord.pos, (int)length(contig));
            // We do not need to adjust the sizes for insertions.
            if (svRecord.kind != StructuralVariantRecord::INDEL || svRecord.size < 0)
                svRecord.size = hostToVirtualPosition(journal, svRecord.pos + svRecord.size) -
                        hostToVirtualPosition(journal, svRecord.pos);
            if (svRecord.targetPos != -1)
                svRecord.targetPos = hostToVirtualPosition(journal, svRecord.targetPos);
            if (options.verbosity >= 2)
                std::cerr << "  => " << svRecord << '\n';

            // Copy from contig to seq with SVs.
            if (options.verbosity >= 3)
                std::cerr << "lastPos == " << lastPos << "\n";
            append(seq, infix(contig, lastPos, svRecord.pos));  // interim chars
            if (options.simulateMethylationLevels)
            {
                append(levelsLargeVariants.forward, infix(levels.forward, lastPos, svRecord.pos));
                append(levelsLargeVariants.reverse, infix(levels.reverse, lastPos, svRecord.pos));
                appendValue(varPoints, std::make_pair(length(seq), false));
            }
            currentPos = length(seq);
            if (options.verbosity >= 3)
                std::cerr << "append(seq, infix(contig, " << lastPos << ", " << svRecord.pos << ") " << __LINE__ << " (interim)\n";
            switch (svRecord.kind)
            {
                case StructuralVariantRecord::INDEL:
                    {
                        if (svRecord.size > 0)  // insertion
                        {
                            SEQAN_ASSERT_EQ((int)length(svRecord.seq), svRecord.size);

                            // Simulate methylation levels for insertion.
                            MethylationLevels lvls;
                            if (options.simulateMethylationLevels)
                            {
                                MethylationLevelSimulator methSim(rng, options);
                                methSim.run(lvls, svRecord.seq);
                            }

                            // Append novel sequence and methylation levels.
                            append(seq, svRecord.seq);
                            if (options.simulateMethylationLevels)
                            {
                                append(levelsLargeVariants.forward, lvls.forward);
                                append(levelsLargeVariants.reverse, lvls.reverse);
                                appendValue(varPoints, std::make_pair(length(seq), false));  // variation point after insertion
                            }
                            if (options.verbosity >= 3)
                                std::cerr << "append(seq, svRecord.seq (length == " << length(svRecord.seq) << ") " << __LINE__ << " (insertion)\n";
                            lastPos = svRecord.pos;
                            SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                            if (!empty(options.outputBreakpointFile))  // write out breakpoints
                                breakpointsOut << ref << "\t" << (currentPos + 1) << "\n"
                                               << ref << "\t" << (length(seq) + 1) << "\n";
                            currentPos = length(seq);
                        }
                        else  // deletion
                        {
                            lastPos = svRecord.pos - svRecord.size;
                            SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                            if (!empty(options.outputBreakpointFile))  // write out breakpoint
                                breakpointsOut << ref << "\t" << currentPos << "\n";
                        }
                    }
                    break;
                case StructuralVariantRecord::INVERSION:
                    {
                        unsigned oldLen = length(seq);
                        append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                        if (options.simulateMethylationLevels)
                        {
                            appendValue(varPoints, std::make_pair(length(seq), false));  // variation point at deletion
                            append(levelsLargeVariants.forward, infix(levels.reverse, svRecord.pos, svRecord.pos + svRecord.size));
                            reverse(infix(levelsLargeVariants.forward, oldLen, length(levelsLargeVariants.forward)));
                            append(levelsLargeVariants.reverse, infix(levels.forward, svRecord.pos, svRecord.pos + svRecord.size));
                            reverse(infix(levelsLargeVariants.reverse, oldLen, length(levelsLargeVariants.reverse)));
                        }

                        if (options.verbosity >= 3)
                            std::cerr << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << " (inversion)\n";
                        reverseComplement(infix(seq, oldLen, length(seq)));
                        lastPos = svRecord.pos + svRecord.size;
                        SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                        if (!empty(options.outputBreakpointFile))  // write out breakpoint
                            breakpointsOut << ref << "\t" << (currentPos + 1) << "\n"
                                           << ref << "\t" << (length(seq) + 1) << "\n";
                        currentPos = length(seq);
                    }
                    break;
                case StructuralVariantRecord::TRANSLOCATION:
                    {
                        SEQAN_ASSERT_GEQ(svRecord.targetPos, svRecord.pos + svRecord.size);
                        append(seq, infix(contig, svRecord.pos + svRecord.size, svRecord.targetPos));
                        if (options.simulateMethylationLevels)
                        {
                            appendValue(varPoints, std::make_pair(length(seq), false));
                            append(levelsLargeVariants.forward, infix(levels.forward, svRecord.pos + svRecord.size, svRecord.targetPos));
                            append(levelsLargeVariants.reverse, infix(levels.reverse, svRecord.pos + svRecord.size, svRecord.targetPos));
                        }
                        append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                        if (options.simulateMethylationLevels)
                        {
                            appendValue(varPoints, std::make_pair(length(seq), false));
                            append(levelsLargeVariants.forward, infix(levels.forward, svRecord.pos, svRecord.pos + svRecord.size));
                            append(levelsLargeVariants.reverse, infix(levels.reverse, svRecord.pos, svRecord.pos + svRecord.size));
                        }
                        if (options.verbosity >= 3)
                            std::cerr << "append(seq, infix(contig, " << svRecord.pos + svRecord.size << ", " << svRecord.targetPos << ") " << __LINE__ << " (translocation)\n"
                                      << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << "\n";
                        lastPos = svRecord.targetPos;
                        SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                        if (!empty(options.outputBreakpointFile))  // write out breakpoint
                            breakpointsOut << ref << "\t" << (currentPos + 1) << "\n"
                                           << ref << "\t" << (currentPos + svRecord.targetPos - svRecord.pos - svRecord.size + 1) << "\n"
                                           << ref << "\t" << (length(seq) + 1) << "\n";
                        currentPos = length(seq);
                    }
                    break;
                case StructuralVariantRecord::DUPLICATION:
                    {
                        append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                        SEQAN_ASSERT_GEQ(svRecord.targetPos, svRecord.pos + svRecord.size);
                        if (options.simulateMethylationLevels)
                        {
                            appendValue(varPoints, std::make_pair(length(seq), false));
                            append(levelsLargeVariants.forward, infix(levels.forward, svRecord.pos, svRecord.pos + svRecord.size));
                            append(levelsLargeVariants.reverse, infix(levels.reverse, svRecord.pos, svRecord.pos + svRecord.size));
                        }
                        append(seq, infix(contig, svRecord.pos + svRecord.size, svRecord.targetPos));
                        if (options.simulateMethylationLevels)
                        {
                            appendValue(varPoints, std::make_pair(length(seq), false));
                            append(levelsLargeVariants.forward, infix(levels.forward, svRecord.pos + svRecord.size, svRecord.targetPos));
                            append(levelsLargeVariants.reverse, infix(levels.reverse, svRecord.pos + svRecord.size, svRecord.targetPos));
                        }
                        append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                        if (options.simulateMethylationLevels)
                        {
                            appendValue(varPoints, std::make_pair(length(seq), false));
                            append(levelsLargeVariants.forward, infix(levels.forward, svRecord.pos, svRecord.pos + svRecord.size));
                            append(levelsLargeVariants.reverse, infix(levels.reverse, svRecord.pos, svRecord.pos + svRecord.size));
                        }
                        if (options.verbosity >= 3)
                            std::cerr << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << " (duplication)\n"
                                      << "append(seq, infix(contig, " << svRecord.pos + svRecord.size << ", " << svRecord.targetPos << ") " << __LINE__ << "\n"
                                      << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << "\n";
                        lastPos = svRecord.targetPos;
                        SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                        if (!empty(options.outputBreakpointFile))  // write out breakpoint
                            breakpointsOut << ref << "\t" << (currentPos + 1) << "\n"
                                           << ref << "\t" << (currentPos + svRecord.pos + svRecord.size - svRecord.pos + 1) << "\n"
                                           << ref << "\t" << (currentPos + svRecord.pos + svRecord.size - svRecord.pos +
                                                              svRecord.targetPos - (svRecord.pos + svRecord.size) + 1) << "\n"
                                           << ref << "\t" << (length(seq) + 1) << "\n";
                        currentPos = length(seq);
                    }
                    break;
                default:
                    return 1;
            }
        }
        if (options.verbosity >= 3)
            std::cerr << "append(seq, infix(contig, " << lastPos << ", " << length(contig) << ") "
                      << __LINE__ << " (last interim)\n";
        append(seq, infix(contig, lastPos, length(contig)));
        if (options.simulateMethylationLevels)
        {
            append(levelsLargeVariants.forward, infix(levels.forward, lastPos, length(contig)));
            append(levelsLargeVariants.reverse, infix(levels.reverse, lastPos, length(contig)));

            SEQAN_ASSERT_EQ(length(seq), length(levelsLargeVariants.forward));
            SEQAN_ASSERT_EQ(length(seq), length(levelsLargeVariants.reverse));

            _fixVariationLevels(levelsLargeVariants, seq, varPoints);
        }
        
        return 0;
    }

    // Apply small indels and SNPs from variants into seq using contig.
    int _materializeSmallVariants(seqan::Dna5String & seq,
                                  TJournalEntries & journal,
                                  MethylationLevels & levelsSmallVariants,
                                  seqan::Dna5String const & contig,
                                  Variants const & variants,
                                  MethylationLevels const & levels,
                                  int /*rId*/,
                                  int hId)
    {
        // Clear journal and output methylation levels.
        reinit(journal, length(contig));
        levelsSmallVariants.clear();
        // Store variation points with a flag whether it is a SNP (true) or a breakpoint (false).
        seqan::String<std::pair<int, bool> > varPoints;

        // Fors this, we have to iterate in parallel over SNP and small indel records.
        //
        // Current index in snp/small indel array.
        unsigned snpsIdx = 0;
        unsigned smallIndelIdx = 0;
        // Current SNP record, default to sentinel.
        SnpRecord snpRecord;
        snpRecord.rId = seqan::maxValue<int>();
        if (snpsIdx < length(variants.snps))
            snpRecord = variants.snps[snpsIdx++];
        // Current small indel record, default to sentinel.
        SmallIndelRecord smallIndelRecord;
        smallIndelRecord.rId = seqan::maxValue<int>();
        if (smallIndelIdx < length(variants.smallIndels))
            smallIndelRecord = variants.smallIndels[smallIndelIdx++];
        // Track last position from contig appended to seq so far.
        int lastPos = 0;
        if (options.verbosity >= 3)
            std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";

        // TODO(holtgrew): Extract contig building into their own functions.
        if (options.verbosity >= 2)
            std::cerr << "building output\n";
        while (snpRecord.rId != seqan::maxValue<int>() || smallIndelRecord.rId != seqan::maxValue<int>())
        {
            // TODO(holtgrew): Extract SNP and small indel handling in functions.
            if (snpRecord.getPos() < smallIndelRecord.getPos())  // process SNP records
            {
                if (snpRecord.haplotype == hId)  // Ignore all but the current contig.
                {
                    if (options.verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << lastPos << ", " << snpRecord.pos << ") " << __LINE__ << "\n";
                    // Append interim sequence and methylation levels.
                    append(seq, infix(contig, lastPos, snpRecord.pos));
                    if (options.simulateMethylationLevels)
                    {
                        append(levelsSmallVariants.forward, infix(levels.forward, lastPos, snpRecord.pos + 1));
                        append(levelsSmallVariants.reverse, infix(levels.reverse, lastPos, snpRecord.pos + 1));
                        appendValue(varPoints, std::make_pair(length(seq), true));      // variation points before/after SNP
                        appendValue(varPoints, std::make_pair(length(seq) + 1, true));
                    }

                    SEQAN_ASSERT_GEQ(snpRecord.pos, lastPos);
                    if (options.verbosity >= 3)
                        std::cerr << "appendValue(seq, " << snpRecord.to << "')\n";
                    appendValue(seq, snpRecord.to);
                    lastPos = snpRecord.pos + 1;
                    if (options.verbosity >= 3)
                        std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";
                }

                if (snpsIdx >= length(variants.snps))
                    snpRecord.rId = seqan::maxValue<int>();
                else
                    snpRecord = variants.snps[snpsIdx++];
            }
            else
            {
                if (smallIndelRecord.haplotype == hId)  // Ignore all but the current contig.
                {
                    if (smallIndelRecord.size > 0)
                    {
                        if (options.verbosity >= 3)
                            std::cerr << "append(seq, infix(contig, " << lastPos << ", " << smallIndelRecord.pos << ") " << __LINE__ << "\n";

                        // Simulate methylation levels for insertion.
                        MethylationLevels lvls;
                        if (options.simulateMethylationLevels)
                        {
                            MethylationLevelSimulator methSim(rng, options);
                            methSim.run(lvls, smallIndelRecord.seq);
                        }

                        // Append interim sequence and methylation levels.
                        append(seq, infix(contig, lastPos, smallIndelRecord.pos));
                        if (options.simulateMethylationLevels)
                        {
                            append(levelsSmallVariants.forward, infix(levels.forward, lastPos, smallIndelRecord.pos));
                            append(levelsSmallVariants.reverse, infix(levels.reverse, lastPos, smallIndelRecord.pos));
                            appendValue(varPoints, std::make_pair(length(seq), false));  // variation point before insertion
                        }

                        SEQAN_ASSERT_GEQ(smallIndelRecord.pos, lastPos);
                        if (options.verbosity >= 3)
                            std::cerr << "append(seq, \"" << smallIndelRecord.seq << "\") " << __LINE__ << "\n";
                        // Append novel sequence and methylation levels.
                        append(seq, smallIndelRecord.seq);
                        if (options.simulateMethylationLevels)
                        {
                            append(levelsSmallVariants.forward, lvls.forward);
                            append(levelsSmallVariants.reverse, lvls.reverse);
                            appendValue(varPoints, std::make_pair(length(seq), false));  // variation point after insertion
                        }
                        lastPos = smallIndelRecord.pos;
                        recordInsertion(journal, hostToVirtualPosition(journal, smallIndelRecord.pos),
                                        0, smallIndelRecord.size);
                        if (options.verbosity >= 3)
                            std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";
                    }
                    else  // deletion
                    {
                        if (options.verbosity >= 3)
                            std::cerr << "append(seq, infix(contig, " << lastPos << ", " << smallIndelRecord.pos << ") " << __LINE__ << "\n";
                        // Append interim sequence and methylation levels.
                        append(seq, infix(contig, lastPos, smallIndelRecord.pos));  // interim chars
                        if (options.simulateMethylationLevels)
                        {
                            appendValue(varPoints, std::make_pair(length(seq), false));  // variation point at deletion
                            append(levelsSmallVariants.forward, infix(levels.forward, lastPos, smallIndelRecord.pos));
                            append(levelsSmallVariants.reverse, infix(levels.reverse, lastPos, smallIndelRecord.pos));
                        }

                        lastPos = smallIndelRecord.pos - smallIndelRecord.size;
                        recordErase(journal,
                                    hostToVirtualPosition(journal, smallIndelRecord.pos),
                                    hostToVirtualPosition(journal, smallIndelRecord.pos - smallIndelRecord.size));
                        if (options.verbosity >= 3)
                            std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";
                    }
                }

                if (smallIndelIdx >= length(variants.smallIndels))
                    smallIndelRecord.rId = seqan::maxValue<int>();
                else
                    smallIndelRecord = variants.smallIndels[smallIndelIdx++];
            }
        }
        // Insert remaining characters.
        if (options.verbosity >= 3)
            std::cerr << "append(seq, infix(contig, " << lastPos << ", " << length(contig) << ")\n";
        append(seq, infix(contig, lastPos, length(contig)));

        if (options.simulateMethylationLevels)
        {
            append(levelsSmallVariants.forward, infix(levels.forward, lastPos, length(contig)));
            append(levelsSmallVariants.reverse, infix(levels.reverse, lastPos, length(contig)));

            SEQAN_ASSERT_EQ(length(seq), length(levelsSmallVariants.forward));
            SEQAN_ASSERT_EQ(length(seq), length(levelsSmallVariants.reverse));

            _fixVariationLevels(levelsSmallVariants, seq, varPoints);
        }

        return 0;
    }

    // Recompute variation levels around the variation points varPoints (position, isSnp) with the given contig and
    // levels.
    void _fixVariationLevels(MethylationLevels & levels,
                             seqan::Dna5String const & contig,
                             seqan::String<std::pair<int, bool> > const & varPoints)
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

    // Write out variants for the given contig to the VCF file.
    int _writeVcf(seqan::Dna5String const & contig, Variants const & variants, int /*rId*/)
    {
        // Current index in snp/small indel and SV array.
        unsigned snpsIdx = 0;
        unsigned smallIndelIdx = 0;
        unsigned svIdx = 0;
        // Current SNP record, default to sentinel.
        SnpRecord snpRecord;
        snpRecord.rId = seqan::maxValue<int>();
        if (snpsIdx < length(variants.snps))
            snpRecord = variants.snps[snpsIdx++];
        // Current small indel record, default to sentinel.
        SmallIndelRecord smallIndelRecord;
        smallIndelRecord.rId = seqan::maxValue<int>();
        if (smallIndelIdx < length(variants.smallIndels))
            smallIndelRecord = variants.smallIndels[smallIndelIdx++];
        // Current SV record, default to sentinel.
        StructuralVariantRecord svRecord;
        svRecord.rId = seqan::maxValue<int>();
        if (svIdx < length(variants.svRecords))
            svRecord = variants.svRecords[svIdx++];
        while (snpRecord.rId != seqan::maxValue<int>() ||
               smallIndelRecord.rId != seqan::maxValue<int>() ||
               svRecord.rId != seqan::maxValue<int>())
        {
            if (snpRecord.rId != seqan::maxValue<int>() && smallIndelRecord.rId != seqan::maxValue<int>())
                SEQAN_ASSERT(snpRecord.getPos() != smallIndelRecord.getPos());  // are generated indendently
            if (snpRecord.rId != seqan::maxValue<int>() && svRecord.rId != seqan::maxValue<int>())
                SEQAN_ASSERT_MSG(snpRecord.getPos() != svRecord.getPos(),
                                 "Should be generated non-overlapping (snp pos = %d, sv pos = %d).",
                                 snpRecord.pos, svRecord.pos);
            if (smallIndelRecord.rId != seqan::maxValue<int>() && svRecord.rId != seqan::maxValue<int>())
                SEQAN_ASSERT(smallIndelRecord.getPos() != svRecord.getPos());  // are generated indendently
            SEQAN_ASSERT_NEQ(snpRecord.pos, 0);   // Not simulated, VCF complexer.
            SEQAN_ASSERT_NEQ(svRecord.pos, 0);   // Not simulated, VCF complexer.
            SEQAN_ASSERT_NEQ(smallIndelRecord.pos, 0);   // Not simulated, VCF complexer.

            // Structure of if/else statement is (1) SNP, (2) small indel, (3) structural variants.
            if (snpRecord.getPos() < smallIndelRecord.getPos() &&
                snpRecord.getPos() < svRecord.getPos())  // process SNP records
            {
                if (_writeVcfSnp(contig, variants, snpRecord, snpsIdx) != 0)
                    return 1;
            }
            else if (smallIndelRecord.getPos() < svRecord.getPos())// process small indel records
            {
                if (_writeVcfSmallIndel(contig, variants, smallIndelRecord, smallIndelIdx) != 0)
                    return 1;
            }
            else  // sv record
            {
                if (options.verbosity >= 2)
                    std::cerr << "  SV record to file.\n";
                SEQAN_ASSERT_GT_MSG(svRecord.pos, 0,
                                    "SV cannot be at genome begin yet and should not be generated as such either.");

                if (svRecord.kind == StructuralVariantRecord::INDEL)
                {
                    if (_writeVcfIndel(contig, svRecord) != 0)
                        return 1;
                }
                else if (svRecord.kind == StructuralVariantRecord::INVERSION)
                {
                    if (_writeVcfInversion(contig, svRecord) != 0)
                        return 1;
                }
                else if (svRecord.kind == StructuralVariantRecord::TRANSLOCATION)
                {
                    if (_writeVcfTranslocation(contig, svRecord) != 0)
                        return 1;
                }
                else if (svRecord.kind == StructuralVariantRecord::DUPLICATION)
                {
                    if (_writeVcfDuplication(contig, svRecord) != 0)
                        return 1;
                }

                if (svIdx >= length(variants.svRecords))
                    svRecord.rId = seqan::maxValue<int>();
                else
                    svRecord = variants.svRecords[svIdx++];
            }
        }

        return 0;
    }

    int _writeVcfSnp(seqan::Dna5String const & contig,
                     Variants const & variants,
                     SnpRecord & snpRecord,
                     unsigned & snpsIdx)
    {
        if (options.verbosity >= 2)
            std::cerr << "  snpRecord record.\n";
        int rId = snpRecord.rId;

        // Store information used below.
        // std::cerr << "from = " << contig[snpRecord.pos] << "\n";
        seqan::Dna5 from = contig[snpRecord.pos];
        std::pair<int, int> pos = snpRecord.getPos();

        // Get the value of each haplotype at the position.
        seqan::String<bool> inTos;
        resize(inTos, 4, false);
        seqan::Dna5String tos;
        resize(tos, options.numHaplotypes, from);
        do
        {
            SEQAN_ASSERT(snpRecord.to != from);
            tos[snpRecord.haplotype] = snpRecord.to;
            inTos[ordValue(seqan::Dna5(snpRecord.to))] = true;

            if (snpsIdx >= length(variants.snps))
                snpRecord.rId = seqan::maxValue<int>();
            else
                snpRecord = variants.snps[snpsIdx++];
        }
        while (snpRecord.rId != seqan::maxValue<int>() &&
               snpsIdx < length(variants.snps) &&
               snpRecord.getPos() == pos);

        // Create VCF vcfRecord.
        seqan::VcfRecord vcfRecord;
        vcfRecord.rID = rId;
        vcfRecord.beginPos = pos.second;
        // TODO(holtgrew): Generate an id?
        appendValue(vcfRecord.ref, from);
        for (unsigned i = 0; i < 4; ++i)
        {
            if (!inTos[i])
                continue;  // no ALT
            if (!empty(vcfRecord.alt))
                appendValue(vcfRecord.alt, ',');
            appendValue(vcfRecord.alt, seqan::Dna5(i));
        }
        vcfRecord.filter = "PASS";
        vcfRecord.info = ".";
        // Build genotype infos.
        appendValue(vcfRecord.genotypeInfos, "");
        for (int hId = 0; hId < options.numHaplotypes; ++hId)
        {
            if (!empty(vcfRecord.genotypeInfos[0]))
                appendValue(vcfRecord.genotypeInfos[0], '|');
            if (tos[hId] == vcfRecord.ref[0])
            {
                appendValue(vcfRecord.genotypeInfos[0], '0');
            }
            else
            {
                char buffer[20];
                for (unsigned i = 0; i < length(vcfRecord.alt); i += 2)
                    if (tos[hId] == vcfRecord.alt[i])
                    {
                        snprintf(buffer, 19, "%d", 1 + i / 2);
                        append(vcfRecord.genotypeInfos[0], buffer);
                    }
            }
        }

        // Write out VCF record.
        if (writeRecord(vcfStream, vcfRecord) != 0)
        {
            std::cerr << "ERROR: Problem writing to " << options.vcfOutFile << "\n";
            return 1;
        }

        return 0;
    }

    int _writeVcfSmallIndel(seqan::Dna5String const & contig,
                            Variants const & variants,
                            SmallIndelRecord & smallIndelRecord,
                            unsigned & smallIndelIdx)
    {
        // Collect small indel records at the same position.
        seqan::String<SmallIndelRecord> records;
        do
        {
            if (options.verbosity >= 3)
                std::cerr << "INDEL\t"
                          << smallIndelRecord.haplotype << "\t"
                          << smallIndelRecord.rId << "\t"
                          << smallIndelRecord.pos << "\t"
                          << smallIndelRecord.size << "\t"
                          << smallIndelRecord.seq << "\n";
            appendValue(records, smallIndelRecord);

            if (smallIndelIdx >= length(variants.smallIndels))
                smallIndelRecord.rId = seqan::maxValue<int>();
            else
                smallIndelRecord = variants.smallIndels[smallIndelIdx++];
        }
        while (smallIndelRecord.rId != seqan::maxValue<int>() &&
               smallIndelIdx < length(variants.smallIndels) &&
               smallIndelRecord.getPos() == variants.smallIndels[smallIndelIdx].getPos());
        SEQAN_ASSERT_NOT(empty(records));

        // Create VCF record.
        seqan::VcfRecord vcfRecord;
        vcfRecord.rID = front(records).rId;
        vcfRecord.beginPos = front(records).pos;
        // TODO(holtgrew): Generate an id?
        vcfRecord.filter = "PASS";
        vcfRecord.info = ".";
        // Build genotype infos.

        // Compute the number of bases in the REF column (1 in case of insertion and (k + 1) in the case of a
        // deletion of length k.
        int numRef = 0;
        for (unsigned i = 0; i < length(records); ++i)
        {
            SEQAN_ASSERT_NEQ(records[i].size, 0);
            if (records[i].size > 0)
                numRef = std::max(numRef, 1);  // assign 1 if 0
            else  // if (records[i].size < 0)
                numRef = std::max(numRef, 1 - records[i].size);
        }
        append(vcfRecord.ref, infix(contig, vcfRecord.beginPos - 1, vcfRecord.beginPos - 1 + numRef));

        // Compute ALT columns and a map to the ALT.
        seqan::String<int> toIds;
        resize(toIds, options.numHaplotypes, 0);
        for (unsigned i = 0; i < length(records); ++i)
        {
            if (i > 0)
                appendValue(vcfRecord.alt, ',');
            toIds[records[i].haplotype] = i + 1;
            if (records[i].size > 0)  // insertion
            {
                appendValue(vcfRecord.alt, vcfRecord.ref[0]);
                append(vcfRecord.alt, records[i].seq);
                append(vcfRecord.alt, suffix(vcfRecord.ref, 1));
            }
            else  // deletion
            {
                appendValue(vcfRecord.alt, vcfRecord.ref[0]);
                append(vcfRecord.alt, suffix(vcfRecord.ref, 1 - records[i].size));
            }
        }

        // Create genotype infos.
        appendValue(vcfRecord.genotypeInfos, "");
        for (int i = 0; i < options.numHaplotypes; ++i)
        {
            if (i > 0)
                appendValue(vcfRecord.genotypeInfos[0], '|');
            char buffer[20];
            snprintf(buffer, 19, "%d", toIds[i]);
            append(vcfRecord.genotypeInfos[0], buffer);
        }

        // Write out VCF record.
        if (writeRecord(vcfStream, vcfRecord) != 0)
        {
            std::cerr << "ERROR: Problem writing to " << options.vcfOutFile << "\n";
            return 1;
        }
        return 0;
    }

    int _writeVcfIndel(seqan::Dna5String const & contig,
                       StructuralVariantRecord const & svRecord)
    {
        // TODO(holtgrew): Large indels can be represented by <INS> and <DEL> and should be.
        if (options.verbosity >= 2)
            std::cerr << "indel\t" << svRecord << "\n";

        // Create VCF record.
        seqan::VcfRecord vcfRecord;
        vcfRecord.rID = svRecord.rId;
        vcfRecord.beginPos = svRecord.pos;
        // TODO(holtgrew): Generate an id?
        vcfRecord.filter = "PASS";
        std::stringstream ss;
        if (svRecord.size > 0)
            ss << "SVTYPE=INS";
        else
            ss << "SVTYPE=DEL";
        ss << ";SVLEN=" << svRecord.size;
        vcfRecord.info = ss.str();

        // Compute the number of bases in the REF column (1 in case of insertion and (k + 1) in the case of a
        // deletion of length k.
        int numRef;
        if (svRecord.size > 0)
            numRef = 1;
        else
            numRef = 1 - svRecord.size;
        append(vcfRecord.ref, infix(contig, vcfRecord.beginPos - 1, vcfRecord.beginPos - 1 + numRef));

        // Compute ALT columns and a map to the ALT.
        if (svRecord.size > 0)  // insertion
        {
            appendValue(vcfRecord.alt, vcfRecord.ref[0]);
            append(vcfRecord.alt, svRecord.seq);
        }
        else
        {
            appendValue(vcfRecord.alt, vcfRecord.ref[0]);
            append(vcfRecord.alt, suffix(vcfRecord.ref, 1 - svRecord.size));
        }

        // Create genotype infos.
        appendValue(vcfRecord.genotypeInfos, "");
        for (int i = 0; i < options.numHaplotypes; ++i)
        {
            if (i > 0)
                appendValue(vcfRecord.genotypeInfos[0], '|');
            if (svRecord.haplotype == i)
                appendValue(vcfRecord.genotypeInfos[0], '1');
            else
                appendValue(vcfRecord.genotypeInfos[0], '0');
        }

        // Write out VCF record.
        if (writeRecord(vcfStream, vcfRecord) != 0)
        {
            std::cerr << "ERROR: Problem writing to " << options.vcfOutFile << "\n";
            return 1;
        }
        return 0;
    }

    int _writeVcfTranslocation(seqan::Dna5String const & contig,
                               StructuralVariantRecord const & svRecord)
    {
        // In this function, we will create VCF records left and right of both cut positions and of the paste position.
        seqan::VcfRecord leftOfCutL, rightOfCutL, leftOfCutR, rightOfCutR, leftOfPaste, rightOfPaste;
        // CHROM ID
        leftOfCutL.rID = svRecord.rId;
        rightOfCutL.rID = svRecord.rId;
        leftOfCutR.rID = svRecord.rId;
        rightOfCutR.rID = svRecord.rId;
        leftOfPaste.rID = svRecord.rId;
        rightOfPaste.rID = svRecord.rId;
        // POS
        leftOfCutL.beginPos = svRecord.pos - 1;
        rightOfCutL.beginPos = svRecord.pos;
        leftOfCutR.beginPos = svRecord.pos - 1 + svRecord.size;
        rightOfCutR.beginPos = svRecord.pos + svRecord.size;
        leftOfPaste.beginPos = svRecord.targetPos - 1;
        rightOfPaste.beginPos = svRecord.targetPos;
        // TODO(holtgrew): INFO entry with type of breakend?
        // ID (none)
        // TODO(holtgrew): Generate an id?
        // REF
        appendValue(leftOfCutL.ref, contig[leftOfCutL.beginPos]);
        appendValue(rightOfCutL.ref, contig[rightOfCutL.beginPos]);
        appendValue(leftOfCutR.ref, contig[leftOfCutR.beginPos]);
        appendValue(rightOfCutR.ref, contig[rightOfCutR.beginPos]);
        appendValue(leftOfPaste.ref, contig[leftOfPaste.beginPos]);
        appendValue(rightOfPaste.ref, contig[rightOfPaste.beginPos]);
        // ALT
        std::stringstream ssLeftOfCutL, ssRightOfCutL, ssLeftOfCutR, ssRightOfCutR,
                ssLeftOfPaste, ssRightOfPaste;
        seqan::CharString refName = vcfStream.header.sequenceNames[svRecord.rId];
        ssLeftOfCutL << leftOfCutL.ref << "]" << refName << ":" << (rightOfCutR.beginPos + 1) << "]";
        leftOfCutL.alt = ssLeftOfCutL.str();
        ssRightOfCutL << "[" << refName << ":" << (leftOfPaste.beginPos + 1) << "[" << rightOfCutL.ref;
        rightOfCutL.alt = ssRightOfCutL.str();
        ssLeftOfCutR << leftOfCutR.ref << "]" << refName << ":" << (rightOfPaste.beginPos + 1) << "]";
        leftOfCutR.alt = ssLeftOfCutR.str();
        ssRightOfCutR << "[" << refName << ":" << (leftOfCutL.beginPos + 1) << "[" << rightOfCutR.ref;
        rightOfCutR.alt = ssRightOfCutR.str();
        ssLeftOfPaste << leftOfPaste.ref << "]" << refName << ":" << (leftOfCutR.beginPos + 1) << "]";
        leftOfPaste.alt = ssLeftOfPaste.str();
        ssRightOfPaste << "[" << refName << ":" << (rightOfCutL.beginPos + 1) << "[" << rightOfPaste.ref;
        rightOfPaste.alt = ssRightOfPaste.str();
        // FILTER
        leftOfCutL.filter = "PASS";
        rightOfCutL.filter = "PASS";
        leftOfCutR.filter = "PASS";
        rightOfCutR.filter = "PASS";
        leftOfPaste.filter = "PASS";
        rightOfPaste.filter = "PASS";
        // INFO
        leftOfCutL.info = "SVTYPE=BND";
        rightOfCutL.info = "SVTYPE=BND";
        leftOfCutR.info = "SVTYPE=BND";
        rightOfCutR.info = "SVTYPE=BND";
        leftOfPaste.info = "SVTYPE=BND";
        rightOfPaste.info = "SVTYPE=BND";

        // Create genotype infos.
        appendValue(leftOfCutL.genotypeInfos, "");
        appendValue(rightOfCutL.genotypeInfos, "");
        appendValue(leftOfCutR.genotypeInfos, "");
        appendValue(rightOfCutR.genotypeInfos, "");
        appendValue(leftOfPaste.genotypeInfos, "");
        appendValue(rightOfPaste.genotypeInfos, "");
        for (int i = 0; i < options.numHaplotypes; ++i)
        {
            if (i > 0)
            {
                appendValue(leftOfCutL.genotypeInfos[0], '|');
                appendValue(rightOfCutL.genotypeInfos[0], '|');
                appendValue(leftOfCutR.genotypeInfos[0], '|');
                appendValue(rightOfCutR.genotypeInfos[0], '|');
                appendValue(leftOfPaste.genotypeInfos[0], '|');
                appendValue(rightOfPaste.genotypeInfos[0], '|');
            }
            if (svRecord.haplotype == i)
            {
                appendValue(leftOfCutL.genotypeInfos[0], '1');
                appendValue(rightOfCutL.genotypeInfos[0], '1');
                appendValue(leftOfCutR.genotypeInfos[0], '1');
                appendValue(rightOfCutR.genotypeInfos[0], '1');
                appendValue(leftOfPaste.genotypeInfos[0], '1');
                appendValue(rightOfPaste.genotypeInfos[0], '1');
            }
            else
            {
                appendValue(leftOfCutL.genotypeInfos[0], '0');
                appendValue(rightOfCutL.genotypeInfos[0], '0');
                appendValue(leftOfCutR.genotypeInfos[0], '0');
                appendValue(rightOfCutR.genotypeInfos[0], '0');
                appendValue(leftOfPaste.genotypeInfos[0], '0');
                appendValue(rightOfPaste.genotypeInfos[0], '0');
            }
        }
                
        // Write out VCF records.
        if (writeRecord(vcfStream, leftOfCutL) != 0 || writeRecord(vcfStream, rightOfCutL) != 0 ||
            writeRecord(vcfStream, leftOfCutR) != 0 || writeRecord(vcfStream, rightOfCutR) != 0 ||
            writeRecord(vcfStream, leftOfPaste) != 0 || writeRecord(vcfStream, rightOfPaste) != 0)
        {
            std::cerr << "ERROR: Problem writing to " << options.vcfOutFile << "\n";
            return 1;
        }

        return 0;
    }

    int _writeVcfInversion(seqan::Dna5String const & contig,
                           StructuralVariantRecord const & svRecord)
    {
        if (options.verbosity >= 2)
            std::cerr << "inversion\t" << svRecord << "\n";
        seqan::VcfRecord vcfRecord;

        vcfRecord.rID = svRecord.rId;
        vcfRecord.beginPos = svRecord.pos;
        // TODO(holtgrew): Generate an id?
        appendValue(vcfRecord.ref, contig[vcfRecord.beginPos]);
        vcfRecord.alt = "<INV>";
        vcfRecord.filter = "PASS";
        std::stringstream ss;
        ss << "SVTYPE=INV;END=" << (svRecord.pos + svRecord.size) << ";SVLEN=" << svRecord.size;
        vcfRecord.info = ss.str();

        // Create genotype infos.
        appendValue(vcfRecord.genotypeInfos, "");
        for (int i = 0; i < options.numHaplotypes; ++i)
        {
            if (i > 0)
                appendValue(vcfRecord.genotypeInfos[0], '|');
            if (svRecord.haplotype == i)
                appendValue(vcfRecord.genotypeInfos[0], '1');
            else
                appendValue(vcfRecord.genotypeInfos[0], '0');
        }

        // Write out VCF record.
        if (writeRecord(vcfStream, vcfRecord) != 0)
        {
            std::cerr << "ERROR: Problem writing to " << options.vcfOutFile << "\n";
            return 1;
        }
        return 0;
    }

    int _writeVcfDuplication(seqan::Dna5String const & contig,
                             StructuralVariantRecord const & svRecord)
    {
        // TODO(holtgrew): Large indels can be represented by <INS> and <DEL> and should be.
        if (options.verbosity >= 2)
            std::cerr << "duplication\t" << svRecord << "\n";

        // Create VCF record.
        seqan::VcfRecord vcfRecord;
        vcfRecord.rID = svRecord.rId;
        vcfRecord.beginPos = svRecord.pos;
        // TODO(holtgrew): Generate an id?
        vcfRecord.filter = "PASS";
        std::stringstream ss;
        ss << "SVTYPE=DUP;SVLEN=" << svRecord.size << ";END=" << svRecord.pos + svRecord.size
           << ";TARGETPOS=" << vcfStream.header.sequenceNames[svRecord.targetRId] << ":" << svRecord.targetPos;
        vcfRecord.info = ss.str();
        appendValue(vcfRecord.ref, contig[vcfRecord.beginPos]);
        vcfRecord.alt = "<DUP>";

        // Create genotype infos.
        appendValue(vcfRecord.genotypeInfos, "");
        for (int i = 0; i < options.numHaplotypes; ++i)
        {
            if (i > 0)
                appendValue(vcfRecord.genotypeInfos[0], '|');
            if (svRecord.haplotype == i)
                appendValue(vcfRecord.genotypeInfos[0], '1');
            else
                appendValue(vcfRecord.genotypeInfos[0], '0');
        }

        // Write out VCF record.
        if (writeRecord(vcfStream, vcfRecord) != 0)
        {
            std::cerr << "ERROR: Problem writing to " << options.vcfOutFile << "\n";
            return 1;
        }
        return 0;
    }

};

// ==========================================================================
// Metafunctions
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(MasonVariatorOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("mason_variator");
    // Set short description, version, and date.
    setShortDescription(parser, "Variation Simulation");
    setVersion(parser, "2.0alpha1");
    setDate(parser, "March 2013");
    setCategory(parser, "Simulators");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-if\\fP \\fIIN.fa\\fP \\fB-ov\\fP \\fIOUT.vcf\\fP [\\fB-of\\fP \\fIOUT.fa\\fP]");
    // addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-if\\fP \\fIIN.fa\\fP \\fB-iv\\fP \\fIIN.vcf\\fP \\fB-of\\fP \\fIOUT.fa\\fP");
    addDescription(parser,
                   "Either simulate variation and write out the result to VCF and optionally FASTA files.");
    // addDescription(parser,
    //                "Either simulate variation and write out the result to VCF and FASTA files "
    //                "or apply the variations from a VCF file and write the results to a FASTA file.");

    // ----------------------------------------------------------------------
    // General Options
    // ----------------------------------------------------------------------

    addSection(parser, "General Options");

    // We require one argument.
    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    addOption(parser, seqan::ArgParseOption("s", "seed", "The seed to use for the random number generator.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "seed", "0");

    // ----------------------------------------------------------------------
    // Input / Output Options
    // ----------------------------------------------------------------------

    addSection(parser, "Input / Output");
    
    // addOption(parser, seqan::ArgParseOption("iv", "in-vcf", "VCF file to load variations from.",
    //                                         seqan::ArgParseOption::INPUTFILE, "VCF"));
    // setValidValues(parser, "in-vcf", "vcf");

    addOption(parser, seqan::ArgParseOption("if", "in-fasta", "FASTA file with reference.",
                                            seqan::ArgParseOption::INPUTFILE, "FASTA"));
    setValidValues(parser, "in-fasta", "fasta fa");
    setRequired(parser, "in-fasta");

    addOption(parser, seqan::ArgParseOption("it", "in-variant-tsv",
                                            "TSV file with variants to simulate.  See Section on the Variant TSV File below.",
                                            seqan::ArgParseOption::INPUTFILE, "VCF"));
    setValidValues(parser, "in-variant-tsv", "tsv txt");

    addOption(parser, seqan::ArgParseOption("ov", "out-vcf", "VCF file to write simulated variations to.",
                                            seqan::ArgParseOption::INPUTFILE, "VCF"));
    setRequired(parser, "out-vcf");
    setValidValues(parser, "out-vcf", "vcf");

    addOption(parser, seqan::ArgParseOption("of", "out-fasta", "FASTA file to write simulated haplotypes to.",
                                            seqan::ArgParseOption::INPUTFILE, "FASTA"));
    setValidValues(parser, "out-fasta", "fasta fa");

    addOption(parser, seqan::ArgParseOption("", "out-breakpoints", "TSV file to write breakpoints in variants to.",
                                            seqan::ArgParseOption::OUTPUTFILE, "TSV"));
    setValidValues(parser, "out-breakpoints", "tsv txt");

    addOption(parser, seqan::ArgParseOption("", "haplotype-name-sep", "Haplotype name separator in output FASTA.",
                                            seqan::ArgParseOption::STRING, "SEP"));
    setDefaultValue(parser, "haplotype-name-sep", "/");

    // ----------------------------------------------------------------------
    // Haplotype / Allele Configuration
    // ----------------------------------------------------------------------

    addSection(parser, "Haplotype / Allele Configuration");

    addOption(parser, seqan::ArgParseOption("n", "num-haplotypes", "The number of haplotypes to simulate.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setMinValue(parser, "num-haplotypes", "1");
    setDefaultValue(parser, "num-haplotypes", "1");

    addOption(parser, seqan::ArgParseOption("", "haplotype-sep",
                                            "The separator between the chromosome and the haplotype name "
                                            "in the output FASTA file.",
                                            seqan::ArgParseOption::STRING, "SEP"));
    setDefaultValue(parser, "haplotype-sep", "/");

    // ----------------------------------------------------------------------
    // Variation Simulation Options
    // ----------------------------------------------------------------------

    addSection(parser, "Variation Simulation");

    addOption(parser, seqan::ArgParseOption("", "snp-rate", "Per-base SNP rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "snp-rate", "0.0");
    setMaxValue(parser, "snp-rate", "1.0");
    setDefaultValue(parser, "snp-rate", "0.0001");

    addOption(parser, seqan::ArgParseOption("", "small-indel-rate", "Small indel rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "small-indel-rate", "0.0");
    setMaxValue(parser, "small-indel-rate", "1.0");
    setDefaultValue(parser, "small-indel-rate", "0.000001");

    addOption(parser, seqan::ArgParseOption("", "min-small-indel-size", "Minimal small indel size to simulate.",
                                            seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "min-small-indel-size", "0");
    setDefaultValue(parser, "min-small-indel-size", "1");

    addOption(parser, seqan::ArgParseOption("", "max-small-indel-size", "Maximal small indel size to simulate.",
                                            seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "max-small-indel-size", "0");
    setDefaultValue(parser, "max-small-indel-size", "6");

    addOption(parser, seqan::ArgParseOption("", "sv-indel-rate", "Per-base SNP rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "sv-indel-rate", "0.0");
    setMaxValue(parser, "sv-indel-rate", "1.0");
    setDefaultValue(parser, "sv-indel-rate", "0.0000001");

    addOption(parser, seqan::ArgParseOption("", "sv-inversion-rate", "Per-base SNP rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "sv-inversion-rate", "0.0");
    setMaxValue(parser, "sv-inversion-rate", "1.0");
    setDefaultValue(parser, "sv-inversion-rate", "0.0000001");

    addOption(parser, seqan::ArgParseOption("", "sv-translocation-rate", "Per-base SNP rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "sv-translocation-rate", "0.0");
    setMaxValue(parser, "sv-translocation-rate", "1.0");
    setDefaultValue(parser, "sv-translocation-rate", "0.0000001");

    addOption(parser, seqan::ArgParseOption("", "sv-duplication-rate", "Per-base SNP rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "sv-duplication-rate", "0.0");
    setMaxValue(parser, "sv-duplication-rate", "1.0");
    setDefaultValue(parser, "sv-duplication-rate", "0.0000001");

    addOption(parser, seqan::ArgParseOption("", "min-sv-size", "Minimal SV size to simulate.",
                                            seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "min-sv-size", "0");
    setDefaultValue(parser, "min-sv-size", "50");

    addOption(parser, seqan::ArgParseOption("", "max-sv-size", "Maximal SV size to simulate.",
                                            seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "max-sv-size", "0");
    setDefaultValue(parser, "max-sv-size", "1000");

    // ----------------------------------------------------------------------
    // Methylation Simulation Options
    // ----------------------------------------------------------------------

    addSection(parser, "Methylation Rates");

    addOption(parser, seqan::ArgParseOption("", "methylation-levels", "Enable simulation of methylation levels."));

    addOption(parser, seqan::ArgParseOption("", "meth-fasta-out", "Path to write methylation levels to as FASTA.  "
                                            "Only written if \\fB-of\\fP/\\fB--out-fasta\\fP is given.",
                                            seqan::ArgParseOption::OUTPUTFILE, "FILE"));
    setValidValues(parser, "meth-fasta-out", "fa fasta");

    addOption(parser, seqan::ArgParseOption("", "meth-c-mu", "Median of beta distribution for methylation "
                                            "level of cytosine.", seqan::ArgParseOption::DOUBLE, "MU"));
    setMinValue(parser, "meth-c-mu", "0");
    setMaxValue(parser, "meth-c-mu", "1");
    setDefaultValue(parser, "meth-c-mu", "0.0001");    

    addOption(parser, seqan::ArgParseOption("", "meth-c-sigma", "Standard deviation of beta distribution for "
                                            "methylation level of cytosine.", seqan::ArgParseOption::DOUBLE, "SIGMA"));
    setMinValue(parser, "meth-c-sigma", "0");
    setMaxValue(parser, "meth-c-sigma", "1");
    setDefaultValue(parser, "meth-c-sigma", "0.00001");

    addOption(parser, seqan::ArgParseOption("", "meth-cg-mu", "Median of beta distribution for methylation "
                                            "level of CpG loci.", seqan::ArgParseOption::DOUBLE, "MU"));
    setMinValue(parser, "meth-cg-mu", "0");
    setMaxValue(parser, "meth-cg-mu", "1");
    setDefaultValue(parser, "meth-cg-mu", "0.6");    

    addOption(parser, seqan::ArgParseOption("", "meth-cg-sigma", "Standard deviation of beta distribution for "
                                            "methylation level of CpG loci.", seqan::ArgParseOption::DOUBLE,
                                            "SIGMA"));
    setMinValue(parser, "meth-cg-sigma", "0");
    setMaxValue(parser, "meth-cg-sigma", "1");
    setDefaultValue(parser, "meth-cg-sigma", "0.03");

    addOption(parser, seqan::ArgParseOption("", "meth-chg-mu", "Median of beta distribution for methylation "
                                            "level of CHG loci.", seqan::ArgParseOption::DOUBLE, "MU"));
    setMinValue(parser, "meth-chg-mu", "0");
    setMaxValue(parser, "meth-chg-mu", "1");
    setDefaultValue(parser, "meth-chg-mu", "0.08");    

    addOption(parser, seqan::ArgParseOption("", "meth-chg-sigma", "Standard deviation of beta distribution for "
                                            "methylation level of CHG loci.", seqan::ArgParseOption::DOUBLE,
                                            "SIGMA"));
    setMinValue(parser, "meth-chg-sigma", "0");
    setMaxValue(parser, "meth-chg-sigma", "1");
    setDefaultValue(parser, "meth-chg-sigma", "0.008");

    addOption(parser, seqan::ArgParseOption("", "meth-chh-mu", "Median of beta distribution for methylation "
                                            "level of CHH loci.", seqan::ArgParseOption::DOUBLE, "MU"));
    setMinValue(parser, "meth-chh-mu", "0");
    setMaxValue(parser, "meth-chh-mu", "1");
    setDefaultValue(parser, "meth-chh-mu", "0.05");    

    addOption(parser, seqan::ArgParseOption("", "meth-chh-sigma", "Standard deviation of beta distribution for "
                                            "methylation level of CHH loci.", seqan::ArgParseOption::DOUBLE,
                                            "SIGMA"));
    setMinValue(parser, "meth-chh-sigma", "0");
    setMaxValue(parser, "meth-chh-sigma", "1");
    setDefaultValue(parser, "meth-chh-sigma", "0.005");

    // ----------------------------------------------------------------------
    // Simulation Details Section
    // ----------------------------------------------------------------------

    addTextSection(parser, "Simulation Details");

    addText(parser,
            "SNPs and small indels are simulated such that at each position, a random experiment is "
            "performed whether to simulate either variation.  In case both variations are to be simulated, "
            "the experiment is repeated.");

    addText(parser, "The indel and SV sizes are picked uniformly at random from the argument size intervals.");

    addText(parser,
            "The simulation of haplotypes works as follows.  For small indels, the indel is placed into "
            "one of the haplotypes that are to be simulated.  The exact haplotype is picked uniformly at "
            "random.  For SNPs, we simulate a random base for each haplotype.  For at least one haplotype, "
            "the base has to be different from the reference or the experiment is repeated.");

    // ----------------------------------------------------------------------
    // Examples Section
    // ----------------------------------------------------------------------

    // TODO(holtgrew): Write me!

    // ----------------------------------------------------------------------
    // Variation TSV File
    // ----------------------------------------------------------------------

    addTextSection(parser, "Variation TSV File");
    addText(parser,
            "Instead of simulating the SVs from per-base rates, the user can specify a TSV (tab separated values) "
            "file to load the variations from with \\fB--in-variant-tsv\\fP/\\fB-it\\fP.  The first two columns of "
            "this TSV file are interpreted as the type of the variation and the size.");
    addText(parser,
            "Indels smaller than 50 bp are considered small indels whereas larger indels are considered structural "
            "variants in the VCF file.");
    addListItem(parser, "INS", "An insertion.");
    addListItem(parser, "DEL", "A deletion.");
    addListItem(parser, "INV", "An inversion.");
    // addListItem(parser, "TRA", "An inter-chromosomal translocation.");  // TODO(holtgrew): Add support.
    addListItem(parser, "CTR", "An intra-chromosomal translocation.");
    addListItem(parser, "DUP", "A duplication");

    // ----------------------------------------------------------------------
    // Methylation Level Simulation
    // ----------------------------------------------------------------------

    addTextSection(parser, "Methylation Level Simulation");
    addText(parser,
            "Simulation of cytosine methylation levels is done using a beta distribution.  There is one distribution "
            "each for cytosines in the context CpG, CHG, and CHH and one distribution for all other cytonsines.  You "
            "can give the parameters mu and sigma of the beta distributions.  The methylation level is determined once "
            "for each base of the reference (0 for all non-cytosines) and stored in a string of levels.  This string "
            "is then modified as small and structural variations are simualted.");
    addText(parser,
            "The simulated methylation levels can then be written out to a FASTA file.  This file will contain two "
            "entries for the original and each haplotype;  the levels for the forward and the reverse strand.  The "
            "sequence will be ASCII characters 0, starting at '!' encoding the level in 1.25% steps.  The character "
            "'>' is ignored and encodes no level.");
    addText(parser,
            "Methylation level simulation increases the memory usage of the program by one byte for each character "
            "in the largest contig.");
    // TODO(holtgrew): Simulate different levels for each haplotype?
            
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    options.verbosity = 1;
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValue(options.seed, parser, "seed");

    // getOptionValue(options.vcfInFile, parser, "in-vcf");
    getOptionValue(options.fastaInFile, parser, "in-fasta");
    getOptionValue(options.vcfOutFile, parser, "out-vcf");
    getOptionValue(options.fastaOutFile, parser, "out-fasta");
    getOptionValue(options.outputBreakpointFile, parser, "out-breakpoints");
    getOptionValue(options.inputSVSizeFile, parser, "in-variant-tsv");

    getOptionValue(options.numHaplotypes, parser, "num-haplotypes");
    getOptionValue(options.haplotypeSep, parser, "haplotype-sep");
    getOptionValue(options.snpRate, parser, "snp-rate");
    getOptionValue(options.smallIndelRate, parser, "small-indel-rate");
    getOptionValue(options.minSmallIndelSize, parser, "min-small-indel-size");
    getOptionValue(options.maxSmallIndelSize, parser, "max-small-indel-size");
    getOptionValue(options.svIndelRate, parser, "sv-indel-rate");
    getOptionValue(options.svInversionRate, parser, "sv-inversion-rate");
    getOptionValue(options.svTranslocationRate, parser, "sv-translocation-rate");
    getOptionValue(options.svDuplicationRate, parser, "sv-duplication-rate");
    getOptionValue(options.minSVSize, parser, "min-sv-size");
    getOptionValue(options.maxSVSize, parser, "max-sv-size");

    getOptionValue(options.simulateMethylationLevels, parser, "methylation-levels");
    getOptionValue(options.methFastaOutFile, parser, "meth-fasta-out");
    getOptionValue(options.methMuC, parser, "meth-c-mu");
    getOptionValue(options.methSigmaC, parser, "meth-c-sigma");
    getOptionValue(options.methMuCG, parser, "meth-cg-mu");
    getOptionValue(options.methSigmaCG, parser, "meth-cg-sigma");
    getOptionValue(options.methMuCHG, parser, "meth-chg-mu");
    getOptionValue(options.methSigmaCHG, parser, "meth-chg-sigma");
    getOptionValue(options.methMuCHH, parser, "meth-chh-mu");
    getOptionValue(options.methSigmaCHH, parser, "meth-chh-sigma");

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    MasonVariatorOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Initialize random number generator.
    TRng rng(options.seed);

    std::cerr << "MASON VARIATOR\n"
              << "==============\n\n";

    print(std::cerr, options);

    std::cerr << "\n__PREPARATION_________________________________________________________________\n"
              << "\n";

    std::cerr << "Loading Reference Index " << options.fastaInFile << " ...";
    seqan::FaiIndex faiIndex;
    if (read(faiIndex, toCString(options.fastaInFile)) != 0)
    {
        std::cerr << " FAILED (not fatal, we can just build it)\n";
        std::cerr << "Building Index        " << options.fastaInFile << ".fai ...";
        if (build(faiIndex, toCString(options.fastaInFile)) != 0)
        {
            std::cerr << "Could not build FAI index.\n";
            return 1;
        }
        std::cerr << " OK\n";
        seqan::CharString faiPath = options.fastaInFile;
        append(faiPath, ".fai");
        std::cerr << "Reference Index       " << faiPath << " ...";
        if (write(faiIndex, toCString(faiPath)) != 0)
        {
            std::cerr << "Could not write FAI index we just built.\n";
            return 1;
        }
        std::cerr << " OK (" << length(faiIndex.indexEntryStore) << " seqs)\n";
    }
    else
    {
        std::cerr << " OK (" << length(faiIndex.indexEntryStore) << " seqs)\n";
    }

    std::cerr << "\n__SIMULATION__________________________________________________________________\n"
              << "\n";

    MasonVariatorApp app(rng, faiIndex, options);
    app.run();

    std::cerr << "\nDONE.\n";
    
    return 0;
}
