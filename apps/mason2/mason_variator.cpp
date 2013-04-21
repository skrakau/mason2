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
// TODO(holtgrew): What about shortcuts for SV duplications with target?
// TODO(holtgrew): Simulate different SNPs/small variations for duplications, input for repeat separation.
// TODO(holtgrew): Does SV rate give the per position rate of the event or the number of bases related to an event?

#include <seqan/arg_parse.h>
#include <seqan/random.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/vcf.h>
#include <seqan/sequence_journaled.h>

#include "variation_size_tsv.h"

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

    MasonVariatorOptions() :
            verbosity(1), seed(0),
            snpRate(0), smallIndelRate(0), minSmallIndelSize(0), maxSmallIndelSize(0), svIndelRate(0),
            svInversionRate(0), svTranslocationRate(0), svDuplicationRate(0), minSVSize(0), maxSVSize(0)
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
        << "\n";
}

// --------------------------------------------------------------------------
// Class SnpRecord
// --------------------------------------------------------------------------

// Represents a SNP in one haplotype.
//
// All coordinates are given as source coordinates.

struct SnpRecord
{
    // Reference id and position on the reference.
    int rId;
    int pos;

    // The haplotype that this variation belongs to.
    int haplotype;

    // The target nucleotide.
    seqan::Dna5 to;

    SnpRecord() : rId(-1), pos(-1), haplotype(-1), to('\0')
    {}

    SnpRecord(int haplotype, int rId, int pos, char to) :
            rId(rId), pos(pos), haplotype(haplotype), to(to)
    {}

    bool operator<(SnpRecord const & other) const
    {
        if (rId < other.rId || (rId == other.rId && pos < other.pos) ||
            (rId == other.rId && pos == other.pos && haplotype < other.haplotype))
            return true;
        return false;
    }

    std::pair<int, int> getPos() const
    {
        return std::make_pair(rId, pos);
    }
};

inline std::ostream & operator<<(std::ostream & out, SnpRecord const & record)
{
    out << "SnpRecord(" << record.haplotype << ", " << record.rId << ", " << record.pos
        << ", " << record.to << ")";
    return out;
}

// --------------------------------------------------------------------------
// Class SmallIndelRecord
// --------------------------------------------------------------------------

// Represents a small indel.
//
// All coordinates are given as source coordinates.

struct SmallIndelRecord
{
    // Reference id and position on the reference.
    int rId;
    int pos;

    // The haplotype that this variation belongs to.
    int haplotype;

    // The size of the indel, negative numbers for deletions, positive numbers for insertions.
    int size;

    // The inserted sequence if any.
    seqan::CharString seq;

    SmallIndelRecord() : rId(-1), pos(-1), haplotype(-1), size(0)
    {}

    SmallIndelRecord(int haplotype, int rId, int pos, int size, seqan::CharString const & seq) :
            rId(rId), pos(pos), haplotype(haplotype), size(size), seq(seq)
    {}

    bool operator<(SmallIndelRecord const & other) const
    {
        return getPos() < other.getPos();
    }

    std::pair<int, int> getPos() const
    {
        return std::make_pair(rId, pos);
    }
};

inline std::ostream & operator<<(std::ostream & out, SmallIndelRecord const & record)
{
    out << "SnpRecord(" << record.haplotype << ", " << record.rId << ", " << record.pos
        << ", " << record.size << ", " << record.seq << ")";
    return out;
}

// --------------------------------------------------------------------------
// Class StructuralVariantRecord
// --------------------------------------------------------------------------

// Store structural variant information.
//
// All coordinates are given as source coordinates.

struct StructuralVariantRecord
{
    enum Kind
    {
        INVALID,
        INDEL,
        INVERSION,
        TRANSLOCATION,
        DUPLICATION
    };

    // The kind of the SV.
    Kind kind;

    // Reference id and position on the reference.
    int rId;
    int pos;

    // The haplotype that this variation belongs to.
    int haplotype;

    // In case of indel, negative numbers give deletions, positionve numbers give insertions.  In case of inversions,
    // the length of the inverted sequence.  In case of translocation the length of the translocated sequence, and in
    // the case the length of the duplicated sequence.
    int size;

    // Target reference id and position on this reference.  Currently, we only have SVs on the same chromosome.  Used in
    // case of translocation and duplication.
    int targetRId;
    int targetPos;

    // The inserted sequence if any.
    seqan::CharString seq;

    StructuralVariantRecord() :
            kind(INVALID), rId(-1), pos(-1), haplotype(-1), size(0), targetRId(-1), targetPos(-1)
    {}

    StructuralVariantRecord(Kind kind, int haplotype, int rId, int pos, int size,
                            int targetRId = -1, int targetPos = -1) :
            kind(kind), rId(rId), pos(pos), haplotype(haplotype), size(size), targetRId(targetRId),
            targetPos(targetPos)
    {}

    bool operator<(StructuralVariantRecord const & other) const
    {
        return getPos() < other.getPos();
    }

    std::pair<int, int> getPos() const
    {
        return std::make_pair(rId, pos);
    }

    // Returns true if query is within one base of the breakend.
    bool nearBreakend(int query) const
    {
        if (pos == -1)
            return false;  // invalid/sentinel has no breakends
        
        switch (kind)
        {
            case INDEL:
                if (size > 0)
                    return (query == pos || query == pos + 1);
                else
                    return (query == pos || query == pos + 1 ||
                            query == pos + size || query == pos + size + 1);
            case INVERSION:
                    return (query == pos || query == pos + 1 ||
                            query == pos + size || query == pos + size + 1);
            case TRANSLOCATION:
                    return (query == pos || query == pos + 1 ||
                            query == targetPos || query == targetPos + 1);
            case DUPLICATION:
                    return (query == pos || query == pos + 1 ||
                            query == targetPos || query == targetPos + 1);
            default:
                return false;
        }
    }

    // Returns end position of SV.
    //
    // Return max value of int in case it is marked as invalid.
    int endPosition() const
    {
        if (pos == -1)
            return seqan::maxValue<int>();
        
        switch (kind)
        {
            case INDEL:
                if (size > 0)
                    return pos;
                else
                    return pos + size;
            case INVERSION:
                return pos + size;
            case TRANSLOCATION:
                return targetPos;
            case DUPLICATION:
                return targetPos;
            default:
                return -1;
        }
    }

    // Return true if this SV overlaps with other within one bp.
    bool overlapsWith(StructuralVariantRecord const & other) const
    {
        return (other.pos <= endPosition()) && (pos <= other.endPosition());
    }
};

inline std::ostream & operator<<(std::ostream & out, StructuralVariantRecord const & record)
{
    char const * kind;
    switch (record.kind)
    {
        case StructuralVariantRecord::TRANSLOCATION:
            kind = "TRANSLOCATION";
            break;
        case StructuralVariantRecord::INDEL:
            kind = "INDEL";
            break;
        case StructuralVariantRecord::INVERSION:
            kind = "INVERSION";
            break;
        case StructuralVariantRecord::DUPLICATION:
            kind = "DUPLICATION";
            break;
        default:
            kind = "INVALID";
            break;
    }
    out << "StructuralVariantRecord(kind=" << kind << ", haplotype=" << record.haplotype
        << ", rId=" << record.rId << ", pos=" << record.pos << ", size=" << record.size
        << ", targetRId=" << record.targetRId << ", targetPos=" << record.targetPos
        << ", seq=\"" << record.seq << "\")";
    return out;
}

// --------------------------------------------------------------------------
// Class Variants
// --------------------------------------------------------------------------

// Contains the simulated variants.

struct Variants
{
    // SNP records.
    seqan::String<SnpRecord> snps;

    // Small indel records.
    seqan::String<SmallIndelRecord> smallIndels;

    // Structural variation record.
    seqan::String<StructuralVariantRecord> svRecords;
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
            rng(rng), faiIndex(faiIndex), variationSizeRecords(variationSizeRecords), options(options)
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

            if (options.verbosity >= 3)
                std::cerr << back(variants.svRecords) << "\n";
        }
    }

    bool simulateSVIndel(Variants & variants, int haploCount, int rId, unsigned pos, int size)
    {
        // Indels are simulated for one haplotype only.
        int hId = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, haploCount - 1));
        seqan::CharString indelSeq;
        reserve(indelSeq, options.maxSVSize);
        bool deletion = (size < 0);
        if (deletion && (pos + size) > (int)sequenceLength(faiIndex, rId))
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
            while (pos > svRecord.endPosition())
            {
                // Skip if near breakend.
                skip = svRecord.nearBreakend(pos);
                
                svIdx += 1;
                if (svIdx < length(variants.svRecords))
                    svRecord = variants.svRecords[svIdx];
                else
                    svRecord.pos = -1;  // mark as sentinel
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
                isSnp = isIndel = false;

            // Simulate either SNP or indel.  In the case of a deletion, advance position such that there
            // is no variation in the deleted sequence.
            if (isSnp)
            {
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
        if (deletion && (pos + indelSize) > (int)sequenceLength(faiIndex, rId))
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

    // FAI Index for loading sequence contig-wise.
    seqan::FaiIndex const & faiIndex;

    // Variation size record.
    seqan::String<VariationSizeRecord> variationSizeRecords;

    MasonVariatorApp(TRng & rng, seqan::FaiIndex const & faiIndex,
                     MasonVariatorOptions const & options) :
            rng(rng), options(options), faiIndex(faiIndex)
    {}

    int run()
    {
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

        // Read in variant size TSV if path is given.
        if (_readVariationSizes() != 0)
            return 1;

        // Actually perform the variant simulation.
        if (options.verbosity >= 1)
            std::cerr << "\nSimulation...\n";
        StructuralVariantSimulator svSim(rng, faiIndex, variationSizeRecords, options);
        SmallVariantSimulator smallSim(rng, faiIndex, options);
        for (int rId = 0; rId < (int)numSeqs(faiIndex); ++rId)  // ref seqs
            _simulateContig(svSim, smallSim, options, rId);
        if (options.verbosity >= 1)
            std::cerr << "OK.\n\n";

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
    int _simulateContig(StructuralVariantSimulator & svSim,
                        SmallVariantSimulator & smallSim,
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
        seqan::Dna5String contig;
        if (readSequence(contig, faiIndex, rId) != 0)
        {
            std::cerr << "Could not read contig seq " << rId << "\n";
            return 1;
        }

        // Write out variants for contig to VCF file.
        if (_writeVcf(contig, variants, rId) != 0)
            return 1;

        // Apply variants to contigs.
        if (!empty(options.fastaOutFile))
            for (int hId = 0; hId < options.numHaplotypes; ++hId)
            {
                if (_writeContigs(contig, variants, rId, hId) != 0)
                    return 1;
            }

        return 0;
    }

    int _writeContigs(seqan::Dna5String const & contig, Variants const & variants, int rId, int hId)
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
        if (_materializeSmallVariants(seqSmallVariants, journal, contig, variants, rId, hId) != 0)
            return 1;

        // Apply structural variants.
        seqan::Dna5String seqLargeVariants;
        if (_materializeLargeVariants(seqLargeVariants, journal, seqSmallVariants, variants, rId, hId) != 0)
            return 1;

        return writeRecord(outSeqStream, id, seqLargeVariants);
    }

    // Apply large structural variants from variants into seq.  The input is given as the sequence including small
    // variants and a journal for translating coordinates in contig into coordinates of the underlying sequence which is
    // used as the coordinate system in variants.
    int _materializeLargeVariants(seqan::Dna5String & seq, TJournalEntries const & journal,
                                  seqan::Dna5String const & contig,
                                  Variants const & variants, int rId, int hId)
    {
        if (options.verbosity >= 2)
            std::cerr << "\nMATERIALIZING LARGE VARIANTS FOR HAPLOTYPE " << hId << "\n\n";
        
        // Track last position from contig appended to seq so far.
        int lastPos = 0;
        if (options.verbosity >= 3)
            std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";

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
            svRecord.size = hostToVirtualPosition(journal, svRecord.pos + svRecord.size) -
                    hostToVirtualPosition(journal, svRecord.pos);
            if (svRecord.targetPos != -1)
                svRecord.targetPos = hostToVirtualPosition(journal, svRecord.targetPos);
            if (options.verbosity >= 2)
                std::cerr << "  => " << svRecord << '\n';

            // Copy from contig to seq with SVs.
            append(seq, infix(contig, lastPos, svRecord.pos));  // interim chars
            if (options.verbosity >= 3)
                std::cerr << "append(seq, infix(contig, " << lastPos << ", " << svRecord.pos << ") " << __LINE__ << " (interim)\n";
            switch (svRecord.kind)
            {
                case StructuralVariantRecord::INDEL:
                    {
                        if (svRecord.size > 0)  // insertion
                        {
                            SEQAN_ASSERT_EQ((int)length(svRecord.seq), svRecord.size);
                            append(seq, svRecord.seq);
                            if (options.verbosity >= 3)
                                std::cerr << "append(seq, svRecord.seq (length == " << length(svRecord.seq) << ") " << __LINE__ << " (insertion)\n";
                            lastPos = svRecord.pos;
                        }
                        else  // deletion
                        {
                            lastPos = svRecord.pos - svRecord.size;
                        }
                    }
                    break;
                case StructuralVariantRecord::INVERSION:
                    {
                        unsigned oldLen = length(seq);
                        append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                        if (options.verbosity >= 3)
                            std::cerr << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << " (inversion)\n";
                        reverseComplement(infix(seq, oldLen, length(seq)));
                        lastPos += svRecord.size;
                    }
                    break;
                case StructuralVariantRecord::TRANSLOCATION:
                    {
                        SEQAN_ASSERT_GEQ(svRecord.targetPos, svRecord.pos + svRecord.size);
                        append(seq, infix(contig, svRecord.pos + svRecord.size, svRecord.targetPos));
                        append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                        if (options.verbosity >= 3)
                            std::cerr << "append(seq, infix(contig, " << svRecord.pos + svRecord.size << ", " << svRecord.targetPos << ") " << __LINE__ << " (translocation)\n"
                                      << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << "\n";
                        lastPos = svRecord.targetPos;
                    }
                    break;
                case StructuralVariantRecord::DUPLICATION:
                    {
                        append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                        SEQAN_ASSERT_GEQ(svRecord.targetPos, svRecord.pos + svRecord.size);
                        append(seq, infix(contig, svRecord.pos + svRecord.size, svRecord.targetPos));
                        append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                        std::cerr << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << " (duplication)\n"
                                  << "append(seq, infix(contig, " << svRecord.pos + svRecord.size << ", " << svRecord.targetPos << ") " << __LINE__ << "\n"
                                  << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << "\n";
                        lastPos = svRecord.targetPos;
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
        
        return 0;
    }

    // Apply small indels and SNPs from variants into seq using contig.
    int _materializeSmallVariants(seqan::Dna5String & seq, TJournalEntries & journal,
                                  seqan::Dna5String const & contig,
                                  Variants const & variants, int rId, int hId)
    {
        reinit(journal, length(contig));

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
                    append(seq, infix(contig, lastPos, snpRecord.pos));  // interim chars
                       
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
                        append(seq, infix(contig, lastPos, smallIndelRecord.pos));  // interim chars

                        SEQAN_ASSERT_GEQ(smallIndelRecord.pos, lastPos);
                        if (options.verbosity >= 3)
                            std::cerr << "append(seq, \"" << smallIndelRecord.seq << "\") " << __LINE__ << "\n";
                        append(seq, smallIndelRecord.seq);
                        lastPos = smallIndelRecord.pos;
                        recordInsertion(journal, hostToVirtualPosition(journal, smallIndelRecord.pos),
                                        0, smallIndelRecord.size);
                        if (options.verbosity >= 3)
                            std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";
                    }
                    else
                    {
                        if (options.verbosity >= 3)
                            std::cerr << "append(seq, infix(contig, " << lastPos << ", " << smallIndelRecord.pos << ") " << __LINE__ << "\n";
                        append(seq, infix(contig, lastPos, smallIndelRecord.pos));  // interim chars
                        lastPos = smallIndelRecord.pos - smallIndelRecord.size;
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
        // Insert interim characters.
        if (options.verbosity >= 3)
            std::cerr << "append(seq, infix(contig, " << lastPos << ", " << length(contig) << ")\n";
        append(seq, infix(contig, lastPos, length(contig)));

        return 0;
    }

    // Write out variants for the given contig to the VCF file.
    int _writeVcf(seqan::Dna5String const & contig, Variants const & variants, int rId)
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
                SEQAN_ASSERT(svRecord.getPos() != svRecord.getPos());  // are generated indendently
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
        vcfRecord.chromId = rId;
        vcfRecord.pos = pos.second;
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
        vcfRecord.chromId = front(records).rId;
        vcfRecord.pos = front(records).pos;
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
        append(vcfRecord.ref, infix(contig, vcfRecord.pos - 1, vcfRecord.pos - 1 + numRef));

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
        vcfRecord.chromId = svRecord.rId;
        vcfRecord.pos = svRecord.pos;
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
        append(vcfRecord.ref, infix(contig, vcfRecord.pos - 1, vcfRecord.pos - 1 + numRef));

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
        leftOfCutL.chromId = svRecord.rId;
        rightOfCutL.chromId = svRecord.rId;
        leftOfCutR.chromId = svRecord.rId;
        rightOfCutR.chromId = svRecord.rId;
        leftOfPaste.chromId = svRecord.rId;
        rightOfPaste.chromId = svRecord.rId;
        // POS
        leftOfCutL.pos = svRecord.pos - 1;
        rightOfCutL.pos = svRecord.pos;
        leftOfCutR.pos = svRecord.pos - 1 + svRecord.size;
        rightOfCutR.pos = svRecord.pos + svRecord.size;
        leftOfPaste.pos = svRecord.targetPos - 1;
        rightOfPaste.pos = svRecord.targetPos;
        // TODO(holtgrew): INFO entry with type of breakend?
        // ID (none)
        // TODO(holtgrew): Generate an id?
        // REF
        appendValue(leftOfCutL.ref, contig[leftOfCutL.pos]);
        appendValue(rightOfCutL.ref, contig[rightOfCutL.pos]);
        appendValue(leftOfCutR.ref, contig[leftOfCutR.pos]);
        appendValue(rightOfCutR.ref, contig[rightOfCutR.pos]);
        appendValue(leftOfPaste.ref, contig[leftOfPaste.pos]);
        appendValue(rightOfPaste.ref, contig[rightOfPaste.pos]);
        // ALT
        std::stringstream ssLeftOfCutL, ssRightOfCutL, ssLeftOfCutR, ssRightOfCutR,
                ssLeftOfPaste, ssRightOfPaste;
        seqan::CharString refName = vcfStream.header.sequenceNames[svRecord.rId];
        ssLeftOfCutL << leftOfCutL.ref << "]" << refName << ":" << (rightOfCutR.pos + 1) << "]";
        leftOfCutL.alt = ssLeftOfCutL.str();
        ssRightOfCutL << "[" << refName << ":" << (leftOfPaste.pos + 1) << "[" << rightOfCutL.ref;
        rightOfCutL.alt = ssRightOfCutL.str();
        ssLeftOfCutR << leftOfCutR.ref << "]" << refName << ":" << (rightOfPaste.pos + 1) << "]";
        leftOfCutR.alt = ssLeftOfCutR.str();
        ssRightOfCutR << "[" << refName << ":" << (leftOfCutL.pos + 1) << "[" << rightOfCutR.ref;
        rightOfCutR.alt = ssRightOfCutR.str();
        ssLeftOfPaste << leftOfPaste.ref << "]" << refName << ":" << (leftOfCutR.pos + 1) << "]";
        leftOfPaste.alt = ssLeftOfPaste.str();
        ssRightOfPaste << "[" << refName << ":" << (rightOfCutL.pos + 1) << "[" << rightOfPaste.ref;
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

        vcfRecord.chromId = svRecord.rId;
        vcfRecord.pos = svRecord.pos;
        // TODO(holtgrew): Generate an id?
        appendValue(vcfRecord.ref, contig[vcfRecord.pos]);
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
        vcfRecord.chromId = svRecord.rId;
        vcfRecord.pos = svRecord.pos;
        // TODO(holtgrew): Generate an id?
        vcfRecord.filter = "PASS";
        std::stringstream ss;
        ss << "SVTYPE=DUP;SVLEN=" << svRecord.size << ";END=" << svRecord.pos + svRecord.size;
        vcfRecord.info = ss.str();
        appendValue(vcfRecord.ref, contig[vcfRecord.pos]);
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
    setVersion(parser, "2.1");
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
    setValidValues(parser, "out-vcf", "vcf");

    addOption(parser, seqan::ArgParseOption("of", "out-fasta", "FASTA file to write simulated haplotypes to.",
                                            seqan::ArgParseOption::INPUTFILE, "FASTA"));
    setValidValues(parser, "out-fasta", "fasta fa");

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
    setDefaultValue(parser, "haplotype-sep", ":");

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

    std::cerr << "__WRITING OUTPUT______________________________________________________________\n"
              << "\n";

    if (!empty(options.vcfOutFile))
        std::cerr << "Writing VCF to " << options.vcfOutFile << "\n";
    
    if (!empty(options.fastaOutFile))
        std::cerr << "Writing FASTA to " << options.fastaOutFile << "\n";

    std::cerr << "\nDONE.\n";
    
    return 0;
}
