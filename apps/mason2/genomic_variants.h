// ==========================================================================
//                         Mason - A Read Simulator
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
// Data structures for storing genomic variants and routines to work on them.
//
// Variants are split into SNPs, small indels, and large structural variants.
// There is a data structure for each class of variant.  Doing such a
// separation allows for almost simple translation between the reference
// coordinate system and the one for the genome including all variants.
//
// There are routines provided for 
// ==========================================================================

#ifndef SANDBOX_MASON2_APPS_MASON2_GENOMIC_VARIANTS_H_
#define SANDBOX_MASON2_APPS_MASON2_GENOMIC_VARIANTS_H_

#include "methylation_levels.h"

// ============================================================================
// Forwards
// ============================================================================

typedef seqan::JournalEntries<seqan::JournalEntry<unsigned, int>, seqan::SortedArray> TJournalEntries;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

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
                            query == pos - size || query == pos - size + 1);
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

    void clear()
    {
        seqan::clear(snps);
        seqan::clear(smallIndels);
        seqan::clear(svRecords);
    }
};

// --------------------------------------------------------------------------
// Class VariantMaterializer
// --------------------------------------------------------------------------

// Materialize variants stored in a Variants object.
//
// Note that the class assumes that all variants come from the same contig and haplotype.

// TODO(holtgrew): Rename to ContigMaterializer.
// TODO(holtgrew): Add translation of coordinates.

class VariantMaterializer
{
public:
    // The random number generator to use.
    TRng * rng;
    // The Variants to materialize for.
    Variants const * variants;
    // Options for the methylation level simulator.  Methylation simulation is required for fixing methylation levels.
    MethylationLevelSimulatorOptions methSimOptions;

    // Verbosity.
    int verbosity;

    VariantMaterializer() : rng(), variants(), verbosity(1)
    {}

    VariantMaterializer(TRng & rng, Variants const & variants) :
            rng(&rng), variants(&variants),  verbosity(1)
    {}

    VariantMaterializer(TRng & rng, Variants const & variants, MethylationLevelSimulatorOptions const & methSimOptions) :
            rng(&rng), variants(&variants),  methSimOptions(methSimOptions), verbosity(1)
    {}

    // Materialize the variants from the haplotype with the given id in *variants to result given the reference sequence refSeq.
    //
    // Breakpoints is a vector of points on the contig.
    int run(seqan::Dna5String & resultSeq,
             std::vector<int> & breakpoints,
             seqan::Dna5String const & refSeq,
             int haplotypeId)
    {
        return _runImpl(&resultSeq, 0, breakpoints, &refSeq, 0, haplotypeId);
    }

    // Same as the run() above, but including reference levels.
    int run(seqan::Dna5String & resultSeq,
            MethylationLevels & resultLvls,
            std::vector<int> & breakpoints,
            seqan::Dna5String const & refSeq,
            MethylationLevels const & refLvls,
            int haplotypeId)
    {
        return _runImpl(&resultSeq, &resultLvls, breakpoints, &refSeq, &refLvls, haplotypeId);
    }

    // Implementation of the materialization, uses pointers instead of references for deciding whether materializing
    // levels or not.
    int _runImpl(seqan::Dna5String * resultSeq,
                 MethylationLevels * resultLvls,
                 std::vector<int> & breakpoints,
                 seqan::Dna5String const * ref,
                 MethylationLevels const * refLvls,
                 int haplotypeId);

    // Materialization of the small variants.
    //
    // Levels passed as NULL if not given.
    int _materializeSmallVariants(seqan::Dna5String & seq,
                                  TJournalEntries & journal,
                                  MethylationLevels * levelsSmallVariants,
                                  seqan::Dna5String const & contig,
                                  Variants const & variants,
                                  MethylationLevels const * levels,
                                  int hId);

    // Materialization of the large variants.
    //
    // Levels passed as NULL if not given.
    int _materializeLargeVariants(seqan::Dna5String & seq,
                                  MethylationLevels * levelsLargeVariants,
                                  std::vector<int> & breakpoints,
                                  TJournalEntries const & journal,
                                  seqan::Dna5String const & contig,
                                  Variants const & variants,
                                  MethylationLevels const * levels,
                                  int hId);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function VariantMaterializer::_runImpl()
// ----------------------------------------------------------------------------

int VariantMaterializer::_runImpl(
        seqan::Dna5String * resultSeq,
        MethylationLevels * resultLvls,
        std::vector<int> & breakpoints,
        seqan::Dna5String const * ref,
        MethylationLevels const * refLvls,
        int haplotypeId)
{
    breakpoints.clear();
    clear(*resultSeq);
    if (resultLvls)
        resultLvls->clear();
    
    // Apply small variants.  We get a sequence with the small variants and a journal of the difference to contig.
    seqan::Dna5String seqSmallVariants;
    TJournalEntries journal;
    MethylationLevels levelsSmallVariants;  // only used if revLevels != 0
    MethylationLevels * smallLvlsPtr = refLvls ? &levelsSmallVariants : 0;
    if (_materializeSmallVariants(seqSmallVariants, journal, smallLvlsPtr, *ref, *variants,
                                  refLvls, haplotypeId) != 0)
        return 1;
    
    // Apply structural variants.
    if (_materializeLargeVariants(*resultSeq, resultLvls, breakpoints, journal, seqSmallVariants, *variants,
                                  smallLvlsPtr, haplotypeId) != 0)
        return 1;

    return 0;
}

// ----------------------------------------------------------------------------
// Function VariantMaterializer::_materializeSmallVariants()
// ----------------------------------------------------------------------------

int VariantMaterializer::_materializeSmallVariants(
        seqan::Dna5String & seq,
        TJournalEntries & journal,
        MethylationLevels * levelsSmallVariants,
        seqan::Dna5String const & contig,
        Variants const & variants,
        MethylationLevels const * levels,
        int hId)
{
    SEQAN_ASSERT_EQ(methSimOptions.simulateMethylationLevels, (levelsSmallVariants != 0));
    SEQAN_ASSERT_EQ(methSimOptions.simulateMethylationLevels, (levels != 0));
    
    // Clear journal and output methylation levels.
    reinit(journal, length(contig));
    if (levelsSmallVariants)
        levelsSmallVariants->clear();
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
    if (verbosity >= 3)
        std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";

    // TODO(holtgrew): Extract contig building into their own functions.
    if (verbosity >= 2)
        std::cerr << "building output\n";
    while (snpRecord.rId != seqan::maxValue<int>() || smallIndelRecord.rId != seqan::maxValue<int>())
    {
        // TODO(holtgrew): Extract SNP and small indel handling in functions.
        if (snpRecord.getPos() < smallIndelRecord.getPos())  // process SNP records
        {
            if (snpRecord.haplotype == hId)  // Ignore all but the current contig.
            {
                if (verbosity >= 3)
                    std::cerr << "append(seq, infix(contig, " << lastPos << ", " << snpRecord.pos << ") " << __LINE__ << "\n";
                // Append interim sequence and methylation levels->
                append(seq, infix(contig, lastPos, snpRecord.pos));
                if (methSimOptions.simulateMethylationLevels)
                {
                    append(levelsSmallVariants->forward, infix(levels->forward, lastPos, snpRecord.pos + 1));
                    append(levelsSmallVariants->reverse, infix(levels->reverse, lastPos, snpRecord.pos + 1));
                    appendValue(varPoints, std::make_pair(length(seq), true));      // variation points before/after SNP
                    appendValue(varPoints, std::make_pair(length(seq) + 1, true));
                }

                SEQAN_ASSERT_GEQ(snpRecord.pos, lastPos);
                if (verbosity >= 3)
                    std::cerr << "appendValue(seq, " << snpRecord.to << "')\n";
                appendValue(seq, snpRecord.to);
                lastPos = snpRecord.pos + 1;
                if (verbosity >= 3)
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
                    if (verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << lastPos << ", " << smallIndelRecord.pos << ") " << __LINE__ << "\n";

                    // Simulate methylation levels for insertion.
                    MethylationLevels lvls;
                    if (methSimOptions.simulateMethylationLevels)
                    {
                        MethylationLevelSimulator methSim(*rng, methSimOptions);
                        methSim.run(lvls, smallIndelRecord.seq);
                    }

                    // Append interim sequence and methylation levels->
                    append(seq, infix(contig, lastPos, smallIndelRecord.pos));
                    if (methSimOptions.simulateMethylationLevels)
                    {
                        append(levelsSmallVariants->forward, infix(levels->forward, lastPos, smallIndelRecord.pos));
                        append(levelsSmallVariants->reverse, infix(levels->reverse, lastPos, smallIndelRecord.pos));
                        appendValue(varPoints, std::make_pair(length(seq), false));  // variation point before insertion
                    }

                    SEQAN_ASSERT_GEQ(smallIndelRecord.pos, lastPos);
                    if (verbosity >= 3)
                        std::cerr << "append(seq, \"" << smallIndelRecord.seq << "\") " << __LINE__ << "\n";
                    // Append novel sequence and methylation levels->
                    append(seq, smallIndelRecord.seq);
                    if (methSimOptions.simulateMethylationLevels)
                    {
                        append(levelsSmallVariants->forward, lvls.forward);
                        append(levelsSmallVariants->reverse, lvls.reverse);
                        appendValue(varPoints, std::make_pair(length(seq), false));  // variation point after insertion
                    }
                    lastPos = smallIndelRecord.pos;
                    recordInsertion(journal, hostToVirtualPosition(journal, smallIndelRecord.pos),
                                    0, smallIndelRecord.size);
                    if (verbosity >= 3)
                        std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";
                }
                else  // deletion
                {
                    if (verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << lastPos << ", " << smallIndelRecord.pos << ") " << __LINE__ << "\n";
                    // Append interim sequence and methylation levels->
                    append(seq, infix(contig, lastPos, smallIndelRecord.pos));  // interim chars
                    if (methSimOptions.simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair(length(seq), false));  // variation point at deletion
                        append(levelsSmallVariants->forward, infix(levels->forward, lastPos, smallIndelRecord.pos));
                        append(levelsSmallVariants->reverse, infix(levels->reverse, lastPos, smallIndelRecord.pos));
                    }

                    lastPos = smallIndelRecord.pos - smallIndelRecord.size;
                    recordErase(journal,
                                hostToVirtualPosition(journal, smallIndelRecord.pos),
                                hostToVirtualPosition(journal, smallIndelRecord.pos - smallIndelRecord.size));
                    if (verbosity >= 3)
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
    if (verbosity >= 3)
        std::cerr << "append(seq, infix(contig, " << lastPos << ", " << length(contig) << ")\n";
    append(seq, infix(contig, lastPos, length(contig)));

    if (methSimOptions.simulateMethylationLevels)
    {
        append(levelsSmallVariants->forward, infix(levels->forward, lastPos, length(contig)));
        append(levelsSmallVariants->reverse, infix(levels->reverse, lastPos, length(contig)));

        SEQAN_ASSERT_EQ(length(seq), length(levelsSmallVariants->forward));
        SEQAN_ASSERT_EQ(length(seq), length(levelsSmallVariants->reverse));

        fixVariationLevels(*levelsSmallVariants, *rng, seq, varPoints, methSimOptions);
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function VariantMaterializer::_materializeLargeVariants()
// ----------------------------------------------------------------------------

int VariantMaterializer::_materializeLargeVariants(
        seqan::Dna5String & seq,
        MethylationLevels * levelsLargeVariants,
        std::vector<int> & breakpoints,
        TJournalEntries const & journal,
        seqan::Dna5String const & contig,
        Variants const & variants,
        MethylationLevels const * levels,
        int hId)
{
    SEQAN_ASSERT_EQ(methSimOptions.simulateMethylationLevels, (levelsLargeVariants != 0));
    SEQAN_ASSERT_EQ(methSimOptions.simulateMethylationLevels, (levels != 0));

    // Clear output methylation levels->
    if (levelsLargeVariants)
        levelsLargeVariants->clear();
    // Store variation points.  We reuse the fixVariationLevels() function from small indel/snp simulation and thus
    // have to store a bool that is always set to false.
    seqan::String<std::pair<int, bool> > varPoints;
        
    // Track last position from contig appended to seq so far.
    int lastPos = 0;
    if (verbosity >= 3)
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
        if (verbosity >= 2)
            std::cerr << "  Translating SvRecord\n  " << svRecord << '\n';
        svRecord.pos = hostToVirtualPosition(journal, svRecord.pos);
        SEQAN_ASSERT_LT(svRecord.pos, (int)length(contig));
        // We do not need to adjust the sizes for insertions.
        if (svRecord.kind != StructuralVariantRecord::INDEL || svRecord.size < 0)
            svRecord.size = hostToVirtualPosition(journal, svRecord.pos + svRecord.size) -
                    hostToVirtualPosition(journal, svRecord.pos);
        if (svRecord.targetPos != -1)
            svRecord.targetPos = hostToVirtualPosition(journal, svRecord.targetPos);
        if (verbosity >= 2)
            std::cerr << "  => " << svRecord << '\n';

        // Copy from contig to seq with SVs.
        if (verbosity >= 3)
            std::cerr << "lastPos == " << lastPos << "\n";
        append(seq, infix(contig, lastPos, svRecord.pos));  // interim chars
        if (methSimOptions.simulateMethylationLevels)
        {
            append(levelsLargeVariants->forward, infix(levels->forward, lastPos, svRecord.pos));
            append(levelsLargeVariants->reverse, infix(levels->reverse, lastPos, svRecord.pos));
            appendValue(varPoints, std::make_pair(length(seq), false));
        }
        currentPos = length(seq);
        if (verbosity >= 3)
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
                        if (methSimOptions.simulateMethylationLevels)
                        {
                            MethylationLevelSimulator methSim(*rng, methSimOptions);
                            methSim.run(lvls, svRecord.seq);
                        }

                        // Append novel sequence and methylation levels->
                        append(seq, svRecord.seq);
                        if (methSimOptions.simulateMethylationLevels)
                        {
                            append(levelsLargeVariants->forward, lvls.forward);
                            append(levelsLargeVariants->reverse, lvls.reverse);
                            appendValue(varPoints, std::make_pair(length(seq), false));  // variation point after insertion
                        }
                        if (verbosity >= 3)
                            std::cerr << "append(seq, svRecord.seq (length == " << length(svRecord.seq) << ") " << __LINE__ << " (insertion)\n";
                        lastPos = svRecord.pos;
                        SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                        // Copy out breakpoints.
                        breakpoints.push_back(currentPos);
                        breakpoints.push_back(length(seq));

                        currentPos = length(seq);
                    }
                    else  // deletion
                    {
                        lastPos = svRecord.pos - svRecord.size;
                        SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                        // Copy out breakpoint.
                        breakpoints.push_back(currentPos);
                    }
                }
                break;
            case StructuralVariantRecord::INVERSION:
                {
                    unsigned oldLen = length(seq);
                    append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                    if (methSimOptions.simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair(length(seq), false));  // variation point at deletion
                        append(levelsLargeVariants->forward, infix(levels->reverse, svRecord.pos, svRecord.pos + svRecord.size));
                        reverse(infix(levelsLargeVariants->forward, oldLen, length(levelsLargeVariants->forward)));
                        append(levelsLargeVariants->reverse, infix(levels->forward, svRecord.pos, svRecord.pos + svRecord.size));
                        reverse(infix(levelsLargeVariants->reverse, oldLen, length(levelsLargeVariants->reverse)));
                    }

                    if (verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << " (inversion)\n";
                    reverseComplement(infix(seq, oldLen, length(seq)));
                    lastPos = svRecord.pos + svRecord.size;
                    SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                    // Copy out breakpoints.
                    breakpoints.push_back(currentPos);
                    breakpoints.push_back(length(seq));

                    currentPos = length(seq);
                }
                break;
            case StructuralVariantRecord::TRANSLOCATION:
                {
                    SEQAN_ASSERT_GEQ(svRecord.targetPos, svRecord.pos + svRecord.size);
                    append(seq, infix(contig, svRecord.pos + svRecord.size, svRecord.targetPos));
                    if (methSimOptions.simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair(length(seq), false));
                        append(levelsLargeVariants->forward, infix(levels->forward, svRecord.pos + svRecord.size, svRecord.targetPos));
                        append(levelsLargeVariants->reverse, infix(levels->reverse, svRecord.pos + svRecord.size, svRecord.targetPos));
                    }
                    append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                    if (methSimOptions.simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair(length(seq), false));
                        append(levelsLargeVariants->forward, infix(levels->forward, svRecord.pos, svRecord.pos + svRecord.size));
                        append(levelsLargeVariants->reverse, infix(levels->reverse, svRecord.pos, svRecord.pos + svRecord.size));
                    }
                    if (verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << svRecord.pos + svRecord.size << ", " << svRecord.targetPos << ") " << __LINE__ << " (translocation)\n"
                                  << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << "\n";
                    lastPos = svRecord.targetPos;
                    SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                    // Copy out breakpoints.
                    breakpoints.push_back(currentPos);
                    breakpoints.push_back(currentPos + svRecord.targetPos - svRecord.pos - svRecord.size);
                    breakpoints.push_back(length(seq));

                    currentPos = length(seq);
                }
                break;
            case StructuralVariantRecord::DUPLICATION:
                {
                    append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                    SEQAN_ASSERT_GEQ(svRecord.targetPos, svRecord.pos + svRecord.size);
                    if (methSimOptions.simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair(length(seq), false));
                        append(levelsLargeVariants->forward, infix(levels->forward, svRecord.pos, svRecord.pos + svRecord.size));
                        append(levelsLargeVariants->reverse, infix(levels->reverse, svRecord.pos, svRecord.pos + svRecord.size));
                    }
                    append(seq, infix(contig, svRecord.pos + svRecord.size, svRecord.targetPos));
                    if (methSimOptions.simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair(length(seq), false));
                        append(levelsLargeVariants->forward, infix(levels->forward, svRecord.pos + svRecord.size, svRecord.targetPos));
                        append(levelsLargeVariants->reverse, infix(levels->reverse, svRecord.pos + svRecord.size, svRecord.targetPos));
                    }
                    append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                    if (methSimOptions.simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair(length(seq), false));
                        append(levelsLargeVariants->forward, infix(levels->forward, svRecord.pos, svRecord.pos + svRecord.size));
                        append(levelsLargeVariants->reverse, infix(levels->reverse, svRecord.pos, svRecord.pos + svRecord.size));
                    }
                    if (verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << " (duplication)\n"
                                  << "append(seq, infix(contig, " << svRecord.pos + svRecord.size << ", " << svRecord.targetPos << ") " << __LINE__ << "\n"
                                  << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << "\n";
                    lastPos = svRecord.targetPos;
                    SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                    // Copy out breakpoints.
                    breakpoints.push_back(currentPos);
                    breakpoints.push_back(currentPos + svRecord.pos + svRecord.size - svRecord.pos);
                    breakpoints.push_back(currentPos + svRecord.pos + svRecord.size - svRecord.pos + svRecord.targetPos - (svRecord.pos + svRecord.size));
                    breakpoints.push_back(length(seq));

                    currentPos = length(seq);
                }
                break;
            default:
                return 1;
        }
    }
    if (verbosity >= 3)
        std::cerr << "append(seq, infix(contig, " << lastPos << ", " << length(contig) << ") "
                  << __LINE__ << " (last interim)\n";
    append(seq, infix(contig, lastPos, length(contig)));
    if (methSimOptions.simulateMethylationLevels)
    {
        append(levelsLargeVariants->forward, infix(levels->forward, lastPos, length(contig)));
        append(levelsLargeVariants->reverse, infix(levels->reverse, lastPos, length(contig)));

        SEQAN_ASSERT_EQ(length(seq), length(levelsLargeVariants->forward));
        SEQAN_ASSERT_EQ(length(seq), length(levelsLargeVariants->reverse));

        fixVariationLevels(*levelsLargeVariants, *rng, seq, varPoints, methSimOptions);
    }
        
    return 0;
}

#endif  // #ifndef SANDBOX_MASON2_APPS_MASON2_GENOMIC_VARIANTS_H_
