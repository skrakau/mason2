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

// ============================================================================
// Forwards
// ============================================================================

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
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_MASON2_APPS_MASON2_GENOMIC_VARIANTS_H_
