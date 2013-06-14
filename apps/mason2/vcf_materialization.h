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

#ifndef SANDBOX_MASON2_APPS_MASON2_VCF_MATERIALIZATION_H_
#define SANDBOX_MASON2_APPS_MASON2_VCF_MATERIALIZATION_H_

// ============================================================================
// Forwards
// ============================================================================

#include <stdexcept>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/vcf_io.h>

#include "genomic_variants.h"

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class VcfMaterializer
// ----------------------------------------------------------------------------

// Allows the contig- and haplotype-wise construction of haplotypes stored in VCF files.

class VcfMaterializer
{
public:
    // ------------------------------------------------------------------------
    // Paths
    // ------------------------------------------------------------------------

    // Path to reference file.
    seqan::CharString fastaFileName;
    // Path to VCF file.
    seqan::CharString vcfFileName;

    // ------------------------------------------------------------------------
    // State for position in reference
    // ------------------------------------------------------------------------

    // The index of the current contig.
    int currRID;
    // The index of the current haplotype of the contig.
    int nextHaplotype;
    // The number of haplotypes (set after call to init()).
    int numHaplotypes;
    // The variants for the current contig.
    Variants contigVariants;
    // Current contig reference sequence.
    seqan::Dna5String contigSeq;

    // ------------------------------------------------------------------------
    // File Input
    // ------------------------------------------------------------------------

    // The FAI Index to load the reference sequence from.
    seqan::FaiIndex faiIndex;
    // The VCF stream to load from.
    seqan::VcfStream vcfStream;
    // The current VCF record.  rID == INVALID_REFID if invalid, used for termination.
    seqan::VcfRecord vcfRecord;

    VcfMaterializer() : currRID(0), nextHaplotype(0), numHaplotypes(0)
    {}

    VcfMaterializer(char const * fastaFileName, char const * vcfFileName) :
            fastaFileName(fastaFileName), vcfFileName(vcfFileName), currRID(-1), nextHaplotype(0), numHaplotypes(0)
    {}

    // Call to open all files.
    //
    // Throws: MasonIOException
    void init();

    // Materialize next contig.
    //
    // Write sequence to seq, reference id to rID, haplotype to haplotype.  Returns true if the materialization could be
    // done and false if there are no more contigs to materialize.
    //
    // Call init() before calling materializeNext().
    //
    // Throws: MasonIOException
    bool materializeNext(seqan::Dna5String & seq, int & rID, int & haplotype);

private:

    // Load variants of next contig into variants.
    int _loadVariantsForContig(Variants & variants, int rID);

    // Append VCF record to variants.
    void _appendToVariants(Variants & variants, seqan::VcfRecord const & vcfRecord);

    // Append chunk of 6 BND records to variants.
    void _appendToVariantsBnd(Variants & variants, std::vector<seqan::VcfRecord> const & vcfRecords);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function contains()
// ----------------------------------------------------------------------------

// Returns true if haystack contains needle.

bool contains(seqan::CharString const & haystack, char const * needle)
{
    return strstr(toCString(haystack), needle);
}

// ----------------------------------------------------------------------------
// Function getSVLen()
// ----------------------------------------------------------------------------

// Get integer value of SVLEN=<number> from VCF INFO column strings.

int getSVLen(seqan::CharString const & str)
{
    seqan::RecordReader<seqan::CharString const, seqan::SinglePass<seqan::StringReader> > reader(str);

    // Parse out key/value pairs and interpret SVLEN.
    seqan::CharString key, val;
    enum { IN_KEY, IN_VALUE } state = IN_KEY;
    for (; !atEnd(reader); goNext(reader))
    {
        if (value(reader) == '=')
        {
            state = IN_VALUE;
            continue;
        }
        else if (value(reader) == ';')
        {
            if (key == "SVLEN")
                return seqan::lexicalCast<int>(val);
            
            clear(val);
            clear(key);
            state = IN_KEY;
            continue;
        }
        else if (state == IN_KEY)
        {
            appendValue(key, value(reader));
        }
        else  // (state == IN_VALUE)
        {
            appendValue(val, value(reader));
        }
    }
    
    if (key == "SVLEN")
        return seqan::lexicalCast<int>(val);

    SEQAN_FAIL("Missing INFO SVLEN %s", toCString(str));
    return 0;
}

// ----------------------------------------------------------------------------
// Function getTargetPos()
// ----------------------------------------------------------------------------

// Returns reference name and position from TARGETPOS=<ref>:<pos> in VCF INFO column strings.

std::pair<seqan::CharString, int> getTargetPos(seqan::CharString const & str)
{
    seqan::RecordReader<seqan::CharString const, seqan::SinglePass<seqan::StringReader> > reader(str);

    // Parse out key/value pairs and interpret SVLEN.
    seqan::CharString key, val;
    enum { IN_KEY, IN_VALUE } state = IN_KEY;
    for (; !atEnd(reader); goNext(reader))
    {
        if (value(reader) == '=')
        {
            state = IN_VALUE;
            continue;
        }
        else if (value(reader) == ';')
        {
            if (key == "TARGETPOS")
            {
                seqan::StringSet<seqan::CharString> xs;
                splitString(xs, val, ':');
                SEQAN_CHECK(length(xs) == 2u, "TARGETPOS has invalid format %s", toCString(val));
                SEQAN_CHECK(!empty(xs[0]) && !empty(xs[1]), "TARGETPOS has invalid format %s", toCString(val));
                return std::make_pair(xs[0], seqan::lexicalCast<int>(xs[1]));
            }
            
            clear(val);
            clear(key);
            state = IN_KEY;
            continue;
        }
        else if (state == IN_KEY)
        {
            appendValue(key, value(reader));
        }
        else  // (state == IN_VALUE)
        {
            appendValue(val, value(reader));
        }
    }
    
    if (key == "TARGETPOS")
    {
        seqan::StringSet<seqan::CharString> xs;
        splitString(xs, val, ':');
        SEQAN_CHECK(length(xs) == 2u, "TARGETPOS has invalid format %s", toCString(val));
        SEQAN_CHECK(!empty(xs[0]) && !empty(xs[1]), "TARGETPOS has invalid format %s", toCString(val));
        return std::make_pair(xs[0], seqan::lexicalCast<int>(xs[1]));
    }

    SEQAN_FAIL("Missing INFO TARGETPOS %s", toCString(str));
    return std::make_pair("", 0);
}

// ----------------------------------------------------------------------------
// Function VcfMaterializer::init()
// ----------------------------------------------------------------------------

void VcfMaterializer::init()
{
    // Open VCF stream.
    open(vcfStream, toCString(vcfFileName));
    if (!isGood(vcfStream))
        throw MasonIOException("Could not open VCF stream.");

    // Read first VCF record.
    if (!atEnd(vcfStream) && readRecord(vcfRecord, vcfStream) != 0)
        throw MasonIOException("Problem reading from VCF file.");

    // Get number of haplotypes in VCF file.
    SEQAN_ASSERT_NOT(empty(vcfRecord.genotypeInfos));
    seqan::StringSet<seqan::CharString> xs;
    seqan::RecordReader<seqan::CharString, seqan::SinglePass<seqan::StringReader> >
            reader(vcfRecord.genotypeInfos[0]);
    numHaplotypes = 1;
    for (; !atEnd(reader); goNext(reader))
        numHaplotypes += (value(reader) == '|' || value(reader) == '/');

    // Open input FASTA file and FAI.
    if (read(faiIndex, toCString(fastaFileName)) != 0)
    {
        if (build(faiIndex, toCString(fastaFileName)) != 0)
            throw MasonIOException("Could not build FAI index.");

        seqan::CharString faiPath = fastaFileName;
        append(faiPath, ".fai");
        if (write(faiIndex, toCString(faiPath)) != 0)
            throw MasonIOException("Could not write FAI index.");
    }
}

// ----------------------------------------------------------------------------
// Function VcfMaterializer::materializeNext()
// ----------------------------------------------------------------------------

bool VcfMaterializer::materializeNext(seqan::Dna5String & seq, int & rID, int & haplotype)
{
    // Number of sequences.
    int numSeqs = length(vcfStream.header.sequenceNames);

    // Stop if there are no more haplotypes to materialize.
    if (currRID >= (numSeqs - 1) && nextHaplotype == numHaplotypes)
        return false;

    // Load variants for next contig if necessary.
    if (currRID == -1 || nextHaplotype == numHaplotypes)
    {
        currRID += 1;
        nextHaplotype = 0;

        _loadVariantsForContig(contigVariants, currRID);
        if (readSequence(contigSeq, faiIndex, currRID) != 0)
            throw MasonIOException("Could not load reference sequence.");
    }

    TRng rng;  // Will not be touched by materializer since methylation simulation is not enabled.
    std::vector<int> breakpoints;  // unused

    // Materialize variants for the current haplotype.
    VariantMaterializer varMat(rng, contigVariants);
    varMat.run(seq, breakpoints, contigSeq, nextHaplotype);

    // Write out rID and haploty
    rID = currRID;
    haplotype = nextHaplotype++;
    return true;
}

// ----------------------------------------------------------------------------
// Function VcfMaterializer::_loadVariantsForContig()
// ----------------------------------------------------------------------------

int VcfMaterializer::_loadVariantsForContig(Variants & variants, int rID)
{
    variants.clear();

    // Compute number of haplotypes.
    SEQAN_ASSERT_NOT(empty(vcfRecord.genotypeInfos));
    seqan::StringSet<seqan::CharString> xs;
    seqan::RecordReader<seqan::CharString, seqan::SinglePass<seqan::StringReader> >
            reader(vcfRecord.genotypeInfos[0]);
    numHaplotypes = 1;
    for (; !atEnd(reader); goNext(reader))
        numHaplotypes += (value(reader) == '|' || value(reader) == '/');

    std::vector<seqan::VcfRecord> chunk;
    while (vcfRecord.rID != -1 && vcfRecord.rID <= rID)
    {
        // Translocations are the only SVs that are stored as breakends (BND).  We collect the BNDs in chunks of 6
        // which then represent a translocation.  Requiring that the BNDs are stored in adjacent chunks of 6 records
        // is a strong limitation but supporting more generic variations would be a bit too much here.
        if (contains(vcfRecord.info, "SVTYPE=BND"))
        {
            chunk.push_back(vcfRecord);
            if (chunk.size() == 6u)
            {
                _appendToVariantsBnd(variants, chunk);
                chunk.clear();
            }
        }
        else
        {
            SEQAN_CHECK(chunk.size() == 0u, "Found chunk of != 6 BND records!");
            _appendToVariants(variants, vcfRecord);
        }

        if (atEnd(vcfStream))
        {
            vcfRecord.rID = -1;
            continue;
        }
        if (readRecord(vcfRecord, vcfStream) != 0)
        {
            std::cerr << "ERROR: Problem reading from VCF\n";
            return 1;
        }
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function VcfMaterializer::_appendToVariants()
// ----------------------------------------------------------------------------

void VcfMaterializer::_appendToVariants(Variants & variants, seqan::VcfRecord const & vcfRecord)
{
    // Compute maximal length of alternative.
    unsigned altLength = 0;
    seqan::StringSet<seqan::CharString> alts;
    splitString(alts, vcfRecord.alt, ',');
    for (unsigned i = 0; i < length(alts); ++i)
        altLength = std::max(altLength, (unsigned)length(alts[i]));
        
    if (contains(vcfRecord.info, "SVTYPE"))  // Structural Variant
    {
        StructuralVariantRecord svRecord;
        svRecord.rId = vcfRecord.rID;
        svRecord.pos = vcfRecord.beginPos;
        svRecord.haplotype = 0;

        SEQAN_ASSERT_EQ(length(alts), 1u);
            
        if (contains(vcfRecord.info, "SVTYPE=INS"))  // Insertion
        {
            svRecord.kind = StructuralVariantRecord::INDEL;
            svRecord.size = getSVLen(vcfRecord.info);
            svRecord.seq = suffix(vcfRecord.alt, 1);
        }
        else if (contains(vcfRecord.info, "SVTYPE=DEL"))  // Deletion
        {
            svRecord.kind = StructuralVariantRecord::INDEL;
            svRecord.size = -getSVLen(vcfRecord.info);
        }
        else if (contains(vcfRecord.info, "SVTYPE=INV"))  // Inversion
        {
            svRecord.kind = StructuralVariantRecord::INVERSION;
            svRecord.size = getSVLen(vcfRecord.info);
        }
        else if (contains(vcfRecord.info, "SVTYPE=DUP"))  // Duplication
        {
            svRecord.kind = StructuralVariantRecord::DUPLICATION;
            svRecord.size = getSVLen(vcfRecord.info);
            std::pair<seqan::CharString, int> pos = getTargetPos(vcfRecord.info);
            unsigned idx = 0;
            if (!getIdByName(vcfStream._context.sequenceNames, pos.first, idx,
                             vcfStream._context.sequenceNamesCache))
                SEQAN_FAIL("Unknown sequence %s", toCString(pos.first));
            svRecord.targetRId = idx;
            svRecord.targetPos = pos.second - 1;
        }
        else if (contains(vcfRecord.info, "SVTYPE=BND"))  // Breakend (Must be Translocation)
        {
            SEQAN_FAIL("Unexpected 'SVTYPE=BND' at this place!");
        }
        else
        {
            SEQAN_FAIL("ERROR: Unknown SVTYPE!\n");
        }

        // Split the target variants.
        SEQAN_ASSERT_NOT(empty(vcfRecord.genotypeInfos));
        seqan::RecordReader<seqan::CharString const, seqan::SinglePass<seqan::StringReader> >
                reader(vcfRecord.genotypeInfos[0]);
        seqan::CharString buffer;
        svRecord.haplotype = 0;
        for (; !atEnd(reader); goNext(reader))
            if ((value(reader) == '|' || value(reader) == '/'))
            {
                if (!empty(buffer))
                {
                    unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer), 1u);
                    if (idx != 0u)  // if not == ref
                        appendValue(variants.svRecords, svRecord);
                }
                svRecord.haplotype++;
                clear(buffer);
            }
            else
            {
                appendValue(buffer, value(reader));
            }
        if (!empty(buffer))
        {
            unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer), 1u);
            if (idx != 0u)  // if not == ref
                appendValue(variants.svRecords, svRecord);
        }
    }
    else if (length(vcfRecord.ref) == 1u && altLength == 1u)  // SNP
    {
        SnpRecord snpRecord;
        snpRecord.rId = vcfRecord.rID;
        snpRecord.pos = vcfRecord.beginPos;

        // Split the alternatives.
        seqan::StringSet<seqan::CharString> alternatives;
        splitString(alternatives, vcfRecord.alt, ',');
            
        // Split the target variants.
        SEQAN_ASSERT_NOT(empty(vcfRecord.genotypeInfos));
        seqan::RecordReader<seqan::CharString const, seqan::SinglePass<seqan::StringReader> >
                reader(vcfRecord.genotypeInfos[0]);
        seqan::CharString buffer;
        snpRecord.haplotype = 0;
        for (; !atEnd(reader); goNext(reader))
            if ((value(reader) == '|' || value(reader) == '/'))
            {
                if (!empty(buffer))
                {
                    unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer),
                                            (unsigned)length(alternatives));
                    if (idx != 0u)  // if not == ref
                    {
                        SEQAN_ASSERT_NOT(empty(alternatives[idx - 1]));
                        snpRecord.to = alternatives[idx - 1][0];
                        appendValue(variants.snps, snpRecord);
                    }
                }
                snpRecord.haplotype++;
                clear(buffer);
            }
            else
            {
                appendValue(buffer, value(reader));
            }
        if (!empty(buffer))
        {
            unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer),
                                    (unsigned)length(alternatives));
            if (idx != 0u)  // if not == ref
            {
                SEQAN_ASSERT_NOT(empty(alternatives[idx - 1]));
                snpRecord.to = alternatives[idx - 1][0];
                appendValue(variants.snps, snpRecord);
            }
        }
    }
    else  // Small Indel
    {
        SmallIndelRecord smallIndel;
        smallIndel.rId = vcfRecord.rID;
        smallIndel.pos = vcfRecord.beginPos;

        SEQAN_ASSERT_NOT(contains(vcfRecord.alt, ","));  // only one alternative
        SEQAN_ASSERT((length(vcfRecord.alt) == 1u) != (length(vcfRecord.ref) == 1u));  // XOR

        smallIndel.haplotype = 0;
        if (length(vcfRecord.ref) == 1u)  // insertion
        {
            smallIndel.seq = suffix(vcfRecord.alt, 1);
            smallIndel.size = length(smallIndel.seq);
        }
        else  // deletion
        {
            smallIndel.size = -(int)(length(vcfRecord.ref) - 1);
        }
            
        // Split the target variants.
        SEQAN_ASSERT_NOT(empty(vcfRecord.genotypeInfos));
        seqan::RecordReader<seqan::CharString const, seqan::SinglePass<seqan::StringReader> >
                reader(vcfRecord.genotypeInfos[0]);
        seqan::CharString buffer;
        smallIndel.haplotype = 0;
        for (; !atEnd(reader); goNext(reader))
            if ((value(reader) == '|' || value(reader) == '/'))
            {
                if (!empty(buffer))
                {
                    unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer), 1u);
                    if (idx != 0u)  // if not == ref
                        appendValue(variants.smallIndels, smallIndel);
                }
                smallIndel.haplotype++;
                clear(buffer);
            }
            else
            {
                appendValue(buffer, value(reader));
            }
        if (!empty(buffer))
        {
            unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer), 1u);
            if (idx != 0u)  // if not == ref
                appendValue(variants.smallIndels, smallIndel);
        }
    }
}

// ----------------------------------------------------------------------------
// Function VcfMaterializer::_appendToVariantsBnd()
// ----------------------------------------------------------------------------

void VcfMaterializer::_appendToVariantsBnd(Variants & variants, std::vector<seqan::VcfRecord> const & vcfRecords)
{
    // Sanity check.
    SEQAN_CHECK(length(vcfRecords) == 6u, "Must be 6-tuple.");
    SEQAN_CHECK(vcfRecords[0].beginPos + 1 == vcfRecords[1].beginPos, "Invalid BND 6-tuple.");
    SEQAN_CHECK(vcfRecords[2].beginPos + 1 == vcfRecords[3].beginPos, "Invalid BND 6-tuple.");
    SEQAN_CHECK(vcfRecords[4].beginPos + 1 == vcfRecords[5].beginPos, "Invalid BND 6-tuple.");
    SEQAN_CHECK(vcfRecords[1].beginPos < vcfRecords[2].beginPos, "Wrong ordering.");
    SEQAN_CHECK(vcfRecords[3].beginPos < vcfRecords[4].beginPos, "Wrong ordering.");

    // Add translocation.
    StructuralVariantRecord svRecord;
    svRecord.kind = StructuralVariantRecord::TRANSLOCATION;
    svRecord.rId = vcfRecords[0].rID;
    svRecord.pos = vcfRecords[1].beginPos;
    svRecord.targetRId = svRecord.rId;
    svRecord.targetPos = vcfRecords[5].beginPos;
    svRecord.size = vcfRecords[3].beginPos - vcfRecords[1].beginPos;
    svRecord.haplotype = 0;

    // Split the target variants.
    SEQAN_ASSERT_NOT(empty(vcfRecord.genotypeInfos));
    seqan::RecordReader<seqan::CharString const, seqan::SinglePass<seqan::StringReader> >
            reader(vcfRecord.genotypeInfos[0]);
    seqan::CharString buffer;
    svRecord.haplotype = 0;
    for (; !atEnd(reader); goNext(reader))
        if ((value(reader) == '|' || value(reader) == '/'))
        {
            if (!empty(buffer))
            {
                unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer), 1u);
                if (idx != 0u)  // if not == ref
                    appendValue(variants.svRecords, svRecord);
            }
            svRecord.haplotype++;
            clear(buffer);
        }
        else
        {
            appendValue(buffer, value(reader));
        }
    if (!empty(buffer))
    {
        unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer), 1u);
        if (idx != 0u)  // if not == ref
            appendValue(variants.svRecords, svRecord);
    }
}

#endif  // #ifndef SANDBOX_MASON2_APPS_MASON2_VCF_MATERIALIZATION_H_
