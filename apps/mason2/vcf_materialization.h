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

    VcfMaterializer() : currRID(-1), nextHaplotype(0), numHaplotypes(0)
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

#endif  // #ifndef SANDBOX_MASON2_APPS_MASON2_VCF_MATERIALIZATION_H_
