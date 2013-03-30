// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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

#ifndef SANDBOX_HERBARIUM_INCLUDE_SEQAN_VCF_WRITE_VCF_H_
#define SANDBOX_HERBARIUM_INCLUDE_SEQAN_VCF_WRITE_VCF_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function writeRecord()                                           [VcfHeader]
// ----------------------------------------------------------------------------

template <typename TStream>
int writeRecord(TStream & stream,
                VcfHeader const & header,
                VcfIOContext const & /*vcfIOContext*/,
                Vcf const & /*tag*/)
{
    for (unsigned i = 0; i < length(header.headerRecords); ++i)
    {
        streamWriteBlock(stream, "##", 2);
        streamWriteBlock(stream, &header.headerRecords[i].key[0], length(header.headerRecords[i].key));
        streamWriteChar(stream, '=');
        streamWriteBlock(stream, &header.headerRecords[i].value[0], length(header.headerRecords[i].value));
        streamWriteChar(stream, '\n');
    }

    streamWriteBlock(stream, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
                     length("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"));
    for (unsigned i = 0; i < length(header.sampleNames); ++i)
    {
        streamWriteChar(stream, '\t');
        streamWriteBlock(stream, &header.sampleNames[i][0], length(header.sampleNames[i]));
    }
    streamWriteChar(stream, '\n');

    return 0;
}

// ----------------------------------------------------------------------------
// Function writeRecord()                                           [VcfRecord]
// ----------------------------------------------------------------------------

template <typename TStream>
int writeRecord(TStream & stream,
                VcfRecord const & record,
                VcfIOContext const & vcfIOContext,
                Vcf const & /*tag*/)
{
    streamWriteBlock(stream, &(*vcfIOContext.sequenceNames)[record.chromId][0],
                     length((*vcfIOContext.sequenceNames)[record.chromId]));
    streamWriteChar(stream, '\t');
    streamPut(stream, record.pos + 1);
    streamWriteChar(stream, '\t');
    if (empty(record.id))
        streamWriteBlock(stream, ".", 1);
    else
        streamWriteBlock(stream, &record.id[0], length(record.id));
    streamWriteChar(stream, '\t');
    if (empty(record.ref))
        streamWriteBlock(stream, ".", 1);
    else
        streamWriteBlock(stream, &record.ref[0], length(record.ref));
    streamWriteChar(stream, '\t');
    if (empty(record.alt))
        streamWriteBlock(stream, ".", 1);
    else
        streamWriteBlock(stream, &record.alt[0], length(record.alt));
    streamWriteChar(stream, '\t');
    if (record.qual != record.qual)  // only way to test for nan
        streamWriteChar(stream, '.');
    else
        streamPut(stream, record.qual);
    streamWriteChar(stream, '\t');
    if (empty(record.filter))
        streamWriteBlock(stream, ".", 1);
    else
        streamWriteBlock(stream, &record.filter[0], length(record.filter));
    streamWriteChar(stream, '\t');
    if (empty(record.info))
        streamWriteBlock(stream, ".", 1);
    else
        streamWriteBlock(stream, &record.info[0], length(record.info));
        streamWriteChar(stream, '\t');
    if (empty(record.format))
        streamWriteBlock(stream, ".", 1);
    else
        streamWriteBlock(stream, &record.format[0], length(record.format));
    for (unsigned i = 0; i < length(record.genotypeInfos); ++i)
    {
        streamWriteChar(stream, '\t');
        if (empty(record.genotypeInfos[i]))
            streamWriteBlock(stream, ".", 1);
        else
            streamWriteBlock(stream, &record.genotypeInfos[i][0], length(record.genotypeInfos[i]));
    }
    streamWriteChar(stream, '\n');

    return 0;
}

}  // namespace seqan

#endif  // #ifndef SANDBOX_HERBARIUM_INCLUDE_SEQAN_VCF_WRITE_VCF_H_
