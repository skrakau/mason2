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
// Management of fragment/read-to-content distribution.
//
// Simulating contig-wise saves memory to a comfortable amount and is easy.
// We first simulate for the fragments 0..(n-1) from which contig they come
// from write their ids to one file per contig.  In a second step, we read
// contig-wise through the id files and simulate the fragments/reads.  In a
// final step, we merge the fragments/reads by their id and write them to
// the final file.
//
// This header provides the data structures and routines to manage this.
// ==========================================================================

#ifndef SANDBOX_MASON2_APPS_MASON2_EXTERNAL_SPLIT_MERGE_H_
#define SANDBOX_MASON2_APPS_MASON2_EXTERNAL_SPLIT_MERGE_H_

#include <vector>
#include <iostream>

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

#include "mason_types.h"

// ============================================================================
// Forwards
// ============================================================================

inline bool ltBamAlignmentRecord(seqan::BamAlignmentRecord const & lhs,
                                 seqan::BamAlignmentRecord const & rhs);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class IdSplitter
// ----------------------------------------------------------------------------

// Allows distributing ids from/to files.
//
// General protocol:
//
// * construct
// * open()
// * write ids, splitting
// * reset()
// * read ids, contig-wise
// * close()

// TODO(holtgrew): Name bogus, FileBundle would be better.

class IdSplitter
{
public:
    // The number of contigs to split to.
    unsigned numContigs;

    // The file pointers for each contig.
    std::vector<FILE *> files;

    IdSplitter() : numContigs(0)
    {}

    IdSplitter(unsigned numContigs) : numContigs(numContigs)
    {}

    ~IdSplitter()
    {
        close();
    }

    // Open files in the splitter.
    void open();

    // Reset all files in the splitter, ready for reading.
    void reset();

    // Close splitter.
    void close();
};

// ----------------------------------------------------------------------------
// Class FastxJoiner
// ----------------------------------------------------------------------------

// Allows joining by id name from FASTA data stored in a IdSplitter.
//
// Construct with IdSplitter after reset() call.

// TODO(holtgrew): Could use a heap/tournament tree.

template <typename TTag>
class FastxJoiner
{
public:
    // The type of the record reader to use.
    typedef seqan::RecordReader<FILE *, seqan::SinglePass<> > TReader;

    // The IdSplitter to use.
    IdSplitter * splitter;
    // Number of active files.
    unsigned numActive;
    // Buffer for id and sequence for each input file.
    seqan::StringSet<seqan::CharString> ids, seqs, quals;
    // Maps files for activeness.
    std::vector<bool> active;
    // Record reads, one for each input file.
    std::vector<TReader *> readers;

    FastxJoiner() : splitter(), numActive(0)
    {}

    FastxJoiner(IdSplitter & splitter) : splitter(&splitter), numActive(0)
    {
        _init();
    }

    ~FastxJoiner()
    {
        for (unsigned i = 0; i < readers.size(); ++i)
            delete readers[i];
        readers.clear();
    }

    void _init()
    {
        resize(ids, splitter->files.size());
        resize(seqs, splitter->files.size());
        resize(quals, splitter->files.size());
        active.resize(splitter->files.size());

        for (unsigned i = 0; i < splitter->files.size(); ++i)
        {
            readers.push_back(new TReader(splitter->files[i]));
            active[i] = _loadNext(ids[i], seqs[i], quals[i], i);
            numActive += (active[i] != false);
        }
    }

    template <typename TSeq>
    bool _loadNext(TSeq & id, TSeq & seq, TSeq & qual, unsigned idx)
    {
        if (seqan::atEnd(*readers[idx]))
            return false;
        if (readRecord(id, seq, qual, *readers[idx], TTag()) != 0)
        {
            std::cerr << "ERROR: Problem reading temporary data.\n";
            exit(1);
        }
        return true;
    }

    bool atEnd() const
    {
        return (numActive == 0);
    }

    int get(seqan::CharString & id, seqan::CharString & seq, seqan::CharString & qual)
    {
        unsigned idx = seqan::maxValue<unsigned>();
        for (unsigned i = 0; i < length(ids); ++i)
        {
            if (!active[i])
                continue;
            if (idx == seqan::maxValue<unsigned>() || ids[i] < ids[idx])
                idx = i;
        }
        if (idx == seqan::maxValue<unsigned>())
            return 1;

        // We use double-buffering and the input parameters as buffers.
        active[idx] = _loadNext(id, seq, qual, idx);
        swap(id, ids[idx]);
        swap(seq, seqs[idx]);
        swap(qual, quals[idx]);
        numActive -= !active[idx];

        return 0;
    }
};

// ----------------------------------------------------------------------------
// Class SamJoiner
// ----------------------------------------------------------------------------

// Allows joining by id name from FASTA data stored in a IdSplitter.
//
// Construct with IdSplitter after reset() call.

// TODO(holtgrew): Could use a heap/tournament tree.

// Compare two BAM alignment records by query name, tie is broken by first/last flag, first < last.

class SamJoiner
{
public:
    // The type of the record reader to use.
    typedef seqan::RecordReader<FILE *, seqan::SinglePass<> > TReader;

    // The IdSplitter to use.
    IdSplitter * splitter;
    // Number of active files.
    unsigned numActive;
    // Buffer for id and sequence for each input file.
    seqan::String<seqan::BamAlignmentRecord> records;
    // Maps files for activeness.
    std::vector<bool> active;
    // Record reads, one for each input file.
    std::vector<TReader *> readers;

    // One of the identical BAM headers.
    seqan::BamHeader header;
    // The reference sequence name store and cache.
    seqan::StringSet<seqan::CharString> nameStore;
    seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > nameStoreCache;
    seqan::BamIOContext<seqan::StringSet<seqan::CharString> > context;

    SamJoiner() : splitter(), numActive(0), nameStoreCache(nameStore), context(nameStore, nameStoreCache)
    {}

    SamJoiner(IdSplitter & splitter) :
            splitter(&splitter), numActive(0), nameStoreCache(nameStore), context(nameStore, nameStoreCache)
    {
        init();
    }

    ~SamJoiner()
    {
        for (unsigned i = 0; i < readers.size(); ++i)
            delete readers[i];
        readers.clear();
    }

    void init();

    bool _loadNext(seqan::BamAlignmentRecord & record, unsigned idx);

    bool atEnd() const
    {
        return (numActive == 0);
    }

    // Get next BAM alignment record to lhs.  If it is paired-end, load the second mate as well.
    int get(seqan::BamAlignmentRecord & record);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function ltBamAlignmentRecord()
// ----------------------------------------------------------------------------

inline bool ltBamAlignmentRecord(seqan::BamAlignmentRecord const & lhs,
                                 seqan::BamAlignmentRecord const & rhs)
{
    seqan::Lexical<> cmp(lhs.qName, rhs.qName);
    if (isLess(cmp) || (isEqual(cmp) && hasFlagFirst(lhs)))
        return true;
    return false;
}

#endif  // #ifndef SANDBOX_MASON2_APPS_MASON2_EXTERNAL_SPLIT_MERGE_H_
