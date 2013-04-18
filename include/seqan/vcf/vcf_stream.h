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

// TODO(holtgrew): Implement writing.

#ifndef SANDBOX_HOLTGREW_INCLUDE_SEQAN_VCF_VCF_STREAM_H_
#define SANDBOX_HOLTGREW_INCLUDE_SEQAN_VCF_VCF_STREAM_H_

#include <memory>
#include <fstream>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

class VcfStream
{
public:
    typedef RecordReader<std::fstream, SinglePass<> > TReader_;

    enum Mode
    {
        INVALID,
        READ,
        WRITE
    };

    std::auto_ptr<std::fstream> _stream;
    CharString _filename;
    std::auto_ptr<TReader_> _reader;
    Mode _mode;
    int _error;
    bool _isGood;
    bool _headerWritten;

    VcfHeader header;
    VcfIOContext _context;

    VcfStream() : _mode(INVALID), _error(0), _isGood(true), _headerWritten(false),
                  _context(header.sequenceNames, header.sampleNames)
    {}

    VcfStream(char const * filename, Mode mode = READ) :
            _filename(filename), _mode(mode), _error(0), _isGood(true), _headerWritten(false),
            _context(header.sequenceNames, header.sampleNames)
    {
        _open(filename, mode);
    }

    bool _open(char const * filename, Mode mode)
    {
        // Reset.
        _filename = filename;
        _mode = mode;
        _error = 0;
        _isGood = true;
        _headerWritten = false;

        if (mode == READ)
        {
            _stream.reset(new std::fstream);
            _stream->open(toCString(_filename), std::ios::binary | std::ios::in);
            if (!_stream->good())
            {
                _isGood = false;
                return false;
            }
            _reader.reset(new TReader_(*_stream));

            int res = read(header, *_reader, _context, Vcf());
            if (res != 0)
            {
                _error = res;
                _isGood = false;
            }
        }
        else if (mode == WRITE)
        {
            _stream.reset(new std::fstream);
            _stream->open(toCString(_filename), std::ios::binary | std::ios::out);
            if (!_stream->good())
            {
                _isGood = false;
                return false;
            }
            _reader.reset();
        }
        return true;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

inline bool open(VcfStream & stream, char const * filename, VcfStream::Mode mode = VcfStream::READ)
{
    return stream._open(filename, mode);
}

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

inline int readRecord(VcfRecord & record,
                      VcfStream & stream)
{
    int res = readRecord(record, *stream._reader, stream._context, Vcf());
    if (res != 0)
        stream._isGood = false;
    return res;
}

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

inline int writeRecord(VcfStream & stream,
                       VcfRecord const & record)
{
    if (!stream._headerWritten)
    {
        int res = writeRecord(*stream._stream, stream.header, stream._context, Vcf());
        if (res != 0)
            stream._isGood = false;
        stream._headerWritten = true;
    }

    int res = writeRecord(*stream._stream, record, stream._context, Vcf());
    if (res != 0)
        stream._isGood = false;
    return res;
}

// ----------------------------------------------------------------------------
// Function flush()
// ----------------------------------------------------------------------------

inline int flush(VcfStream & stream)
{
    stream._stream->flush();
    return 0;
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

inline int close(VcfStream & stream)
{
    stream._stream->close();
    return 0;
}

// ----------------------------------------------------------------------------
// Function isGood()
// ----------------------------------------------------------------------------

inline bool isGood(VcfStream const & stream)
{
    return stream._isGood;
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

inline bool atEnd(VcfStream const & stream)
{
    return atEnd(*stream._reader);
}

inline bool atEnd(VcfStream & stream)
{
    return atEnd(*stream._reader);
}

}  // namespace seqan

#endif  // #ifndef SANDBOX_HOLTGREW_INCLUDE_SEQAN_VCF_VCF_STREAM_H_
