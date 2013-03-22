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
// Temporary fstream generation.
// ==========================================================================

#ifndef SANDBOX_MASON2_APPS_MASON2_FSTREAM_TEMP_H_
#define SANDBOX_MASON2_APPS_MASON2_FSTREAM_TEMP_H_

#include <fstream>
#include <stdlib.h>

#include <seqan/sequence.h>

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

void openTemp(std::fstream & stream)
{
    seqan::CharString tmpFilename;

    // Try to get the temporary directory from the environment variables TMPDIR or TMP.
    if ((getuid() == geteuid()) && (getgid() == getegid()))
    {
        char * res;
        if ((res = getenv("TMPDIR")) != NULL)
            tmpFilename = res;
        else if ((res = getenv("TMP")) != NULL)
            tmpFilename = res;        
    }
    
    // Fall back to "/tmp" if this fails.
    if (empty(tmpFilename))
        tmpFilename = "/tmp";

    // Create template for the temporary file.
    append(tmpFilename, "/SQNXXXXXX");
    
    // Open temporary file and unlink it immediately afterwards so the memory is released when the program exits.
    int oldMode = umask(077);  // Create with restrictive permissions.
    int fd = -1;  // File descriptor.
    if ((fd = mkstemp(toCString(tmpFilename))) == -1)
    {
        umask(oldMode);  // Reset umask mode.
        throw std::string("Could not open temporary file!");
    }

    // Open the file with stream.
    stream.open(toCString(tmpFilename), std::ios::binary | std::ios::in | std::ios::out | std::ios::trunc);
    if (!stream.good())
        throw std::string("Could not re-open temporary file!");

    // Close and delete the file, we have it open in stream now, and after stream is destructed, the file will be removed.
    close(fd);
    // if (unlink(toCString(tmpFilename)) != 0)
    //     throw std::string("Could not remove temporary file!");
}

#endif  // #ifndef SANDBOX_MASON2_APPS_MASON2_FSTREAM_TEMP_H_
