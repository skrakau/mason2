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
// Compute transcripts from a genome and a GFF/GTF file.  Optionally, you
// can apply a VCF file to the genome before splicing.
//
// Transcripts must not span structural variants.
// ==========================================================================

// Note: We treat all given variants as phased.

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/vcf_io.h>

#include "vcf_materialization.h"
#include "mason_options.h"
#include "mason_types.h"

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class MasonSplicingApp
// --------------------------------------------------------------------------

class MasonSplicingApp
{
public:
    // The configuration to use.
    MasonSplicingOptions const & options;

    // The random number generation.
    TRng rng;

    // Materialization of VCF.
    VcfMaterializer vcfMat;

    // Output sequence stream.
    seqan::SequenceStream outStream;

    MasonSplicingApp(MasonSplicingOptions const & _options) :
            options(_options), rng(options.seed),
            vcfMat(rng, toCString(options.matOptions.fastaFileName), toCString(options.matOptions.vcfFileName))
    {}

    int run()
    {
        // Intialization
        std::cerr << "__INITIALIZATION_____________________________________________________________\n"
                  << "\n";

        std::cerr << "Opening files...";
        try
        {
            vcfMat.init();

            open(outStream, toCString(options.outputFileName), seqan::SequenceStream::WRITE);
            if (!isGood(outStream))
                throw MasonIOException("Could not open output file.");
        }
        catch (MasonIOException e)
        {
            std::cerr << "\nERROR: " << e.what() << "\n";
            return 1;
        }
        std::cerr << " OK\n";

        // Perform genome simulation.
        std::cerr << "\n__COMPUTING TRANSCRIPTS______________________________________________________\n"
                  << "\n";

        // The identifiers of the just materialized data.
        int rID = 0, hID = 0;
        seqan::Dna5String seq;
        std::cerr << "Splicing";
        while (vcfMat.materializeNext(seq, rID, hID))
        {
            std::stringstream ssName;
            std::cerr << "  " << sequenceName(vcfMat.faiIndex, rID) << " (allele " << (hID + 1) << ") ";
            
            if (writeRecord(outStream, ssName.str(), seq) != 0)
            {
                std::cerr << "ERROR: Could not write materialized sequence to output.\n";
                return 1;
            }
        }
        std::cerr << " DONE\n";

        std::cerr << "\nDone splicing FASTA.\n";

        return 0;
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(MasonSplicingOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("mason_splicing");
    // Set short description, version, and date.
    setShortDescription(parser, "Generating Transcripts");
    setVersion(parser, "2.0");
    setDate(parser, "July 2012");
    setCategory(parser, "Simulators");

    // Define usage line and long description.
    addUsageLine(parser,
                 "[OPTIONS] \\fB-ir\\fP \\fIIN.fa\\fP \\fB-ig\\fP \\fIIN.gff\\fP [\\fB-iv\\fP \\fIIN.vcf\\fP] \\fB-o\\fP \\fIOUT.fa\\fP");
    addDescription(parser,
                   "Create transcripts from \\fIIN.fa\\fP using the annotations from \\fIIN.gff\\fP.  The resulting "
                   "transcripts are written to \\fIOUT.fa\\fP.");
    addDescription(parser,
                   "You can pass an optional VCF file \\fIIN.vcf\\fP and the transcripts will be created from the "
                   "haplotypes stored in the VCF file.");

    // Add option and text sections.
    options.addOptions(parser);
    options.addTextSections(parser);
    
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    options.getOptionValues(parser);

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    MasonSplicingOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cerr << "MASON SPLICING\n"
              << "==============\n\n";
    
    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cerr << "__OPTIONS____________________________________________________________________\n"
                  << "\n";
        options.print(std::cerr);
    }

    MasonSplicingApp app(options);
    return app.run();
}
