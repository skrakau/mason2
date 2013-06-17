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
// Apply variants from a VCF file to a genomic sequence.
//
// The variants must be equivalent to the variants written by mason_variator.
// See the documentation of mason_materializer and mason_variator for details
// on this.
// ==========================================================================

// Note: We treat all given variants as phased.

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/vcf_io.h>

#include "vcf_materialization.h"
#include "mason_types.h"

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class MasonMaterializerOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct MasonMaterializerOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The path to the input FASTA file.
    seqan::CharString inFastaFile;
    // The path to the input VCF file.
    seqan::CharString inVcfFile;
    // The output file name.
    seqan::CharString outputFilename;

    // Separator for haplotype numbers.
    seqan::CharString haplotypeNameSep;

    // The seed to use for the RNG.
    int seed;

    MasonMaterializerOptions() : verbosity(1), seed(0)
    {}
};

// --------------------------------------------------------------------------
// Class MasonMaterializerApp
// --------------------------------------------------------------------------

class MasonMaterializerApp
{
public:
    // The configuration to use.
    MasonMaterializerOptions const & options;

    // Materialization of VCF.
    VcfMaterializer vcfMat;
    
    // Output sequence stream.
    seqan::SequenceStream outStream;

    MasonMaterializerApp(MasonMaterializerOptions const & options) :
            options(options), vcfMat(toCString(options.inFastaFile), toCString(options.inVcfFile))
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
            open(outStream, toCString(options.outputFilename), seqan::SequenceStream::WRITE);
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
        std::cerr << "\n__MATERIALIZING______________________________________________________________\n"
                  << "\n";

        // The identifiers of the just materialized data.
        int rID = 0, hID = 0;
        seqan::Dna5String seq;
        std::cerr << "Materializing...";
        while (vcfMat.materializeNext(seq, rID, hID))
        {
            std::stringstream ssName;
            ssName << vcfMat.vcfStream.header.sequenceNames[rID] << options.haplotypeNameSep << hID;
            std::cerr << " " << ssName.str();
            
            if (writeRecord(outStream, ssName.str(), seq) != 0)
            {
                std::cerr << "ERROR: Could not write materialized sequence to output.\n";
                return 1;
            }
        }
        std::cerr << " DONE\n";

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
parseCommandLine(MasonMaterializerOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("mason_materializers");
    // Set short description, version, and date.
    setShortDescription(parser, "Materialize VCF into FASTA");
    setVersion(parser, "2.1");
    setDate(parser, "June 2013");
    setCategory(parser, "Simulators");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-iv\\fP \\fIIN.vcf\\fP \\fB-if\\fP \\fIIN.fa\\fP "
                 "\\fB-of\\fP \\fIOUT.fa\\fP");
    addDescription(parser,
                   "Apply the variants from \\fIIN.vcf\\fP to the reference in \\fIIN.fa\\fP.  The resulting "
                   "sequence is written to \\fIOUT.fa\\fP");

    // General Options
    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    addOption(parser, seqan::ArgParseOption("", "haplotype-name-sep", "Separator between contig and haplotype name.",
                                            seqan::ArgParseOption::STRING, "STR"));
    setDefaultValue(parser, "haplotype-name-sep", "/");
    
    // Input Output Options
    addSection(parser, "Input / Output");
    addOption(parser, seqan::ArgParseOption("if", "in-fasta", "Reference input file",
                                            seqan::ArgParseOption::INPUTFILE, "IN.fa"));
    setValidValues(parser, "in-fasta", "fa fasta");
    setRequired(parser, "in-fasta");

    addOption(parser, seqan::ArgParseOption("iv", "in-vcf", "Variants input file",
                                            seqan::ArgParseOption::INPUTFILE, "IN.vcf"));
    setValidValues(parser, "in-vcf", "vcf");
    setRequired(parser, "in-vcf");

    addOption(parser, seqan::ArgParseOption("of", "out-fasta", "Sequence output file",
                                            seqan::ArgParseOption::OUTPUTFILE, "OUT.fa"));
    setValidValues(parser, "out-fasta", "fa fasta");
    setRequired(parser, "out-fasta");

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBmason_materializer\\fP \\fB-iv\\fP \\fIone.vcf\\fP \\fB-if\\fP \\fIhg19.fa\\fP "
                "\\fB-of\\fP \\fIone.fa\\fP",
                "Apply the variants in \\fIone.vcf\\fP to the genome in \\fIhg19.fa\\fP and write the results to "
                "\\fIone.fa\\fP");

    // Add Text Section.
    addTextSection(parser, "Notes");
    addText(parser,
            "All haplotypes of the first individual in the VCF file will be materialized.  All others "
            "will be ignored.");
    
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValue(options.inFastaFile, parser, "in-fasta");
    getOptionValue(options.inVcfFile, parser, "in-vcf");
    getOptionValue(options.outputFilename, parser, "out-fasta");

    getOptionValue(options.haplotypeNameSep, parser, "haplotype-name-sep");

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    MasonMaterializerOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cout << "MASON VARIANT MATERIALIZER\n"
              << "==========================\n\n";
    
    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY    \t" << options.verbosity << '\n'
                  << "\n"
                  << "HAPLOTYPE SEP\t" << options.haplotypeNameSep << "\n"
                  << "\n"
                  << "INPUT FASTA  \t" << options.inFastaFile << "\n"
                  << "INPUT VCF    \t" << options.inVcfFile << "\n"
                  << "OUTPUT FILE  \t" << options.outputFilename << "\n"
                  << "\n\n";
    }

    MasonMaterializerApp app(options);
    return app.run();
}
