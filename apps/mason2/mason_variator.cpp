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
// Given a genome, create variations thereof and export them VCF and
// materialized into a FASTA file.  Variations can also be imported from a VCF
// file and materialized into a FASTA file.
// ==========================================================================

#include  <seqan/arg_parse.h>
#include  <seqan/sequence.h>
#include  <seqan/seq_io.h>
#include  <seqan/vcf.h>

// ==========================================================================
// Forwards
// ==========================================================================

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class MasonVariatorOptions
// --------------------------------------------------------------------------

struct MasonVariatorOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // ----------------------------------------------------------------------
    // Input / Output Options
    // ----------------------------------------------------------------------

    // VCF file to import.
    seqan::CharString vcfInFile;
    // VCF file to write out.
    seqan::CharString vcfOutFile;
    // FASTA file to write out with variations.
    seqan::CharString fastaOutFile;

    // Path to a TSV file where the first two columns giving the type of the SV to simulate and the size of the SV.
    // This overrides the simulation of SV from the sv*Rate parameters.
    seqan::CharString inputSVSizeFile;

    // ----------------------------------------------------------------------
    // Haplotype / Allele Configuration
    // ----------------------------------------------------------------------

    // The number of haplotypes to simulate.
    int numHaplotypes;

    // The string to use to separate the haplotype identifier from the chromosome name in the output FASTA ifle.
    seqan::CharString haplotypeSep;

    // ----------------------------------------------------------------------
    // Variation Simulation
    // ----------------------------------------------------------------------

    // Per-base probability for SNPs and small-scale indels.
    double snpRate;
    double smallIndelRate;

    // Minimal and maximal size for small indels.  Indels will be simulated uniformly in this range.  The range is
    // stored internally as [min, max) but given as [min, max] from the command line.
    int minSmallIndelSize;
    int maxSmallIndelSize;
    
    // Per-base probability for having a structural variation.
    double svIndelRate;
    double svInversionRate;
    double svTranslocationRate;
    double svDuplicationRate;

    // Minimal and maximal size for structural variations.  SVs will be simulated uniformly in this range.  The range is
    // stored internally as [min, max) but given as [min, max] from the command line.
    int minSVSize;
    int maxSVSize;

    MasonVariatorOptions() :
            snpRate(0), smallIndelRate(0), minSmallIndelSize(0), maxSmallIndelSize(0), svIndelRate(0),
            svInversionRate(0), svTranslocationRate(0), svDuplicationRate(0), minSVSize(0), maxSVSize(0)
    {}
};

// ==========================================================================
// Metafunctions
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(MasonVariatorOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("mason_variator");
    // Set short description, version, and date.
    setShortDescription(parser, "Variation Simulation");
    setVersion(parser, "2.1");
    setDate(parser, "March 2013");
    setCategory(parser, "Simulators");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-ov\\fP \\fIOUT.vcf\\fP [\\fB-of\\fP \\fIOUT.fa\\fP]");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-iv\\fP \\fIIN.vcf\\fP \\fB-of\\fP \\fIOUT.fa\\fP");
    addDescription(parser,
                   "Either simulate variation and write out the result to VCF and FASTA files "
                   "or apply the variations from a VCF file and write the results to a FASTA file.");

    // ----------------------------------------------------------------------
    // General Options
    // ----------------------------------------------------------------------

    addSection(parser, "General Options");

    // We require one argument.
    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    addOption(parser, seqan::ArgParseOption("s", "seed", "The seed to use for the random number generator.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "seed", "0");

    // ----------------------------------------------------------------------
    // Input / Output Options
    // ----------------------------------------------------------------------

    addSection(parser, "Input / Output");
    
    addOption(parser, seqan::ArgParseOption("iv", "in-vcf", "VCF file to load variations from.",
                                            seqan::ArgParseOption::INPUTFILE, "VCF"));
    setValidValues(parser, "in-vcf", "vcf");

    addOption(parser, seqan::ArgParseOption("it", "in-variant-tsv",
                                            "TSV file with variants to simulate.  See Section on the Variant TSV File below.",
                                            seqan::ArgParseOption::INPUTFILE, "VCF"));
    setValidValues(parser, "in-vcf", "vcf");

    addOption(parser, seqan::ArgParseOption("ov", "out-vcf", "VCF file to write simulated variations to.",
                                            seqan::ArgParseOption::INPUTFILE, "VCF"));
    setValidValues(parser, "out-vcf", "vcf");

    addOption(parser, seqan::ArgParseOption("of", "out-fasta", "FASTA file to write simulated haplotypes to.",
                                            seqan::ArgParseOption::INPUTFILE, "FASTA"));
    setValidValues(parser, "out-fasta", "fasta");

    // ----------------------------------------------------------------------
    // Haplotype / Allele Configuration
    // ----------------------------------------------------------------------

    addSection(parser, "Haplotype / Allele Configuration");

    addOption(parser, seqan::ArgParseOption("n", "num-haplotypes", "The number of haplotypes to simulate.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setMinValue(parser, "num-haplotypes", "1");
    setDefaultValue(parser, "num-haplotypes", "1");

    addOption(parser, seqan::ArgParseOption("", "haplotype-sep",
                                            "The separator between the chromosome and the haplotype name "
                                            "in the output FASTA file.",
                                            seqan::ArgParseOption::STRING, "SEP"));
    setDefaultValue(parser, "haplotype-sep", ":");

    // ----------------------------------------------------------------------
    // Variation Simulation Options
    // ----------------------------------------------------------------------

    addSection(parser, "Variation Simulation");

    addOption(parser, seqan::ArgParseOption("", "snp-rate", "Per-base SNP rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "snp-rate", "0.0");
    setMaxValue(parser, "snp-rate", "1.0");

    addOption(parser, seqan::ArgParseOption("", "small-indel-rate", "Small indel rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "small-indel-rate", "0.0");
    setMaxValue(parser, "small-indel-rate", "1.0");

    addOption(parser, seqan::ArgParseOption("", "min-small-indel-size", "Minimal small indel size to simulate.",
                                            seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "min-small-indel-size", "0");
    setDefaultValue(parser, "min-small-indel-size", "50");

    addOption(parser, seqan::ArgParseOption("", "max-small-indel-size", "Maximal small indel size to simulate.",
                                            seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "max-small-indel-size", "0");
    setDefaultValue(parser, "max-small-indel-size", "1000");

    addOption(parser, seqan::ArgParseOption("", "sv-indel-rate", "Per-base SNP rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "sv-indel-rate", "0.0");
    setMaxValue(parser, "sv-indel-rate", "1.0");

    addOption(parser, seqan::ArgParseOption("", "sv-inversion-rate", "Per-base SNP rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "sv-inversion-rate", "0.0");
    setMaxValue(parser, "sv-inversion-rate", "1.0");

    addOption(parser, seqan::ArgParseOption("", "sv-translocation-rate", "Per-base SNP rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "sv-translocation-rate", "0.0");
    setMaxValue(parser, "sv-translocation-rate", "1.0");

    addOption(parser, seqan::ArgParseOption("", "sv-duplication-rate", "Per-base SNP rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "sv-duplication-rate", "0.0");
    setMaxValue(parser, "sv-duplication-rate", "1.0");

    addOption(parser, seqan::ArgParseOption("", "min-sv-size", "Minimal SV size to simulate.",
                                            seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "min-sv-size", "0");
    setDefaultValue(parser, "min-sv-size", "50");

    addOption(parser, seqan::ArgParseOption("", "max-sv-size", "Maximal SV size to simulate.",
                                            seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "max-sv-size", "0");
    setDefaultValue(parser, "max-sv-size", "1000");

    // ----------------------------------------------------------------------
    // Simulation Details Section
    // ----------------------------------------------------------------------

    addTextSection(parser, "Simulation Details");
    addText(parser, "The indel and SV sizes are picked uniformly at random from the size intervals.");

    // ----------------------------------------------------------------------
    // Examples Section
    // ----------------------------------------------------------------------

    // TODO(holtgrew): Write me!

    // ----------------------------------------------------------------------
    // Variation TSV File
    // ----------------------------------------------------------------------

    addTextSection(parser, "Variation TSV File");
    addText(parser,
            "Instead of simulating the variants from per-base rates, the user can specify a TSV (tab separated values) "
            "file to load the variations from with \\fB--in-variant-tsv\\fP/\\fB-it\\fP.  The first two columns of "
            "this TSV file are interpreted as the type of the variation and the size");
    addText(parser, "The following types are defined");
    addListItem(parser, "SNP", "A SNP, the second column is ignored");
    addListItem(parser, "INS", "An insertion.");
    addListItem(parser, "DEL", "A deletion.");
    addListItem(parser, "INV", "An inversion.");
    addListItem(parser, "TRA", "An inter-chromosomal translocation.");
    addListItem(parser, "CTR", "An intra-chromosomal translocation.");
    addListItem(parser, "DUP", "A duplication");

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

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    MasonVariatorOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cout << "MASON VARIATOR\n"
              << "==============\n";

    return 0;
}
