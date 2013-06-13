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

#include "genomic_variants.h"

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

bool contains(seqan::CharString const & haystack, char const * needle)
{
    return strstr(toCString(haystack), needle);
}

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

class MasonMaterializerApp
{
public:
    // The configuration to use.
    MasonMaterializerOptions const & options;

    // -----------------------------------------------------------------------
    // VCF loading.
    //
    // We will load from vcfStream into vcfRecord.
    // -----------------------------------------------------------------------
    // The index of the current contig.
    int rID;
    // The number of haplotypes in vcfStream.
    unsigned numHaplotypes;
    // The VCF stream to load from.
    seqan::VcfStream vcfStream;
    // The current VCF record.  rID == INVALID_REFID if invalid.
    seqan::VcfRecord vcfRecord;

    // FAI index for reading input sequence file.
    seqan::FaiIndex faiIndex;
    // Output sequence stream.
    seqan::SequenceStream outStream;

    MasonMaterializerApp(MasonMaterializerOptions const & options) :
            options(options), rID(-1), numHaplotypes(0)
    {}

    int run()
    {
        if (_init() != 0)
            return 1;

        // Perform genome simulation.
        std::cout << "__MATERIALIZING______________________________________________________________\n"
                  << "\n";

        // We do not use the methylation simulation in the materializer so the RNG will not be touched.
        TRng rng;
        
        // Load variants for each contig and materialize them.
        Variants variants;
        std::vector<int> breakpoints;  // unused
        seqan::Dna5String refSeq, seq;  // reference sequence and materialized sequence
        for (unsigned rID = 0; rID < length(vcfStream.header.sequenceNames); ++rID)
        {
            // Load reference sequence.
            unsigned idx = 0;
            if (!getIdByName(faiIndex, vcfStream.header.sequenceNames[rID], idx))
            {
                std::cerr << "ERROR: Sequence " << vcfStream.header.sequenceNames[rID]
                          << " from VCF not known in FAI.\n";
                return 1;
            }
            if (readSequence(refSeq, faiIndex, idx) != 0)
            {
                std::cerr << "ERROR: Could not load sequence from FAI.\n";
                return 1;
            }

            // Load variants.
            if (_loadVariantsForContig(variants, rID) != 0)
                return 1;

            // Materialize variants and write them out.
            for (unsigned hID = 0; hID < numHaplotypes; ++hID)
            {
                VariantMaterializer varMat(rng, variants);
                varMat.run(seq, breakpoints, refSeq, hID);

                std::stringstream ssName;
                ssName << vcfStream.header.sequenceNames[rID] << options.haplotypeNameSep
                       << hID;

                if (writeRecord(outStream, ssName.str(), seq) != 0)
                {
                    std::cerr << "ERROR: Could not write materialized sequence to output.\n";
                    return 1;
                }
            }
        }

        return 0;
    }

    int _init()
    {
        // Perform genome simulation.
        std::cout << "__OPENING FILES______________________________________________________________\n"
                  << "\n";

        // Open VCF stream.
        if (options.verbosity >= 1)
            std::cerr << "Opening VCF           \t" << options.inVcfFile << "...";
        open(vcfStream, toCString(options.inVcfFile));
        if (!isGood(vcfStream))
        {
            std::cerr << "ERROR: Could not open input VCF stream.\n";
            return 1;
        }
        if (options.verbosity >= 1)
            std::cerr << " OK\n";

        // Read first VCF record.
        if (!atEnd(vcfStream) && readRecord(vcfRecord, vcfStream) != 0)
        {
            std::cerr << "ERROR: Problem reading from VCF file.\n";
            return 1;
        }

        // Open input FASTA file and FAI.
        std::cerr << "Loading Reference Index " << options.inFastaFile << " ...";
        if (read(faiIndex, toCString(options.inFastaFile)) != 0)
        {
            std::cerr << " FAILED (not fatal, we can just build it)\n";
            std::cerr << "Building Index        " << options.inFastaFile << ".fai ...";
            if (build(faiIndex, toCString(options.inFastaFile)) != 0)
            {
                std::cerr << "Could not build FAI index.\n";
                return 1;
            }
            std::cerr << " OK\n";
            seqan::CharString faiPath = options.inFastaFile;
            append(faiPath, ".fai");
            std::cerr << "Reference Index       " << faiPath << " ...";
            if (write(faiIndex, toCString(faiPath)) != 0)
            {
                std::cerr << "Could not write FAI index we just built.\n";
                return 1;
            }
            std::cerr << " OK (" << length(faiIndex.indexEntryStore) << " seqs)\n";
        }
        else
        {
            std::cerr << " OK (" << length(faiIndex.indexEntryStore) << " seqs)\n";
        }

        // Open output FASTA file.
        if (options.verbosity >= 1)
            std::cerr << "Opening Output File     " << options.outputFilename << "...";
        open(outStream, toCString(options.outputFilename), seqan::SequenceStream::WRITE);
        if (!isGood(outStream))
        {
            std::cerr << "ERROR: Problem opening output FASTA file for writing.\n";
            return 1;
        }
        if (options.verbosity >= 1)
            std::cerr << " OK\n";

        return 0;
    }

    // Load variants of next contig into variants.
    int _loadVariantsForContig(Variants & variants, int rID)
    {
        variants.clear();

        // Compute number of haplotypes.
        SEQAN_ASSERT_NOT(empty(vcfRecord.genotypeInfos));
        seqan::StringSet<seqan::CharString> xs;
        splitString(xs, vcfRecord.genotypeInfos[0]);
        numHaplotypes = length(xs);

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
                    _appendToVariantsBnd(variants, chunk);
                chunk.clear();
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

    // Append VCF record to variants.
    void _appendToVariants(Variants & variants, seqan::VcfRecord const & vcfRecord)
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

    // Append chunk of 6 BND records to variants.
    void _appendToVariantsBnd(Variants & variants, std::vector<seqan::VcfRecord> const & vcfRecords)
    {
        // SEQAN_FAIL("Check chunk for being a valid translocation!");
        SEQAN_FAIL("Implement me!");
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
