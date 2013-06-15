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
// Simulate sequencing process from a genome.
// ==========================================================================

#include <vector>
#include <utility>

#include "mason_options.h"
#include "mason_types.h"
#include "vcf_materialization.h"
#include "external_split_merge.h"

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class ContigPicker
// --------------------------------------------------------------------------

// Distribute to contig and haplotypes.
//
// Contigs are picked with a probability proportional to their length and haplotypes are picked uniformly at random.

class ContigPicker
{
public:
    // The random number generator to use.
    TRng & rng;

    // The length of the contigs.
    std::vector<int> lengthSums;
    // The number of haplotypes.
    int numHaplotypes;
    
    ContigPicker(TRng & rng) : rng(rng)
    {}

    // Return a position (contig, haplotype) to distribute the read to.
    std::pair<int, int> pick()
    {
        // Pick reference id.
        int rID = 0;
        if (lengthSums.size() > 1u)
        {
            seqan::Pdf<seqan::Uniform<int> > pdf(0, lengthSums.back() - 1);
            int x = pickRandomNumber(rng, pdf);
            for (unsigned i = 0; i < lengthSums.size(); ++i)
            {
                if (x >= lengthSums[i])
                    rID = i + 1;
                if (x < lengthSums[i])
                    break;
            }
        }

        // Pick haplotype id.
        int hID = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, numHaplotypes - 1));

        return std::make_pair(rID, hID);
    }

    // Convert a position (contig, haplotype) to an integer.
    int toId(std::pair<int, int> pos) const
    {
        return pos.first * numHaplotypes + pos.second;
    }            
};

// --------------------------------------------------------------------------
// Class MasonSimulatorApp
// --------------------------------------------------------------------------

class MasonSimulatorApp
{
public:
    // The configuration to use for the simulation.
    MasonSimulatorOptions options;

    // The random number generator to use for the simulation.
    TRng rng;

    // ----------------------------------------------------------------------
    // VCF Materialization
    // ----------------------------------------------------------------------

    // Materialization of the contigs from a VCF file.
    VcfMaterializer vcfMat;

    // ----------------------------------------------------------------------
    // Sample Source Distribution
    // ----------------------------------------------------------------------

    // Helper for distributing reads/pairs to contigs/haplotypes.
    ContigPicker contigPicker;
    // Helper for storing the read ids for each contig/haplotype pair.
    IdSplitter fragmentIdSplitter;
    // Helper for storing the simulated reads for each contig/haplotype pair.  We will write out SAM files with the
    // alignment information relative to the materialized sequence.
    IdSplitter fragmentSplitter;
    // Helper for joining the SAM files.
    // SamJoiner fragmentJoiner;

    // ----------------------------------------------------------------------
    // Fragment and Read Simulation
    // ----------------------------------------------------------------------

    // Generation of infixes on the genome.
    // FragmentGenerator fragGenerator;
    // Sequencing of reads from given fragments.
    // std::auto_ptr<SequencingSimulator> seqSimulator;

    // ----------------------------------------------------------------------
    // File Output
    // ----------------------------------------------------------------------

    // For writing left/right reads.
    // seqan::SequenceStream outSeqsLeft, outSeqsRight;
    // For writing the final SAM/BAM file.
    // seqan::BamStream outBamStream;

    MasonSimulatorApp(MasonSimulatorOptions const & options) :
            options(options), rng(options.seed), contigPicker(rng)
    {}

    int run()
    {
        // Initialize.
        _init();

        // Print the header and the options.
        _printHeader();
        return 0;
    }

    void _init()
    {
        /*
        // Initialize VCF materialization (reference FASTA and input VCF).
        vcfMat.init();

        // Configure contigPicker.
        contigPicker.numHaplotypes = vcfMat.numHaplotypes;
        contigPicker.lengthSums.clear();
        for (unsigned i = 0; i < numSeqs(vcfMat.faiIndex); ++i)
        {
            contigPicker.contigLengths.push_back(sequenceLength(vcfMat.faiIndex, i));
            if (i > 0u)
                contigPicker[i] += contigPicker[i - 1];
        }

        // Configure fragment splitter and joiner.
        fragmentIdSplitter.numContigs = numSeqs(vcfMat.faiIndex);
        fragmentSplitter.numContigs = fragmentIdSplitter.numContigs;
        fragmentJoiner.splitter = &fragmentSplitter;
        fragmentJoiner.init();

        // Initialize seqSimulator.
        */
    }

    void _printHeader()
    {
        std::cerr << "MASON SIMULATOR\n"
                  << "===============\n"
                  << "\n";
        options.print(std::cerr);
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(MasonSimulatorOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("mason_simulator");
    // Set short description, version, and date.
    setShortDescription(parser, "Read Simulation");
    setVersion(parser, "2.0");
    setDate(parser, "July 2012");

    // Define usage line and long description.
    addUsageLine(parser,
                 "[OPTIONS] \\fB-ir\\fP \\fIIN.fa\\fP [\\fB-iv\\fP \\fIIN.vcf\\fP] \\fB-o\\fP \\fILEFT.fq\\fP "
                 "[\\fB-or\\fP \\fIRIGHT.fq\\fP]");
    addDescription(parser,
                   "Simulate reads from the reference sequence \\fIIN.fa\\fP, potentially with variants "
                   "from \\fIIN.vcf\\fP.  In case that both \\fB-o\\fP and \\fB-or\\fP are given, write out "
                   "paired-end data, if only \\fB-io\\fP is given, only single-end reads are simulated.");

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

int main(int argc, char const ** argv)
{
    // Parse options.
    MasonSimulatorOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Initialize Global State
    //
    // Random number generator to use throughout mason.
    TRng rng(options.seed);

    // Run the application.
    MasonSimulatorApp app(options);
    return app.run();
}
