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
    std::vector<__int64> lengthSums;
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
            seqan::Pdf<seqan::Uniform<__int64> > pdf(0, lengthSums.back() - 1);
            __int64 x = pickRandomNumber(rng, pdf);
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
// Class ReadSimulatorThread
// --------------------------------------------------------------------------

// State for one thread for simulation of reads.

class ReadSimulatorThread
{
public:
    // Options for the read simulation.
    MasonSimulatorOptions const * options;

    // The random number generator to use for this thread.
    TRng rng;

    // The ids of the fragments.
    std::vector<int> fragmentIds;

    // Buffer with ids and sequence of reads simulated in this thread.
    seqan::StringSet<seqan::CharString, seqan::Owner<seqan::ConcatDirect<> > > ids;
    seqan::StringSet<seqan::Dna5String, seqan::Owner<seqan::ConcatDirect<> > > seqs;
    seqan::StringSet<seqan::CharString, seqan::Owner<seqan::ConcatDirect<> > > quals;

    ReadSimulatorThread() : options()
    {}

    void init(int seed, MasonSimulatorOptions const & newOptions)
    {
        reSeed(rng, seed);
        options = &newOptions;
    }

    // Simulate next chunk.
    void run()
    {
        clear(ids);
        clear(seqs);
        clear(quals);
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

    // Threads used for simulation.
    std::vector<ReadSimulatorThread> threads;
    
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
    seqan::SequenceStream outSeqsLeft, outSeqsRight;
    // For writing the final SAM/BAM file.
    seqan::BamStream outBamStream;

    MasonSimulatorApp(MasonSimulatorOptions const & options) :
            options(options), rng(options.seed),
            vcfMat(toCString(options.matOptions.fastaFileName), toCString(options.matOptions.vcfFileName)),
            contigPicker(rng)
    {}

    int run()
    {
        // Print the header and the options.
        _printHeader();
        // Initialize.
        _init();
        // Simulate reads.
        _simulateReads();

        return 0;
    }

    void _simulateReads()
    {
        std::cerr << "\n____READ SIMULATION___________________________________________________________\n"
                  << "\n";

        // (1) Distribute read ids to the contigs/haplotypes.
        //
        // We will simulate the reads in the order of contigs/haplotypes and in a final join script generate output file
        // that are sorted by read id.
        int seqCount = numSeqs(vcfMat.faiIndex);
        int haplotypeCount = vcfMat.numHaplotypes;
        std::cerr << "Distributing fragments to " << seqCount << " contigs (" << haplotypeCount
                  << " haplotypes each) ...";
        for (int i = 0; i < options.numFragments; ++i)
            fwrite(&i, sizeof(int), 1, fragmentIdSplitter.files[contigPicker.toId(contigPicker.pick())]);
        fragmentIdSplitter.reset();
        std::cerr << " OK\n";

        // (2) Simulate the reads in the order of contigs/haplotypes.
        std::cerr << "\nSimulating Reads:\n\n";
        seqan::Dna5String contigSeq;  // materialized contig
        int rID = 0;  // current reference id
        int hID = 0;  // current haplotype id
        int contigFragmentCount = 0;  // number of reads on the contig
        // Note that all shared variables are correctly synchronized by implicit flushes at the critical sections below.
        while (vcfMat.materializeNext(contigSeq, rID, hID))
        {
            std::cerr << sequenceName(vcfMat.faiIndex, rID) << " (allele " << (hID + 1) << ") ";
            contigFragmentCount = 0;

            SEQAN_OMP_PRAGMA(num_threads(options.numThreads))
            {
                int tID = omp_get_thread_num();

                while (true)  // Execute as long as there are fragments left.
                {
                    // Read in the ids of the fragments to simulate.
                    threads[tID].fragmentIds.resize(options.chunkSize);  // make space
                    SEQAN_OMP_PRAGMA(critical(fragment_ids))
                    {
                        // Load the fragment ids to simulate for.
                        int numRead = fread(&threads[tID].fragmentIds[0], sizeof(int), options.chunkSize,
                                            fragmentIdSplitter.files[rID * haplotypeCount + hID]);
                        contigFragmentCount += numRead;
                        if (numRead == 0)
                            break;  // No more work left.
                        threads[tID].fragmentIds.resize(numRead);
                    }
                    
                    // Perform the simulation.
                    threads[tID].run();
                    
                    // Write out the temporary sequence.
                    SEQAN_OMP_PRAGMA(critical(seq_io))
                    {
                        if (write2(fragmentSplitter.files[rID * haplotypeCount + hID],
                                   threads[tID].ids, threads[tID].seqs, threads[tID].quals, seqan::Fastq()))
                            throw MasonIOException("Could not write out temporary sequence.");
                    }
                    
                    SEQAN_OMP_PRAGMA(critical(io_log))
                    {
                        std::cerr << '.' << std::flush;
                    }
                }
            }
            
            std::cerr << " (" << contigFragmentCount << " fragments) OK\n";
        }
        std::cerr << "\nDone simulating reads.\n";
    }

    void _init()
    {
        std::cerr << "\n____INITIALIZING______________________________________________________________\n"
                  << "\n";
        
        // Initialize VCF materialization (reference FASTA and input VCF).
        std::cerr << "Opening reference and variants file ...";
        vcfMat.init();
        std::cerr << " OK\n";

        // Configure contigPicker and fragment id splitter.
        std::cerr << "Initializing fragment-to-contig distribution ...";
        // Contig picker.
        contigPicker.numHaplotypes = vcfMat.numHaplotypes;
        contigPicker.lengthSums.clear();
        for (unsigned i = 0; i < numSeqs(vcfMat.faiIndex); ++i)
        {
            contigPicker.lengthSums.push_back(sequenceLength(vcfMat.faiIndex, i));
            if (i > 0u)
                contigPicker.lengthSums[i] += contigPicker.lengthSums[i - 1];
        }
        // Fragment id splitter.
        fragmentIdSplitter.numContigs = numSeqs(vcfMat.faiIndex) * vcfMat.numHaplotypes;
        fragmentIdSplitter.open();
        // Splitter for sequence.
        fragmentSplitter.numContigs = fragmentIdSplitter.numContigs;
        fragmentSplitter.open();
        std::cerr << " OK\n";

        // Configure fragment splitter and joiner.
        // fragmentJoiner.splitter = &fragmentSplitter;
        // fragmentJoiner.init();

        // Initialize simulation threads.
        std::cerr << "Initializing simulation threads ...";
        threads.resize(options.numThreads);
        for (int i = 0; i < options.numThreads; ++i)
            threads[i].init(options.seed + i * options.seedSpacing, options);
        std::cerr << " OK\n";

        // Open output files.
        std::cerr << "Opening output file " << options.outFileNameLeft << " ...";
        open(outSeqsLeft, toCString(options.outFileNameLeft), seqan::SequenceStream::WRITE);
        if (!isGood(outSeqsLeft))
            throw MasonIOException("Could not open left/single-end output file.");
        std::cerr << " OK\n";

        if (!empty(options.outFileNameRight))
        {
            std::cerr << "Opening output file " << options.outFileNameRight << " ...";
            open(outSeqsRight, toCString(options.outFileNameRight), seqan::SequenceStream::WRITE);
            if (!isGood(outSeqsRight))
                throw MasonIOException("Could not open right/single-end output file.");
            std::cerr << " OK\n";
        }

        if (!empty(options.outFileNameSam))
        {
            std::cerr << "Opening output file " << options.outFileNameSam << "...";
            open(outBamStream, toCString(options.outFileNameSam), seqan::BamStream::WRITE);
            if (!isGood(outBamStream))
                throw MasonIOException("Could not open SAM/BAM output file.");
            std::cerr << " OK\n";
        }
    }

    void _printHeader()
    {
        std::cerr << "MASON SIMULATOR\n"
                  << "===============\n";
        if (options.verbosity >= 2)
        {
            std::cerr << "\n";
            options.print(std::cerr);
        }
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
                 "[OPTIONS] \\fB-ir\\fP \\fIIN.fa\\fP \\fB-n\\fP \\fINUM\\fP [\\fB-iv\\fP \\fIIN.vcf\\fP] \\fB-o\\fP \\fILEFT.fq\\fP "
                 "[\\fB-or\\fP \\fIRIGHT.fq\\fP]");
    addDescription(parser,
                   "Simulate \\fINUM\\fP reads/pairs from the reference sequence \\fIIN.fa\\fP, potentially with "
                   "variants from \\fIIN.vcf\\fP.  In case that both \\fB-o\\fP and \\fB-or\\fP are given, write out "
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
