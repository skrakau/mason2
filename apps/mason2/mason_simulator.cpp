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

#include "fragment_generation.h"
#include "sequencing.h"
#include "mason_options.h"
#include "mason_types.h"
#include "vcf_materialization.h"
#include "external_split_merge.h"

// ==========================================================================
// Classes
// ==========================================================================

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

    // The fragment generator and fragment buffer.
    std::vector<Fragment> fragments;
    FragmentSampler * fragSampler;

    // The sequencing simulator to use.
    SequencingSimulator * seqSimulator;

    // Buffer with ids and sequence of reads simulated in this thread.
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;
    seqan::StringSet<seqan::CharString> quals;
    std::vector<SequencingSimulationInfo> infos;
    // Buffer for the BAM alignment records.
    bool buildAlignments;  // Whether or not compute the BAM alignment records.
    std::vector<seqan::BamAlignmentRecord> alignmentRecords;

    ReadSimulatorThread() : options(), fragSampler(), seqSimulator(), buildAlignments(false)
    {}

    ~ReadSimulatorThread()
    {
        delete fragSampler;
        delete seqSimulator;
    }

    void init(int seed, MasonSimulatorOptions const & newOptions)
    {
        reSeed(rng, seed);
        options = &newOptions;
        buildAlignments = !empty(options->outFileNameSam);

        // Initialize fragment generator here with reference to RNG and options.
        fragSampler = new FragmentSampler(rng, options->fragSamplerOptions);

        // Create sequencing simulator.
        SequencingSimulatorFactory simFactory(rng, options->seqOptions, options->illuminaOptions,
                                              options->rocheOptions, options->sangerOptions);
        std::auto_ptr<SequencingSimulator> ptr = simFactory.make();
        seqSimulator = ptr.release();
    }

    void _setId(seqan::CharString & str, std::stringstream & ss, int fragId, int num,
                SequencingSimulationInfo const & info, bool forceNoEmbed = false)
    {
        ss.clear();
        ss.str("");
        ss << options->seqOptions.readNamePrefix;
        if (num == 0 || forceNoEmbed)
            ss << (fragId + 1);
        else if (num == 1)
            ss << (fragId + 1) << "/1";
        else  // num == 2
            ss << (fragId + 1) << "/2";
        if (options->seqOptions.embedReadInfo && !forceNoEmbed)
        {
            ss << ' ';
            info.serialize(ss);
        }
        str = ss.str();
    }

    void _simulatePairedEnd(seqan::Dna5String const & seq, int rID, int hID)
    {
        std::stringstream ss;

        for (unsigned i = 0; i < 2 * fragmentIds.size(); i += 2)
        {
            TFragment frag(seq, fragments[i / 2].beginPos, fragments[i / 2].endPos);
            seqSimulator->simulatePairedEnd(seqs[i], quals[i], infos[i],
                                            seqs[i + 1], quals[i + 1], infos[i + 1],
                                            frag);
            infos[i].rID = infos[i + 1].rID = rID;
            infos[i].hID = infos[i + 1].hID = hID;
            _setId(ids[i], ss, fragmentIds[i / 2], 1, infos[i]);
            _setId(ids[i + 1], ss, fragmentIds[i / 2], 2, infos[i + 1]);
            if (buildAlignments)
            {
                alignmentRecords[i].flag = 0;
                alignmentRecords[i + 1].flag = 0;

                if (!infos[i].isForward)
                {
                    alignmentRecords[i].flag |= seqan::BAM_FLAG_RC;
                    alignmentRecords[i + 1].flag |= seqan::BAM_FLAG_NEXT_RC;
                }
                if (!infos[i + 1].isForward)
                {
                    alignmentRecords[i].flag |= seqan::BAM_FLAG_NEXT_RC;
                    alignmentRecords[i + 1].flag |= seqan::BAM_FLAG_RC;
                }

                // Build alignment records, filling qName, flag, rID, and pos only to save disk space.  We will
                // compute the whole record and reference coordinates when writing out.
                _setId(alignmentRecords[i].qName, ss, fragmentIds[i / 2], 1, infos[i], true);
                alignmentRecords[i].flag |= seqan::BAM_FLAG_ALL_PROPER | seqan::BAM_FLAG_MULTIPLE |
                        seqan::BAM_FLAG_FIRST;
                alignmentRecords[i].rID = rID;
                alignmentRecords[i].beginPos = infos[i].beginPos;
                if (!infos[i].isForward)
                {
                    reverseComplement(seqs[i]);
                    reverse(quals[i]);
                    reverse(infos[i].cigar);
                }
                alignmentRecords[i].cigar = infos[i].cigar;
                alignmentRecords[i].seq = seqs[i];
                alignmentRecords[i].qual = quals[i];
                if (!infos[i].isForward)
                {
                    reverseComplement(seqs[i]);
                    reverse(quals[i]);
                    reverse(infos[i].cigar);
                }

                _setId(alignmentRecords[i + 1].qName, ss, fragmentIds[i / 2], 1, infos[i + 1], true);
                alignmentRecords[i + 1].flag |= seqan::BAM_FLAG_ALL_PROPER | seqan::BAM_FLAG_MULTIPLE |
                        seqan::BAM_FLAG_LAST;
                alignmentRecords[i + 1].rID = rID;
                alignmentRecords[i + 1].beginPos = infos[i + 1].beginPos;
                if (!infos[i + 1].isForward)
                {
                    reverseComplement(seqs[i + 1]);
                    reverse(quals[i + 1]);
                    reverse(infos[i + 1].cigar);
                }
                alignmentRecords[i + 1].cigar = infos[i + 1].cigar;
                alignmentRecords[i + 1].seq = seqs[i + 1];
                alignmentRecords[i + 1].qual = quals[i + 1];
                if (!infos[i + 1].isForward)
                {
                    reverseComplement(seqs[i + 1]);
                    reverse(quals[i + 1]);
                    reverse(infos[i + 1].cigar);
                }

                alignmentRecords[i].rNextId = alignmentRecords[i + 1].rID;
                alignmentRecords[i + 1].rNextId = alignmentRecords[i].rID;
                alignmentRecords[i].pNext = alignmentRecords[i + 1].beginPos;
                alignmentRecords[i + 1].pNext = alignmentRecords[i].beginPos;
                if (alignmentRecords[i].beginPos < alignmentRecords[i + 1].beginPos)
                {
                    alignmentRecords[i].tLen = alignmentRecords[i + 1].beginPos + (int)length(seqs[i]) - alignmentRecords[i].beginPos;
                    alignmentRecords[i + 1].tLen = -alignmentRecords[i].tLen;
                }
                else
                {
                    alignmentRecords[i + 1].tLen = alignmentRecords[i].beginPos + (int)length(seqs[i + 1]) - alignmentRecords[i + 1].beginPos;
                    alignmentRecords[i].tLen = -alignmentRecords[i + 1].tLen;
                }
            }
        }
    }

    void _simulateSingleEnd(seqan::Dna5String const & seq, int rID, int hID)
    {
        std::stringstream ss;

        for (unsigned i = 0; i < fragmentIds.size(); ++i)
        {
            TFragment frag(seq, fragments[i].beginPos, fragments[i].endPos);
            seqSimulator->simulateSingleEnd(seqs[i], quals[i], infos[i], frag);
            _setId(ids[i], ss, fragmentIds[i], 0, infos[i]);
            if (buildAlignments)
            {
                // Build alignment records, filling qName, flag, rID, and pos only to save disk space.  We will
                // compute the whole record and reference coordinates when writing out.
                infos[i].rID = rID;
                infos[i].hID = hID;
                _setId(alignmentRecords[i].qName, ss, fragmentIds[i], 1, infos[i], true);
                alignmentRecords[i].flag = 0;
                alignmentRecords[i].rID = rID;
                alignmentRecords[i].beginPos = infos[i].beginPos;
                if (!infos[i].isForward)
                {
                    alignmentRecords[i].flag = seqan::BAM_FLAG_RC;
                    reverseComplement(seqs[i]);
                    reverse(quals[i]);
                    reverse(infos[i].cigar);
                }
                alignmentRecords[i].cigar = infos[i].cigar;
                alignmentRecords[i].seq = seqs[i];
                alignmentRecords[i].qual = quals[i];
                if (!infos[i].isForward)
                {
                    reverseComplement(seqs[i]);
                    reverse(quals[i]);
                    reverse(infos[i].cigar);
                }
            }
        }
    }

    // Simulate next chunk.
    void run(seqan::Dna5String const & seq, int rID, int hID)
    {
        // Sample fragments.
        fragSampler->generateMany(fragments, rID, length(seq), fragmentIds.size());

        // Simulate reads.
        int seqCount = (options->seqOptions.simulateMatePairs ? 2 : 1) * fragmentIds.size();
        resize(ids, seqCount);
        resize(seqs, seqCount);
        resize(quals, seqCount);
        infos.resize(seqCount);
        if (buildAlignments)
            alignmentRecords.resize(seqCount);
        // TODO(holtgrew): Optimize number of virtual function calls.
        if (options->seqOptions.simulateMatePairs)
            _simulatePairedEnd(seq, rID, hID);
        else
            _simulateSingleEnd(seq, rID, hID);
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
    // Helper for joining the FASTQ files.
    std::auto_ptr<FastxJoiner<seqan::Fastq> > fastxJoiner;
    // Helper for storing SAM records for each contig/haplotype pair.  In the end, we will join this again.
    IdSplitter alignmentSplitter;
    // Helper for joining the SAM files.
    std::auto_ptr<SamJoiner> alignmentJoiner;

    // ----------------------------------------------------------------------
    // Header used for writing temporary SAM.
    // ----------------------------------------------------------------------

    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
    typedef seqan::BamIOContext<TNameStore>     TBamIOContext;
    TNameStore nameStore;
    seqan::BamHeader header;
    TNameStoreCache nameStoreCache;
    TBamIOContext   bamIOContext;

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
            contigPicker(rng), nameStoreCache(nameStore), bamIOContext(nameStore, nameStoreCache)
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

    void _simulateReadsDoSimulation()
    {
        std::cerr << "\nSimulating Reads:\n";
        int haplotypeCount = vcfMat.numHaplotypes;
        seqan::Dna5String contigSeq;  // materialized contig
        int rID = 0;  // current reference id
        int hID = 0;  // current haplotype id
        int contigFragmentCount = 0;  // number of reads on the contig
        // Note that all shared variables are correctly synchronized by implicit flushes at the critical sections below.
        while (vcfMat.materializeNext(contigSeq, rID, hID))
        {
            std::cerr << "  " << sequenceName(vcfMat.faiIndex, rID) << " (allele " << (hID + 1) << ") ";
            contigFragmentCount = 0;

            SEQAN_OMP_PRAGMA(parallel num_threads(options.numThreads))
            {
                int tID = omp_get_thread_num();

                while (true)  // Execute as long as there are fragments left.
                {
                    // Read in the ids of the fragments to simulate.
                    threads[tID].fragmentIds.resize(options.chunkSize);  // make space
                    bool doBreak = false;
                    SEQAN_OMP_PRAGMA(critical(fragment_ids))
                    {
                        // Load the fragment ids to simulate for.
                        int numRead = fread(&threads[tID].fragmentIds[0], sizeof(int), options.chunkSize,
                                            fragmentIdSplitter.files[rID * haplotypeCount + hID]);
                        contigFragmentCount += numRead;
                        if (numRead == 0)
                            doBreak = true;
                        threads[tID].fragmentIds.resize(numRead);
                    }
                    if (doBreak)
                        break;  // No more work left.
                    
                    // Perform the simulation.
                    threads[tID].run(contigSeq, rID, hID);
                    
                    // Write out the temporary sequence.
                    SEQAN_OMP_PRAGMA(critical(seq_io))
                    {
                        if (write2(fragmentSplitter.files[rID * haplotypeCount + hID],
                                   threads[tID].ids, threads[tID].seqs, threads[tID].quals, seqan::Fastq()))
                            throw MasonIOException("Could not write out temporary sequence.");
                        if (!empty(options.outFileNameSam))
                            for (unsigned i = 0; i < length(threads[tID].alignmentRecords); ++i)
                            {
                                if (write2(alignmentSplitter.files[rID * haplotypeCount + hID],
                                           threads[tID].alignmentRecords[i], bamIOContext, seqan::Sam()) != 0)
                                    throw MasonIOException("Could not write out temporary alignment record.");
                            }
                    }
                    
                    SEQAN_OMP_PRAGMA(critical(io_log))
                    {
                        std::cerr << '.' << std::flush;
                    }
                }
            }
            
            std::cerr << " (" << contigFragmentCount << " fragments) OK\n";
        }
        std::cerr << "  Done simulating reads.\n";
    }

    void _simulateReadsJoin()
    {
        std::cerr << "\nJoining temporary files ...";
        fragmentSplitter.reset();
        fastxJoiner.reset(new FastxJoiner<seqan::Fastq>(fragmentSplitter));
        FastxJoiner<seqan::Fastq> & joiner = *fastxJoiner.get();  // Shortcut
        // TODO(holtgrew): Use bulk-reading calls.
        seqan::CharString id, seq, qual;
        if (options.seqOptions.simulateMatePairs)
            while (!joiner.atEnd())
            {
                joiner.get(id, seq, qual);
                if (writeRecord(outSeqsLeft, id, seq, qual) != 0)
                    throw MasonIOException("Problem joining sequences.");
                joiner.get(id, seq, qual);
                if (writeRecord(outSeqsRight, id, seq, qual) != 0)
                    throw MasonIOException("Problem joining sequences.");
            }
        else
            while (!joiner.atEnd())
            {
                joiner.get(id, seq, qual);
                if (writeRecord(outSeqsLeft, id, seq, qual) != 0)
                    throw MasonIOException("Problem joining sequences.");
            }
        if (!empty(options.outFileNameSam))
        {
            alignmentSplitter.reset();
            alignmentJoiner.reset(new SamJoiner(alignmentSplitter));

            outBamStream.header = alignmentJoiner->header;

            SamJoiner & joiner = *alignmentJoiner.get();  // Shortcut
            seqan::BamAlignmentRecord record;
            while (!joiner.atEnd())
            {
                joiner.get(record);
                if (writeRecord(outBamStream, record) != 0)
                    throw MasonIOException("Problem writing to alignment out file.");
            }
        }
        std::cerr << " OK\n";
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
        _simulateReadsDoSimulation();

        // (3) Merge the sequences from external files into the output stream.
        _simulateReadsJoin();
    }

    // Initialize the alignment splitter data structure.
    void _initAlignmentSplitter()
    {
        // Open alignment splitters.
        alignmentSplitter.numContigs = fragmentIdSplitter.numContigs;
        alignmentSplitter.open();
        // Build and write out header, fill ref name store.
        seqan::BamHeaderRecord vnHeaderRecord;
        vnHeaderRecord.type = seqan::BAM_HEADER_FIRST;
        appendValue(vnHeaderRecord.tags, seqan::Pair<seqan::CharString>("VN", "1.4"));
        appendValue(header.records, vnHeaderRecord);
        resize(header.sequenceInfos, numSeqs(vcfMat.faiIndex));
        for (unsigned i = 0; i < numSeqs(vcfMat.faiIndex); ++i)
        {
            if (!empty(options.matOptions.vcfFileName))
                header.sequenceInfos[i].i1 = vcfMat.vcfStream.header.sequenceNames[i];
            else
                header.sequenceInfos[i].i1 = sequenceName(vcfMat.faiIndex, i);
            unsigned idx = 0;
            if (!getIdByName(vcfMat.faiIndex, header.sequenceInfos[i].i1, idx))
            {
                std::stringstream ss;
                ss << "Could not find " << header.sequenceInfos[i].i1 << " from VCF file in FAI index.";
                throw MasonIOException(ss.str());
            }
            header.sequenceInfos[i].i2 = sequenceLength(vcfMat.faiIndex, idx);
            appendValue(nameStore, sequenceName(vcfMat.faiIndex, idx));
            seqan::BamHeaderRecord seqHeaderRecord;
            seqHeaderRecord.type = seqan::BAM_HEADER_REFERENCE;
            appendValue(seqHeaderRecord.tags, seqan::Pair<seqan::CharString>("SN", header.sequenceInfos[i].i1));
            std::stringstream ss;
            ss << header.sequenceInfos[i].i2;
            appendValue(seqHeaderRecord.tags, seqan::Pair<seqan::CharString>("LN", ss.str().c_str()));
            appendValue(header.records, seqHeaderRecord);
        }
        refresh(nameStoreCache);
        for (unsigned i = 0; i < alignmentSplitter.files.size(); ++i)
            if (write2(alignmentSplitter.files[i], header, bamIOContext, seqan::Sam()) != 0)
                throw MasonIOException("Could not write out SAM header to temporary file.");
    }

    // Configure contigPicker.
    void _initContigPicker()
    {
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
        // Splitter for alignments, only required when writing out SAM/BAM.
        if (!empty(options.outFileNameSam))
            _initAlignmentSplitter();
        std::cerr << " OK\n";
    }

    // Open the output files.
    void _initOpenOutputFiles()
    {
        std::cerr << "Opening output file " << options.outFileNameLeft << " ...";
        open(outSeqsLeft, toCString(options.outFileNameLeft), seqan::SequenceStream::WRITE);
        outSeqsLeft.outputOptions = seqan::SequenceOutputOptions(0);  // also FASTA in one line
        if (!isGood(outSeqsLeft))
            throw MasonIOException("Could not open left/single-end output file.");
        std::cerr << " OK\n";

        if (!empty(options.outFileNameRight))
        {
            std::cerr << "Opening output file " << options.outFileNameRight << " ...";
            open(outSeqsRight, toCString(options.outFileNameRight), seqan::SequenceStream::WRITE);
            outSeqsRight.outputOptions = seqan::SequenceOutputOptions(0);  // also FASTA in one line
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

    void _init()
    {
        std::cerr << "\n____INITIALIZING______________________________________________________________\n"
                  << "\n";
        
        // Initialize VCF materialization (reference FASTA and input VCF).
        std::cerr << "Opening reference and variants file ...";
        vcfMat.init();
        std::cerr << " OK\n";

        // Configure contigPicker and fragment id splitter.
        _initContigPicker();

        // Initialize simulation threads.
        std::cerr << "Initializing simulation threads ...";
        threads.resize(options.numThreads);
        for (int i = 0; i < options.numThreads; ++i)
            threads[i].init(options.seed + i * options.seedSpacing, options);
        std::cerr << " OK\n";

        // Open output files.
        _initOpenOutputFiles();
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
