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
#include  <seqan/random.h>
#include  <seqan/sequence.h>
#include  <seqan/seq_io.h>
#include  <seqan/vcf.h>

// ==========================================================================
// Forwards
// ==========================================================================

// ==========================================================================
// Classes
// ==========================================================================

typedef seqan::Rng<> TRng;

// --------------------------------------------------------------------------
// Class MasonVariatorOptions
// --------------------------------------------------------------------------

struct MasonVariatorOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Seed for RNG.
    int seed;

    // ----------------------------------------------------------------------
    // Input / Output Options
    // ----------------------------------------------------------------------

    // VCF file to import.
    seqan::CharString vcfInFile;
    // FASTA file to import.
    seqan::CharString fastaInFile;
    // VCF file to write out.
    seqan::CharString vcfOutFile;
    // FASTA file to write out with variations.
    seqan::CharString fastaOutFile;

    // Path to a TSV file where the first two columns giving the type of the SV to simulate and the size of the SV.
    // This overrides the simulation of SV from the sv*Rate parameters.
    seqan::CharString inputSVSizeFile;

    // Separator for haplotype contigs in resulting FASTA file.
    seqan::CharString haplotypeNameSep;

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
            verbosity(1), seed(0),
            snpRate(0), smallIndelRate(0), minSmallIndelSize(0), maxSmallIndelSize(0), svIndelRate(0),
            svInversionRate(0), svTranslocationRate(0), svDuplicationRate(0), minSVSize(0), maxSVSize(0)
    {}
};

void print(std::ostream & out, MasonVariatorOptions const & options)
{
    out << "__OPTIONS_____________________________________________________________________\n"
        << "\n"
        << "VCF IN               \t" << options.vcfInFile << "\n"
        << "FASTA IN             \t" << options.fastaInFile << "\n"
        << "SV SIZE TSV IN       \t" << options.inputSVSizeFile << "\n"
        << "VCF OUT              \t" << options.vcfOutFile << "\n"
        << "FASTA OUT            \t" << options.fastaOutFile << "\n"
        << "HAPLOTYPE NAME SEP   \t" << options.haplotypeNameSep << "\n"
        << "\n"
        << "NUM HAPLOTYPES       \t" << options.numHaplotypes << "\n"
        << "HAPLOTYPE SEP        \t\"" << options.haplotypeSep << "\"\n"
        << "\n"
        << "SNP RATE             \t" << options.snpRate << "\n"
        << "SMALL INDEL RATE     \t" << options.smallIndelRate << "\n"
        << "\n"
        << "MIN SMALL INDEL SIZE \t" << options.minSmallIndelSize << "\n"
        << "MAX SMALL INDEL SIZE \t" << options.maxSmallIndelSize << "\n"
        << "\n"
        << "SV INDEL RATE        \t" << options.svIndelRate << "\n"
        << "SV INVERSION RATE    \t" << options.svInversionRate << "\n"
        << "SV TRANSLOCATION RATE\t" << options.svTranslocationRate << "\n"
        << "SV DUPLICATION RATE  \t" << options.svDuplicationRate << "\n"
        << "\n"
        << "MIN SV SIZE          \t" << options.minSVSize << "\n"
        << "MAX SV SIZE          \t" << options.maxSVSize << "\n"
        << "\n";
}

// --------------------------------------------------------------------------
// Class SnpRecord
// --------------------------------------------------------------------------

// Represents a SNP in one haplotype.

struct SnpRecord
{
    // Reference id and position on the reference.
    int rId;
    int pos;

    // The haplotype that this variation belongs to.
    int haplotype;

    // The target nucleotide.
    seqan::Dna5 to;

    SnpRecord() : rId(-1), pos(-1), haplotype(-1), to('\0')
    {}

    SnpRecord(int haplotype, int rId, int pos, char to) :
            rId(rId), pos(pos), haplotype(haplotype), to(to)
    {}

    bool operator<(SnpRecord const & other) const
    {
        if (rId < other.rId || (rId == other.rId && pos < other.pos) ||
            (rId == other.rId && pos == other.pos && haplotype < other.haplotype))
            return true;
        return false;
    }

    std::pair<int, int> getPos() const
    {
        return std::make_pair(rId, pos);
    }
};

// --------------------------------------------------------------------------
// Class SmallIndelRecord
// --------------------------------------------------------------------------

// Represents a small indel.

struct SmallIndelRecord
{
    // Reference id and position on the reference.
    int rId;
    int pos;

    // The haplotype that this variation belongs to.
    int haplotype;

    // The size of the indel, negative numbers for deletions, positive numbers for insertions.
    int size;

    // The inserted sequence if any.
    seqan::CharString seq;

    SmallIndelRecord() : rId(-1), pos(-1), haplotype(-1), size(0)
    {}

    SmallIndelRecord(int haplotype, int rId, int pos, int size, seqan::CharString const & seq) :
            rId(rId), pos(pos), haplotype(haplotype), size(size), seq(seq)
    {}

    bool operator<(SmallIndelRecord const & other) const
    {
        if (rId < other.rId || (rId == other.rId && pos < other.pos) ||
            (rId == other.rId && pos == other.pos && haplotype < other.haplotype))
            return true;
        return false;
    }

    std::pair<int, int> getPos() const
    {
        return std::make_pair(rId, pos);
    }
};

// --------------------------------------------------------------------------
// Class Variants
// --------------------------------------------------------------------------

// Contains the simulated variants.

struct Variants
{
    seqan::String<SnpRecord> snps;
    seqan::String<SmallIndelRecord> smallIndels;
};

// --------------------------------------------------------------------------
// Class SmallVariantSimulator
// --------------------------------------------------------------------------

// Simulation of small variants.

class SmallVariantSimulator
{
public:
    // Random number generator.
    TRng & rng;

    // FAI Index for loading sequence contig-wise.
    seqan::FaiIndex const & faiIndex;

    // The variator options.
    MasonVariatorOptions options;

    SmallVariantSimulator(TRng & rng, seqan::FaiIndex const & faiIndex, MasonVariatorOptions const & options) :
            rng(rng), faiIndex(faiIndex), options(options)
    {}

    // Perform simulation for one contig.
    void simulateContig(Variants & variants, unsigned rId, int haploCount)
    {
        seqan::CharString seq;
        if (readSequence(seq, faiIndex, rId) != 0)
        {
            std::cerr << "ERROR: Could not read sequence " << rId << " from FASTA file!\n";
            return;
        }

        // For each base, compute the whether to simulate a SNP and/or small indel.
        for (unsigned pos = 0; pos < length(seq); ++pos)
        {
            seqan::Pdf<seqan::Uniform<double> > pdf(0, 1);

            // Perform experiment for SNP and small indel.
            bool isSnp = (pickRandomNumber(rng, pdf) < options.snpRate);
            double isIndel = (pickRandomNumber(rng, pdf) < options.smallIndelRate);
            while (isSnp && isIndel)  // TODO(holtgrew): Add limit?
            {
                isSnp = (pickRandomNumber(rng, pdf) < options.snpRate);
                isIndel = (pickRandomNumber(rng, pdf) < options.smallIndelRate);
                if (pos == 0)
                    isIndel = false;  // No indel at beginning, complex VCF case.
            }

            // Simulate either SNP or indel.  In the case of a deletion, advance position such that there
            // is no variation in the deleted sequence.
            if (isSnp)
            {
                simulateSnp(variants, seq, haploCount, rId, pos);
            }
            else if (isIndel)
            {
                simulateSmallIndel(variants, seq, haploCount, rId, pos);
                if (back(variants.smallIndels).size < 0)
                    pos += -back(variants.smallIndels).size;
            }
        }
    }

    void simulateSnp(Variants & variants, seqan::CharString & seq, int haploCount, int rId, unsigned pos)
    {
        // We simulate an alternative base for each haplotype.

        seqan::Dna5 from = seq[pos];
        for (int hId = 0; hId < haploCount; ++hId)
        {
            int toInt = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, 2));
            if (ordValue(from) <= toInt)
                toInt += 1;
            // std::cerr << hId << "\t" << rId << "\t" << pos << "\t" << from << "\t" << seqan::Dna5(toInt) << "\n";
            SEQAN_ASSERT_NEQ((int)ordValue(from), toInt);
            seqan::Dna5 to(toInt);
            appendValue(variants.snps, SnpRecord(hId, rId, pos, to));
        }
    }

    void simulateSmallIndel(Variants & variants, seqan::CharString & /*seq*/, int haploCount, int rId, unsigned pos)
    {
        // Indels are simulated for one haplotype only.
        int hId = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, haploCount - 1));
        seqan::CharString indelSeq;
        reserve(indelSeq, options.maxSmallIndelSize);
        int indelSize = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(options.minSmallIndelSize,
                                                                               options.maxSmallIndelSize));
        bool deletion = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, 1));
        indelSize = deletion ? -indelSize : indelSize;
        seqan::Pdf<seqan::Uniform<int> > pdf(0, 3);
        for (int i = 0; i < indelSize; ++i)  // not executed in case of deleted sequence
            appendValue(indelSeq, seqan::Dna5(pickRandomNumber(rng, pdf)));
        appendValue(variants.smallIndels, SmallIndelRecord(hId, rId, pos, indelSize, indelSeq));
    }
};

// --------------------------------------------------------------------------
// Class MasonVariatorApp
// --------------------------------------------------------------------------

class MasonVariatorApp
{
public:
    TRng & rng;

    MasonVariatorOptions options;

    seqan::VcfStream vcfStream;
    seqan::SequenceStream outSeqStream;

    // FAI Index for loading sequence contig-wise.
    seqan::FaiIndex const & faiIndex;

    MasonVariatorApp(TRng & rng, seqan::FaiIndex const & faiIndex,
                     MasonVariatorOptions const & options) :
            rng(rng), options(options), faiIndex(faiIndex)
    {}

    int run()
    {
        // Open VCF stream to write to.
        open(vcfStream, toCString(options.vcfOutFile), seqan::VcfStream::WRITE);
        if (!isGood(vcfStream))
        {
            std::cerr << "Could not open " << options.vcfOutFile << " for writing.\n";
            return 1;
        }
        // Copy over sequence names.
        for (unsigned i = 0; i < numSeqs(faiIndex); ++i)
            appendName(*vcfStream._context.sequenceNames,
                       sequenceName(faiIndex, i),
                       vcfStream._context.sequenceNamesCache);
        // Copy over sample names.
        appendName(*vcfStream._context.sampleNames, "simulated", vcfStream._context.sampleNamesCache);

        // Open output FASTA file if necessary.
        if (!empty(options.fastaOutFile))
        {
            open(outSeqStream, toCString(options.fastaOutFile), seqan::SequenceStream::WRITE,
                 seqan::SequenceStream::FASTA);
            if (!isGood(outSeqStream))
            {
                std::cerr << "ERROR: Could not open " << options.fastaOutFile << " for writing!\n";
                return 1;
            }
        }

        // TODO(holtgrew): Write header.

        // Actually perform the variant simulation.
        SmallVariantSimulator smallSim(rng, faiIndex, options);
        for (int rId = 0; rId < (int)numSeqs(faiIndex); ++rId)  // ref seqs
            _simulateContig(smallSim, options, rId);

        return 0;
    }

    // Perform simulation of one contig.
    int _simulateContig(SmallVariantSimulator & sim,
                        MasonVariatorOptions const & options,
                        int rId)
    {
        if (options.verbosity >= 1)
            std::cerr << "  " << sequenceName(faiIndex, rId) << "\n";

        // Simulate variants.
        Variants variants;
        sim.simulateContig(variants, rId, options.numHaplotypes);
        if (options.verbosity >= 1)
            std::cerr << "  snps:         " << length(variants.snps) << "\n"
                      << "  small indels: " << length(variants.smallIndels) << "\n";

        // Load contig seq.
        seqan::Dna5String contig;
        if (readSequence(contig, faiIndex, rId) != 0)
        {
            std::cerr << "Could not read contig seq " << rId << "\n";
            return 1;
        }

        // Write out variants for contig to VCF file.
        if (_writeVcf(contig, variants, rId) != 0)
            return 1;

        // Apply variants to contigs.
        if (!empty(options.fastaOutFile))
            for (int hId = 0; hId < options.numHaplotypes; ++hId)
            {
                if (_writeContigs(contig, variants, rId, hId) != 0)
                    return 1;
            }

        return 0;
    }

    int _writeContigs(seqan::Dna5String const & contig, Variants const & variants, int rId, int hId)
    {
        // Build sequence id.
        seqan::CharString id = sequenceName(faiIndex, rId);
        append(id, options.haplotypeSep);
        char buffer[20];
        snprintf(buffer, 19, "%d", hId);
        append(id, buffer);

        // Build contig.
        seqan::CharString seq;
        // Fors this, we have to iterate in parallel over SNP and small indel records.
        //
        // Current index in snp/small indel array.
        unsigned snpsIdx = 0;
        unsigned smallIndelIdx = 0;
        // Current SNP record, default to sentinel.
        SnpRecord snpRecord;
        snpRecord.rId = seqan::maxValue<int>();
        if (snpsIdx < length(variants.snps))
            snpRecord = variants.snps[snpsIdx++];
        // Current small indel record, default to sentinel.
        SmallIndelRecord smallIndelRecord;
        smallIndelRecord.rId = seqan::maxValue<int>();
        if (smallIndelIdx < length(variants.smallIndels))
            smallIndelRecord = variants.smallIndels[smallIndelIdx++];
        int lastPos = 0;
        if (options.verbosity >= 2)
            std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";

        if (options.verbosity >= 2)
            std::cerr << "building output\n";
        while (snpRecord.rId != seqan::maxValue<int>() || smallIndelRecord.rId != seqan::maxValue<int>())
        {
            // TODO(holtgrew): Extract SNP and small indel handling in functions.
            if (snpRecord.getPos() < smallIndelRecord.getPos())  // process SNP records
            {
                if (snpRecord.haplotype == hId)  // Ignore all but the current contig.
                {
                    if (options.verbosity >= 2)
                        std::cerr << "append(seq, infix(contig, " << lastPos << ", " << snpRecord.pos << ") " << __LINE__ << "\n";
                    append(seq, infix(contig, lastPos, snpRecord.pos));  // interim chars
                       
                    SEQAN_ASSERT_GEQ(snpRecord.pos, lastPos);
                    if (options.verbosity >= 2)
                        std::cerr << "appendValue(seq, " << snpRecord.to << "')\n";
                    appendValue(seq, snpRecord.to);
                    lastPos = snpRecord.pos + 1;
                    if (options.verbosity >= 2)
                        std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";
                }

                if (snpsIdx >= length(variants.snps))
                    snpRecord.rId = seqan::maxValue<int>();
                else
                    snpRecord = variants.snps[snpsIdx++];
            }
            else
            {
                if (smallIndelRecord.haplotype == hId)  // Ignore all but the current contig.
                {
                    if (smallIndelRecord.size > 0)
                    {
                        if (options.verbosity >= 2)
                            std::cerr << "append(seq, infix(contig, " << lastPos << ", " << smallIndelRecord.pos << ") " << __LINE__ << "\n";
                        append(seq, infix(contig, lastPos, smallIndelRecord.pos));  // interim chars

                        SEQAN_ASSERT_GEQ(smallIndelRecord.pos, lastPos);
                        if (options.verbosity >= 2)
                            std::cerr << "append(seq, \"" << smallIndelRecord.seq << "\") " << __LINE__ << "\n";
                        append(seq, smallIndelRecord.seq);
                        lastPos = smallIndelRecord.pos;
                        if (options.verbosity >= 2)
                            std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";
                    }
                    else
                    {
                        if (options.verbosity >= 2)
                            std::cerr << "append(seq, infix(contig, " << lastPos << ", " << smallIndelRecord.pos << ") " << __LINE__ << "\n";
                        append(seq, infix(contig, lastPos, smallIndelRecord.pos));  // interim chars
                        lastPos = smallIndelRecord.pos - smallIndelRecord.size;
                        if (options.verbosity >= 2)
                            std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";
                    }
                }

                if (smallIndelIdx >= length(variants.smallIndels))
                    smallIndelRecord.rId = seqan::maxValue<int>();
                else
                    smallIndelRecord = variants.smallIndels[smallIndelIdx++];
            }
        }
        // Insert interim characters.
        if (options.verbosity >= 2)
            std::cerr << "append(seq, infix(contig, infix(contig, " << lastPos << ", " << length(contig) << ")\n";
        append(seq, infix(contig, lastPos, length(contig)));

        return writeRecord(outSeqStream, id, seq);
    }

    // Write out variants for the given contig to the VCF file.
    int _writeVcf(seqan::Dna5String const & contig, Variants const & variants, int rId)
    {
        // Current index in snp/small indel array.
        unsigned snpsIdx = 0;
        unsigned smallIndelIdx = 0;
        // Current SNP record, default to sentinel.
        SnpRecord snpRecord;
        snpRecord.rId = seqan::maxValue<int>();
        if (snpsIdx < length(variants.snps))
            snpRecord = variants.snps[snpsIdx++];
        // Current small indel record, default to sentinel.
        SmallIndelRecord smallIndelRecord;
        smallIndelRecord.rId = seqan::maxValue<int>();
        if (smallIndelIdx < length(variants.smallIndels))
            smallIndelRecord = variants.smallIndels[smallIndelIdx++];
        while (snpRecord.rId != seqan::maxValue<int>() || smallIndelRecord.rId != seqan::maxValue<int>())
        {
            SEQAN_ASSERT(snpRecord.getPos() != smallIndelRecord.getPos());  // are generated indendently

            // TODO(holtgrew): Extract SNP and small indel generation in their own methods.
            if (snpRecord.getPos() < smallIndelRecord.getPos())  // process SNP records
            {
                if (options.verbosity >= 2)
                    std::cerr << "  snpRecord record.\n";

                // Store information used below.
                // std::cerr << "from = " << contig[snpRecord.pos] << "\n";
                seqan::Dna5 from = contig[snpRecord.pos];
                std::pair<int, int> pos = snpRecord.getPos();

                // Get the value of each haplotype at the position.
                seqan::String<bool> inTos;
                resize(inTos, 4, false);
                seqan::Dna5String tos;
                resize(tos, options.numHaplotypes, from);
                do
                {
                    SEQAN_ASSERT(snpRecord.to != from);
                    tos[snpRecord.haplotype] = snpRecord.to;
                    inTos[ordValue(seqan::Dna5(snpRecord.to))] = true;

                    if (snpsIdx >= length(variants.snps))
                        snpRecord.rId = seqan::maxValue<int>();
                    else
                        snpRecord = variants.snps[snpsIdx++];
                }
                while (snpRecord.rId != seqan::maxValue<int>() &&
                       snpsIdx < length(variants.snps) &&
                       snpRecord.getPos() == pos);

                // Create VCF vcfRecord.
                seqan::VcfRecord vcfRecord;
                vcfRecord.chromId = rId;
                vcfRecord.pos = pos.second;
                // TODO(holtgrew): Generate an id?
                appendValue(vcfRecord.ref, from);
                for (unsigned i = 0; i < 4; ++i)
                {
                    if (!inTos[i])
                        continue;  // no ALT
                    if (!empty(vcfRecord.alt))
                        appendValue(vcfRecord.alt, ',');
                    appendValue(vcfRecord.alt, seqan::Dna5(i));
                }
                vcfRecord.filter = "PASS";
                vcfRecord.info = ".";
                // Build genotype infos.
                appendValue(vcfRecord.genotypeInfos, "");
                for (int hId = 0; hId < options.numHaplotypes; ++hId)
                {
                    if (!empty(vcfRecord.genotypeInfos[0]))
                        appendValue(vcfRecord.genotypeInfos[0], '|');
                    if (tos[hId] == vcfRecord.ref[0])
                    {
                        appendValue(vcfRecord.genotypeInfos[0], '0');
                    }
                    else
                    {
                        char buffer[20];
                        for (unsigned i = 0; i < length(vcfRecord.alt); i += 2)
                            if (tos[hId] == vcfRecord.alt[i])
                            {
                                snprintf(buffer, 19, "%d", 1 + i / 2);
                                append(vcfRecord.genotypeInfos[0], buffer);
                            }
                    }
                }

                // Write out VCF record.
                if (writeRecord(vcfStream, vcfRecord) != 0)
                {
                    std::cerr << "ERROR: Problem writing to " << options.vcfOutFile << "\n";
                    return 1;
                }
            }
            else  // process small indel records
            {
                SEQAN_ASSERT_NEQ(smallIndelRecord.pos, 0);   // Not simulated, VCF complexer.
            
                // Collect small indel records at the same position.
                seqan::String<SmallIndelRecord> records;
                do
                {
                    if (options.verbosity >= 2)
                        std::cerr << "INDEL\t"
                                  << smallIndelRecord.haplotype << "\t"
                                  << smallIndelRecord.rId << "\t"
                                  << smallIndelRecord.pos << "\t"
                                  << smallIndelRecord.size << "\t"
                                  << smallIndelRecord.seq << "\n";
                    appendValue(records, smallIndelRecord);

                    if (smallIndelIdx >= length(variants.smallIndels))
                        smallIndelRecord.rId = seqan::maxValue<int>();
                    else
                        smallIndelRecord = variants.smallIndels[smallIndelIdx++];
                }
                while (smallIndelRecord.rId != seqan::maxValue<int>() &&
                       smallIndelIdx < length(variants.smallIndels) &&
                       smallIndelRecord.getPos() == variants.smallIndels[smallIndelIdx].getPos());
                SEQAN_ASSERT_NOT(empty(records));


                // Create VCF record.
                seqan::VcfRecord vcfRecord;
                vcfRecord.chromId = rId;
                vcfRecord.pos = front(records).pos;
                // TODO(holtgrew): Generate an id?
                vcfRecord.filter = "PASS";
                vcfRecord.info = ".";
                // Build genotype infos.

                // Compute the number of bases in the REF column (1 in case of insertion and (k + 1) in the case of a
                // deletion of length k.
                int numRef = 0;
                for (unsigned i = 0; i < length(records); ++i)
                {
                    SEQAN_ASSERT_NEQ(records[i].size, 0);
                    if (records[i].size > 0)
                        numRef = std::max(numRef, 1);  // assign 1 if 0
                    else  // if (records[i].size < 0)
                        numRef = std::max(numRef, 1 - records[i].size);
                }
                append(vcfRecord.ref, infix(contig, vcfRecord.pos - 1, vcfRecord.pos - 1 + numRef));

                // Compute ALT columns and a map to the ALT.
                seqan::String<int> toIds;
                resize(toIds, options.numHaplotypes, 0);
                for (unsigned i = 0; i < length(records); ++i)
                {
                    if (i > 0)
                        appendValue(vcfRecord.alt, ',');
                    toIds[records[i].haplotype] = i + 1;
                    if (records[i].size > 0)  // insertion
                    {
                        appendValue(vcfRecord.alt, vcfRecord.ref[0]);
                        append(vcfRecord.alt, records[i].seq);
                        append(vcfRecord.alt, suffix(vcfRecord.ref, 1));
                    }
                    else  // deletion
                    {
                        appendValue(vcfRecord.alt, vcfRecord.ref[0]);
                        append(vcfRecord.alt, suffix(vcfRecord.ref, 1 - records[i].size));
                    }
                }

                // Create genotype infos.
                appendValue(vcfRecord.genotypeInfos, "");
                for (int i = 0; i < options.numHaplotypes; ++i)
                {
                    if (i > 0)
                        appendValue(vcfRecord.genotypeInfos[0], '|');
                    char buffer[20];
                    snprintf(buffer, 19, "%d", toIds[i]);
                    append(vcfRecord.genotypeInfos[0], buffer);
                }

                // Write out VCF record.
                if (writeRecord(vcfStream, vcfRecord) != 0)
                {
                    std::cerr << "ERROR: Problem writing to " << options.vcfOutFile << "\n";
                    return 1;
                }
            }
        }

        return 0;
    }
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
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-if\\fP \\fIIN.fa\\fP \\fB-ov\\fP \\fIOUT.vcf\\fP [\\fB-of\\fP \\fIOUT.fa\\fP]");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-if\\fP \\fIIN.fa\\fP \\fB-iv\\fP \\fIIN.vcf\\fP \\fB-of\\fP \\fIOUT.fa\\fP");
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

    addOption(parser, seqan::ArgParseOption("if", "in-fasta", "FASTA file with reference.",
                                            seqan::ArgParseOption::INPUTFILE, "FASTA"));
    setValidValues(parser, "in-fasta", "fasta fa");
    setRequired(parser, "in-fasta");

    addOption(parser, seqan::ArgParseOption("it", "in-variant-tsv",
                                            "TSV file with variants to simulate.  See Section on the Variant TSV File below.",
                                            seqan::ArgParseOption::INPUTFILE, "VCF"));
    setValidValues(parser, "in-vcf", "vcf");

    addOption(parser, seqan::ArgParseOption("ov", "out-vcf", "VCF file to write simulated variations to.",
                                            seqan::ArgParseOption::INPUTFILE, "VCF"));
    setValidValues(parser, "out-vcf", "vcf");

    addOption(parser, seqan::ArgParseOption("of", "out-fasta", "FASTA file to write simulated haplotypes to.",
                                            seqan::ArgParseOption::INPUTFILE, "FASTA"));
    setValidValues(parser, "out-fasta", "fasta fa");

    addOption(parser, seqan::ArgParseOption("", "haplotype-name-sep", "Haplotype name separator in output FASTA.",
                                            seqan::ArgParseOption::STRING, "SEP"));
    setDefaultValue(parser, "haplotype-name-sep", "/");

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
    setDefaultValue(parser, "snp-rate", "0.0001");

    addOption(parser, seqan::ArgParseOption("", "small-indel-rate", "Small indel rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "small-indel-rate", "0.0");
    setMaxValue(parser, "small-indel-rate", "1.0");
    setDefaultValue(parser, "small-indel-rate", "0.00005");

    addOption(parser, seqan::ArgParseOption("", "min-small-indel-size", "Minimal small indel size to simulate.",
                                            seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "min-small-indel-size", "0");
    setDefaultValue(parser, "min-small-indel-size", "1");

    addOption(parser, seqan::ArgParseOption("", "max-small-indel-size", "Maximal small indel size to simulate.",
                                            seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "max-small-indel-size", "0");
    setDefaultValue(parser, "max-small-indel-size", "6");

    addOption(parser, seqan::ArgParseOption("", "sv-indel-rate", "Per-base SNP rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "sv-indel-rate", "0.0");
    setMaxValue(parser, "sv-indel-rate", "1.0");
    setDefaultValue(parser, "sv-indel-rate", "0.0");

    addOption(parser, seqan::ArgParseOption("", "sv-inversion-rate", "Per-base SNP rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "sv-inversion-rate", "0.0");
    setMaxValue(parser, "sv-inversion-rate", "1.0");
    setDefaultValue(parser, "sv-inversion-rate", "0.0");

    addOption(parser, seqan::ArgParseOption("", "sv-translocation-rate", "Per-base SNP rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "sv-translocation-rate", "0.0");
    setMaxValue(parser, "sv-translocation-rate", "1.0");
    setDefaultValue(parser, "sv-translocation-rate", "0.0");

    addOption(parser, seqan::ArgParseOption("", "sv-duplication-rate", "Per-base SNP rate.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "sv-duplication-rate", "0.0");
    setMaxValue(parser, "sv-duplication-rate", "1.0");
    setDefaultValue(parser, "sv-duplication-rate", "0.0");

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

    addText(parser,
            "SNPs and small indels are simulated such that at each position, a random experiment is "
            "performed whether to simulate either variation.  In case both variations are to be simulated, "
            "the experiment is repeated.");

    addText(parser, "The indel and SV sizes are picked uniformly at random from the argument size intervals.");

    addText(parser,
            "The simulation of haplotypes works as follows.  For small indels, the indel is placed into "
            "one of the haplotypes that are to be simulated.  The exact haplotype is picked uniformly at "
            "random.  For SNPs, we simulate a random base for each haplotype.  For at least one haplotype, "
            "the base has to be different from the reference or the experiment is repeated.");

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
    options.verbosity = 1;
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValue(options.seed, parser, "seed");

    getOptionValue(options.vcfInFile, parser, "in-vcf");
    getOptionValue(options.fastaInFile, parser, "in-fasta");
    getOptionValue(options.vcfOutFile, parser, "out-vcf");
    getOptionValue(options.fastaOutFile, parser, "out-fasta");
    getOptionValue(options.inputSVSizeFile, parser, "in-variant-tsv");
    getOptionValue(options.haplotypeNameSep, parser, "haplotype-name-sep");

    getOptionValue(options.numHaplotypes, parser, "num-haplotypes");
    getOptionValue(options.haplotypeSep, parser, "haplotype-sep");
    getOptionValue(options.snpRate, parser, "snp-rate");
    getOptionValue(options.smallIndelRate, parser, "small-indel-rate");
    getOptionValue(options.minSmallIndelSize, parser, "min-small-indel-size");
    getOptionValue(options.maxSmallIndelSize, parser, "max-small-indel-size");
    getOptionValue(options.svIndelRate, parser, "sv-indel-rate");
    getOptionValue(options.svInversionRate, parser, "sv-inversion-rate");
    getOptionValue(options.svTranslocationRate, parser, "sv-translocation-rate");
    getOptionValue(options.svDuplicationRate, parser, "sv-duplication-rate");
    getOptionValue(options.minSVSize, parser, "min-sv-size");
    getOptionValue(options.maxSVSize, parser, "max-sv-size");

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

    // Initialize random number generator.
    TRng rng(options.seed);

    std::cerr << "MASON VARIATOR\n"
              << "==============\n\n";

    print(std::cerr, options);

    std::cerr << "\n__PREPARATION_________________________________________________________________\n"
              << "\n";

    std::cerr << "Loading Reference Index " << options.fastaInFile << " ...";
    seqan::FaiIndex faiIndex;
    if (read(faiIndex, toCString(options.fastaInFile)) != 0)
    {
        std::cerr << " FAILED (not fatal, we can just build it)\n";
        std::cerr << "Building Index        " << options.fastaInFile << ".fai ...";
        if (build(faiIndex, toCString(options.fastaInFile)) != 0)
        {
            std::cerr << "Could not build FAI index.\n";
            return 1;
        }
        std::cerr << " OK\n";
        seqan::CharString faiPath = options.fastaInFile;
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

    std::cerr << "\n__SIMULATION__________________________________________________________________\n"
              << "\n";

    std::cerr << "SNP/Small Indel Simulation:\n";

    MasonVariatorApp app(rng, faiIndex, options);
    app.run();
    
    std::cerr << "Structural Variation Simulation ... NOT IMPLEMENTED\n\n";

    std::cerr << "__WRITING OUTPUT______________________________________________________________\n"
              << "\n";

    if (!empty(options.vcfOutFile))
        std::cerr << "Writing VCF to " << options.vcfOutFile << "\n";
    
    if (!empty(options.fastaOutFile))
        std::cerr << "Writing FASTA to " << options.fastaOutFile << "\n";

    std::cerr << "\nDONE.\n";
    
    return 0;
}
