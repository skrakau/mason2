#include "genomic_variants.h"

std::ostream & operator<<(std::ostream & out, SnpRecord const & record)
{
    out << "SnpRecord(" << record.haplotype << ", " << record.rId << ", " << record.pos
        << ", " << record.to << ")";
    return out;
}

std::ostream & operator<<(std::ostream & out, SmallIndelRecord const & record)
{
    out << "SnpRecord(" << record.haplotype << ", " << record.rId << ", " << record.pos
        << ", " << record.size << ", " << record.seq << ")";
    return out;
}

int StructuralVariantRecord::endPosition() const
{
    if (pos == -1)
        return seqan::maxValue<int>();
        
    switch (kind)
    {
        case INDEL:
            if (size > 0)
                return pos;
            else
                return pos + size;
        case INVERSION:
            return pos + size;
        case TRANSLOCATION:
            return targetPos;
        case DUPLICATION:
            return targetPos;
        default:
            return -1;
    }
}

bool StructuralVariantRecord::nearBreakend(int query) const
{
    if (pos == -1)
        return false;  // invalid/sentinel has no breakends
        
    switch (kind)
    {
        case INDEL:
            if (size > 0)
                return (query == pos || query == pos + 1);
            else
                return (query == pos || query == pos + 1 ||
                        query == pos - size || query == pos - size + 1);
        case INVERSION:
            return (query == pos || query == pos + 1 ||
                    query == pos + size || query == pos + size + 1);
        case TRANSLOCATION:
            return (query == pos || query == pos + 1 ||
                    query == targetPos || query == targetPos + 1);
        case DUPLICATION:
            return (query == pos || query == pos + 1 ||
                    query == targetPos || query == targetPos + 1);
        default:
            return false;
    }
}

std::ostream & operator<<(std::ostream & out, StructuralVariantRecord const & record)
{
    char const * kind;
    switch (record.kind)
    {
        case StructuralVariantRecord::TRANSLOCATION:
            kind = "TRANSLOCATION";
            break;
        case StructuralVariantRecord::INDEL:
            kind = "INDEL";
            break;
        case StructuralVariantRecord::INVERSION:
            kind = "INVERSION";
            break;
        case StructuralVariantRecord::DUPLICATION:
            kind = "DUPLICATION";
            break;
        default:
            kind = "INVALID";
            break;
    }
    out << "StructuralVariantRecord(kind=" << kind << ", haplotype=" << record.haplotype
        << ", rId=" << record.rId << ", pos=" << record.pos << ", size=" << record.size
        << ", targetRId=" << record.targetRId << ", targetPos=" << record.targetPos
        << ", seq=\"" << record.seq << "\")";
    return out;
}

// ----------------------------------------------------------------------------
// Function VariantMaterializer::_runImpl()
// ----------------------------------------------------------------------------

int VariantMaterializer::_runImpl(
        seqan::Dna5String * resultSeq,
        MethylationLevels * resultLvls,
        std::vector<int> & breakpoints,
        seqan::Dna5String const * ref,
        MethylationLevels const * refLvls,
        int haplotypeId)
{
    breakpoints.clear();
    clear(*resultSeq);
    if (resultLvls)
        resultLvls->clear();
    
    // Apply small variants.  We get a sequence with the small variants and a journal of the difference to contig.
    seqan::Dna5String seqSmallVariants;
    TJournalEntries journal;
    MethylationLevels levelsSmallVariants;  // only used if revLevels != 0
    MethylationLevels * smallLvlsPtr = refLvls ? &levelsSmallVariants : 0;
    if (_materializeSmallVariants(seqSmallVariants, journal, smallLvlsPtr, *ref, *variants,
                                  refLvls, haplotypeId) != 0)
        return 1;
    
    // Apply structural variants.
    if (_materializeLargeVariants(*resultSeq, resultLvls, breakpoints, journal, seqSmallVariants, *variants,
                                  smallLvlsPtr, haplotypeId) != 0)
        return 1;

    return 0;
}

// ----------------------------------------------------------------------------
// Function VariantMaterializer::_materializeSmallVariants()
// ----------------------------------------------------------------------------

int VariantMaterializer::_materializeSmallVariants(
        seqan::Dna5String & seq,
        TJournalEntries & journal,
        MethylationLevels * levelsSmallVariants,
        seqan::Dna5String const & contig,
        Variants const & variants,
        MethylationLevels const * levels,
        int hId)
{
    if (methSimOptions)
    {
        SEQAN_ASSERT_EQ(methSimOptions->simulateMethylationLevels, (levelsSmallVariants != 0));
        SEQAN_ASSERT_EQ(methSimOptions->simulateMethylationLevels, (levels != 0));
    }
    
    // Clear journal and output methylation levels.
    reinit(journal, length(contig));
    if (levelsSmallVariants)
        levelsSmallVariants->clear();
    // Store variation points with a flag whether it is a SNP (true) or a breakpoint (false).
    seqan::String<std::pair<int, bool> > varPoints;

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
    // Track last position from contig appended to seq so far.
    int lastPos = 0;
    if (verbosity >= 3)
        std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";

    // TODO(holtgrew): Extract contig building into their own functions.
    if (verbosity >= 2)
        std::cerr << "building output\n";
    while (snpRecord.rId != seqan::maxValue<int>() || smallIndelRecord.rId != seqan::maxValue<int>())
    {
        // TODO(holtgrew): Extract SNP and small indel handling in functions.
        if (snpRecord.getPos() < smallIndelRecord.getPos())  // process SNP records
        {
            if (snpRecord.haplotype == hId)  // Ignore all but the current contig.
            {
                if (verbosity >= 3)
                    std::cerr << "append(seq, infix(contig, " << lastPos << ", " << snpRecord.pos << ") " << __LINE__ << "\n";
                // Append interim sequence and methylation levels->
                append(seq, infix(contig, lastPos, snpRecord.pos));
                if (methSimOptions && methSimOptions->simulateMethylationLevels)
                {
                    append(levelsSmallVariants->forward, infix(levels->forward, lastPos, snpRecord.pos + 1));
                    append(levelsSmallVariants->reverse, infix(levels->reverse, lastPos, snpRecord.pos + 1));
                    appendValue(varPoints, std::make_pair(length(seq), true));      // variation points before/after SNP
                    appendValue(varPoints, std::make_pair(length(seq) + 1, true));
                }

                SEQAN_ASSERT_GEQ(snpRecord.pos, lastPos);
                if (verbosity >= 3)
                    std::cerr << "appendValue(seq, " << snpRecord.to << "')\n";
                appendValue(seq, snpRecord.to);
                lastPos = snpRecord.pos + 1;
                if (verbosity >= 3)
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
                    if (verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << lastPos << ", " << smallIndelRecord.pos << ") " << __LINE__ << "\n";

                    // Simulate methylation levels for insertion.
                    MethylationLevels lvls;
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        MethylationLevelSimulator methSim(*rng, *methSimOptions);
                        methSim.run(lvls, smallIndelRecord.seq);
                    }

                    // Append interim sequence and methylation levels->
                    append(seq, infix(contig, lastPos, smallIndelRecord.pos));
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        append(levelsSmallVariants->forward, infix(levels->forward, lastPos, smallIndelRecord.pos));
                        append(levelsSmallVariants->reverse, infix(levels->reverse, lastPos, smallIndelRecord.pos));
                        appendValue(varPoints, std::make_pair(length(seq), false));  // variation point before insertion
                    }

                    SEQAN_ASSERT_GEQ(smallIndelRecord.pos, lastPos);
                    if (verbosity >= 3)
                        std::cerr << "append(seq, \"" << smallIndelRecord.seq << "\") " << __LINE__ << "\n";
                    // Append novel sequence and methylation levels->
                    append(seq, smallIndelRecord.seq);
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        append(levelsSmallVariants->forward, lvls.forward);
                        append(levelsSmallVariants->reverse, lvls.reverse);
                        appendValue(varPoints, std::make_pair(length(seq), false));  // variation point after insertion
                    }
                    lastPos = smallIndelRecord.pos;
                    recordInsertion(journal, hostToVirtualPosition(journal, smallIndelRecord.pos),
                                    0, smallIndelRecord.size);
                    if (verbosity >= 3)
                        std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";
                }
                else  // deletion
                {
                    if (verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << lastPos << ", " << smallIndelRecord.pos << ") " << __LINE__ << "\n";
                    // Append interim sequence and methylation levels->
                    append(seq, infix(contig, lastPos, smallIndelRecord.pos));  // interim chars
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair(length(seq), false));  // variation point at deletion
                        append(levelsSmallVariants->forward, infix(levels->forward, lastPos, smallIndelRecord.pos));
                        append(levelsSmallVariants->reverse, infix(levels->reverse, lastPos, smallIndelRecord.pos));
                    }

                    lastPos = smallIndelRecord.pos - smallIndelRecord.size;
                    recordErase(journal,
                                hostToVirtualPosition(journal, smallIndelRecord.pos),
                                hostToVirtualPosition(journal, smallIndelRecord.pos - smallIndelRecord.size));
                    if (verbosity >= 3)
                        std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";
                }
            }

            if (smallIndelIdx >= length(variants.smallIndels))
                smallIndelRecord.rId = seqan::maxValue<int>();
            else
                smallIndelRecord = variants.smallIndels[smallIndelIdx++];
        }
    }
    // Insert remaining characters.
    if (verbosity >= 3)
        std::cerr << "append(seq, infix(contig, " << lastPos << ", " << length(contig) << ")\n";
    append(seq, infix(contig, lastPos, length(contig)));

    if (methSimOptions && methSimOptions->simulateMethylationLevels)
    {
        append(levelsSmallVariants->forward, infix(levels->forward, lastPos, length(contig)));
        append(levelsSmallVariants->reverse, infix(levels->reverse, lastPos, length(contig)));

        SEQAN_ASSERT_EQ(length(seq), length(levelsSmallVariants->forward));
        SEQAN_ASSERT_EQ(length(seq), length(levelsSmallVariants->reverse));

        fixVariationLevels(*levelsSmallVariants, *rng, seq, varPoints, *methSimOptions);
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function VariantMaterializer::_materializeLargeVariants()
// ----------------------------------------------------------------------------

int VariantMaterializer::_materializeLargeVariants(
        seqan::Dna5String & seq,
        MethylationLevels * levelsLargeVariants,
        std::vector<int> & breakpoints,
        TJournalEntries const & journal,
        seqan::Dna5String const & contig,
        Variants const & variants,
        MethylationLevels const * levels,
        int hId)
{
    if (methSimOptions)
    {
        SEQAN_ASSERT_EQ(methSimOptions->simulateMethylationLevels, (levelsLargeVariants != 0));
        SEQAN_ASSERT_EQ(methSimOptions->simulateMethylationLevels, (levels != 0));
    }

    // Clear output methylation levels->
    if (levelsLargeVariants)
        levelsLargeVariants->clear();
    // Store variation points.  We reuse the fixVariationLevels() function from small indel/snp simulation and thus
    // have to store a bool that is always set to false.
    seqan::String<std::pair<int, bool> > varPoints;
        
    // Track last position from contig appended to seq so far.
    int lastPos = 0;
    if (verbosity >= 3)
        std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";

    // Number of bytes written out so far/current position in variant.
    int currentPos = 0;

    for (unsigned i = 0; i < length(variants.svRecords); ++i)
    {
        if (variants.svRecords[i].haplotype != hId)  // Ignore all but the current contig.
            continue;
        // We obtain a copy of the current SV record since we translate its positions below.
        StructuralVariantRecord svRecord = variants.svRecords[i];

        // Translate positions and lengths of SV record.
        if (verbosity >= 2)
            std::cerr << "  Translating SvRecord\n  " << svRecord << '\n';
        svRecord.pos = hostToVirtualPosition(journal, svRecord.pos);
        SEQAN_ASSERT_LT(svRecord.pos, (int)length(contig));
        // We do not need to adjust the sizes for insertions.
        if (svRecord.kind != StructuralVariantRecord::INDEL || svRecord.size < 0)
            svRecord.size = hostToVirtualPosition(journal, svRecord.pos + svRecord.size) -
                    hostToVirtualPosition(journal, svRecord.pos);
        if (svRecord.targetPos != -1)
            svRecord.targetPos = hostToVirtualPosition(journal, svRecord.targetPos);
        if (verbosity >= 2)
            std::cerr << "  => " << svRecord << '\n';

        // Copy from contig to seq with SVs.
        if (verbosity >= 3)
            std::cerr << "lastPos == " << lastPos << "\n";
        append(seq, infix(contig, lastPos, svRecord.pos));  // interim chars
        if (methSimOptions && methSimOptions->simulateMethylationLevels)
        {
            append(levelsLargeVariants->forward, infix(levels->forward, lastPos, svRecord.pos));
            append(levelsLargeVariants->reverse, infix(levels->reverse, lastPos, svRecord.pos));
            appendValue(varPoints, std::make_pair(length(seq), false));
        }
        currentPos = length(seq);
        if (verbosity >= 3)
            std::cerr << "append(seq, infix(contig, " << lastPos << ", " << svRecord.pos << ") " << __LINE__ << " (interim)\n";
        switch (svRecord.kind)
        {
            case StructuralVariantRecord::INDEL:
                {
                    if (svRecord.size > 0)  // insertion
                    {
                        SEQAN_ASSERT_EQ((int)length(svRecord.seq), svRecord.size);

                        // Simulate methylation levels for insertion.
                        MethylationLevels lvls;
                        if (methSimOptions && methSimOptions->simulateMethylationLevels)
                        {
                            MethylationLevelSimulator methSim(*rng, *methSimOptions);
                            methSim.run(lvls, svRecord.seq);
                        }

                        // Append novel sequence and methylation levels->
                        append(seq, svRecord.seq);
                        if (methSimOptions && methSimOptions->simulateMethylationLevels)
                        {
                            append(levelsLargeVariants->forward, lvls.forward);
                            append(levelsLargeVariants->reverse, lvls.reverse);
                            appendValue(varPoints, std::make_pair(length(seq), false));  // variation point after insertion
                        }
                        if (verbosity >= 3)
                            std::cerr << "append(seq, svRecord.seq (length == " << length(svRecord.seq) << ") " << __LINE__ << " (insertion)\n";
                        lastPos = svRecord.pos;
                        SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                        // Copy out breakpoints.
                        breakpoints.push_back(currentPos);
                        breakpoints.push_back(length(seq));

                        currentPos = length(seq);
                    }
                    else  // deletion
                    {
                        lastPos = svRecord.pos - svRecord.size;
                        SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                        // Copy out breakpoint.
                        breakpoints.push_back(currentPos);
                    }
                }
                break;
            case StructuralVariantRecord::INVERSION:
                {
                    unsigned oldLen = length(seq);
                    append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair(length(seq), false));  // variation point at deletion
                        append(levelsLargeVariants->forward, infix(levels->reverse, svRecord.pos, svRecord.pos + svRecord.size));
                        reverse(infix(levelsLargeVariants->forward, oldLen, length(levelsLargeVariants->forward)));
                        append(levelsLargeVariants->reverse, infix(levels->forward, svRecord.pos, svRecord.pos + svRecord.size));
                        reverse(infix(levelsLargeVariants->reverse, oldLen, length(levelsLargeVariants->reverse)));
                    }

                    if (verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << " (inversion)\n";
                    reverseComplement(infix(seq, oldLen, length(seq)));
                    lastPos = svRecord.pos + svRecord.size;
                    SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                    // Copy out breakpoints.
                    breakpoints.push_back(currentPos);
                    breakpoints.push_back(length(seq));

                    currentPos = length(seq);
                }
                break;
            case StructuralVariantRecord::TRANSLOCATION:
                {
                    SEQAN_ASSERT_GEQ(svRecord.targetPos, svRecord.pos + svRecord.size);
                    append(seq, infix(contig, svRecord.pos + svRecord.size, svRecord.targetPos));
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair(length(seq), false));
                        append(levelsLargeVariants->forward, infix(levels->forward, svRecord.pos + svRecord.size, svRecord.targetPos));
                        append(levelsLargeVariants->reverse, infix(levels->reverse, svRecord.pos + svRecord.size, svRecord.targetPos));
                    }
                    append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair(length(seq), false));
                        append(levelsLargeVariants->forward, infix(levels->forward, svRecord.pos, svRecord.pos + svRecord.size));
                        append(levelsLargeVariants->reverse, infix(levels->reverse, svRecord.pos, svRecord.pos + svRecord.size));
                    }
                    if (verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << svRecord.pos + svRecord.size << ", " << svRecord.targetPos << ") " << __LINE__ << " (translocation)\n"
                                  << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << "\n";
                    lastPos = svRecord.targetPos;
                    SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                    // Copy out breakpoints.
                    breakpoints.push_back(currentPos);
                    breakpoints.push_back(currentPos + svRecord.targetPos - svRecord.pos - svRecord.size);
                    breakpoints.push_back(length(seq));

                    currentPos = length(seq);
                }
                break;
            case StructuralVariantRecord::DUPLICATION:
                {
                    append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                    SEQAN_ASSERT_GEQ(svRecord.targetPos, svRecord.pos + svRecord.size);
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair(length(seq), false));
                        append(levelsLargeVariants->forward, infix(levels->forward, svRecord.pos, svRecord.pos + svRecord.size));
                        append(levelsLargeVariants->reverse, infix(levels->reverse, svRecord.pos, svRecord.pos + svRecord.size));
                    }
                    append(seq, infix(contig, svRecord.pos + svRecord.size, svRecord.targetPos));
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair(length(seq), false));
                        append(levelsLargeVariants->forward, infix(levels->forward, svRecord.pos + svRecord.size, svRecord.targetPos));
                        append(levelsLargeVariants->reverse, infix(levels->reverse, svRecord.pos + svRecord.size, svRecord.targetPos));
                    }
                    append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair(length(seq), false));
                        append(levelsLargeVariants->forward, infix(levels->forward, svRecord.pos, svRecord.pos + svRecord.size));
                        append(levelsLargeVariants->reverse, infix(levels->reverse, svRecord.pos, svRecord.pos + svRecord.size));
                    }
                    if (verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << " (duplication)\n"
                                  << "append(seq, infix(contig, " << svRecord.pos + svRecord.size << ", " << svRecord.targetPos << ") " << __LINE__ << "\n"
                                  << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << "\n";
                    lastPos = svRecord.targetPos;
                    SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                    // Copy out breakpoints.
                    breakpoints.push_back(currentPos);
                    breakpoints.push_back(currentPos + svRecord.pos + svRecord.size - svRecord.pos);
                    breakpoints.push_back(currentPos + svRecord.pos + svRecord.size - svRecord.pos + svRecord.targetPos - (svRecord.pos + svRecord.size));
                    breakpoints.push_back(length(seq));

                    currentPos = length(seq);
                }
                break;
            default:
                return 1;
        }
    }
    if (verbosity >= 3)
        std::cerr << "append(seq, infix(contig, " << lastPos << ", " << length(contig) << ") "
                  << __LINE__ << " (last interim)\n";
    append(seq, infix(contig, lastPos, length(contig)));
    if (methSimOptions && methSimOptions->simulateMethylationLevels)
    {
        append(levelsLargeVariants->forward, infix(levels->forward, lastPos, length(contig)));
        append(levelsLargeVariants->reverse, infix(levels->reverse, lastPos, length(contig)));

        SEQAN_ASSERT_EQ(length(seq), length(levelsLargeVariants->forward));
        SEQAN_ASSERT_EQ(length(seq), length(levelsLargeVariants->reverse));

        fixVariationLevels(*levelsLargeVariants, *rng, seq, varPoints, *methSimOptions);
    }
        
    return 0;
}

// --------------------------------------------------------------------------
// Function PositionMap::overlapsWithVariant()
// --------------------------------------------------------------------------

bool PositionMap::overlapsWithVariant(int svBeginPos, int svEndPos) const
{
    seqan::String<GenomicInterval> intervals;
    findIntervals(svIntervalTree, svBeginPos, svEndPos, intervals);
    return !empty(intervals);
}

// --------------------------------------------------------------------------
// Function PositionMap::getGenomicInterval()
// --------------------------------------------------------------------------

GenomicInterval PositionMap::getGenomicInterval(int svPos) const
{
    seqan::String<GenomicInterval> intervals;
    findIntervals(svIntervalTree, svPos, intervals);
    SEQAN_ASSERT_EQ(length(intervals), 1u);
    return intervals[0];
}

// --------------------------------------------------------------------------
// Function PositionMap::toSmallVarInterval()
// --------------------------------------------------------------------------

std::pair<int, int> PositionMap::toSmallVarInterval(int svBeginPos, int svEndPos) const
{
    SEQAN_ASSERT(!overlapsWithVariant(svBeginPos, svEndPos));
    GenomicInterval gi = getGenomicInterval(svBeginPos);
    if (gi.kind == GenomicInterval::INSERTED)
    {
        // novel sequence, cannot be projected
        return std::make_pair(-1, -1);
    }
    if (gi.kind != GenomicInterval::INVERTED)
    {
        // forward
        return std::make_pair(gi.smallVarBeginPos + (svBeginPos - gi.svBeginPos),
                              gi.smallVarBeginPos + (svEndPos - gi.svBeginPos));
    }
    else
    {
        // reverse
        return std::make_pair(gi.smallVarBeginPos + (gi.svEndPos - svBeginPos),
                              gi.smallVarBeginPos + (gi.svEndPos - svEndPos));
    }

    return std::make_pair(-1, -1);  // cannot reach here
}

// --------------------------------------------------------------------------
// Function PositionMap::toOriginalInterval()
// --------------------------------------------------------------------------

std::pair<int, int> PositionMap::toOriginalInterval(int smallVarBeginPos, int smallVarEndPos) const
{
    // TODO(holtgrew): Project to the left of gaps as documented.
    int refBeginPos = hostToVirtualPosition(smallVariantJournal, smallVarBeginPos);
    int refEndPos = hostToVirtualPosition(smallVariantJournal, smallVarEndPos);
    return std::make_pair(refBeginPos, refEndPos);
}

// --------------------------------------------------------------------------
// Function PositionMap::reinit()
// --------------------------------------------------------------------------

void PositionMap::reinit(unsigned contigLength)
{
    // Reset the journal.
    seqan::reinit(smallVariantJournal, contigLength);
    // Reset the interval tree.
    // TODO(holtgrew): Better API support for IntervalTree?
    svIntervalTree = TIntervalTree();
}
