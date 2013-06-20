#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>

#include "sequencing.h"
#include "genomic_variants.h"

SEQAN_DEFINE_TEST(mason_tests_append_orientation_elementary_operations)
{
    // Below, we test for match, mismatch, insertion, deletion and insertion.
    {
        TCigarString cigar;
        
        std::pair<int, int> v = appendOperation(cigar, 'M');
        
        SEQAN_ASSERT_EQ(v.first, 1);
        SEQAN_ASSERT_EQ(v.second, 1);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'M');
        SEQAN_ASSERT_EQ(cigar[0].count, 1u);
    }
    {
        TCigarString cigar;
        
        std::pair<int, int> v = appendOperation(cigar, 'X');
        
        SEQAN_ASSERT_EQ(v.first, 1);
        SEQAN_ASSERT_EQ(v.second, 1);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'X');
        SEQAN_ASSERT_EQ(cigar[0].count, 1u);
    }
    {
        TCigarString cigar;
        
        std::pair<int, int> v = appendOperation(cigar, 'I');
        
        SEQAN_ASSERT_EQ(v.first, 1);
        SEQAN_ASSERT_EQ(v.second, 0);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'I');
        SEQAN_ASSERT_EQ(cigar[0].count, 1u);
    }
    {
        TCigarString cigar;
        
        std::pair<int, int> v = appendOperation(cigar, 'D');
        
        SEQAN_ASSERT_EQ(v.first, 0);
        SEQAN_ASSERT_EQ(v.second, 1);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'D');
        SEQAN_ASSERT_EQ(cigar[0].count, 1u);
    }
}

SEQAN_DEFINE_TEST(mason_tests_append_orientation_combination)
{
    // Test with the combination of equal operations.
    {
        TCigarString cigar;

        appendValue(cigar, seqan::CigarElement<>('M', 1));
        std::pair<int, int> v = appendOperation(cigar, 'M');
        
        SEQAN_ASSERT_EQ(v.first, 1);
        SEQAN_ASSERT_EQ(v.second, 1);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'M');
        SEQAN_ASSERT_EQ(cigar[0].count, 2u);
    }
    {
        TCigarString cigar;
        
        appendValue(cigar, seqan::CigarElement<>('X', 1));
        std::pair<int, int> v = appendOperation(cigar, 'X');
        
        SEQAN_ASSERT_EQ(v.first, 1);
        SEQAN_ASSERT_EQ(v.second, 1);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'X');
        SEQAN_ASSERT_EQ(cigar[0].count, 2u);
    }
    {
        TCigarString cigar;
        
        appendValue(cigar, seqan::CigarElement<>('I', 1));
        std::pair<int, int> v = appendOperation(cigar, 'I');
        
        SEQAN_ASSERT_EQ(v.first, 1);
        SEQAN_ASSERT_EQ(v.second, 0);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'I');
        SEQAN_ASSERT_EQ(cigar[0].count, 2u);
    }
    {
        TCigarString cigar;
        
        appendValue(cigar, seqan::CigarElement<>('D', 1));
        std::pair<int, int> v = appendOperation(cigar, 'D');
        
        SEQAN_ASSERT_EQ(v.first, 0);
        SEQAN_ASSERT_EQ(v.second, 1);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'D');
        SEQAN_ASSERT_EQ(cigar[0].count, 2u);
    }
}

SEQAN_DEFINE_TEST(mason_tests_append_orientation_canceling_out)
{
    // Test with the combination of operations that cancel each other out (I/D, D/I)
    {
        TCigarString cigar;

        appendValue(cigar, seqan::CigarElement<>('I', 1));
        std::pair<int, int> v = appendOperation(cigar, 'D');
        
        SEQAN_ASSERT_EQ(v.first, -1);
        SEQAN_ASSERT_EQ(v.second, 0);
        SEQAN_ASSERT_EQ(length(cigar), 0u);
    }
    {
        TCigarString cigar;
        
        appendValue(cigar, seqan::CigarElement<>('D', 1));
        std::pair<int, int> v = appendOperation(cigar, 'I');
        
        SEQAN_ASSERT_EQ(v.first, 0);
        SEQAN_ASSERT_EQ(v.second, -1);
        SEQAN_ASSERT_EQ(length(cigar), 0u);
    }
}

SEQAN_DEFINE_TEST(mason_tests_position_map_inversion)
{
    typedef PositionMap::TInterval TInterval;

    PositionMap positionMap;

    // Inversion: --1000-->|<--1000--|--1000-->
    GenomicInterval gi1(   0, 1000,    0, 1000, '+', GenomicInterval::NORMAL);
    GenomicInterval gi2(1000, 2000, 1000, 2000, '-', GenomicInterval::INVERTED);
    GenomicInterval gi3(2000, 3000, 2000, 3000, '+', GenomicInterval::NORMAL);
    
    // Build interval tree.
    seqan::String<TInterval> intervals;
    appendValue(intervals, TInterval(gi1.svBeginPos, gi1.svEndPos, gi1));
    appendValue(intervals, TInterval(gi2.svBeginPos, gi2.svEndPos, gi2));
    appendValue(intervals, TInterval(gi3.svBeginPos, gi3.svEndPos, gi3));
    createIntervalTree(positionMap.svIntervalTree, intervals);

    // Add breakpoints.
    positionMap.svBreakpoints.insert(0);
    positionMap.svBreakpoints.insert(gi1.svEndPos);
    positionMap.svBreakpoints.insert(gi2.svEndPos);
    positionMap.svBreakpoints.insert(gi3.svEndPos);

    // Tests for overlapsWithBreakpoint()
    SEQAN_ASSERT_NOT(positionMap.overlapsWithBreakpoint(   0, 1000));
    SEQAN_ASSERT_NOT(positionMap.overlapsWithBreakpoint(1000, 2000));
    SEQAN_ASSERT_NOT(positionMap.overlapsWithBreakpoint(2000, 3000));

    SEQAN_ASSERT(positionMap.overlapsWithBreakpoint( 999, 1001));
    SEQAN_ASSERT(positionMap.overlapsWithBreakpoint(1999, 2001));
    SEQAN_ASSERT(positionMap.overlapsWithBreakpoint(2999, 3001));

    // Tests for getGenomicInterval()
    SEQAN_ASSERT(gi1 == positionMap.getGenomicInterval(   0));
    SEQAN_ASSERT(gi1 == positionMap.getGenomicInterval( 999));
    SEQAN_ASSERT(gi2 == positionMap.getGenomicInterval(1000));
    SEQAN_ASSERT(gi2 == positionMap.getGenomicInterval(1999));
    SEQAN_ASSERT(gi3 == positionMap.getGenomicInterval(2000));
    SEQAN_ASSERT(gi3 == positionMap.getGenomicInterval(2999));

    // Tests for toSmallVarInterval().
    typedef std::pair<int, int> TPair;

    TPair i1 = positionMap.toSmallVarInterval(0, 100);
    SEQAN_ASSERT_EQ(i1.first, 0);
    SEQAN_ASSERT_EQ(i1.second, 100);

    TPair i2 = positionMap.toSmallVarInterval(900, 1000);
    SEQAN_ASSERT_EQ(i2.first, 900);
    SEQAN_ASSERT_EQ(i2.second, 1000);

    TPair i3 = positionMap.toSmallVarInterval(1000, 1100);
    SEQAN_ASSERT_EQ(i3.first, 2000);
    SEQAN_ASSERT_EQ(i3.second, 1900);

    TPair i4 = positionMap.toSmallVarInterval(1900, 2000);
    SEQAN_ASSERT_EQ(i4.first, 1100);
    SEQAN_ASSERT_EQ(i4.second, 1000);

    TPair i5 = positionMap.toSmallVarInterval(2000, 2100);
    SEQAN_ASSERT_EQ(i5.first, 2000);
    SEQAN_ASSERT_EQ(i5.second, 2100);

    TPair i6 = positionMap.toSmallVarInterval(2900, 3000);
    SEQAN_ASSERT_EQ(i6.first, 2900);
    SEQAN_ASSERT_EQ(i6.second, 3000);
}

SEQAN_DEFINE_TEST(mason_tests_position_map_translocation)
{
    typedef PositionMap::TInterval TInterval;

    PositionMap positionMap;

    // Translocation: --A--> --B-->
    //                --B--> --A-->
    GenomicInterval gi1(   0, 1000, 1000, 2000, '+', GenomicInterval::NORMAL);
    GenomicInterval gi2(1000, 2000,    0, 1000, '-', GenomicInterval::NORMAL);
    
    // Build interval tree.
    seqan::String<TInterval> intervals;
    appendValue(intervals, TInterval(gi1.svBeginPos, gi1.svEndPos, gi1));
    appendValue(intervals, TInterval(gi2.svBeginPos, gi2.svEndPos, gi2));
    createIntervalTree(positionMap.svIntervalTree, intervals);

    // Add breakpoints.
    positionMap.svBreakpoints.insert(0);
    positionMap.svBreakpoints.insert(gi1.svEndPos);
    positionMap.svBreakpoints.insert(gi2.svEndPos);

    // Tests for overlapsWithBreakpoint()
    SEQAN_ASSERT_NOT(positionMap.overlapsWithBreakpoint(   0, 1000));
    SEQAN_ASSERT_NOT(positionMap.overlapsWithBreakpoint(1000, 2000));

    SEQAN_ASSERT(positionMap.overlapsWithBreakpoint( 999, 1001));
    SEQAN_ASSERT(positionMap.overlapsWithBreakpoint(1999, 2001));

    // Tests for getGenomicInterval()
    SEQAN_ASSERT(gi1 == positionMap.getGenomicInterval(   0));
    SEQAN_ASSERT(gi1 == positionMap.getGenomicInterval( 999));
    SEQAN_ASSERT(gi2 == positionMap.getGenomicInterval(1000));
    SEQAN_ASSERT(gi2 == positionMap.getGenomicInterval(1999));

    // Tests for toSmallVarInterval().
    typedef std::pair<int, int> TPair;

    TPair i1 = positionMap.toSmallVarInterval(0, 100);
    SEQAN_ASSERT_EQ(i1.first, 1000);
    SEQAN_ASSERT_EQ(i1.second, 1100);

    TPair i2 = positionMap.toSmallVarInterval(900, 1000);
    SEQAN_ASSERT_EQ(i2.first, 1900);
    SEQAN_ASSERT_EQ(i2.second, 2000);

    TPair i3 = positionMap.toSmallVarInterval(1000, 1100);
    SEQAN_ASSERT_EQ(i3.first, 0);
    SEQAN_ASSERT_EQ(i3.second, 100);

    TPair i4 = positionMap.toSmallVarInterval(1900, 2000);
    SEQAN_ASSERT_EQ(i4.first, 900);
    SEQAN_ASSERT_EQ(i4.second, 1000);
}

SEQAN_BEGIN_TESTSUITE(mason_tests)
{
    SEQAN_CALL_TEST(mason_tests_append_orientation_elementary_operations);
    SEQAN_CALL_TEST(mason_tests_append_orientation_combination);
    SEQAN_CALL_TEST(mason_tests_append_orientation_canceling_out);

    SEQAN_CALL_TEST(mason_tests_position_map_inversion);
    SEQAN_CALL_TEST(mason_tests_position_map_translocation);
}
SEQAN_END_TESTSUITE
