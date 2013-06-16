#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>

#include "sequencing.h"

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

SEQAN_BEGIN_TESTSUITE(mason_tests)
{
    SEQAN_CALL_TEST(mason_tests_append_orientation_elementary_operations);
    SEQAN_CALL_TEST(mason_tests_append_orientation_combination);
    SEQAN_CALL_TEST(mason_tests_append_orientation_canceling_out);
}
SEQAN_END_TESTSUITE
