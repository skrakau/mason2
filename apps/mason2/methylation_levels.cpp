#include "methylation_levels.h"

// ----------------------------------------------------------------------------
// Function VariantMaterializer::_fixVariationLevels()
// ----------------------------------------------------------------------------

void fixVariationLevels(MethylationLevels & levels,
                        TRng & rng,
                        seqan::Dna5String const & contig,
                        seqan::String<std::pair<int, bool> > const & varPoints,
                        MethylationLevelSimulatorOptions const & options)
{
    MethylationLevelSimulator methSim(rng, options);
    seqan::Shape<seqan::Dna5> shape2, shape3;
    resize(shape2, 2);
    resize(shape3, 3);

    for (unsigned i = 0; i < length(varPoints); ++i)
    {
        int pos = varPoints[i].first;
        if (varPoints[i].second)  // is SNP
        {
            if (pos > 2)
            {
                levels.forward[pos - 2] = levels.reverse[pos - 2] = '!';
                methSim.handleOneMer(levels, pos - 2, ordValue(contig[pos - 2]));
                methSim.handleTwoMer(levels, pos - 2, hash(shape2, iter(contig, pos - 2, seqan::Standard())));
                methSim.handleThreeMer(levels, pos - 2, hash(shape3, iter(contig, pos - 2, seqan::Standard())));
            }
            if (pos > 1)
            {
                levels.forward[pos - 1] = levels.reverse[pos - 1] = '!';
                methSim.handleOneMer(levels, pos - 1, ordValue(contig[pos - 1]));
                methSim.handleTwoMer(levels, pos - 1, hash(shape2, iter(contig, pos - 1, seqan::Standard())));
            }
            levels.forward[pos] = levels.reverse[pos] = '!';
            methSim.handleOneMer(levels, pos, ordValue(contig[pos]));
            if (pos + 1 < (int)length(contig))
            {
                levels.forward[pos + 1] = levels.reverse[pos + 1] = '!';
                methSim.handleOneMer(levels, pos + 1, ordValue(contig[pos + 1]));
                methSim.handleTwoMer(levels, pos, hash(shape2, iter(contig, pos, seqan::Standard())));
            }
            if (pos + 2 < (int)length(contig))
            {
                levels.forward[pos + 2] = levels.reverse[pos + 2] = '!';
                methSim.handleOneMer(levels, pos + 2, ordValue(contig[pos + 2]));
                methSim.handleTwoMer(levels, pos + 1, hash(shape2, iter(contig, pos + 1, seqan::Standard())));
                methSim.handleThreeMer(levels, pos, hash(shape3, iter(contig, pos, seqan::Standard())));
            }
        }
        else  // is no SNP but breakpoint
        {
            // TODO(holtgrew): Double-check for correctness, might recompute too much around breakpoints.
            if (pos > 2)
            {
                levels.forward[pos - 2] = levels.reverse[pos - 2] = '!';
                methSim.handleOneMer(levels, pos - 2, ordValue(contig[pos - 2]));
                methSim.handleTwoMer(levels, pos - 2, hash(shape2, iter(contig, pos - 2, seqan::Standard())));
                methSim.handleThreeMer(levels, pos - 2, hash(shape3, iter(contig, pos - 2, seqan::Standard())));
            }
            if (pos > 1)
            {
                levels.forward[pos - 1] = levels.reverse[pos - 1] = '!';
                methSim.handleOneMer(levels, pos - 1, ordValue(contig[pos - 1]));
                methSim.handleTwoMer(levels, pos - 1, hash(shape2, iter(contig, pos - 1, seqan::Standard())));
            }
            levels.forward[pos] = levels.reverse[pos] = '!';
            methSim.handleOneMer(levels, pos, ordValue(contig[pos]));
            if (pos + 1 < (int)length(contig))
            {
                methSim.handleTwoMer(levels, pos, hash(shape2, iter(contig, pos, seqan::Standard())));
                levels.forward[pos + 1] = levels.reverse[pos + 1] = '!';
                methSim.handleOneMer(levels, pos + 1, ordValue(contig[pos + 2]));
            }
            if (pos + 2 < (int)length(contig))
            {
                methSim.handleTwoMer(levels, pos + 1, hash(shape2, iter(contig, pos + 1, seqan::Standard())));
                methSim.handleThreeMer(levels, pos, hash(shape3, iter(contig, pos, seqan::Standard())));
            }
        }
    }
}
