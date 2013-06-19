#include "sequencing.h"

// ===========================================================================
// Class SequencingSimulationInfo
// ===========================================================================

// ---------------------------------------------------------------------------
// Function SequencingSimulationInfo::()
// ---------------------------------------------------------------------------

// ===========================================================================
// Class SequencingSimulator
// ===========================================================================

// ---------------------------------------------------------------------------
// Function SequencingSimulator::simulatePairedEnds()
// ---------------------------------------------------------------------------

// Simulate paired-end sequencing from a fragment.
void SequencingSimulator::simulatePairedEnd(TRead & seqL, TQualities & qualsL, SequencingSimulationInfo & infoL,
                                            TRead & seqR, TQualities & qualsR, SequencingSimulationInfo & infoR,
                                            TFragment const & frag,
                                            MethylationLevels const * levels)
{
    bool isForward = (pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, 1)) == 1);
    if (!seqOptions->bsSeqOptions.bsSimEnabled)
    {
        _simulatePairedEnd(seqL, qualsL, infoL, seqR, qualsR, infoR, frag, isForward);
    }
    else
    {
        SEQAN_ASSERT(levels);
        bool bsForward = isForward;
        // Re-pick strandedness of the BS-treated fragment.
        if (seqOptions->bsSeqOptions.bsProtocol != BSSeqOptions::DIRECTIONAL)
            bsForward = (pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, 1)) == 1);
        _simulateBSTreatment(methFrag, frag, *levels, !bsForward);
        _simulatePairedEnd(seqL, qualsL, infoL, seqR, qualsR, infoR,
                           infix(methFrag, 0, length(methFrag)), isForward);
    }
}

// ---------------------------------------------------------------------------
// Function SequencingSimulator::simulateSingleEnd()
// ---------------------------------------------------------------------------

// Simulate single-end sequencing from a fragment.
void SequencingSimulator::simulateSingleEnd(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                                            TFragment const & frag,
                                            MethylationLevels const * levels)
{
    bool isForward = (pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, 1)) == 1);
    if (!seqOptions->bsSeqOptions.bsSimEnabled)
    {
        _simulateSingleEnd(seq, quals, info, frag, isForward);
    }
    else
    {
        SEQAN_ASSERT(levels);
        bool bsForward = isForward;
        // Re-pick strandedness of the BS-treated fragment.
        if (seqOptions->bsSeqOptions.bsProtocol != BSSeqOptions::DIRECTIONAL)
            bsForward = (pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, 1)) == 1);
        _simulateBSTreatment(methFrag, frag, *levels, !bsForward);
        _simulateSingleEnd(seq, quals, info, infix(methFrag, 0, length(methFrag)), isForward);
    }
}

// ---------------------------------------------------------------------------
// Function SequencingSimulator::_simulateBSSeq()
// ---------------------------------------------------------------------------

// Simulate single-end sequencing from a fragment.
void SequencingSimulator::_simulateBSTreatment(seqan::Dna5String & methFragment,
                                               TFragment const & frag,
                                               MethylationLevels const & levels,
                                               bool reverse)
{
    methFragment = frag;
    for (unsigned pos = 0; pos != length(frag); ++pos)
    {
        double level = reverse ? levels.levelR(pos + beginPosition(frag)) : levels.levelF(pos + beginPosition(frag));
        if ((!reverse && methFragment[pos] != 'C') || (reverse && methFragment[pos] != 'G'))  // Skip all non-cyteline chars
        {
            SEQAN_ASSERT_EQ_MSG(level, 0.0,
                                "Methylation for non-C should be 0 (pos+beginPosition(frag)=%d, reverse=%u",
                                pos + beginPosition(frag), reverse);
            continue;
        }

        // Decide whether methFragment[pos] is methylated.  If this is the case then we leave it untouched.
        seqan::Pdf<seqan::Uniform<double> > pdf(0, 1);
        if (pickRandomNumber(rng, pdf) < level)
            continue;

        // Otherwise, pick whether we will convert.
        if (pickRandomNumber(rng, pdf) < seqOptions->bsSeqOptions.bsConversionRate)
            methFragment[pos] = reverse ? 'A' : 'T';
    }
}

// ---------------------------------------------------------------------------
// Function SequencingSimulator::simulatePairedEnds()
// ---------------------------------------------------------------------------

// Simulate paired-end sequencing from a fragment.
void SequencingSimulator::_simulatePairedEnd(TRead & seqL, TQualities & qualsL, SequencingSimulationInfo & infoL,
                                             TRead & seqR, TQualities & qualsR, SequencingSimulationInfo & infoR,
                                             TFragment const & frag,
                                             bool isForward)
{
    // TODO(holtgrew): Use a table for this to simplify things?
    Strand leftStrand = isForward ? FORWARD : REVERSE;
    switch (seqOptions->mateOrientation)
    {
        case SequencingOptions::FORWARD_REVERSE:
            if (leftStrand == FORWARD)
            {
                this->simulateRead(seqL, qualsL, infoL, frag, LEFT, FORWARD);
                this->simulateRead(seqR, qualsR, infoR, frag, RIGHT, REVERSE);
            }
            else
            {
                this->simulateRead(seqL, qualsL, infoL, frag, RIGHT, REVERSE);
                this->simulateRead(seqR, qualsR, infoR, frag, LEFT, FORWARD);
            }
            break;
        case SequencingOptions::REVERSE_FORWARD:
            if (leftStrand == FORWARD)
            {
                this->simulateRead(seqL, qualsL, infoL, frag, LEFT, REVERSE);
                this->simulateRead(seqR, qualsR, infoR, frag, RIGHT, FORWARD);
            }
            else
            {
                this->simulateRead(seqL, qualsL, infoL, frag, RIGHT, FORWARD);
                this->simulateRead(seqR, qualsR, infoR, frag, LEFT, REVERSE);
            }
            break;
        case SequencingOptions::FORWARD_FORWARD:
            if (leftStrand == FORWARD)
            {
                this->simulateRead(seqL, qualsL, infoL, frag, LEFT, FORWARD);
                this->simulateRead(seqR, qualsR, infoR, frag, RIGHT, FORWARD);
            }
            else
            {
                this->simulateRead(seqL, qualsL, infoL, frag, RIGHT, REVERSE);
                this->simulateRead(seqR, qualsR, infoR, frag, LEFT, REVERSE);
            }
            break;
        case SequencingOptions::FORWARD_FORWARD2:
            if (leftStrand == FORWARD)
            {
                this->simulateRead(seqL, qualsL, infoL, frag, RIGHT, FORWARD);
                this->simulateRead(seqR, qualsR, infoR, frag, LEFT, FORWARD);
            }
            else
            {
                this->simulateRead(seqL, qualsL, infoL, frag, LEFT, REVERSE);
                this->simulateRead(seqR, qualsR, infoR, frag, RIGHT, REVERSE);
            }
            break;
    }
    // std::cerr << "\n";
}

// ---------------------------------------------------------------------------
// Function SequencingSimulator::_simulateSingleEnd()
// ---------------------------------------------------------------------------

// Simulate single-end sequencing from a fragment.
void SequencingSimulator::_simulateSingleEnd(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                                             TFragment const & frag,
                                             bool isForward)
{
    Strand strand = isForward ? FORWARD : REVERSE;
    this->simulateRead(seq, quals, info, frag, LEFT, strand);
}
