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
                                            TFragment const & frag)
{
    bool isForward = (pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, 1)) == 1);
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
// Function SequencingSimulator::simulateSingleEnd()
// ---------------------------------------------------------------------------

// Simulate single-end sequencing from a fragment.
void SequencingSimulator::simulateSingleEnd(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                                            TFragment const & frag)
{
    bool isForward = (pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, 1)) == 1);
    Strand strand = isForward ? FORWARD : REVERSE;
    this->simulateRead(seq, quals, info, frag, LEFT, strand);
}
