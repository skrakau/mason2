#include "sequencing.h"

// ===========================================================================
// Class SequencingSimulationInfo
// ===========================================================================

// ---------------------------------------------------------------------------
// Function SequencingSimulationInfo::()
// ---------------------------------------------------------------------------

// ===========================================================================
// Class SequencingOptions
// ===========================================================================

// ---------------------------------------------------------------------------
// Function SequencingOptions::yesNo()
// ---------------------------------------------------------------------------

const char * SequencingOptions::yesNo(bool b)
{
    return b ? "YES" : "NO";
}

// ---------------------------------------------------------------------------
// Function SequencingOptions::mateOrientationStr()
// ---------------------------------------------------------------------------

const char * SequencingOptions::mateOrientationStr(MateOrientation o)
{
    switch (o)
    {
        case FORWARD_REVERSE:
            return "FR (R1 --> <-- R2)";
        case REVERSE_FORWARD:
            return "RF (R1 <-- --> R2)";
        case FORWARD_FORWARD:
            return "FF (R1 --> --> R2)";
        case FORWARD_FORWARD2:
            return "FF2 (R2 --> --> R1)";
        default:
            return "INVALID";
    }
}            

// ---------------------------------------------------------------------------
// Function SequencingOptions::strandsStr()
// ---------------------------------------------------------------------------

const char * SequencingOptions::strandsStr(SourceStrands s)
{
    switch (s)
    {
        case BOTH:
            return "BOTH";
        case FORWARD:
            return "FORWARD";
        case REVERSE:
            return "REVERSE";
        default:
            return "INVALID";
    }
}

// ---------------------------------------------------------------------------
// Function SequencingOptions::print()
// ---------------------------------------------------------------------------

void SequencingOptions::print(std::ostream & out)
{
    out << "SIMULATION OPTIONS\n"
        << "\n";
    if (numReads > 0)  // is 0 if doing simulatin from fragment (induces numReads)
        out << "  NUMBER OF READS               \t" << numReads << "\n";
    out << "  SIMULATE QUALITIES            \t" << yesNo(simulateQualities) << "\n"
        << "  SIMULATE MATE PAIRS           \t" << yesNo(simulateMatePairs) << "\n"
        << "  MATE ORIENTATION              \t" << mateOrientationStr(mateOrientation) << "\n"
        << "  STRANDS                       \t" << strandsStr(strands) << "\n";
}

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
    switch (options->mateOrientation)
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
