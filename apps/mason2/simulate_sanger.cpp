#include "sequencing.h"

// ===========================================================================
// Class SangerSequencingSimulator
// ===========================================================================

// ---------------------------------------------------------------------------
// Helper Function _simulateSequence()
// ---------------------------------------------------------------------------

namespace {

// TODO(holtgrew): Copy-and-paste from IlluminaSequencingSimulator.
//
// Simulate the characters that polymorphisms turn into and inserted characters.
//
// Through the usage of ModifiedString, we will always go from the left to the right end.
template <typename TFrag>
void _simulateSequence(TRead & read, TRng & rng, TFrag const & frag,
                       TCigarString const & cigar)
{
    clear(read);

    typedef typename seqan::Iterator<TFrag>::Type TFragIter;
    TFragIter it = begin(frag, seqan::Standard());

    for (unsigned i = 0; i < length(cigar); ++i)
    {
        //unsigned numSimulate = 0;
        if (cigar[i].operation == 'M')
        {
            for (unsigned j = 0; j < cigar[i].count; ++j, ++it)
                appendValue(read, *it);
            continue;
        }
        else if (cigar[i].operation == 'D')
        {
            it += cigar[i].count;
            continue;
        }

        // Otherwise, we have insertions or mismatches.
        for (unsigned j = 0; j < cigar[i].count; ++j)
        {
            // Pick a value between 0 and 1.
            double x = 1.0;
            while (x == 1.0)
                x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double> >(0, 1));
            int num = x / 0.25;

            // NOTE: We can only insert CGAT, but we can have a polymorphism to N.

            if (cigar[i].operation == 'I')
                appendValue(read, seqan::Dna5(num));
            else
                appendValue(read, seqan::Dna5(num + (num == ordValue(*it))));
        }
    }
}

}  // namespace (anonymous)

// ---------------------------------------------------------------------------
// Function SangerSequencingSimulator::readLength()
// ---------------------------------------------------------------------------

unsigned SangerSequencingSimulator::readLength()
{
    if (sangerOptions.readLengthIsUniform)
    {
        // Pick uniformly.
        double minLen = sangerOptions.readLengthMin;
        double maxLen = sangerOptions.readLengthMax;
        double len = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double> >(minLen, maxLen));
        return static_cast<unsigned>(round(len));
    }
    else
    {
        // Pick normally distributed.
        double len = pickRandomNumber(rng, seqan::Pdf<seqan::Normal>(sangerOptions.readLengthMean, sangerOptions.readLengthError));
        return static_cast<unsigned>(round(len));
    }
}

// ---------------------------------------------------------------------------
// Function SangerSequencingSimulator::simulateRead()
// ---------------------------------------------------------------------------

// Actually simulate read and qualities from fragment and direction forward/reverse strand.
void SangerSequencingSimulator::simulateRead(
        TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
        TFragment const & frag, Direction dir, Strand strand)
{
    // TODO(holtgrew): Pick better names for read length here.
    // Pick sampled length.
    unsigned readLength = this->readLength();

    if (readLength > length(frag))
    {
        throw std::runtime_error("Sanger read is too long, increase fragment length");
    }

    // Simulate CIGAR string.
    TCigarString cigar;
    this->_simulateCigar(cigar, readLength);

    // Simulate sequence (materialize mismatches and insertions).
    typedef seqan::ModifiedString<seqan::ModifiedString<TFragment, seqan::ModView<seqan::FunctorComplement<seqan::Dna5> > >, seqan::ModReverse> TRevCompFrag;
    if ((dir == LEFT) && (strand == FORWARD))
        _simulateSequence(seq, rng, prefix(frag, readLength), cigar);
    else if ((dir == LEFT) && (strand == REVERSE))
        _simulateSequence(seq, rng, TRevCompFrag(prefix(frag, readLength)), cigar);
    else if ((dir == RIGHT) && (strand == FORWARD))
        _simulateSequence(seq, rng, suffix(frag, length(frag) - readLength), cigar);
    else  // ((dir == RIGHT) && (strand == REVERSE))
        _simulateSequence(seq, rng, TRevCompFrag(suffix(frag, length(frag) - readLength)), cigar);

    // Simulate Qualities.
    this->_simulateQualities(quals, cigar, readLength);

    // Reverse qualities if necessary.
    if (strand == REVERSE)
        reverse(quals);

    // Write out sequencing information info if configured to do so.
    if (seqOptions->embedReadInfo)
    {
        info.cigar = cigar;
        unsigned len = 0;
        _getLengthInRef(cigar, len);
        if (dir == LEFT)
            info.sampleSequence = prefix(frag, len);
        else
            info.sampleSequence = suffix(frag, length(frag) - len);
        info.isForward = (strand == FORWARD);
        if (strand == REVERSE)
            reverseComplement(info.sampleSequence);
    }
}

// ---------------------------------------------------------------------------
// Function SangerSequencingSimulator::_simulateCigar()
// ---------------------------------------------------------------------------

// Simulate CIGAR string.  We can do this with position specific parameters only and thus independent of any
// context.
void SangerSequencingSimulator::_simulateCigar(TCigarString & cigar, unsigned readLength)
{
    clear(cigar);

    for (unsigned i = 0; i < readLength;)
    {
        double x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double> >(0, 1));
        double pos = 1.0 * i / (readLength - 1);
        double pMismatch = sangerOptions.probabilityMismatchBegin + pos * (sangerOptions.probabilityMismatchEnd - sangerOptions.probabilityMismatchBegin);
        double pInsert   = sangerOptions.probabilityInsertBegin + pos * (sangerOptions.probabilityInsertEnd - sangerOptions.probabilityInsertBegin);
        double pDelete   = sangerOptions.probabilityDeleteBegin + pos * (sangerOptions.probabilityDeleteEnd - sangerOptions.probabilityDeleteBegin);
        double pMatch    = 1.0 - pMismatch - pInsert - pDelete;

        // Simulate mutation/insertion/deletion events.  If possible we reuse the last CIGAR entry.  Adjacent
        // insertion/deletion pairs cancel each other out.  We count i up to the input read length, thus using the
        // "second" member of appendOperation()'s result.

        // TODO(holtgrew): No indels at beginning or end of read.

        if (x < pMatch)  // match
            i += appendOperation(cigar, 'M').second;
        else if (x < pMatch + pMismatch)  // point polymorphism
            i += appendOperation(cigar, 'X').second;
        else if (x < pMatch + pMismatch + pInsert) // insertion
            i += appendOperation(cigar, 'I').second;
        else  // deletion
            i += appendOperation(cigar, 'D').second;
    }
}

// ---------------------------------------------------------------------------
// Function SangerSequencingSimulator::_simulateQualities()
// ---------------------------------------------------------------------------

void SangerSequencingSimulator::_simulateQualities(
        TQualities & quals, TCigarString const & cigar, unsigned sampleLength)
{
    clear(quals);

    unsigned pos = 0;   // Position in result.
    unsigned rPos = 0;  // Position in fragment.
    for (unsigned i = 0; i < length(cigar); ++i)
    {
        for (unsigned j = 0; j < cigar[i].count; ++j, ++pos)
        {
            double mean = 0.0, stdDev = 0.0;
            double relPos = 1.0 * rPos / (1.0 * sampleLength);

            rPos += (cigar[i].operation != 'D');

            if (cigar[i].operation == 'D')
            {
                continue;  // No quality to give out.
            }
            else if (cigar[i].operation == 'I' || cigar[i].operation == 'X')
            {
                mean = sangerOptions.qualityMatchStartMean + relPos * (sangerOptions.qualityMatchEndMean - sangerOptions.qualityMatchStartMean);
                stdDev = sangerOptions.qualityMatchStartStdDev + relPos * (sangerOptions.qualityMatchEndStdDev - sangerOptions.qualityMatchStartStdDev);
            }
            else  // cigar[i].operation == 'M'
            {
                mean = sangerOptions.qualityErrorStartMean + relPos * (sangerOptions.qualityErrorEndMean - sangerOptions.qualityErrorStartMean);
                stdDev = sangerOptions.qualityErrorStartStdDev + relPos * (sangerOptions.qualityErrorEndStdDev - sangerOptions.qualityErrorStartStdDev);
            }

            seqan::Pdf<seqan::Normal> pdf(mean, stdDev);
            int q = static_cast<int>(pickRandomNumber(rng, pdf));
            q = std::max(0, std::min(40, q));
            appendValue(quals, (char)('!' + q));
        }
    }
}
