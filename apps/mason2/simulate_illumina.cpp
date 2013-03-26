#include "sequencing.h"

// ===========================================================================
// Class IlluminaSequencingOptions
// ===========================================================================

class IlluminaModel
{
public:
    // Probabilities for a mismatch at a given position.
    seqan::String<double> mismatchProbabilities;

    // Standard deviations for the normal distributions of base
    // qualities for the mismatch case.
    seqan::String<double> mismatchQualityMeans;
    // Standard deviations for the normal distributions of base
    // qualities for the mismatch case.
    seqan::String<double> mismatchQualityStdDevs;

    // Standard deviations for the normal distributions of base
    // qualities for the non-mismatch case.
    seqan::String<double> qualityMeans;
    // Standard deviations for the normal distributions of base
    // qualities for the non-mismatch case.
    seqan::String<double> qualityStdDevs;

    IlluminaModel()
    {}
};

// ===========================================================================
// Class IlluminaSequencingOptions
// ===========================================================================

// ---------------------------------------------------------------------------
// Function IlluminaSequencingOptions::print()
// ---------------------------------------------------------------------------

void IlluminaSequencingOptions::print(std::ostream & out)
{
    SequencingOptions::print(out);

    out << "\n"
        << "  SIMULATED TECHNOLOGY          \tillumina\n"
        << "\n"
        << "ILLUMINA OPTIONS\n"
        << "\n"
        << "  PROB INSERTION                \t" << probabilityInsert << "\n"
        << "  PROB DELETION                 \t" << probabilityDelete << "\n"
        << "  PROB MISMATCH SCALE           \t" << probabilityMismatchScale << "\n"
        << "  PROB MISMATCH                 \t" << probabilityMismatch << "\n"
        << "  PROB MISMATCH BEGIN           \t" << probabilityMismatchBegin << "\n"
        << "  PROB MISMATCH END             \t" << probabilityMismatchEnd << "\n"
        << "  POSITION RAISE                \t" << positionRaise << "\n"
        << "  QUALITY MEAN BEGIN            \t" << meanQualityBegin << "\n"
        << "  QUALITY MEAN END              \t" << meanQualityEnd << "\n"
        << "  QUALITY STD DEV BEGIN         \t" << stdDevQualityBegin << "\n"
        << "  QUALITY STD DEV END           \t" << stdDevQualityEnd << "\n"
        << "  MISMATCH QUALITY MEAN BEGIN   \t" << meanMismatchQualityBegin << "\n"
        << "  MISMATCH QUALITY MEAN END     \t" << meanMismatchQualityEnd << "\n"
        << "  MISMATCH QUALITY STD DEV BEGIN\t" << stdDevMismatchQualityBegin << "\n"
        << "  MISMATCH QUALITY STD DEV END  \t" << stdDevMismatchQualityEnd << "\n";
}

// ===========================================================================
// Class IlluminaSequencingSimulator
// ===========================================================================

// ---------------------------------------------------------------------------
// Constructor IlluminaSequencingSimulator::IlluminaSequencingSimulator
// ---------------------------------------------------------------------------

IlluminaSequencingSimulator::IlluminaSequencingSimulator(
        TRng & rng, IlluminaSequencingOptions const & illuminaOptions) :
        SequencingSimulator(rng, illuminaOptions), illuminaOptions(illuminaOptions),
        model(new IlluminaModel())
{
    this->_initModel();
}

// ---------------------------------------------------------------------------
// Function IlluminaSequencingSimulator::_initModel()
// ---------------------------------------------------------------------------

void IlluminaSequencingSimulator::_initModel()
{
    // Compute mismatch probabilities, piecewise linear function.
    resize(model->mismatchProbabilities, illuminaOptions.readLength);
    // Compute probability at raise point.
    double y_r = 2 * illuminaOptions.probabilityMismatch - illuminaOptions.positionRaise * illuminaOptions.probabilityMismatchBegin - illuminaOptions.probabilityMismatchEnd + illuminaOptions.probabilityMismatchEnd * illuminaOptions.positionRaise;
    if (illuminaOptions.verbosity >= 2)
    {
        std::cerr << "Illumina error curve:\n"
                  << "  (0, " << illuminaOptions.probabilityMismatchBegin << ") -- (" << illuminaOptions.positionRaise << ", " << y_r << ") -- (1, " << illuminaOptions.probabilityMismatchEnd << ")\n";
    }
    // std::cout << "y_r = " << y_r << std::endl;
    // Compute mismatch probability at each base.
    if (!empty(illuminaOptions.probabilityMismatchFile))
    {
        // Open file.
        std::fstream file;
        file.open(toCString(illuminaOptions.probabilityMismatchFile), std::ios_base::in);
        if (!file.is_open())
        {
            std::cerr << "Failed to load mismatch probabilities from " << illuminaOptions.probabilityMismatchFile << std::endl;
            // return 1;
        }
        // Load probabilities.
        double x;
        file >> x;
        unsigned i;
        for (i = 0; i < illuminaOptions.readLength && !file.eof(); ++i) {
            model->mismatchProbabilities[i] = x;
            file >> x;
        }
        if (i != illuminaOptions.readLength)
        {
            std::cerr << "Not enough mismatch probabilites in " << illuminaOptions.probabilityMismatchFile << " (" << i << " < " << illuminaOptions.readLength << ")!" << std::endl;
            // return 1;
        }
    } else {
        // Use piecewise linear function for mismatch probability simulation.
        for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
            double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
            if (x < illuminaOptions.positionRaise) {
                double b = illuminaOptions.probabilityMismatchBegin;
                double m = (y_r - illuminaOptions.probabilityMismatchBegin) / illuminaOptions.positionRaise;
                model->mismatchProbabilities[i] = m * x + b;
                // std::cout << "model->mismatchProbabilities[" << i << "] = " << model->mismatchProbabilities[i] << std::endl;
            } else {
                double b = y_r;
                double m = (illuminaOptions.probabilityMismatchEnd - y_r) / (1 - illuminaOptions.positionRaise);
                x -= illuminaOptions.positionRaise;
                model->mismatchProbabilities[i] = m * x + b;
                // std::cout << "model->mismatchProbabilities[" << i << "] = " << model->mismatchProbabilities[i] << std::endl;
            }
        }
    }
    if (illuminaOptions.probabilityMismatchScale != 1.0) {
        for (unsigned i = 0; i < illuminaOptions.readLength; ++i)
            model->mismatchProbabilities[i] *= illuminaOptions.probabilityMismatchScale;
    }

    // Compute match/mismatch means and standard deviations.
    resize(model->mismatchQualityMeans, illuminaOptions.readLength);
    for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
        double b = illuminaOptions.meanMismatchQualityBegin;
        double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
        double m = (illuminaOptions.meanMismatchQualityEnd - illuminaOptions.meanMismatchQualityBegin);
        model->mismatchQualityMeans[i] = m * x + b;
        // std::cout << "model->mismatchQualityMeans[" << i << "] = " << model->mismatchQualityMeans[i] << std::endl;
    }
    resize(model->mismatchQualityStdDevs, illuminaOptions.readLength);
    for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
        double b = illuminaOptions.stdDevMismatchQualityBegin;
        double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
        double m = (illuminaOptions.stdDevMismatchQualityEnd - illuminaOptions.stdDevMismatchQualityBegin);
        model->mismatchQualityStdDevs[i] = m * x + b;
        // std::cout << "model->mismatchQualityStdDevs[" << i << "] = " << model->mismatchQualityStdDevs[i] << std::endl;
    }
    resize(model->qualityMeans, illuminaOptions.readLength);
    for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
        double b = illuminaOptions.meanQualityBegin;
        double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
        double m = (illuminaOptions.meanQualityEnd - illuminaOptions.meanQualityBegin);
        model->qualityMeans[i] = m * x + b;
        // std::cout << "model->qualityMeans[" << i << "] = " << model->qualityMeans[i] << std::endl;
    }
    resize(model->qualityStdDevs, illuminaOptions.readLength);
    for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
        double b = illuminaOptions.stdDevQualityBegin;
        double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
        double m = (illuminaOptions.stdDevQualityEnd - illuminaOptions.stdDevQualityBegin);
        model->qualityStdDevs[i] = m * x + b;
        // std::cout << "model->qualityStdDevs[" << i << "] = " << model->qualityStdDevs[i] << std::endl;
    }
}

// ---------------------------------------------------------------------------
// Function _simulateRead()
// ---------------------------------------------------------------------------

namespace {

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
        unsigned numSimulate = 0;
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
// Function IlluminaSequencingSimulator::simulateRead()
// ---------------------------------------------------------------------------

// Actually simulate read and qualities from fragment and direction forward/reverse strand.
void IlluminaSequencingSimulator::simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                                               TFragment const & frag, Direction dir, Strand strand)
{
    // std::cerr << "simulateRead(" << (char const *)(dir == LEFT ? "L" : "R") << ", " << (char const *)(strand == FORWARD ? "-->" : "<--") << ")\n";
    // Simulate sequencing operations.
    TCigarString cigar;
    _simulateCigar(cigar);
    unsigned lenInRef = 0;
    _getLengthInRef(cigar, lenInRef);

    // TODO(holtgrew): Check that the CIGAR string does not need more sequence than we have in frag.

    // Simulate sequence (materialize mismatches and insertions).
    typedef seqan::ModifiedString<seqan::ModifiedString<TFragment, seqan::ModView<seqan::FunctorComplement<seqan::Dna5> > >, seqan::ModReverse> TRevCompFrag;
    if ((dir == LEFT) && (strand == FORWARD))
        _simulateSequence(seq, rng, prefix(frag, lenInRef), cigar);
    else if ((dir == LEFT) && (strand == REVERSE))
        _simulateSequence(seq, rng, TRevCompFrag(prefix(frag, lenInRef)), cigar);
    else if ((dir == RIGHT) && (strand == FORWARD))
        _simulateSequence(seq, rng, suffix(frag, length(frag) - lenInRef), cigar);
    else  // ((dir == RIGHT) && (strand == REVERSE))
        _simulateSequence(seq, rng, TRevCompFrag(suffix(frag, length(frag) - lenInRef)), cigar);

    // Simulate qualities.
    _simulateQualities(quals, cigar);
    SEQAN_ASSERT_EQ(length(seq), length(quals));

    // Reverse qualities if necessary.
    if (strand == REVERSE)
        reverse(quals);

    // Write out sequencing information info if configured to do so.
    if (illuminaOptions.embedReadInfo)
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
    // std::cerr << "  frag=" << frag << "\n";
    // std::cerr << "  " << seq << "\t" << quals << "\n";
}

// ---------------------------------------------------------------------------
// Function IlluminaSequencingSimulator::_simulateQualities()
// ---------------------------------------------------------------------------

// Simulate PHRED qualities from the CIGAR string.
void IlluminaSequencingSimulator::_simulateQualities(TQualities & quals, TCigarString const & cigar)
{
    clear(quals);

    unsigned pos = 0;
    for (unsigned i = 0; i < length(cigar); ++i)
    {
        for (unsigned j = 0; j < cigar[i].count; ++j)
        {
            int q = 0;
            if (cigar[i].operation == 'M')
            {
                seqan::Pdf<seqan::Normal> pdf(model->qualityMeans[pos], model->qualityStdDevs[pos]);
                q = static_cast<int>(pickRandomNumber(rng, pdf));
                ++pos;
            }
            else if (cigar[i].operation == 'I' || cigar[i].operation == 'X')
            {
                seqan::Pdf<seqan::Normal> pdf(model->mismatchQualityMeans[pos], model->mismatchQualityStdDevs[pos]);
                q = static_cast<int>(pickRandomNumber(rng, pdf));
                ++pos;
            }
            else
            {
                // Deletion/padding, no quality required.
                continue;
            }
            q = std::max(0, std::min(q, 40));  // limit quality to 0..40
            appendValue(quals, (char)('!' + q));
        }
    }
}

// ---------------------------------------------------------------------------
// Function IlluminaSequencingSimulator::_simulateCigar()
// ---------------------------------------------------------------------------

// Simulate CIGAR string.  We can do this with position specific parameters only and thus independent of any
// context.
void IlluminaSequencingSimulator::_simulateCigar(TCigarString & cigar)
{
    clear(cigar);
    unsigned len = this->readLength();

    for (int i = 0; i < (int)len;)
    {
        double x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double> >(0, 1));
        double pMismatch = model->mismatchProbabilities[i];
        double pInsert   = illuminaOptions.probabilityInsert;
        double pDelete   = illuminaOptions.probabilityDelete;
        double pMatch    = 1.0 - pMismatch - pInsert - pDelete;

        // Simulate mutation/insertion/deletion events.  If possible we reuse the last CIGAR entry.  Adjacent
        // insertion/deletion pairs cancel each other out.

        // TODO(holtgrew): No indel at beginning or ending! Same for other simulators!

        if (x < pMatch)  // match
            i += appendOperation(cigar, 'M').first;
        else if (x < pMatch + pMismatch)  // point polymorphism
            i += appendOperation(cigar, 'X').first;
        else if (x < pMatch + pMismatch + pInsert) // insertion
            i += appendOperation(cigar, 'I').first;
        else  // deletion
            i += appendOperation(cigar, 'D').first;
    }
}

// ============================================================================
// Class SequencingSimulatorFactory
// ============================================================================

// ----------------------------------------------------------------------------
// Function SequencingSimulatorFactory::make()
// ----------------------------------------------------------------------------

std::auto_ptr<SequencingSimulator> SequencingSimulatorFactory::make()
{
    std::auto_ptr<SequencingSimulator> res;
        
    switch (tech)
    {
        case ILLUMINA:
            res.reset(new IlluminaSequencingSimulator(rng, *dynamic_cast<IlluminaSequencingOptions const *>(options)));
            break;
        case SANGER:
            res.reset(new SangerSequencingSimulator(rng, *dynamic_cast<SangerSequencingOptions const *>(options)));
            break;
        case ROCHE_454:
            res.reset(new Roche454SequencingSimulator(rng, *dynamic_cast<Roche454SequencingOptions const *>(options)));
            break;
    }

    return res;
}
