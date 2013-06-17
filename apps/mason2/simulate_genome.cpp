#include "simulate_genome.h"

// ----------------------------------------------------------------------------
// Function simulateGenome()
// ----------------------------------------------------------------------------

// Simulate a genome given the simulation options.
//
// The resulting sequence is written to stream.

int simulateGenome(seqan::SequenceStream & stream, MasonSimulateGenomeOptions const & options)
{
    // Initialize RNG and PDF.
    seqan::Rng<seqan::MersenneTwister>  rng(options.seed);
    seqan::Pdf<seqan::Uniform<double> > pdf(0, 1);

    seqan::CharString id;
    seqan::Dna5String contig;
    
    for (unsigned i = 0; i < length(options.contigLengths); ++i)
    {
        clear(id);
        clear(contig);

        std::stringstream ss;
        ss << (i + 1);
        id = ss.str();
        
        std::cerr << "contig " << id << " ...";

        for (int j = 0; j < options.contigLengths[i];)
        {
            double x = pickRandomNumber(rng, pdf);
            if (x < 0.25)
                appendValue(contig, 'A');
            else if (x < 0.5)
                appendValue(contig, 'C');
            else if (x < 0.75)
                appendValue(contig, 'G');
            else if (x < 1.0)
                appendValue(contig, 'T');
            else
                continue;  // Redraw.
            ++j;
        }

        if (writeRecord(stream, id, contig) != 0)
        {
            std::cerr << "\nERROR: Could not write contig " << id << " to output file.\n";
            return 1;
        }

        std::cerr << " DONE\n";
    }

    return 0;
} 

// ----------------------------------------------------------------------------
// Function simulateGenome()
// ----------------------------------------------------------------------------

int simulateGenome(char const * filename, MasonSimulateGenomeOptions const & options)
{
    seqan::SequenceStream stream;
    open(stream, filename, seqan::SequenceStream::WRITE, seqan::SequenceStream::FASTA);
    if (!isGood(stream))
    {
        std::cerr << "ERROR: Could not open " << filename << "for writing!\n";
        return 1;
    }

    return simulateGenome(stream, options);
} 
