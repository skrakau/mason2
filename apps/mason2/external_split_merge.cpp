#include "external_split_merge.h"

// ---------------------------------------------------------------------------
// Function IdSplitter::open()
// ---------------------------------------------------------------------------

void IdSplitter::open()
{
    close();
    for (unsigned i = 0; i < numContigs; ++i)
    {
        char const * tmpdir = 0;
        if ((tmpdir = getenv("TMPDIR")) == NULL)
            tmpdir = "/tmp";
        std::string pathTpl = tmpdir;
        pathTpl += "/MASON_XXXXXX";

        int fd = mkstemp(&pathTpl[0]);
        files.push_back(fopen(pathTpl.c_str(), "w+b"));
        remove(pathTpl.c_str());
        ::close(fd);

        if (!files.back())
        {
            std::cerr << "ERROR: Could not open temporary file!\n";
            exit(1);
        }
    }
}

// ---------------------------------------------------------------------------
// Function IdSplitter::reset()
// ---------------------------------------------------------------------------

void IdSplitter::reset()
{
    for (unsigned i = 0; i < files.size(); ++i)
        if (files[i] != 0)
        {
            SEQAN_ASSERT(!ferror(files[i]));
            fflush(files[i]);
            int res = fseek(files[i], 0L, SEEK_SET);
            (void)res;
            SEQAN_ASSERT_EQ(res, 0);
            SEQAN_ASSERT(!ferror(files[i]));
        }
}

// ---------------------------------------------------------------------------
// Function IdSplitter::close()
// ---------------------------------------------------------------------------

void IdSplitter::close()
{
    for (unsigned i = 0; i < files.size(); ++i)
        if (files[i])
        {
            fclose(files[i]);
            files[i] = 0;
        }
    files.clear();
}

// ---------------------------------------------------------------------------
// Function SamJoiner::init()
// ---------------------------------------------------------------------------

void SamJoiner::init()
{
    resize(records, splitter->files.size());
    active.resize(splitter->files.size());

    for (unsigned i = 0; i < splitter->files.size(); ++i)
    {
        readers.push_back(new TReader(splitter->files[i]));

        // We use a separate header structure and name stores and caches.  Since the headers of all files are equal, we
        // will write out the first one only.
        seqan::BamHeader tmpHeader;
        seqan::StringSet<seqan::CharString> tmpNameStore;
        seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > tmpNameStoreCache(tmpNameStore);
        seqan::BamIOContext<seqan::StringSet<seqan::CharString> > tmpContext(tmpNameStore, tmpNameStoreCache);
        if (readRecord(tmpHeader, tmpContext, *readers[i], seqan::Sam()) != 0)
            throw MasonIOException("Could not load SAM header from temporary file.");
        if (i == 0u)
        {
            header = tmpHeader;
            nameStore = tmpNameStore;
        }

        active[i] = _loadNext(records[i], i);
        numActive += (active[i] != false);
    }

    // Refresh name store cache.
    refresh(nameStoreCache);
}

// ---------------------------------------------------------------------------
// Function SamJoiner::_loadNext()
// ---------------------------------------------------------------------------

bool SamJoiner::_loadNext(seqan::BamAlignmentRecord & record, unsigned idx)
{
    if (seqan::atEnd(*readers[idx]))
        return false;
    if (readRecord(record, context, *readers[idx], seqan::Sam()) != 0)
    {
        std::cerr << "ERROR: Problem reading temporary data.\n";
        exit(1);
    }
    return true;
}

// ---------------------------------------------------------------------------
// Function SamJoiner::get()
// ---------------------------------------------------------------------------

int SamJoiner::get(seqan::BamAlignmentRecord & record)
{
    unsigned idx = seqan::maxValue<unsigned>();
    for (unsigned i = 0; i < length(records); ++i)
    {
        if (!active[i])
            continue;
        if (idx == seqan::maxValue<unsigned>() || ltBamAlignmentRecord(records[i], records[idx]))
            idx = i;
    }
    if (idx == seqan::maxValue<unsigned>())
        return 1;

    // We use double-buffering and the input parameters as buffers.
    using std::swap;
    active[idx] = _loadNext(record, idx);
    swap(record, records[idx]);
    numActive -= !active[idx];

    return 0;
}

// ---------------------------------------------------------------------------
// Function ContigPicker::pick()
// ---------------------------------------------------------------------------

std::pair<int, int> ContigPicker::pick()
{
    // Pick reference id.
    int rID = 0;
    if (lengthSums.size() > 1u)
    {
        seqan::Pdf<seqan::Uniform<__int64> > pdf(0, lengthSums.back() - 1);
        __int64 x = pickRandomNumber(rng, pdf);
        for (unsigned i = 0; i < lengthSums.size(); ++i)
        {
            if (x >= lengthSums[i])
                rID = i + 1;
            if (x < lengthSums[i])
                break;
        }
    }

    // Pick haplotype id.
    int hID = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, numHaplotypes - 1));

    return std::make_pair(rID, hID);
}
