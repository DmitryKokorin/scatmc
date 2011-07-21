#include <algorithm>
#include <cstdio>

#include "partition.h"
#include "partchunk.h"

Partition::Partition() :
    m_chunks(),
    m_maxRectsCount(0),
    m_maxKnotsCount(0)
{
}


Partition::~Partition()
{
    ChunkList::iterator i;

    for (i = m_chunks.begin(); i != m_chunks.end(); ++i) {

        delete (*i);
    }
}

bool Partition::load(const std::string& name)
{
    FILE *file = fopen(name.c_str(), "r");
    if (NULL == file)
        return false;

    int res = 0;

    unsigned long int chunksNum, i;
    res = fscanf(file, "%lu", &chunksNum);

    if (0 == res)
        return false;

    for (i = 0; i < chunksNum; ++i) {

        PartitionChunk *chunk = new PartitionChunk();

        chunk->load(file);
        addChunk(chunk);
    }

    fclose(file);

    return true;
}

bool Partition::save(const std::string& name)
{
    FILE *file = fopen(name.c_str(), "w");
    if (NULL == file)
        return false;

    fprintf(file, "%lu\n", (unsigned long int)m_chunks.size());

    ChunkList::iterator i;

    for (i = m_chunks.begin(); i != m_chunks.end(); ++i) {

        (*i)->save(file);
    }

    fclose(file);

    return true;
}

PartitionChunk* Partition::addChunk(const Float kMinAngle, const Float kMaxAngle, const int kIterations)
{
    PartitionChunk *chunk = new PartitionChunk();

    fprintf(stderr, "creating partition chunk:\nangle: %f - %f\titerations: %d...\n", kMinAngle, kMaxAngle, kIterations);
    chunk->create(kMinAngle, kMaxAngle, kIterations);
    fprintf(stderr, "done, %lu rects\n", (unsigned long int)chunk->getRectsCount());
    
    addChunk(chunk);

    return chunk;
}

void Partition::addChunk(PartitionChunk *chunk)
{
    m_chunks.push_back(chunk);
    m_maxRectsCount = std::max(m_maxRectsCount, chunk->getRectsCount());
    m_maxKnotsCount = std::max(m_maxKnotsCount, chunk->getKnotsCount());
}

PartitionChunk* Partition::getChunk(const Float angle)
{
    ChunkList::iterator i = m_chunks.begin();

    for (;i != m_chunks.end(); ++i)
        if ((*i)->isAngleInRange(angle))
            return *i;

    return NULL;
}
