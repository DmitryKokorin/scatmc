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

    fclose(file);

    return true;
}

bool Partition::save(const std::string& name)
{
    FILE *file = fopen(name.c_str(), "w");
    if (NULL == file)
        return false;

    fclose(file);

    return true;
}

PartitionChunk* Partition::addChunk(const Float kMinAngle, const Float kMaxAngle, const int kIterations)
{
    PartitionChunk *chunk = new PartitionChunk();

    fprintf(stderr, "creating partition chunk:\nangle: %f - %f\titerations: %d...\n", kMinAngle, kMaxAngle, kIterations);
    chunk->create(kMinAngle, kMaxAngle, kIterations);
    fprintf(stderr, "done, %lu rects\n", (unsigned long int)chunk->getRectsCount());
    
    m_chunks.push_back(chunk);
    m_maxRectsCount = std::max(m_maxRectsCount, chunk->getRectsCount());
    m_maxKnotsCount = std::max(m_maxKnotsCount, chunk->getKnotsCount());

    return chunk;
}

PartitionChunk* Partition::getChunk(const Float angle)
{
    ChunkList::iterator i = m_chunks.begin();

    for (;i != m_chunks.end(); ++i)
        if ((*i)->isAngleInRange(angle))
            return *i;

    return NULL;

}
