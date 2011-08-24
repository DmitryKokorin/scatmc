#ifndef _PARTITION_H_

#define _PARTITION_H_

#include <list>
#include <string>
#include <cstdio>

#include "common.h"
#include "partchunk.h"


typedef std::list<PartitionChunk*> ChunkList;


class Partition
{
public:

    Partition();
    ~Partition();

    bool load(const std::string& name);
    bool save(const std::string& name);

    template <class T>
    PartitionChunk* addChunk(const Float kMinAngle, const Float kMaxAngle, const int kIterations);

    void addChunk(PartitionChunk *chunk);
    PartitionChunk* getChunk(const Float kAngle);

    inline size_t getMaxRectsCount() const { return m_maxRectsCount; }
    inline size_t getMaxKnotsCount() const { return m_maxKnotsCount; }

private:

    ChunkList m_chunks;

    size_t    m_maxRectsCount;
    size_t    m_maxKnotsCount;

    Partition(const Partition&);
    Partition& operator=(const Partition&);
};

template <class T>
PartitionChunk* Partition::addChunk(const Float kMinAngle, const Float kMaxAngle, const int kIterations)
{
    PartitionChunk *chunk = new PartitionChunk();

    fprintf(stderr, "creating partition chunk:\nangle: %f - %f\titerations: %d...\n", kMinAngle, kMaxAngle, kIterations);
    chunk->create<T>(kMinAngle, kMaxAngle, kIterations);
    fprintf(stderr, "done, %lu rects\n", (unsigned long int)chunk->getRectsCount());
    
    addChunk(chunk);

    return chunk;
}



#endif /* end of include guard: _PARTITION_H_ */
