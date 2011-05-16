#ifndef PARTITION_P92OD9XA

#define PARTITION_P92OD9XA

#include <list>
#include <string>

#include "common.h"

class PartitionChunk;
typedef std::list<PartitionChunk*> ChunkList;


class Partition
{
public:

    Partition();
    ~Partition();

    bool load(const std::string& name);
    bool save(const std::string& name);

    PartitionChunk* addChunk(const Float kMinAngle, const Float kMaxAngle, const int kIterations);
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


#endif /* end of include guard: PARTITION_P92OD9XA */
