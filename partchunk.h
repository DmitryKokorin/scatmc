#ifndef _PARTCHUNK_H_
#define _PARTCHUNK_H_

#include <vector>
#include <string>
#include <map>

#include "common.h"
#include "node.h"
#include "rect.h"

typedef std::vector<Rect> RectsVector;
typedef std::vector<Knot> KnotsVector;

class PartitionChunk
{
public:

    PartitionChunk();
    virtual ~PartitionChunk();
    bool create(const Float minAngle, const Float maxAngle, const int iterations);
   
    bool load(FILE *file);
    bool save(FILE *file);

	void setData(Float** const data, const Float& cellSquare);
	void refine();

	inline size_t getRectsCount() const { return m_rects.size(); }
	inline size_t getKnotsCount() const { return m_knots.size(); }
    inline bool  isAngleInRange(const Float angle) const { return angle >= m_minAngle && angle <= m_maxAngle;}

	static const Float kEpsilon;
	static const int   kThetaDegree = 10;           
	static const int   kThetaSize   = (1 << kThetaDegree) + 1;    //1025
	static const int   kPhiDegree   = 10;                     
	static const int   kPhiSize     = (1 << kPhiDegree) + 1;      //257


	static const Float kThetaResolution;
	static const Float kPhiResolution;

	
	
	RectsVector m_rects;
	KnotsVector m_knots;


private:
	
	Float integral(const GreedRect& rect);
	Float approxIntegral(const GreedRect& rect);

	void createPartitionTree();
	void createRectsList();
	
	Float m_minAngle;
	Float m_maxAngle;
	int   m_iterations;

	Float m_iterationStep; 

	
	Node*    m_root;
	Float**  m_data;             //knots
	Float**  m_cellIntegrals;    //integrals of elementary cells

//	int m_rectCount;

	std::map<int, int> m_knotsMap;

	//used on tree creation
	void refineNode(Node* node);

	//used on list creation
	void processTreeNode(Node* node);

	Float m_cellSquare;
	Float m_fullIntegral;

	//disable copying
	PartitionChunk(const PartitionChunk&);
	PartitionChunk& operator=(const PartitionChunk&);
};


#endif /* _PARTCHUNK_H_ */
