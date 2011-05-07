#ifndef _PARTITION_H_
#define _PARTITION_H_

#include <vector>
#include <string>
#include <map>

#include "common.h"
#include "node.h"
#include "rect.h"

typedef std::vector<Rect> RectsVector;
typedef std::vector<Knot> KnotsVector;

class Partition
{
public:

    Partition();
    virtual ~Partition();
    bool create();
    
    bool load(const std::string& name);
    bool save(const std::string& name);

	void setData(Float** const data, const Float& cellSquare);
	void refine();

	static const Float kEpsilon;
	static const int   kThetaDegree = 10;           
	static const int   kThetaSize   = (1 << kThetaDegree) + 1;    //1025
	static const int   kPhiDegree   = 8;                     
	static const int   kPhiSize     = (1 << kPhiDegree) + 1;      //257


	static const Float kThetaResolution;
	static const Float kPhiResolution;

	static const int   kIterations = 2000;
	static const Float kIterationStep; 
	
	int m_rectCount;

	RectsVector m_rects;
	KnotsVector m_knots;


private:
	
	Float integral(const GreedRect& rect);
	Float approxIntegral(const GreedRect& rect);

	void createPartitionTree();
	void createRectsList();
	
	
	Node*    m_root;
	Float**  m_data;             //knots
	Float**  m_cellIntegrals;    //integrals of elementary cells

	std::map<int, int> m_knotsMap;

	//used on tree creation
	void refineNode(Node* node);

	//used on list creation
	void processTreeNode(Node* node);

	Float m_cellSquare;
	Float m_fullIntegral;

	//disable copying
	Partition(const Partition&);
	Partition& operator=(const Partition&);
};




#endif /* _PARTITION_H_ */
