#ifndef _PARTITION_H_
#define _PARTITION_H_

#include <list>
#include <string>

#include "common.h"
#include "node.h"
#include "rect.h"


class Partition
{
public:

    Partition();
    virtual ~Partition();
    bool create();
    
    bool load(const std::string& name);
    bool save(const std::string& name);

	void setData(Float** const data_, const Float& cellSquare_);
	void refine();

	static const Float epsilon;
	static const int   degree = 10;           
	static const int   size   = (1 << degree) + 1;    //1025x1025
	
	static int rectCount;

	std::list<Rect> rects;


private:
	
	Float integral(const GreedRect& rect);
	Float approxIntegral(const GreedRect& rect);

	void preparePartitionTree();
	
	void processTreeNode(Node* pNode);

	Node* pRoot;
	Float**  data;             //knots
	Float**  cellIntegrals;    //integrals of elementary cells

	void refineNode(Node* pNode);

	Float cellSquare;
	Float fullIntegral;

	//disable copying
	Partition(const Partition&);
	Partition& operator=(const Partition&);
};




#endif /* _PARTITION_H_ */
