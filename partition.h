#ifndef _PARTITION_H_
#define _PARTITION_H_

#include <list>

#include "common.h"
#include "node.h"

struct Knot
{
	Knot(const Float x_, const Float y_, const int id_) :
		x(x_),
		y(y_),
		id(id_)
        {}

	Float x, y;
	int id;
};

class Partition
{
public:

    Partition();
    virtual ~Partition();

	void setData(Float** const data_, const Float& cellSquare_);
	void refine();

	static const Float epsilon;
	static const int   degree = 10;           
	static const int   size   = (1 << degree) + 1;    //1025x1025
	
	static int rectCount;


private:
	
	Float integral(const GreedRect& rect);
	Float approxIntegral(const GreedRect& rect);

	void preparePartitionTree();
	void processTree();	

	Node* pRoot;
	Float**  data;             //knots
	Float**  cellIntegrals;    //integrals of elementary cells

	void refineNode(Node* pNode);

	Float cellSquare;
	Float fullIntegral;

	std::list<Knot> knots;

	//disable copying
	Partition(const Partition&);
	Partition& operator=(const Partition&);
};




#endif /* _PARTITION_H_ */
