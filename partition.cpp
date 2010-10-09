#include <cmath>
#include <cstdio>

#include "common.h"

#include "partition.h"

#if !defined NULL
#define NULL 0
#endif


///////////////////////////////////////////////////////////////////////////////


Node::Node(const GreedRect& rect_) :
	rect(rect_),
	pParent(NULL),
	pChild1(NULL),
	pChild2(NULL)
{
}

Node::~Node()
{ 
	delete pChild1;
	delete pChild2;
}

bool Node::isLeaf()
{
	return (NULL == pChild1) && (NULL == pChild2); 
}

bool Node::splitX() 
{ 
	if (!rect.canSplitX())
		return false;

	pChild1 = new Node(rect.leftHalf());
	pChild2 = new Node(rect.rightHalf());
	pChild1->pParent = this;
	pChild2->pParent = this;

	return true;
}

bool Node::splitY() 
{ 
	if (!rect.canSplitY())
		return false;

	pChild1 = new Node(rect.topHalf());
	pChild2 = new Node(rect.bottomHalf());
	pChild1->pParent = this;
	pChild2->pParent = this;

	return true;
}




//////////////////////////////////////////////////////////////////////////////////////

const Float Partition::epsilon = 0.01;
int Partition::rectCount = 0;

Partition::Partition() :
	pRoot(NULL),
	data(NULL),
	cellIntegrals(NULL),
	cellSquare(0.),
	fullIntegral(0.)
{
	pRoot = new Node(GreedRect(0,0, size-1, size-1));
	rectCount = 1;

	cellIntegrals = allocate2dArray<Float>(size-1, size-1);
}

Partition::~Partition()
{
	delete pRoot;
	free2dArray<Float>(cellIntegrals);
}

void Partition::setData(Float** const data_, const Float& cellSquare_)
{
	Partition::data = data_;
	cellSquare = cellSquare_;
	fullIntegral = 0.;

	for (int i = 0; i < size-1; ++i)
		for (int j = 0; j < size-1; ++j) {

			cellIntegrals[j][i] = approxIntegral(GreedRect(i, j, i+1, j+1));
			fullIntegral += cellIntegrals[j][i];
		}
			
}

void Partition::refine()
{
	refineNode(pRoot);
}

void Partition::refineNode(Node* pNode)
{
	if (pNode->isLeaf()) {

		Float nodeIntegral    = integral(pNode->rect);
		Float rectMaxError    = nodeIntegral*epsilon;

		if (pNode->rect.canSplitX() && pNode->rect.canSplitY()) {

			Float xSplitApproxIntegral = approxIntegral(pNode->rect.leftHalf()) + approxIntegral(pNode->rect.rightHalf());
			Float ySplitApproxIntegral = approxIntegral(pNode->rect.topHalf()) + approxIntegral(pNode->rect.bottomHalf());
 
			Float xSplitError  = fabs(nodeIntegral - xSplitApproxIntegral);
			Float ySplitError  = fabs(nodeIntegral - ySplitApproxIntegral);


			if ((xSplitError > rectMaxError) && (xSplitError > ySplitError)) {

				pNode->splitX();
				++rectCount;
			}
			else if ((ySplitError > rectMaxError) && (ySplitError > xSplitError)) {

				pNode->splitY();
				++rectCount;
			}
		}
		else if (pNode->rect.canSplitX()) {

			Float xSplitApproxIntegral = approxIntegral(pNode->rect.leftHalf()) + approxIntegral(pNode->rect.rightHalf());
			Float xSplitError  = fabs(nodeIntegral - xSplitApproxIntegral);

			if (xSplitError > rectMaxError) {

				pNode->splitX();
				++rectCount;
			}

		}
		else if (pNode->rect.canSplitY()) {

			Float ySplitApproxIntegral = approxIntegral(pNode->rect.topHalf()) + approxIntegral(pNode->rect.bottomHalf());
			Float ySplitError  = fabs(nodeIntegral - ySplitApproxIntegral);

			if (ySplitError > rectMaxError) {

				pNode->splitY();
				++rectCount;
			}
		}
	}

	if (!pNode->isLeaf()) {

		refineNode(pNode->pChild1);
		refineNode(pNode->pChild2);
	}
}

Float Partition::integral(const GreedRect& rect)
{
	/*
	//rough estimation
	Float res;
	
	if (rect.canSplitX() && rect.canSplitY()) {

		res = approxIntegral(rect.leftHalf().topHalf()) + 
		      approxIntegral(rect.leftHalf().bottomHalf()) + 
			  approxIntegral(rect.rightHalf().topHalf()) + 
              approxIntegral(rect.rightHalf().bottomHalf());
	}
	else {
		
		res = approxIntegral(rect);
	}
	*/
	
	Float res = 0.;
	
	for (int i = rect.x1; i < rect.x2; ++i)
		for (int j = rect.y1; j < rect.y2; ++j)
			res += cellIntegrals[j][i];

	return res;
}

Float Partition::approxIntegral(const GreedRect& rect)
{
	return 0.25 * (data[rect.y1][rect.x1] + data[rect.y2][rect.x1] + data[rect.y1][rect.x2] + data[rect.y2][rect.x2]) * (rect.square * cellSquare);
}
