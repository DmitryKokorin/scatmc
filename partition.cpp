#include <cmath>
#include <cstdio>
#include <memory.h>

#include "common.h"
#include "indicatrix.h"
#include "partition.h"



//////////////////////////////////////////////////////////////////////////////////////

const Float Partition::epsilon = 0.01;
int Partition::rectCount = 0;

Partition::Partition() :
	pRoot(NULL),
	data(NULL),
	cellIntegrals(NULL),
	cellSquare(0.),
	fullIntegral(0.),
	knots()
{
	pRoot = new Node(GreedRect(0,0, size-1, size-1));
	rectCount = 1;

	cellIntegrals = allocate2dArray<Float>(size-1, size-1);
	
	preparePartitionTree();
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


void Partition::preparePartitionTree()
{
	Float ** data = allocate2dArray<Float>(Partition::size, Partition::size);

	const Float phiStep   = 2*M_PI/Partition::size;
	const Float thetaStep =   M_PI/Partition::size;

	const Float idxIterations = 100;
	const Float idxThetaStep  = M_PI / idxIterations;

	for (int k = 0; k < idxIterations; ++k) {

		Direction d_i = Direction(k*idxThetaStep, 0.);
		Indicatrix ind(d_i);

		Vector3 oX = Vector3(1, 0, 0);
		Vector3 iVector = d_i.toVector3();
				
		Float   angle = acos(oX*iVector);

		Float m[3][3];  //rotation matrix

		if (angle) {

			Vector3 axis  = crossProduct(iVector, oX).normalize(); 

			{
				Float rcos = cos(angle);
				Float rsin = sin(angle);

				Float u = axis.x(); 
				Float v = axis.y(); 
				Float w = axis.z(); 

				m[0][0] =      rcos + u*u*(1-rcos);
				m[1][0] =  w * rsin + v*u*(1-rcos);
				m[2][0] = -v * rsin + w*u*(1-rcos);
				m[0][1] = -w * rsin + u*v*(1-rcos);
				m[1][1] =      rcos + v*v*(1-rcos);
				m[2][1] =  u * rsin + w*v*(1-rcos);
				m[0][2] =  v * rsin + u*w*(1-rcos);
				m[1][2] = -u * rsin + v*w*(1-rcos);
				m[2][2] =      rcos + w*w*(1-rcos);
			}

		}
		else {

			memset(&m, 0, sizeof(m));
			m[0][0] = m[1][1] = m[2][2] = 1.;
		}

		for (int i = 0; i < Partition::size; ++i)
			for (int j = 0; j < Partition::size; ++j) {

				/*
				Direction d_s = Direction(i*thetaStep, j*phiStep);
				data[j][i] = ind(d_s)*d_s.sintheta;
				*/

				Float t = i*thetaStep;
				Float p = j*phiStep;

				Vector3 rel = Direction(t,p).toVector3();
				Vector3 abs(m[0][0]*rel.x() + m[0][1]*rel.y() + m[0][2]*rel.z(),
							m[1][0]*rel.x() + m[1][1]*rel.y() + m[1][2]*rel.z(),
							m[2][0]*rel.x() + m[2][1]*rel.y() + m[2][2]*rel.z());

				Direction d_s = Direction(abs);
				data[j][i] = ind(d_s)*d_s.sintheta;
			}

		setData(data, (M_PI/Partition::size) * (2*M_PI/Partition::size));
		refine();

		//printf("%d nodes...\n", p.rectCount);
	}

	free2dArray(data);
}

void Partition::processTree()
{
}

