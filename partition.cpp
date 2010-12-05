#include <cmath>
#include <cstdio>
#include <memory.h>

#include "common.h"
#include "optics.h"
#include "indicatrix.h"
#include "matrix3.h"
#include "coords.h"
#include "rect.h"
#include "partition.h"



////////////////////////////////////////////////////////////////////////////////

const Float Partition::epsilon = 0.01;
int Partition::rectCount = 0;

Partition::Partition() :
	rects(),
	pRoot(NULL),
	data(NULL),
	cellIntegrals(NULL),
	cellSquare(0.),
	fullIntegral(0.)
{
}

bool Partition::create()
{
	pRoot = new Node(GreedRect(0,0, size-1, size-1));
	rectCount = 1;

	cellIntegrals = allocate2dArray<Float>(size-1, size-1);
	
	preparePartitionTree();

	return true;
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

			Float xSplitApproxIntegral = approxIntegral(pNode->rect.leftHalf()) +
											approxIntegral(pNode->rect.rightHalf());
											
			Float ySplitApproxIntegral = approxIntegral(pNode->rect.topHalf()) +
											approxIntegral(pNode->rect.bottomHalf());
 
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

			Float xSplitApproxIntegral = approxIntegral(pNode->rect.leftHalf()) +
											approxIntegral(pNode->rect.rightHalf());
											
			Float xSplitError  = fabs(nodeIntegral - xSplitApproxIntegral);

			if (xSplitError > rectMaxError) {

				pNode->splitX();
				++rectCount;
			}
		}
		else if (pNode->rect.canSplitY()) {

			Float ySplitApproxIntegral = approxIntegral(pNode->rect.topHalf()) +
											approxIntegral(pNode->rect.bottomHalf());
											
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
	return 0.25 * (	data[rect.y1][rect.x1] +
					data[rect.y2][rect.x1] +
					data[rect.y1][rect.x2] +
					data[rect.y2][rect.x2]) *
		(rect.square * cellSquare);
}

void Partition::preparePartitionTree()
{
	Float ** data = allocate2dArray<Float>(Partition::size, Partition::size);

	const Float phiStep   = M_PI/Partition::size;
	const Float thetaStep = M_PI/Partition::size;

	const Float idxIterations = 1000;

/*	FILE* f1 = fopen("idx1.txt", "w");
	FILE* f2 = fopen("idx2.txt", "w");
*/
	for (int k = 1; k < idxIterations; ++k) { //don't refine partition for "k_i || n" case

		//angle between k_i and director
		Angle   a_i = Angle(k*thetaStep);

		//here we construct some k_i, that has a_i angle to director
		Vector3 s_i = createSomeDeviantVector(Optics::n, a_i).normalize();

		//now we can create coordinate system
		Vector3 v2 = crossProduct(s_i, Optics::n).normalize();
		Vector3 v3 = crossProduct(s_i, v2).normalize();

		
		//create matrix
		Matrix3 mtx = createTransformMatrix(s_i, v2, v3);

		Vector3 n_i = mtx*Optics::n;
		Vector3 k_i = Vector3(1., 0., 0.)*Optics::ne(a_i);

		Indicatrix ind = Indicatrix(k_i, n_i);

		//calculate array values
		for (int i = 0; i < Partition::size; ++i) {

			Float t = i*thetaStep;

			for (int j = 0; j < Partition::size; ++j) {

				Float p = j*phiStep;

				Vector3 k_s = Vector3(cos(t), sin(t)*sin(p), sin(t)*cos(p));
				Float val = ind(k_s)*sin(t);
				data[j][i]  = val*sin(t);

				/*if (k == 25) {
					fprintf(f1, "%f %f %f\n", t, p, val);
				}

				if (k == 50) {
					fprintf(f2, "%f %f %f\n", t, p, val);
				}*/

			}
		}

		setData(data, (M_PI/Partition::size) * (M_PI/Partition::size));
		refine();
	}

	free2dArray(data);

/*	fclose(f1);
	fclose(f2);
*/
	fprintf(stderr, "Partitioning done, %d rects...\n", rectCount);
}

void Partition::processTreeNode(Node* pNode)
{
	if (pNode->isLeaf()) {

	}
	else {

		processTreeNode(pNode->pChild1);
		processTreeNode(pNode->pChild2);
	}
}

bool Partition::load(const std::string& name)
{
	FILE* file = fopen(name.c_str(), "r");

	if (!file) {

		return false;
	}

	fclose(file);

	return true;
}

bool Partition::save(const std::string& name)
{
	FILE* file = fopen(name.c_str(), "w");

	if (!file) {

		return false;
	}

	fclose(file);

	return true;
}

