#include <cmath>
#include <cstdio>
#include <memory.h>
//#include <map>

#include "common.h"
#include "optics.h"
#include "indicatrix.h"
#include "matrix3.h"
#include "coords.h"
#include "rect.h"
#include "partition.h"



////////////////////////////////////////////////////////////////////////////////
//
// theta - x
// phi   - y
//

const Float Partition::kEpsilon = 0.03;
const Float Partition::kThetaResolution = M_PI / Partition::kThetaSize;
const Float Partition::kPhiResolution   = M_PI / Partition::kPhiSize;
const Float Partition::kIterationStep   = 0.5*M_PI / Partition::kIterations;


Partition::Partition() :
	m_rectCount(0),
	m_rects(),
	m_knots(),
	m_root(NULL),
	m_data(NULL),
	m_cellIntegrals(NULL),
	m_knotsMap(),
	m_cellSquare(0.),
	m_fullIntegral(0.)
{
}

bool Partition::create()
{
	m_root = new Node(GreedRect(0,0, kThetaSize-1, kPhiSize-1));
	m_rectCount = 1;

	m_cellIntegrals = allocate2dArray<Float>(kThetaSize-1, kPhiSize-1);

	createPartitionTree();
	createRectsList();

	return true;
}

Partition::~Partition()
{
	delete m_root;
}

void Partition::setData(Float** const data, const Float& cellSquare)
{
	m_data = data;
	m_cellSquare = cellSquare;
	m_fullIntegral = 0.;

	for (int i = 0; i < kThetaSize-1; ++i)
		for (int j = 0; j < kPhiSize-1; ++j) {

			m_cellIntegrals[j][i] = approxIntegral(GreedRect(i, j, i+1, j+1));
			m_fullIntegral += m_cellIntegrals[j][i];
		}
}

void Partition::refine()
{
	refineNode(m_root);
}

void Partition::refineNode(Node* node)
{
	if (node->isLeaf()) {

		Float nodeIntegral    = integral(node->rect);
		//Float rectMaxError    = nodeIntegral*kEpsilon;
		Float rectMaxError    = std::max(nodeIntegral*kEpsilon, kEpsilon*m_fullIntegral*((Float)(node->rect.square)/m_root->rect.square));

		if (node->rect.canSplitX() && node->rect.canSplitY()) {

			Float xSplitApproxIntegral = approxIntegral(node->rect.leftHalf()) +
			                             approxIntegral(node->rect.rightHalf());

			Float ySplitApproxIntegral = approxIntegral(node->rect.topHalf()) +
			                             approxIntegral(node->rect.bottomHalf());

			Float xSplitError  = fabs(nodeIntegral - xSplitApproxIntegral);
			Float ySplitError  = fabs(nodeIntegral - ySplitApproxIntegral);


			if ((xSplitError > rectMaxError) && (xSplitError >= ySplitError)) {

				node->splitX();
				++m_rectCount;
			}
			else if ((ySplitError > rectMaxError) && (ySplitError > xSplitError)) {

				node->splitY();
				++m_rectCount;
			}
		}
		else if (node->rect.canSplitX()) {

			Float xSplitApproxIntegral = approxIntegral(node->rect.leftHalf()) +
			                             approxIntegral(node->rect.rightHalf());

			Float xSplitError  = fabs(nodeIntegral - xSplitApproxIntegral);

			if (xSplitError > rectMaxError) {

				node->splitX();
				++m_rectCount;
			}
		}
		else if (node->rect.canSplitY()) {

			Float ySplitApproxIntegral = approxIntegral(node->rect.topHalf()) +
			                             approxIntegral(node->rect.bottomHalf());

			Float ySplitError  = fabs(nodeIntegral - ySplitApproxIntegral);

			if (ySplitError > rectMaxError) {

				node->splitY();
				++m_rectCount;
			}
		}
	}

	if (!node->isLeaf()) {

		refineNode(node->pChild1);
		refineNode(node->pChild2);
	}
}

Float Partition::integral(const GreedRect& rect)
{
	Float res = 0.;

	for (int i = rect.x1; i < rect.x2; ++i)
		for (int j = rect.y1; j < rect.y2; ++j)
			res += m_cellIntegrals[j][i];

	return res;
}

Float Partition::approxIntegral(const GreedRect& rect)
{
	return 0.25 * (	m_data[rect.y1][rect.x1] +
					m_data[rect.y2][rect.x1] +
					m_data[rect.y1][rect.x2] +
					m_data[rect.y2][rect.x2]) *
		(rect.square * m_cellSquare);
}

void Partition::createPartitionTree()
{
	Float ** data = allocate2dArray<Float>(Partition::kThetaSize, Partition::kPhiSize);

	const Float phiStep   = Partition::kPhiResolution;
	const Float thetaStep = Partition::kThetaResolution;

	for (int k = 1; k < kIterations; ++k) { //don't refine partition in "k_i || n" case

		//angle between k_i and director
		Angle   a_i = Angle(k*kIterationStep);

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

        Float t, p;

		for (int i = 0; i < Partition::kThetaSize; ++i) {

			t = i*thetaStep;

			for (int j = 0; j < Partition::kPhiSize; ++j) {

			    p = j*phiStep;

				Vector3 s_s = Vector3(cos(t), sin(t)*sin(p), sin(t)*cos(p));
				Float val = ind(s_s)*sin(t);
				data[j][i]  = val*sin(t);
			}
		}

		setData(data, Partition::kThetaResolution * Partition::kPhiResolution);
		refine();
	}

	free2dArray(data);

	fprintf(stderr, "Partitioning done, %d rects...\n", m_rectCount);
}

void Partition::createRectsList()
{
	Rect::s_knots = &m_knots;
	processTreeNode(m_root);
}

void Partition::processTreeNode(Node* node)
{
	if (node->isLeaf()) {

		int  keys[4];
		int  indeces[4];
		Knot knots[4];

        keys[0] = node->rect.x1 + node->rect.y1*kThetaSize; //tl
        keys[1] = node->rect.x2 + node->rect.y1*kThetaSize; //tr
        keys[2] = node->rect.x1 + node->rect.y2*kThetaSize; //bl
        keys[3] = node->rect.x2 + node->rect.y2*kThetaSize; //br

        Float x1 = node->rect.x1 * kThetaResolution;
		Float x2 = node->rect.x2 * kThetaResolution;
		Float y1 = node->rect.y1 * kPhiResolution;
		Float y2 = node->rect.y2 * kPhiResolution;

		knots[0] = Knot(x1, y1);
		knots[1] = Knot(x2, y1);
		knots[2] = Knot(x1, y2);
		knots[3] = Knot(x2, y2);

		for (int i = 0; i < 4; ++i) {

			if (m_knotsMap.find(keys[i]) == m_knotsMap.end()) {

				indeces[i] = m_knots.size();
				m_knots.push_back(knots[i]);
				m_knotsMap[keys[i]] = indeces[i];
			}
			else {
				indeces[i] = m_knotsMap[keys[i]];
			}
		}

		m_rects.push_back(Rect(indeces[0], indeces[1], indeces[2], indeces[3]));	
	}
	else {

		processTreeNode(node->pChild1);
		processTreeNode(node->pChild2);
	}
}

bool Partition::load(const std::string& name)
{
	FILE* file = fopen(name.c_str(), "r");

	if (!file)
		return false;

	int res = 0;

	int knotsNumber;
	res = fscanf(file, "%d", &knotsNumber);

	for (int i = 0; i < knotsNumber; ++i) {

		Float x, y;
		res = fscanf(file, "%le\t%le", &x, &y);
		m_knots.push_back(Knot(x, y));
	}


	Rect::s_knots = &m_knots;

	unsigned long int rectsNumber;

	res = fscanf(file, "%lu", &rectsNumber);

	for (unsigned long int i = 0; i < rectsNumber; ++i) {

		int tl, tr, bl, br;
		res = fscanf(file, "%d\t%d\t%d\t%d", &tl, &tr, &bl, &br);
		m_rects.push_back(Rect(tl, tr, bl, br));
	}


	fclose(file);

	return true;
}

bool Partition::save(const std::string& name)
{
	FILE* file = fopen(name.c_str(), "w");

	if (!file)
		return false;

	//knots
	{
		KnotsVector::iterator i;

		fprintf(file, "%lu\n", (long unsigned int)m_knots.size());

		for (i = m_knots.begin(); i != m_knots.end(); ++i) {

			fprintf(file, "%.17e\t%.17e\n", (*i).x, (*i).y);
		}
	}

	//rects
	{
		RectsVector::iterator i;

		fprintf(file, "%lu\n", (long unsigned int)m_rects.size());

		for (i = m_rects.begin(); i != m_rects.end(); ++i) {

			fprintf(file, "%d\t%d\t%d\t%d\n", (*i).tl, (*i).tr, (*i).bl, (*i).br);
		}
	}


	fclose(file);

	return true;
}
