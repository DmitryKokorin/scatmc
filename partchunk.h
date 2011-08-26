#ifndef _PARTCHUNK_H_
#define _PARTCHUNK_H_

#include <vector>
#include <string>
#include <map>

#include "common.h"
#include "node.h"
#include "rect.h"
#include "optics.h"
#include "indicatrix.h"
#include "matrix3.h"
#include "coords.h"


typedef std::vector<Rect> RectsVector;
typedef std::vector<Knot> KnotsVector;

class PartitionChunk
{
public:

    PartitionChunk();
    virtual ~PartitionChunk();

    template <class T>
    bool create(const Float minAngle, const Float maxAngle, const int iterations);
   
    bool load(FILE *file);
    bool save(FILE *file);

	void setData(Float** const data, const Float& cellSquare);
	void refine();

	inline size_t getRectsCount() const { return m_rects.size(); }
	inline size_t getKnotsCount() const { return m_knots.size(); }
    inline bool  isAngleInRange(const Float angle) const { return angle > m_minAngle && angle <= m_maxAngle;}

	static const Float kEpsilon;
	static const int   kThetaDegree = 10;           
	static const int   kThetaSize   = (1 << kThetaDegree) + 1;    //1025
	static const int   kPhiDegree   = 10;                     
	static const int   kPhiSize     = (1 << kPhiDegree) + 1;      //1025


	static const Float kThetaResolution;
	static const Float kPhiResolution;

	
	
	RectsVector m_rects;
	KnotsVector m_knots;


private:
	
	Float integral(const GreedRect& rect);
	Float approxIntegral(const GreedRect& rect);
	Float rectError(const GreedRect& rect);

	template <class T>
	void createPartitionTree();

	void createRectsList();

	void cleanUp();
	
	Float m_minAngle;
	Float m_maxAngle;
	int   m_iterations;

	Float m_iterationStep; 

	
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
	PartitionChunk(const PartitionChunk&);
	PartitionChunk& operator=(const PartitionChunk&);
};

template <class T>
bool PartitionChunk::create(const Float minAngle, const Float maxAngle, const int iterations)
{
    m_minAngle   = minAngle;
    m_maxAngle   = maxAngle;
    m_iterations = iterations;

    m_iterationStep = (m_maxAngle - m_minAngle) / m_iterations;


	m_root = new Node(GreedRect(0,0, kThetaSize-1, kPhiSize-1));

	m_cellIntegrals = allocate2dArray<Float>(kThetaSize-1, kPhiSize-1);

	createPartitionTree<T>();
	createRectsList();

	cleanUp();

	return true;
}

template <class T>
void PartitionChunk::createPartitionTree()
{
	Float ** data = allocate2dArray<Float>(PartitionChunk::kThetaSize, PartitionChunk::kPhiSize);

	for (int k = 0; k < m_iterations; ++k) {

		//angle between k_i and director
		Angle   a_i = Angle(m_minAngle + k*m_iterationStep);

		//here we construct some k_i, that has a_i angle to director
		Vector3 s_i = createSomeDeviantVector(Optics::director, a_i).normalize();

		//now we can create coordinate system
		Vector3 v2 = crossProduct(s_i, Optics::director).normalize();
		Vector3 v3 = crossProduct(s_i, v2).normalize();


		//create matrix
		Matrix3 mtx = createTransformMatrix(v2, v3, s_i);

		Vector3 nn = mtx*Optics::director;
		Vector3 ss_i = Vector3(0., 0., 1.);

		T ind(ss_i, nn);

		//calculate array values

        Float t, p;

		for (int i = 0; i < kThetaSize; ++i) {

			t = i*kThetaResolution;

			for (int j = 0; j < kPhiSize; ++j) {

			    p = j*kPhiResolution;

				Vector3 ss_s = Vector3(sin(t)*cos(p), sin(t)*sin(p), cos(t));
				data[j][i]  = ind(ss_s)*sin(t);
			}
		}

		setData(data, kThetaResolution * kPhiResolution);
		refine();
	}

	free2dArray(data);
}





#endif /* _PARTCHUNK_H_ */
