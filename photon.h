#ifndef _PHOTON_H_
#define _PHOTON_H_

#ifdef WIN32
#include <random>
#else
#include <tr1/random>
#endif

#include <limits.h>

#include "vector3.h"
#include "matrix3.h"

class ExtLength;
class Partition;
class PartitionChunk;
class EscFunction;
class Angle;

typedef std::vector<Float> ValuesVector;

typedef std::tr1::mt19937 RngEngine;


class Photon

{
public:

	Photon();
    static void init(FreePathEE* length, Partition* partition, EscFunction* escFunction,
            unsigned long seed = 1000);

	void move();
	void scatter();

	Vector3   pos;
	Vector3   s_i;

	int scatterings;
	Float weight;
	Float fullIntegral;

protected:

    void createTransformToPartitionCoords(Matrix3& mtx, Vector3& nn, Angle& a_i);
    void selectPartitionChunk(const Float theta);
    void calcPartitionValues(const Vector3& nn);
	void choosePointInRect(Float& x, Float& y, const int rectIdx, const Float randX, const Float randY);


	static RngEngine	rng_engine;

	static inline Float random() { return (Float)Photon::rng_engine() / UINT_MAX; }


	static FreePathEE* s_length;
	static Partition* s_partition;
	static EscFunction* s_escFunction;

    //these are to simulate static behavior for a reference (without ugly pointer syntax)
	FreePathEE& length;
	Partition& partition;
	EscFunction& escFunction;


	PartitionChunk *m_chunk;

	ValuesVector    m_knotValues;
	ValuesVector    m_rectValues;

private:

    Photon(const Photon&);
    Photon& operator=(const Photon&);
};


#endif /* _PHOTON_H_ */
