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


class LinearInterpolation;
class Partition;
class PartitionChunk;
//class EscFunctionEE;
class Angle;


typedef std::vector<Float> ValuesVector;

typedef std::tr1::mt19937 RngEngine;


class Photon

{
public:

	Photon(const Vector3& s = Vector3(0., 0., 1.), const int channel_ = Optics::ECHANNEL);
    static void init(LinearInterpolation* length, Partition* partition, EscFunctionEE* escFunction,
            unsigned long seed = 1000);

	void move();
	void scatter();

	Vector3   pos;  //position
	Vector3   s_i;  //normalized wave vector

	int   scatterings;
	Float weight;
	Float fullIntegral;
	int   channel;

protected:

    void createTransformToPartitionCoords(Matrix3& mtx, Vector3& nn, Angle& a_i);
    void selectPartitionChunk(const Float theta);
    void calcPartitionValues(const Vector3& nn);
	void choosePointInRect(Float& x, Float& y, const int rectIdx, const Float randX, const Float randY);


	static RngEngine	rng_engine;

	static inline Float random() { return (Float)Photon::rng_engine() / UINT_MAX; }


	static LinearInterpolation* s_length;
	static Partition*           s_partition;
	static EscFunctionEE*         s_escFunction;

    //these are to simulate static behavior for a reference (without ugly pointer syntax)
	LinearInterpolation&    length;
	Partition&              partition;
	EscFunctionEE&          escFunction;


	PartitionChunk *m_chunk;

	ValuesVector    m_knotValues;
	ValuesVector    m_rectValues;

private:

    Photon(const Photon&);
    Photon& operator=(const Photon&);
};


#endif /* _PHOTON_H_ */
