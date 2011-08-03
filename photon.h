#ifndef _PHOTON_H_
#define _PHOTON_H_

#ifdef WIN32
#include <random>
#else
#include <tr1/random>
#endif

#include <limits.h>

#include "vector3.h"
#include "angle.h"

class ExtLength;
class Partition;
class PartitionChunk;
class EscFunction;

typedef std::vector<Float> ValuesVector;

typedef std::tr1::mt19937 RngEngine;


class Photon

{
public:

	Photon();
    static void init(FreePath* length, Partition* partition, EscFunction* escFunction,
            unsigned long seed = 1000);

	void move();
	void scatter();

	Vector3   pos;
	Vector3   s_i;
//	Angle     a_i;

	int scatterings;
	Float weight;
	Float fullIntegral;

protected:

	void choosePointInRect(Float& x, Float& y, const int rectIdx, const Float randX, const Float randY);


	static RngEngine	rng_engine;

	static inline Float random() { return (Float)Photon::rng_engine() / /*ULONG_MAX*/UINT_MAX; }


	static FreePath* s_length;
	static Partition* s_partition;
	static EscFunction* s_escFunction;

    //these are to simulate static behavior for a reference (without ugly pointer syntax)
	FreePath& length;
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
