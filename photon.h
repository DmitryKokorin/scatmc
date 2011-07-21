#ifndef _PHOTON_H_
#define _PHOTON_H_

#ifdef WIN32
#include <random>
#else
#include <tr1/random>
#endif

#include "vector3.h"
#include "angle.h"

class ExtLength;
class Partition;
class PartitionChunk;
class EscFunction;

typedef std::vector<Float> ValuesVector;

#define MAX_LONG	0xffffffff

//typedef std::tr1::mersenne_twister< long long, 32, 624, 397, 31, 0x9908b0df, 11, 7, 0x9d2c5680, 15, 0xefc60000, 18> mt19937_;
typedef std::tr1::mt19937 RngEngine;
//typedef mt19937_ RngEngine;
typedef std::tr1::uniform_real<Float> RngDistrib;
typedef std::tr1::variate_generator<RngEngine, RngDistrib> Rng;


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
	Angle     a_i;

	int scatterings;
	Float weight;
	Float fullIntegral;

protected:

	void choosePointInRect(Float& x, Float& y, const int rectIdx, const Float randX, const Float randY);


	static RngEngine	rng_engine;
	static RngDistrib	rng_distrib;

	// have to use this due to 40263 gcc bug in uniform_real realization
	static Rng* rng;

	static inline Float random() { return (Float)Photon::rng_engine() / MAX_LONG; }


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
