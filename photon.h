#ifndef _PHOTON_H_
#define _PHOTON_H_

#include <tr1/random>

#include "vector3.h"
#include "angle.h"

class ExtLength;
class Partition;


class Photon

{
public:
	
	Photon();
    static void init(ExtLength* length, Partition* partition, unsigned long seed = 1000);


	void move();
	void scatter();

	Vector3   pos;
	Vector3   s_i;
	Angle     a_i;

	int scatterings;
	Float weight;
	Float fullIntegral;

protected:

	static std::tr1::mt19937 rng_core;  
	static std::tr1::uniform_real<Float> dist; 

	// have to use this due to 40263 gcc bug in uniform_real realization
	static std::tr1::variate_generator<std::tr1::mt19937, std::tr1::uniform_real<Float> > rng;


	//these two are to simulate static behaviour for a reference (without ugly pointer syntax)
	static ExtLength* s_length;  
	ExtLength& length;   

	static Partition* s_partition;  
	Partition& partition;   


	static const int kThetaIterations = 1000;
	static const int kPhiIterations   = 1000;
	
	static Float probs[kThetaIterations*kPhiIterations];
};


#endif /* _PHOTON_H_ */
