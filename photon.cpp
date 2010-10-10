#include <omp.h>

#include "extlength.h"
#include "indicatrix.h"
#include "photon.h"


const Float Photon::kThetaStep = M_PI / kThetaIterations;
const Float Photon::kPhiStep   = 2. * M_PI / kPhiIterations;

Float Photon::probs[kThetaIterations*kPhiIterations] = {0};
ExtLength* Photon::s_length       = NULL;
Partition* Photon::s_partition    = NULL;


using namespace std::tr1;

mt19937 Photon::rng_core = mt19937();
uniform_real<Float> Photon::dist = uniform_real<Float> (0., 1.); 

variate_generator<mt19937, uniform_real<Float> > Photon::rng = variate_generator<mt19937, uniform_real<Float> > (rng_core, dist);




void Photon::init(ExtLength* length_, Partition* partition_, unsigned long seed_)
{
	s_length    = length_;
	s_partition = partition_;

	Photon::rng_core.seed(seed_);
}


Photon::Photon() :
	pos(0.,0.,0.),
	d_i(M_PI*0.5, 0.),
	scatterings(0),
	weight(1.),
	length(*s_length),
	partition(*s_partition)
{
}


void Photon::move()
{
	Float rnd;
	#pragma omp critical
	{
		rnd = 1. - rng();
	}

	Float d   = -log(rnd) * length(d_i);

	pos += d*d_i.toVector3();
}

void Photon::scatter()
{
	//prepare probs array

	Indicatrix ind = Indicatrix(d_i);

	Float sum   = 0.;
	int    index = 0;
 
	for (int i = 0; i < kThetaIterations; ++i) {
			
		Float theta_s = i*kThetaStep;
		Float phi_s = 0.;

		Direction d = Direction(theta_s, phi_s);

		for (int j = 0; j < kPhiIterations; ++j, phi_s += kPhiStep) {
				
			d.phi    = phi_s;
			d.sinphi = sin(phi_s);
			d.cosphi = cos(phi_s);

			sum += d.sintheta * ind(d) * kPhiStep*kThetaStep;
			probs[index++] = sum;
		}
	}
	
	//choose scattering direction
	std::tr1::uniform_real<Float> probs_dist(0, sum);
	std::tr1::variate_generator<std::tr1::mt19937, std::tr1::uniform_real<Float> > rng_direction(rng_core, probs_dist);

	Float rnd = rng_direction();
	
	index = 0;
	while (probs[index] < rnd)
		index++;

	Float theta = (Float)(index / kPhiIterations)*kThetaStep;
	Float phi   = (Float)(index % kPhiIterations)*kPhiStep;
	
	d_i = Direction(theta, phi);

	//	printf("%d\n", scatterings);

	scatterings++;
}
