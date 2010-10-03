#include <cmath>
#include <omp.h>

#include <cstdio>


#include "direction.h"
#include "optics.h"
#include "indicatrix.h"
#include "extlength.h"


using namespace Optics;

const Float ExtLength::kResolution = 0.5 * M_PI / (kPoints-1);



ExtLength::ExtLength(const int kThetaIterations /*= 1000*/, const int kPhiIterations /*= 1000*/)
{

   	const Float kThetaStep = M_PI / kThetaIterations;
	const Float kPhiStep   = 2. * M_PI / kPhiIterations;

	
	Float theta_i = 0.;


	for (int i = 0; i < kPoints; ++i, theta_i += kResolution) {

		Indicatrix ind = Indicatrix(Direction(theta_i, 0.));
		Float integral = 0.;

	
		#pragma omp parallel
		{
			Float t_integral = 0;
 
			#pragma omp for
			for (int j = 0; j < kThetaIterations; ++j) {
			
				Float theta_s = j*kThetaStep;
				Float phi_s   = 0.;

				Direction d = Direction(theta_s, phi_s);

				for (int k = 0; k < kPhiIterations; ++k, phi_s += kPhiStep) {
				
					d.phi    = phi_s;
					d.sinphi = sin(phi_s);
					d.cosphi = cos(phi_s);

					t_integral += d.sintheta * ind(d)*cosde(d)/f2(d) ;
				}
			}

            #pragma omp critical
			integral += t_integral;
		}
		
		lengths[i] =/* 1. /*/ (integral*kThetaStep*kPhiStep);
	}

	lengths[0] = lengths[1];  //dirty fix of NaN
}

ExtLength::~ExtLength()
{
}

Float ExtLength::operator()(const Direction& d) const//TODO: rewrite
{
	if (d.theta < 0 || d.theta >= M_PI)
		return 0.;

	Float theta = (d.theta < M_PI * 0.5) ? d.theta : (M_PI - d.theta);

	//locate index
	Float mu = theta / kResolution;
	int  idx = (int)(mu);

	mu  = mu - idx;  

	return lengths[idx]*(1. - mu) + lengths[idx+1] * mu;
}
