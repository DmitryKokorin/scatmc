#include <cmath>
#include <omp.h>

#include <cstdio>


//#include "direction.h"
#include "optics.h"
#include "indicatrix.h"
#include "extlength.h"


using namespace Optics;

const Float ExtLength::kResolution = 0.5 * M_PI / (kPoints-1);



ExtLength::ExtLength(const int kThetaIterations /*= 1000*/,
                     const int kPhiIterations /*= 1000*/)
{

   	const Float kThetaStep = M_PI / kThetaIterations;
	const Float kPhiStep   = 2. * M_PI / kPhiIterations;

	Float theta_i = 0.;
	Float phi_i = 0.;


	for (int i = 0; i < kPoints; ++i, theta_i += kResolution) {

		Angle   a_i = Angle(theta_i);
		Vector3 k_i = Vector3(a_i.costheta, a_i.sintheta*sin(phi_i), a_i.sintheta*cos(phi_i))*ne(a_i);
		Indicatrix ind = Indicatrix(k_i, n);

		Float integral = 0.;

	
		#pragma omp parallel
		{
			Float t_integral = 0;
 
			#pragma omp for
			for (int j = 0; j < kThetaIterations; ++j) {
			
				Float theta_s = j*kThetaStep;
				Float phi_s   = 0.;
			
				for (int k = 0; k < kPhiIterations; ++k, phi_s += kPhiStep) {
				
					Vector3 k_s = Vector3(cos(theta_s), sin(theta_s)*sin(phi_s), sin(theta_s)*cos(phi_s))*ne(a_i);
					Angle a_s   = Angle(k_s, n);

					t_integral += sin(theta_s) * ind(k_s)*cosde(a_s)/f2(a_s) ;
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

Float ExtLength::operator()(const Angle& a) const//TODO: rewrite
{
	if (a.theta < 0 || a.theta >= M_PI)
		return 0.;

	Float theta = (a.theta < M_PI * 0.5) ? a.theta : (M_PI - a.theta);

	//locate index
	Float mu = theta / kResolution;
	int  idx = (int)(mu);

	mu  = mu - idx;  

	return lengths[idx]*(1. - mu) + lengths[idx+1] * mu;
}
