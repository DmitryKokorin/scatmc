#ifndef _FREEPATH_H 
#define _FREEPATH_H_

#include <omp.h>

#include "common.h"
#include "mathcompat.h"
#include "optics.h"
#include "indicatrix.h"
#include "coords.h"
#include "linterpol.h"


Float symmetrizeTheta(const Float theta);


template <class T>
void createFreePath(LinearInterpolation& li, const int kPoints = 400, const int kThetaIterations = 1000, const int kPhiIterations = 1000)
{

    li.resize(kPoints);
    li.setRange(0., 0.5*M_PI);


#ifdef EXPERIMENTAL
   	for (int i = 0; i < kPoints; ++i) {

   	    li[i] = Optics::l;
    }
#else

    const Float kThetaStep = M_PI / kThetaIterations;
	const Float kPhiStep   = 2. * M_PI / kPhiIterations;
	const Float kResolution = li.resolution();

	Float theta_i = 0.;

	Vector3 v2, v3;
	Vector3 nn;

	for (int i = 1; i < kPoints; ++i) {

		theta_i     = i*kResolution;
		Angle   a_i = Angle(theta_i);
		Vector3 s_i = createSomeDeviantVector(Optics::director, a_i).normalize();

		//create coordinate system

		v2 = crossProduct(s_i, Optics::director).normalize();
		v3 = crossProduct(v2, s_i);

		Matrix3 mtx = invert(createTransformMatrix(s_i, v2, v3));
		nn = mtx*Optics::director;

        Indicatrix<T, Optics::OBeam> indO = Indicatrix<T, Optics::OBeam>(Vector3(1., 0., 0.), nn);
		Indicatrix<T, Optics::EBeam> indE = Indicatrix<T, Optics::EBeam>(Vector3(1., 0., 0.), nn);

		Float integral = 0.;

		#pragma omp parallel
		{
			Float t_integral = 0;

			#pragma omp for
			for (int j = 0; j < kThetaIterations; ++j) {

				Float theta_s = j*kThetaStep;
				Float phi_s   = 0.;

				for (int k = 0; k < kPhiIterations; ++k, phi_s += kPhiStep) {

					Vector3 s_s = Vector3(  cos(theta_s),
											sin(theta_s)*sin(phi_s),
											sin(theta_s)*cos(phi_s));

					Angle a_s   = Angle(s_s, nn);

					t_integral += sin(theta_s) * (  indO(s_s)*Optics::OBeam::cosd(a_s)/Optics::OBeam::f2(a_s) +
					                                indE(s_s)*Optics::EBeam::cosd(a_s)/Optics::EBeam::f2(a_s) ) /T::cosd(a_i);
				}
			}

            #pragma omp critical
			integral += t_integral;
		}

		li[i] = 1./(integral*kThetaStep*kPhiStep);
	}

	li[0] = li[1];

#endif
}


#endif /* _FREEPATH_H_ */
