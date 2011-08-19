#ifndef _CHANNEL_H_
#define _CHANNEL_H_

#include "common.h"
#include "linterpol.h"
#include "vector3.h"
#include "matrix3.h"
#include "optics.h"
#include "indicatrix.h"

template <class T>
void createEChannelProb(LinearInterpolation& li, const int kPoints = 400, const int kThetaIterations = 1000, const int kPhiIterations = 1000)
{

    li.resize(kPoints);
    li.setRange(0., 0.5*M_PI);


#ifdef EXPERIMENTAL
   	for (int i = 0; i < kPoints; ++i) {

   	    li[i] = 1.0;
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

		Float oIntegral = 0.;
		Float eIntegral = 0.;

		#pragma omp parallel
		{
			Float t_oIntegral = 0.;
			Float t_eIntegral = 0.;

			#pragma omp for
			for (int j = 0; j < kThetaIterations; ++j) {

				Float theta_s = j*kThetaStep;
				Float phi_s   = 0.;
				Float sintheta_s = sin(theta_s);
				Float costheta_s = cos(theta_s);


				for (int k = 0; k < kPhiIterations; ++k, phi_s += kPhiStep) {

					Vector3 s_s = Vector3(  costheta_s,
											sintheta_s*sin(phi_s),
											sintheta_s*cos(phi_s));

					Angle a_s   = Angle(s_s, nn);

					t_oIntegral += sintheta_s * indO(s_s);
					t_eIntegral += sintheta_s * indE(s_s);
				}
			}

            #pragma omp critical
            {
			    oIntegral += t_oIntegral;
			    eIntegral += t_eIntegral;
            }
		}

		li[i] = eIntegral /(eIntegral + oIntegral);
	}

	li[0] = li[1];

#endif
}


#endif /* end of include guard: _CHANNEL_H_ */

