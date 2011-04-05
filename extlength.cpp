#include <cmath>
#include <omp.h>

#include <cstdio>

#include "optics.h"
#include "indicatrix.h"
#include "extlength.h"
#include "coords.h"


const Float ExtLength::kResolution = 0.5 * M_PI / (kPoints-1);



bool ExtLength::create(const int kThetaIterations /*= 1000*/,
                       const int kPhiIterations /*= 1000*/)
{
   	const Float kThetaStep = M_PI / kThetaIterations;
	const Float kPhiStep   = 2. * M_PI / kPhiIterations;

	Float theta_i = 0.;

	Vector3 v2, v3;
	Vector3 nn;
	Vector3 k_i;

	for (int i = 1; i < kPoints; ++i) {

		theta_i = i*kResolution;
		Angle   a_i = Angle(theta_i);
		Vector3 s_i = createSomeDeviantVector(Optics::n, a_i).normalize();

		//create coordinate system
		
		v2 = crossProduct(s_i, Optics::n).normalize();
		v3 = crossProduct(v2, s_i);

		Matrix3 mtx = invert(createTransformMatrix(s_i, v2, v3));
		nn = mtx*Optics::n;

		k_i = Vector3(1., 0., 0.)*Optics::ne(a_i);
		
		Indicatrix ind = Indicatrix(k_i, nn);

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

					t_integral += sin(theta_s) * ind(s_s)*
								Optics::cosde(a_s)/Optics::cosde(a_i)/Optics::f2(a_s) ;
				}
			}

            #pragma omp critical
			integral += t_integral;
		}
		
		lengths[i] = 1./(integral*kThetaStep*kPhiStep);
	}

	lengths[0] = lengths[1];

	return true;
}

Float ExtLength::operator()(const Angle& a) const
{
//	if (a.theta < 0 || a.theta >= M_PI)
//		return 0.;
    Float theta = a.theta;

    theta = fabs(theta);
    theta = (theta >= M_PI) ? theta - M_PI : theta;


	theta = (theta < M_PI * 0.5) ? theta : (M_PI - theta);

	//locate index
	Float mu = theta / kResolution;
	int  idx = (int)(mu);

	mu  = mu - idx;  

	return lengths[idx]*(1. - mu) + lengths[idx+1] * mu;
}

bool ExtLength::load(const std::string& name)
{
	FILE* file = fopen(name.c_str(), "r");

	if (!file)
		return false;

	int maxPoints;
	int res = 0;

	res = fscanf(file, "%d", &maxPoints);

	if (EOF == res || maxPoints != kPoints) {

		fclose(file);
		return false;
	}

	Float theta;
	
	for (int i = 0; i < kPoints; ++i) {

		res = fscanf(file, "%le\t%le", &theta, &lengths[i]);

		if (EOF == res && i != kPoints) {

			fclose(file);
			return false;
		}
	}


	fclose(file);
	return true;
}

bool ExtLength::save(const std::string& name)
{
	FILE* file = fopen(name.c_str(), "w");

	if (!file)
		return false;


	fprintf(file, "%d\n", kPoints);

	for (int i = 0; i < kPoints; ++i) {

		fprintf(file, "%.17e\t%.17e\n", i*kResolution, lengths[i]);
	}

	fclose(file);

	return true;
}

