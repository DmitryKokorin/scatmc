#ifndef _ANGLE_H_
#define _ANGLE_H_

#include <cmath>

#include "vector3.h"

//helper struct 

struct Angle
{
	Angle() :
		theta(0.),
		sintheta(0.),
		costheta(0.),
		sin2theta(0.),
		cos2theta(0.)
	{
	}


    Angle(const Float theta_) :
		theta(theta_),
		sintheta(0.),
		costheta(0.),
		sin2theta(0.),
		cos2theta(0.)
	{
		theta     = theta_;
		sintheta  = sin(theta);
		costheta  = cos(theta);
		sin2theta = sintheta*sintheta;
		cos2theta = 1 - sin2theta;
	}

	Angle(const Float costheta_, const Float sintheta_) :
		theta(0.),
		sintheta(sintheta_),
		costheta(costheta_),
		sin2theta(0.),
		cos2theta(0.)
	{
		theta     = acos(costheta_);
		sin2theta = sintheta*sintheta;
		cos2theta = 1 - sin2theta;
	}

	Angle(const Vector3& v1, const Vector3& v2) :
		theta(0.),
		sintheta(0.),
		costheta(0.),
		sin2theta(0.),
		cos2theta(0.)
	{
		Float n1 = v1.norm();
		Float n2 = v2.norm();

		if (n1 > kMachineEpsilon && n2 > kMachineEpsilon) {

			costheta = v1*v2/(n1*n2);
			theta = acos(costheta);
			sintheta = sin(theta);
			cos2theta = costheta*costheta;
			sin2theta = 1 - cos2theta;
		}
	}

	Float theta;

	Float sintheta;
	Float costheta;
	
	Float sin2theta;
	Float cos2theta;
};


#endif /* _ANGLE_H_ */
