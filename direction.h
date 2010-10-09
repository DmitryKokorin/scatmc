#ifndef _DIRECTION_H_
#define _DIRECTION_H_

#include <cmath>

#include "vector3.h"

//helper struct 

struct Direction
{
	Direction() :
		theta(0.),
		phi(0.),
		sintheta(0.),
		costheta(0.),
		sinphi(0.),
		cosphi(0.),
		sin2theta(0.),
		cos2theta(0.)
	{
	}


    Direction(const Float theta_, const Float phi_) :
		theta(theta_),
		phi(phi_),
		sintheta(0.),
		costheta(0.),
		sinphi(0.),
		cosphi(0.),
		sin2theta(0.),
		cos2theta(0.)
	{
		theta     = theta_;
		phi       = phi_;
		sintheta  = sin(theta);
		costheta  = cos(theta);
		sinphi    = sin(phi);
		cosphi    = cos(phi);
		sin2theta = sintheta*sintheta;
		cos2theta = 1 - sin2theta;
	}

    Direction(const Float costheta_, const Float sintheta_, const Float cosphi_, const Float sinphi_) :
		theta(0.),
		phi(0.),
		sintheta(sintheta_),
		costheta(costheta_),
		sinphi(sinphi_),
		cosphi(cosphi_),
		sin2theta(0.),
		cos2theta(0.)
	{
		theta     = acos(costheta_);
		phi       = sinphi_ ? acos(costheta_) : 2.*M_PI - acos(costheta_);
		sin2theta = sintheta*sintheta;
		cos2theta = 1 - sin2theta;
	}

	Direction(const Vector3& vec) :
		theta(0.),
		phi(0.),
		sintheta(0.),
		costheta(0.),
		sinphi(0.),
		cosphi(0.),
		sin2theta(0.),
		cos2theta(0.)
	{
		theta = acos(vec.x());
		phi   = atan2(vec.y(), vec.z());
		sintheta  = sin(theta);
		costheta  = cos(theta);
		sinphi    = sin(phi);
		cosphi    = cos(phi);
		sin2theta = sintheta*sintheta;
		cos2theta = 1 - sin2theta;

	}

	
	Vector3 toVector3() const
	{
		return Vector3(costheta, sintheta*sinphi, sintheta*cosphi);
	}

	Float theta;
	Float phi;

	Float sintheta;
	Float costheta;
	Float sinphi;
	Float cosphi;
	
	Float sin2theta;
	Float cos2theta;
};


#endif /* _DIRECTION_H_ */
