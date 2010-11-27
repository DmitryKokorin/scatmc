#include <cmath>
#include <cstdio>

#include "optics.h"
#include "indicatrix.h"


using namespace Optics;


const Float Indicatrix::kCalculationEpsilon = kMachineEpsilon;



Indicatrix::Indicatrix(const Vector3& k_i_, const Vector3& n_) :
	k_i(k_i_),
	n(n_),
	a_i(k_i, n),
	ee_i(),
	een_i(),
	factor1()
{
	ee_i  = ee(k_i, n, a_i);
	een_i = een(a_i); 
	factor1 = s0/(ne(a_i)*cosde(a_i));
}

Indicatrix::~Indicatrix()
{
}

Float Indicatrix::operator()(const Vector3& k_s)
{
	Angle a_s = Angle(k_s, n);
	Float cosde_s = cosde(a_s);
	Float res = factor1 * ne(a_s)/(cosde_s*cosde_s*cosde_s);

	Vector3 qee     = k_s- k_i;
	Vector3 qeeperp = qee - n*(qee*n);

	Vector3 a1 = qeeperp; 

	if (a1.norm() < qee.norm()*kCalculationEpsilon) { //along the optical axis
		
		//we need some a1 vector here, any unit vector that is perpendicular to n

		if (fabs(n.y()) > kCalculationEpsilon || fabs(n.z()) > kCalculationEpsilon)
			a1 = Vector3(1., 0., 0.);
		else
			a1 = Vector3(0., 1., 0.);

		a1 -= n*(a1*n); 
	}

	a1.normalize();

	Vector3 a2 = crossProduct(n, a1);
	a2.normalize();   //to be sure
	
	Vector3 ee_s  = ee(k_s, n, a_s);
	Float   een_s = een(a_s);

	Float eea1, eea2, eeia1, eeia2;

	eea1  = ee_s*a1;
	eea2  = ee_s*a2;
	eeia1 = ee_i*a1;
	eeia2 = ee_i*a2;

	res *= (eeia1*eeia1*een_s*een_s + 2.*eeia1*een_s*een_i*eea1 + eea1*eea1*een_i*een_i) /
		       (t1*qeeperp*qeeperp +	 qee.x()*qee.x() + add) + 
	       (eeia2*eeia2*een_s*een_s + 2.*eeia2*een_s*een_i*eea2 + eea2*eea2*een_i*een_i) / 
	           (t2*qeeperp*qeeperp +	 qee.x()*qee.x() + add);

	return res;
}
