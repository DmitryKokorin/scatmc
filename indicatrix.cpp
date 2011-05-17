#include <cmath>
#include <cstdio>

#include "optics.h"
#include "coords.h"
#include "indicatrix.h"


const Float Indicatrix::kCalculationEpsilon = kMachineEpsilon;


Indicatrix::Indicatrix(const Vector3& s_i_, const Vector3& n_i_) :
	s_i(s_i_),
	k_i(),
	n_i(n_i_),
	a_i(s_i, n_i_),
	ee_i(),
	een_i(),
	factor1()
{
	k_i   = Optics::ke(s_i, a_i);
	ee_i  = Optics::ee(k_i, n_i, a_i);
	een_i = Optics::een(a_i); 
	factor1 = Optics::s0/(Optics::ne(a_i)*Optics::cosde(a_i));
}

Indicatrix::~Indicatrix()
{
}

Float Indicatrix::operator()(const Vector3& s_s)
{
	Angle a_s     = Angle(s_s, n_i);
	Vector3 k_s   = Optics::ke(s_s, a_s);
	Float cosde_s = Optics::cosde(a_s);

	Float res = factor1 * Optics::f2(a_s) * Optics::ne(a_s)/(cosde_s*cosde_s*cosde_s);

	Vector3 qee     = k_s - k_i;
	Vector3 qeepar  = n_i*(qee*n_i);
	Vector3 qeeperp = qee - qeepar;

	Vector3 a1 = qeeperp; 

	if (a1.norm() < qee.norm()*kCalculationEpsilon) { //along the optical axis

		//we need some a1 vector here, any unit vector that is perpendicular to n

		a1 = createSomePerpendicular(n_i);
	}

	a1.normalize();

	Vector3 a2 = crossProduct(n_i, a1);
	a2.normalize();   //to be sure
	
	Vector3 ee_s  = Optics::ee(k_s, n_i, a_s);
	Float   een_s = Optics::een(a_s);

	Float eea1, eea2, eeia1, eeia2;

	eea1  = ee_s*a1;
	eea2  = ee_s*a2;
	eeia1 = ee_i*a1;
	eeia2 = ee_i*a2;

	res *= (eeia1*eeia1*een_s*een_s + 2.*eeia1*een_s*een_i*eea1 + eea1*eea1*een_i*een_i) /
		       (Optics::t1*qeeperp*qeeperp + qeepar*qeepar + Optics::add) + 
	       (eeia2*eeia2*een_s*een_s + 2.*eeia2*een_s*een_i*eea2 + eea2*eea2*een_i*een_i) / 
	           (Optics::t2*qeeperp*qeeperp + qeepar*qeepar + Optics::add);

	return res;
}
