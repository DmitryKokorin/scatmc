#include <cmath>
#include <cstdio>

#include "optics.h"
#include "indicatrix.h"


using namespace Optics;


const Float Indicatrix::kCalculationEpsilon = kMachineEpsilon;



Indicatrix::Indicatrix(const Direction& d_i_) :
	d_i(d_i_),
	k_i(),
	ee_i(),
	een_i(),
	factor1()
{
	ee_i  = ee(d_i);
	een_i = een(d_i); 
	k_i   = ke(d_i);
	factor1 = s0/(ne(d_i)*cosde(d_i));
}

Indicatrix::~Indicatrix()
{
}

Float Indicatrix::operator()(const Direction& d_s)
{
	Float cosde_s = cosde(d_s);
	Float res = factor1 * ne(d_s)/(cosde_s*cosde_s*cosde_s);

	Vector3 qee = ke(d_s)- k_i;
	Float   qeeperp = sqrt(qee.y()*qee.y() + qee.z()*qee.z());
	
	Vector3 ee_s  = ee(d_s);
	Float   een_s = een(d_s);

	Float eea1, eea2, eeia1, eeia2;

	if (qeeperp > qee.norm()*kCalculationEpsilon) {

		Float  iqeeperp = 1. / qeeperp;

		eea1  = (-ee_s.y()*qee.y() - ee_s.z()*qee.z()) * iqeeperp;
		eea2  = ( ee_s.y()*qee.z() - ee_s.z()*qee.y()) * iqeeperp;
		eeia1 = (-ee_i.y()*qee.y() - ee_i.z()*qee.z()) * iqeeperp;
		eeia2 = ( ee_i.y()*qee.z() - ee_i.z()*qee.y()) * iqeeperp;
	}
	else {  

		// 0/0 indeterminateness, we can choose a1 and a2. 
		// a1 = (0,1,0), a2 = (0,0,1)

		eea1  = ee_s.y();
		eea2  = ee_s.z();
		eeia1 = ee_i.y();
		eeia2 = ee_i.z();
 	}
	
	res *= (eeia1*eeia1*een_s*een_s + 2.*eeia1*een_s*een_i*eea1 + eea1*eea1*een_i*een_i) /
		       (t1*qeeperp*qeeperp +	 qee.x()*qee.x() + add) + 
	       (eeia2*eeia2*een_s*een_s + 2.*eeia2*een_s*een_i*eea2 + eea2*eea2*een_i*een_i) / 
	           (t2*qeeperp*qeeperp +	 qee.x()*qee.x() + add);

	return res;
}
