#ifndef _INDICATRIX_H_
#define _INDICATRIX_H_

#include "mathcompat.h"
#include "vector3.h"
#include "angle.h"
#include "coords.h"

#include "optics.h"


template <class T1, class T2>
class Indicatrix
{
public:

    Indicatrix(const Vector3& s_i, const Vector3& director_i);
    Float operator()(const Vector3& s_s);

protected:
	
	Vector3   s_i;
	Vector3   k_i;
	Vector3   director_i;   //director in k_i based coordinate system
	Angle     a_i;          //angle between k_i and n vectors 
	Vector3   e_i;

	Float     ei_n;
	Float     factor1;

	static const Float kCalculationEpsilon;
};

template <class T1, class T2>
const Float Indicatrix<T1, T2>::kCalculationEpsilon = kMachineEpsilon;


template <class T1, class T2>
Indicatrix<T1, T2>::Indicatrix(const Vector3& s_i_, const Vector3& director_i_) :
   	s_i(s_i_),
	k_i(),
	director_i(director_i_),
	a_i(s_i, director_i_),
	e_i(),
	ei_n(),
	factor1()
{
	k_i  = T1::k(s_i, a_i);
	e_i  = T1::e(k_i, director_i, a_i);

	ei_n = e_i*director_i;
	factor1 = Optics::s0/(T1::n(a_i)*T1::cosd(a_i));
}

template <class T1, class T2>
Float Indicatrix<T1, T2>::operator()(const Vector3& s_s)
{
	Angle a_s     = Angle(s_s, director_i);
	Vector3 k_s   = T2::k(s_s, a_s);
	Float cosd_s  = T2::cosd(a_s);

	Float res = factor1 * T2::f2(a_s) * T2::n(a_s)/(cosd_s*cosd_s*cosd_s);

	Vector3 q       = k_s - k_i;
	Vector3 q_par   = director_i*(q*director_i);
	Vector3 q_perp  = q - q_par;

	Vector3 a1 = q_perp; 

	if (a1.norm() < q.norm()*kCalculationEpsilon) { //along the optical axis

		//we need some a1 vector here, any unit vector that is perpendicular to n
		a1 = createSomePerpendicular(director_i);
	}

	a1.normalize();

	Vector3 a2 = crossProduct(director_i, a1);
	a2.normalize();   //to be sure
	
	Vector3 e_s  = T2::e(k_s, director_i, a_s);
	Float   es_n = e_s*director_i;

	Float es_a1, es_a2, ei_a1, ei_a2;

	es_a1  = e_s*a1;
	es_a2  = e_s*a2;
	ei_a1  = e_i*a1;
	ei_a2  = e_i*a2;

	res *= (ei_a1*ei_a1*es_n*es_n + 2.*ei_a1*es_n*ei_n*es_a1 + es_a1*es_a1*ei_n*ei_n) /
		       (Optics::t1*q_perp*q_perp + q_par*q_par + Optics::add) + 
	       (ei_a2*ei_a2*es_n*es_n + 2.*ei_a2*es_n*ei_n*es_a2 + es_a2*es_a2*ei_n*ei_n) / 
	           (Optics::t2*q_perp*q_perp + q_par*q_par + Optics::add);

	return res;
}

typedef Indicatrix<Optics::EBeam, Optics::EBeam> IndicatrixEE;


#endif /* _INDICATRIX_H_ */
