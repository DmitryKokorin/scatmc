#include <cmath>
#include <cstdio>

#include "optics.h"
#include "coords.h"
#include "indicatrix.h"


const Float IndicatrixBase::kCalculationEpsilon = kMachineEpsilon;



IndicatrixBase::IndicatrixBase(const Vector3& s_i_, const Vector3& director_i_) :
	s_i(s_i_),
	k_i(),
	director_i(director_i_),
	a_i(s_i, director_i_),
	e_i(),
	ei_n(),
	factor1()
{
}

/////////////////////////////////////////////////////////

IndicatrixEE::IndicatrixEE(const Vector3& s_i_, const Vector3& director_i_) :
    IndicatrixBase(s_i_, director_i_)
{
	k_i  = Optics::EBeam::k(s_i, a_i);
	e_i  = Optics::EBeam::e(k_i, director_i, a_i);

	ei_n = e_i*director_i;
	factor1 = Optics::s0/(Optics::EBeam::n(a_i)*Optics::EBeam::cosd(a_i));
}


Float IndicatrixEE::operator()(const Vector3& s_s)
{
	Angle a_s     = Angle(s_s, director_i);
	Vector3 k_s   = Optics::EBeam::k(s_s, a_s);
	Float cosd_s  = Optics::EBeam::cosd(a_s);

	Float res = factor1 * Optics::EBeam::f2(a_s) * Optics::EBeam::n(a_s)/(cosd_s*cosd_s*cosd_s);

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
	
	Vector3 e_s  = Optics::EBeam::e(k_s, director_i, a_s);
	Float es_n = e_s*director_i;

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


/////////////////////////////////////////////////////////

IndicatrixEO::IndicatrixEO(const Vector3& s_i_, const Vector3& director_i_) :
    IndicatrixBase(s_i_, director_i_)
{
	k_i  = Optics::EBeam::k(s_i, a_i);
	e_i  = Optics::EBeam::e(k_i, director_i, a_i);

	ei_n = e_i*director_i; 
	factor1 = Optics::s0/(Optics::EBeam::n(a_i)*Optics::EBeam::cosd(a_i));
}

Float IndicatrixEO::operator()(const Vector3& s_s)
{
	Angle a_s     = Angle(s_s, director_i);
	Vector3 k_s   = Optics::OBeam::k(s_s, a_s);
	Float cosd_s  = Optics::OBeam::cosd(a_s);

	Float res = factor1 * Optics::OBeam::f2(a_s) * Optics::OBeam::n(a_s)/(cosd_s*cosd_s*cosd_s);

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
	
	Vector3 e_s = Optics::OBeam::e(k_s, director_i, a_s);
    Float es_n  = e_s*director_i;

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


/////////////////////////////////////////////////////////

IndicatrixOE::IndicatrixOE(const Vector3& s_i_, const Vector3& director_i_) :
    IndicatrixBase(s_i_, director_i_)
{
	k_i  = Optics::OBeam::k(s_i, a_i);
	e_i  = Optics::OBeam::e(k_i, director_i, a_i);

	ei_n = e_i*director_i; 
	factor1 = Optics::s0/(Optics::OBeam::n(a_i)*Optics::OBeam::cosd(a_i));
}


Float IndicatrixOE::operator()(const Vector3& s_s)
{
	Angle a_s     = Angle(s_s, director_i);
	Vector3 k_s   = Optics::EBeam::k(s_s, a_s);
	Float cosd_s  = Optics::EBeam::cosd(a_s);

	Float res = factor1 * Optics::EBeam::f2(a_s) * Optics::EBeam::n(a_s)/(cosd_s*cosd_s*cosd_s);

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
	
	Vector3 e_s  = Optics::EBeam::e(k_s, director_i, a_s);
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



///////////////////////////////////////////////////////

#ifdef EXPERIMENTAL

IndicatrixHenyey::IndicatrixHenyey(const Vector3& s_i_, const Vector3& director_i_) :
    IndicatrixBase(s_i_, director_i_)
{
}

Float IndicatrixHenyey::operator()(const Vector3& s_s)
{
    return Optics::HenyeyGreenstein(s_i*s_s);
}

#endif
