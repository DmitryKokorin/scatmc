#ifndef _OPTICS_H_
#define _OPTICS_H_

#include "mathcompat.h"

#include "angle.h"
#include "vector3.h"


namespace Optics {

//physical constants

extern const Vector3 n;  //director

extern const Float eps_par;
extern const Float eps_perp;
extern const Float eps_a;

extern const Float K3;
extern const Float t1;
extern const Float t2;
extern const Float K1;
extern const Float K2;

extern const Float lambda;
extern const Float k0;
extern const Float xi;

extern const Float T;
extern const Float Kb;

extern const Float H;
extern const Float Hi_alpha;

//precalculated constants

extern const Float s0;
extern const Float add;


#ifdef EXPERIMENTAL

extern const Float g;

extern const Float l;				//free path in cm
extern const Float ls;

#endif




inline Float ne(const Angle& a)
{
	return sqrt(eps_perp*eps_par / (eps_perp + eps_a*a.cos2theta));
}

inline Float een(const Angle& a)
{
	return eps_perp*a.sintheta/sqrt(eps_perp*eps_perp*a.sin2theta + eps_par*eps_par*a.cos2theta);
}


inline Float cosde(const Angle& a)
{
	return (eps_perp*a.sin2theta + eps_par*a.cos2theta) /                               \
		sqrt(eps_perp*eps_perp*a.sin2theta + eps_par*eps_par*a.cos2theta);
}

inline Float f2(const Angle& a)
{
	return (eps_perp*a.sin2theta + eps_par*a.cos2theta) *                               \
		   (eps_perp*eps_perp*a.sin2theta + eps_par*eps_par*a.cos2theta) /              \
		   (eps_par*eps_perp*eps_perp);
}

inline Vector3 ke(const Vector3& direction, const Angle& a)
{
	Float nn = ne(a);
	return Vector3(direction).normalize()*nn;
//    return Vector3(direction).normalize()*nn*k0;
}

inline Vector3 ke(const Vector3& direction, const Vector3& nn)
{
	return ke(direction, Angle(direction, nn));
}

//s is a unit vector, s = k/|k|
inline Vector3 ee(const Vector3& k, const Vector3& nn, const Angle& a)
{
	if (fabs(a.sintheta) > kMachineEpsilon) {

		Vector3 s = Vector3(k).normalize();
//		Float iNo = 1./ (a.sintheta*sqrt(eps_perp*eps_perp*a.sin2theta + eps_par*eps_par*a.cos2theta));
//		return (s*eps_par*a.costheta - n*(eps_par*a.cos2theta + eps_perp*a.sin2theta)) * iNo; 
		return (s*eps_par*a.costheta - nn*(eps_par*a.cos2theta + eps_perp*a.sin2theta)).normalize();
	}
	else { //along the optical axis, so we use the expression for the ordinary beam polarization

		return crossProduct(k, nn).normalize();
	}
}

inline Vector3 ee(const Vector3& s, const Vector3& n)
{
	return ee(s, n, Angle(s, n));
}

}  //namespace Optics

#endif /* _OPTICS_H_ */
