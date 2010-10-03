#ifndef _OPTICS_H_
#define _OPTICS_H_

#include <cmath>

#include "direction.h"
#include "vector3.h"


namespace Optics {

//hardware constants
	
extern const Float kMachineEpsilon;

//physical constants

extern const Float eps_par;
extern const Float eps_perp;
extern const Float eps_a;

extern const Float K3;
extern const Float t1;
extern const Float t2;
extern const Float K1;
extern const Float K2;

extern const Float lambda;
extern const Float xi;

extern const Float T;
extern const Float Kb;

extern const Float H;
extern const Float Hi_alpha;

//precalculated constants

extern const Float s0;
extern const Float add;



inline Float ne(const Direction& d)
{
	return sqrt(eps_perp*eps_par / (eps_perp + eps_a*d.cos2theta));
}

inline Float een(const Direction& d)
{
	return eps_perp*d.sintheta/sqrt(eps_perp*eps_perp*d.sin2theta + eps_par*eps_par*d.cos2theta);
}


inline Float cosde(const Direction& d)
{
	return (eps_perp*d.sin2theta + eps_par*d.cos2theta) /                               \
		sqrt(eps_perp*eps_perp*d.sin2theta + eps_par*eps_par*d.cos2theta);
}

inline Float f2(const Direction& d)
{
	return (eps_perp*d.sin2theta + eps_par*d.cos2theta) *                               \
		   (eps_perp*eps_perp*d.sin2theta + eps_par*eps_par*d.cos2theta) /              \
		   (eps_par*eps_perp*eps_perp);
}


inline Vector3 ke(const Direction& d)
{
	Float  n = ne(d);
	return  (d.toVector3())*=n;
}


__inline__ __attribute__((always_inline)) Vector3 ee(const Direction& d) 
{
	if (fabs(d.sintheta) > kMachineEpsilon) {
	
		Float iNo = (eps_perp + eps_a*d.cos2theta) / 
			(eps_perp*d.sintheta*sqrt(eps_perp*eps_perp*d.sin2theta + eps_par*eps_par*d.cos2theta));

		return (Vector3(eps_perp, 0., 0.) - ne(d)*d.costheta*ke(d)) * iNo;
	}
	else {  // along the optical axis, so we use an expression of the ordinary beam polarization
		
        return crossProduct(d.toVector3(), Vector3(1., 0., 0.)).normalize();
	} 
} 

}  //namespace Optics

#endif /* _OPTICS_H_ */
