#ifndef _OPTICS_H_
#define _OPTICS_H_

#include "mathcompat.h"

#include "angle.h"
#include "vector3.h"


namespace Optics {

//physical constants

extern const Vector3 director;         //vector of director

extern const Float eps_par;
extern const Float eps_perp;
extern const Float eps_a;
extern const Float _no;

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


inline Float HenyeyGreenstein(const Float ct)
{
    Float t = 1./(1 + Optics::g*Optics::g - 2.*Optics::g*ct);
    return 0.5*(1 - Optics::g*Optics::g)*t*sqrt(t);
}

#endif

/////////////////////////////////////////////////////

enum Channel
{
    OCHANNEL,
    ECHANNEL
};


class OBeam {
public:    

    static const int channel = OCHANNEL;

    static inline Float n(const Angle& /*a*/)
    {
        return _no;
    }

    static inline Vector3 k(const Vector3& direction, const Angle& a)
    {
        return Vector3(direction).normalize()*n(a);
    }

    static inline Float cosd(const Angle& /*a*/)
    {
	    return 1.0;
    }

    static inline Float f2(const Angle& /*a*/)
    {
	    return 1.0;
    }

    static inline Vector3 e(const Vector3& s, const Vector3& director_, const Angle& /*a*/)
    {
	    return crossProduct(s, director_).normalize();
    }
};


class EBeam {
public:

    static const int channel = ECHANNEL;

    //refraction index
    static inline Float n(const Angle& a)
    {
    	return sqrt(eps_perp*eps_par / (eps_perp + eps_a*a.cos2theta));
    }

    static inline Vector3 k(const Vector3& direction, const Angle& a)
    {
	    Float nn = n(a);
    	return Vector3(direction).normalize()*nn;
    }

    static inline Float cosd(const Angle& a)
    {
	    return (eps_perp*a.sin2theta + eps_par*a.cos2theta) /                               \
		    sqrt(eps_perp*eps_perp*a.sin2theta + eps_par*eps_par*a.cos2theta);
    }

    static inline Float f2(const Angle& a)
    {
	    return (eps_perp*a.sin2theta + eps_par*a.cos2theta) *                               \
		       (eps_perp*eps_perp*a.sin2theta + eps_par*eps_par*a.cos2theta) /              \
		    (eps_par*eps_perp*eps_perp);
    }

    //polarization vector
//    static inline Vector3 e(const Vector3& k, const Vector3& nn, const Angle& a)
    static inline Vector3 e(const Vector3& s, const Vector3& director_, const Angle& a)
    {
	    if (fabs(a.sintheta) > kMachineEpsilon) {

		    //Vector3 s = Vector3(k).normalize();
		    return (s*eps_par*a.costheta - director_*(eps_par*a.cos2theta + eps_perp*a.sin2theta)).normalize();
	    }
    	else { //along the optical axis, so we use the expression for the ordinary beam polarization

		    return Optics::OBeam::e(s, director_, a);
	    }
    }
};


}  //namespace Optics

#endif /* _OPTICS_H_ */
