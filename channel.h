#ifndef _CHANNEL_H_
#define _CHANNEL_H_

#include "common.h"
#include "linterpol.h"
#include "vector3.h"
#include "matrix3.h"
#include "optics.h"
#include "indicatrix.h"
#include "integrate.h"

namespace channel {

template <class T1, class T2>
class phiFunctor
{
public:

    phiFunctor(Indicatrix<T1, T2>& ind_,
                const Float sint_s_, const Float cost_s_) :
        ind(ind_),
        sint_s(sint_s_),
        cost_s(cost_s_)
    {}

    inline Float operator()(const Float phi)
    {  
        Vector3 s_s = Vector3(sint_s*cos(phi), sint_s*sin(phi), cost_s);

		return  sint_s * ind(s_s);
    }

protected:

    Indicatrix<T1, T2>& ind;

    Float sint_s;
    Float cost_s;
};


template <class T1, class T2>
class thetaFunctor
{
public:

    thetaFunctor(Indicatrix<T1, T2>& ind_) :
        ind(ind_)
    {}

    inline Float operator()(const Float theta)
    {  
        Float cost_s = cos(theta);
        Float sint_s = sin(theta);
     
        phiFunctor<T1, T2> functor(ind, sint_s, cost_s);
        Adapt s(1.0e-4);
        return s.integrate(functor, 0., 2*M_PI);
    }

protected:

    Indicatrix<T1, T2>& ind;
};
 
}

template <class T>
void createEChannelProb(LinearInterpolation& li, const int kPoints = 400)
{

    li.resize(kPoints);
    li.setRange(0., 0.5*M_PI);


#ifdef EXPERIMENTAL
   	for (int i = 0; i < kPoints; ++i) {

   	    li[i] = 1.0;
    }
#else

	const Float kResolution = li.resolution();

    #pragma omp parallel for
	for (int i = 1; i < kPoints; ++i) {

		Float theta_i = i*kResolution;

		Angle   a_i = Angle(theta_i);
		Vector3 s_i = createSomeDeviantVector(Optics::director, a_i).normalize();

		//create coordinate system

		Vector3 v2 = crossProduct(s_i, Optics::director).normalize();
		Vector3 v3 = crossProduct(v2, s_i);

		Matrix3 mtx = invert(createTransformMatrix(s_i, v2, v3));
		Vector3 nn = mtx*Optics::director;

        Indicatrix<T, Optics::OBeam> indO = Indicatrix<T, Optics::OBeam>(Vector3(1., 0., 0.), nn);
		Indicatrix<T, Optics::EBeam> indE = Indicatrix<T, Optics::EBeam>(Vector3(1., 0., 0.), nn);

        channel::thetaFunctor<T, Optics::OBeam> oFunctor(indO);
        channel::thetaFunctor<T, Optics::EBeam> eFunctor(indE);

        Adapt s(1.0e-4);
		Float oIntegral = s.integrate(oFunctor, 0., M_PI);
		Float eIntegral = s.integrate(eFunctor, 0., M_PI);

		li[i] = eIntegral /(eIntegral + oIntegral);
	}

	li[0] = li[1];

#endif
}


#endif /* end of include guard: _CHANNEL_H_ */

