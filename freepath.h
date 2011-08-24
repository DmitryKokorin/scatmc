#ifndef _FREEPATH_H 
#define _FREEPATH_H_

#include <omp.h>

#include "common.h"
#include "mathcompat.h"
#include "optics.h"
#include "indicatrix.h"
#include "coords.h"
#include "linterpol.h"
#include "integrate.h"


Float symmetrizeTheta(const Float theta);


namespace freepath {

template <class T>
class phiFunctor
{
public:

    phiFunctor(Indicatrix<T, Optics::OBeam>& indO_, Indicatrix<T, Optics::EBeam>& indE_,
                const Vector3& nn_, const Angle& a_i_,
                const Float sint_s_, const Float cost_s_) :
        indO(indO_),
        indE(indE_),
        nn(nn_),
        a_i(a_i_),
        sint_s(sint_s_),
        cost_s(cost_s_)
    {}

    inline Float operator()(const Float phi)
    {  
        Vector3 s_s = Vector3(sint_s*cos(phi), sint_s*sin(phi), cost_s);

		Angle a_s   = Angle(s_s, nn);

		return  sint_s * (  indO(s_s)*Optics::OBeam::cosd(a_s)/Optics::OBeam::f2(a_s) +
		                    indE(s_s)*Optics::EBeam::cosd(a_s)/Optics::EBeam::f2(a_s) ) /T::cosd(a_i);
    }

protected:

    Indicatrix<T, Optics::OBeam>& indO;
    Indicatrix<T, Optics::EBeam>& indE;

    const Vector3& nn;
    const Angle& a_i;

    Float sint_s;
    Float cost_s;
};


template <class T>
class thetaFunctor
{
public:

    thetaFunctor(Indicatrix<T, Optics::OBeam>& indO_, Indicatrix<T, Optics::EBeam>& indE_,
                    const Vector3& nn_, const Angle& a_i_) :
        indO(indO_),
        indE(indE_),
        nn(nn_),
        a_i(a_i_)
    {}

    inline Float operator()(const Float theta)
    {  
        Float cost_s = cos(theta);
        Float sint_s = sin(theta);
     
        phiFunctor<T> functor = phiFunctor<T>(indO, indE, nn, a_i, sint_s, cost_s);
        Adapt s(1.0e-8);
        return s.integrate(functor, 0., 2*M_PI);
    }

protected:

    Indicatrix<T, Optics::OBeam>& indO;
    Indicatrix<T, Optics::EBeam>& indE;

    const Vector3& nn;
    const Angle& a_i;
};
 
}


template <class T>
void createFreePath(LinearInterpolation& li, const int kPoints = 400)
{

    li.resize(kPoints);
    li.setRange(0., 0.5*M_PI);


#ifdef EXPERIMENTAL
   	for (int i = 0; i < kPoints; ++i) {

   	    li[i] = Optics::l;
    }
#else

	const Float kResolution = li.resolution();

    #pragma omp parallel for
	for (int i = 1; i < kPoints; ++i) {

		Float theta_i     = i*kResolution;
		Angle   a_i = Angle(theta_i);
		Vector3 s_i = createSomeDeviantVector(Optics::director, a_i).normalize();

		//create coordinate system

		Vector3 v2 = crossProduct(s_i, Optics::director).normalize();
		Vector3 v3 = crossProduct(v2, s_i);

		Matrix3 mtx = invert(createTransformMatrix(s_i, v2, v3));
		Vector3 nn = mtx*Optics::director;

        Indicatrix<T, Optics::OBeam> indO = Indicatrix<T, Optics::OBeam>(Vector3(1., 0., 0.), nn);
		Indicatrix<T, Optics::EBeam> indE = Indicatrix<T, Optics::EBeam>(Vector3(1., 0., 0.), nn);

        Adapt s(1.0e-8);
        freepath::thetaFunctor<T> functor = freepath::thetaFunctor<T> (indO, indE, nn, a_i);
        Float integral = s.integrate(functor, 0., M_PI);

		li[i] = 1./(integral);
	}

	li[0] = li[1];

#endif
}


#endif /* _FREEPATH_H_ */
