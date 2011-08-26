#ifndef _ESCFUNCTION_H_
#define _ESCFUNCTION_H_

#include <string>
#include <cstdio>
#include "mathcompat.h"

#include "optics.h"
#include "indicatrix.h"
#include "freepath.h"
#include "integrate.h"


class EscFunction
{
public:

    EscFunction();
    ~EscFunction();

    template <class T>
    bool create(const LinearInterpolation& FreePathO,
                const LinearInterpolation& FreePAthE,
                const ULong thetaSize,
                const ULong phiSize,
                const ULong zSize,
                const Float maxZ);

    Float operator()(const Float theta, const Float phi, const Float z) const;

    bool load(const std::string& name);
    bool save(const std::string& name);


private:



    void recalcSteps();

    Float m_maxZ;

    ULong m_thetaSize;
    ULong m_phiSize;
    ULong m_zSize;

    Float m_thetaStep;
    Float m_phiStep;
    Float m_zStep;

    Float ***m_array;

    EscFunction(const EscFunction&);
    EscFunction& operator=(const EscFunction&);
};

namespace escfunction {

template <class T1, class T2>
class phiFunctor
{
public:

    phiFunctor(T1& indO_, T2& indE_, const Float sint_s_, const Float cost_s_, const Float z_,
            const LinearInterpolation& FreePathO_,
            const LinearInterpolation& FreePathE_) :
        indO(indO_),
        indE(indE_),
        sint_s(sint_s_),
        cost_s(cost_s_),
        z(z_),
        FreePathO(FreePathO_),
        FreePathE(FreePathE_)
    {}

    inline Float operator()(const Float phi)
    {  
            Vector3 s_s = Vector3(sint_s*cos(phi), sint_s*sin(phi), cost_s);

            Float dist = z / cost_s;
            Float theta = symmetrizeTheta(Angle(s_s, Optics::director).theta);

            return (indO(s_s)*exp(dist/FreePathO(theta)) +
                    indE(s_s)*exp(dist/FreePathE(theta)) ) *  sint_s;

        return 0.;
    }

protected:

    T1& indO;
    T2& indE;

    Float sint_s;
    Float cost_s;
    Float z;

    const LinearInterpolation& FreePathO;
    const LinearInterpolation& FreePathE;
};

template <class T1, class T2>
class phiFunctorNorm
{
public:

    phiFunctorNorm(T1& indO_, T2& indE_, const Float sint_s_, const Float cost_s_) :
        indO(indO_),
        indE(indE_),
        sint_s(sint_s_),
        cost_s(cost_s_)
    {}


    inline Float operator()(const Float phi)
    {
        Vector3 s_s = Vector3(sint_s*cos(phi), sint_s*sin(phi), cost_s);

        return (indO(s_s) + indE(s_s))*sint_s;
    }

protected:

    T1& indO;
    T2& indE;

    Float sint_s;
    Float cost_s;
};


template <class T1, class T2>
class thetaFunctor
{
public:

    thetaFunctor(T1& indO_, T2& indE_, const Float z_,
            const LinearInterpolation& FreePathO_, 
            const LinearInterpolation& FreePathE_) :
        indO(indO_),
        indE(indE_),
        z(z_),
        FreePathO(FreePathO_),
        FreePathE(FreePathE_)
    {}

    inline Float operator()(const Float theta)
    {  
        Float cost_s = cos(theta);
        Float sint_s = sin(theta);
     
        if (cost_s < -kMachineEpsilon) {

            phiFunctor<T1, T2> functor(indO, indE, sint_s, cost_s, z, FreePathO, FreePathE);
            Adapt s(1.0e-6);
            return s.integrate(functor, 0., 2*M_PI);
        }

        return 0.;
    }

protected:

    T1& indO;
    T2& indE;

    Float z;

    const LinearInterpolation& FreePathO;
    const LinearInterpolation& FreePathE;
};


template <class T1, class T2>
class thetaFunctorNorm
{
public:

    thetaFunctorNorm(T1& indO_, T2& indE_) :
        indO(indO_),
        indE(indE_)
    {}

    inline Float operator()(const Float theta)
    {  
        Float cost_s = cos(theta);
        Float sint_s = sin(theta);
     
        phiFunctorNorm<T1, T2>  functorNorm(indO, indE, sint_s, cost_s);
        
        Adapt s(1.0e-6);
        return s.integrate(functorNorm, 0., 2*M_PI);
    }

protected:

    T1& indO;
    T2& indE;
};


} // namespace 





template <class T>
bool EscFunction::create(   const LinearInterpolation& oFreePath,
                            const LinearInterpolation& eFreePath,
                            const ULong thetaSize,
                            const ULong phiSize,
                            const ULong zSize,
                            const Float maxZ)
{

    typedef Indicatrix<T, Optics::OBeam> IndicatrixO;
    typedef Indicatrix<T, Optics::EBeam> IndicatrixE;


    m_array = allocate3dArray<Float>(zSize, phiSize, thetaSize);

    fprintf(stderr, "calculating escape function\n");

    m_maxZ      = maxZ;
    m_thetaSize = thetaSize;
    m_phiSize   = phiSize;
    m_zSize     = zSize;

    recalcSteps();

    #pragma omp parallel for
    for (int i = 0; i < (int) thetaSize; ++i) {

        Float t_i = i*m_thetaStep;
        Float cost_i = cos(t_i);
        Float sint_i = sin(t_i);

#ifndef EXPERIMENTAL
        for (ULong j = 0; j < phiSize; ++j) {

            Float p_i = j*m_phiStep;
#else
        {
            Float p_i = 0.;
#endif

            
            Vector3 s_i = Vector3(sint_i*cos(p_i), sint_i*sin(p_i), cost_i);
            IndicatrixO indO(s_i, Optics::director);
            IndicatrixE indE(s_i, Optics::director);


            for (ULong k = 0; k < zSize; ++k) {

                Float z = k*m_zStep;

                //double integration
                escfunction::thetaFunctor<IndicatrixO, IndicatrixE>     functor(indO, indE, z, oFreePath, eFreePath);
                escfunction::thetaFunctorNorm<IndicatrixO, IndicatrixE> functorNorm(indO, indE);
                
                Adapt s(1.0e-4);

                Float res = s.integrate(functor, 0.5*M_PI, M_PI);
                Float norm = s.integrate(functorNorm, 0., M_PI);


#ifndef EXPERIMENTAL
                m_array[k][j][i] = res / norm;
#else
                for (ULong j = 0; j < phiSize; ++j)
                    m_array[k][j][i] = res / norm;
#endif
            }
#ifndef EXPERIMENTAL
            fprintf(stderr, "%lu\t%lu\n", (ULong)i, j);
#else
            fprintf(stderr, "%lu\n", (ULong)i);
#endif
        }
    }

    return true;
}




#endif /* end of include guard: _ESCFUNCTION_H_ */
