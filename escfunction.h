#ifndef _ESCFUNCTION_H_
#define _ESCFUNCTION_H_

#include <string>
#include <cstdio>
#include "mathcompat.h"

#include "optics.h"
#include "indicatrix.h"
#include "freepath.h"
#include "integrate.h"


template <class T>
class EscFunction
{
public:

    EscFunction();
    ~EscFunction();

    bool create(const LinearInterpolation& lengths,
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


template <class T>
class phiFunctor
{
public:

    phiFunctor(T& indicatrix_, const Float sint_s_, const Float cost_s_, const Float z_, const LinearInterpolation& length_) :
        indicatrix(indicatrix_),
        sint_s(sint_s_),
        cost_s(cost_s_),
        z(z_),
        length(length_)
    {}

    inline Float operator()(const Float phi)
    {  
//        if (cost_s < -kMachineEpsilon) {

            Vector3 s_s = Vector3(sint_s*cos(phi), sint_s*sin(phi), cost_s);

            Float dist = z / cost_s;

            //fprintf(stderr, "dist = %.17e\n", dist);
            //fprintf(stderr, "dist = %.17e\n", dist/length(symmetrizeTheta(Angle(s_s, Optics::director).theta)));


            return indicatrix(s_s)*sint_s*exp(dist/length(symmetrizeTheta(Angle(s_s, Optics::director).theta)));
//        }

        return 0.;
    }

protected:

    T& indicatrix;
    Float sint_s;
    Float cost_s;
    Float z;

    const LinearInterpolation& length;
};

template <class T>
class phiFunctorNorm
{
public:

    phiFunctorNorm(T& indicatrix_, const Float sint_s_, const Float cost_s_) :
        indicatrix(indicatrix_),
        sint_s(sint_s_),
        cost_s(cost_s_)
    {}


    inline Float operator()(const Float phi)
    {
        Vector3 s_s = Vector3(sint_s*cos(phi), sint_s*sin(phi), cost_s);

        return indicatrix(s_s)*sint_s;
    }

protected:

    T& indicatrix;
    Float sint_s;
    Float cost_s;
};


template <class T>
class thetaFunctor
{
public:

    thetaFunctor(T& indicatrix_, const Float z_, const LinearInterpolation& length_) :
        indicatrix(indicatrix_),
        z(z_),
        length(length_)
    {}

    inline Float operator()(const Float theta)
    {  
        Float cost_s = cos(theta);
        Float sint_s = sin(theta);
     
        if (cost_s < -kMachineEpsilon) {

            phiFunctor<T> functor = phiFunctor<T>(indicatrix, sint_s, cost_s, z, length);
            //return  qsimp(functor, 0., 2*M_PI, 0.00001);
            Adapt s(1.0e-6);
            return s.integrate(functor, 0., 2*M_PI);
        }

        return 0.;
    }

protected:

    T& indicatrix;
    Float z;

    const LinearInterpolation& length;
};


template <class T>
class thetaFunctorNorm
{
public:

    thetaFunctorNorm(T& indicatrix_) :
        indicatrix(indicatrix_)
    {}

    inline Float operator()(const Float theta)
    {  
        Float cost_s = cos(theta);
        Float sint_s = sin(theta);
     
        phiFunctorNorm<T>  functorNorm = phiFunctorNorm<T>(indicatrix, sint_s, cost_s);
        //return  qsimp(functorNorm, 0., 2*M_PI, 0.00001);
        Adapt s(1.0e-6);
        return s.integrate(functorNorm, 0., 2*M_PI);
    }

protected:

    T& indicatrix;
};






template <class T>
EscFunction<T>::EscFunction() :
    m_maxZ(0.),
    m_thetaSize(0),
    m_phiSize(0),
    m_zSize(0),
    m_thetaStep(0.),
    m_phiStep(0.),
    m_zStep(0.),
    m_array(NULL)
{
}

template <class T>
EscFunction<T>::~EscFunction()
{
//    free3dArray(m_array);
}


template <class T>
bool EscFunction<T>::create(const LinearInterpolation& length,
                            const ULong thetaSize,
                            const ULong phiSize,
                            const ULong zSize,
                            const Float maxZ)
{
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
            T ind = T(s_i, Optics::director);

            for (ULong k = 0; k < zSize; ++k) {

                Float z = k*m_zStep;

                //double integration
                thetaFunctor<T> functor = thetaFunctor<T>(ind, z, length);
                //Float res = qsimp(functor, 0.5*M_PI + 10*kMachineEpsilon, M_PI, 0.00001);

                thetaFunctorNorm<T> functorNorm = thetaFunctorNorm<T>(ind);
                //Float norm = qsimp(functorNorm, 0., M_PI, 0.00001);
                

                Adapt s(1.0e-4);

                Float res = s.integrate(functor, 0.5*M_PI, M_PI);
                Float norm = s.integrate(functorNorm, 0., M_PI);


                //fprintf(stderr, "%.17e %.17e\n", res, norm);

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

template <class T>
Float EscFunction<T>::operator()(const Float theta, const Float phi, const Float z) const
{
    if (z >= m_maxZ)
        return 0;

    Float phi_ = phi >= 0 ? phi : 2.*M_PI + phi;
    phi_ = phi_ < M_PI ? phi_ : phi_ - M_PI;	// n and -n are equivalent

    int zIdx = (int)(z/m_zStep);
    int phiIdx = (int)(phi_/m_phiStep);
    int thetaIdx = (int)(theta/m_thetaStep);

    if (zIdx == (int)(m_zSize - 1))
        return m_array[zIdx][phiIdx][thetaIdx];

    Float mu = z - zIdx*m_zStep;

    return (1. - mu)*m_array[zIdx][phiIdx][thetaIdx] +
                 mu* m_array[zIdx+1][phiIdx][thetaIdx];
}

template <class T>
bool EscFunction<T>::load(const std::string& name)
{
    FILE *file = fopen(name.c_str(), "r");
    if (!file)
        return false;

    int res = 0;

    res = fscanf(file, "%lu %lu %lu %lf", &m_thetaSize, &m_phiSize, &m_zSize, &m_maxZ);
    recalcSteps();

    m_array = allocate3dArray<Float>(m_zSize, m_phiSize, m_thetaSize);

    for (size_t i = 0; i < m_thetaSize; ++i)
        for (size_t j = 0; j < m_phiSize; ++j)
            for (size_t k = 0; k < m_zSize; ++k)
                res = fscanf(file, "%lf", &(m_array[k][j][i]));

    fclose(file);

    return true;
}

template <class T>
bool EscFunction<T>::save(const std::string& name)
{
    FILE *file = fopen(name.c_str(), "w");
    if (!file)
        return false;

    fprintf(file, "%lu %lu %lu %.17e\n", m_thetaSize,
             m_phiSize, m_zSize, m_maxZ);


    for (size_t i = 0; i < m_thetaSize; ++i)
        for (size_t j = 0; j < m_phiSize; ++j)
            for (size_t k = 0; k < m_zSize; ++k)
                fprintf(file, "%.17e\n", m_array[k][j][i]);


    fclose(file);

    return true;
}

template <class T>
void EscFunction<T>::recalcSteps()
{
    m_thetaStep = M_PI / m_thetaSize;
    m_phiStep = M_PI / m_phiSize;
    m_zStep = m_maxZ / m_zSize;
}

typedef EscFunction<Indicatrix<Optics::EBeam, Optics::EBeam> > EscFunctionEE;


#endif /* end of include guard: _ESCFUNCTION_H_ */
