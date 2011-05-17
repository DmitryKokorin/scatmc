#include <cstdio>
#include <cmath>

#include "optics.h"
#include "indicatrix.h"
#include "freepath.h"
#include "escfunction.h"

// [z][phi][theta]

EscFunction::EscFunction() :
    m_maxZ(0.),
    m_thetaStep(0.),
    m_phiStep(0.),
    m_zStep(0.),
    m_array(NULL)
{
}

EscFunction::~EscFunction()
{
    free3dArray(m_array);
}

bool EscFunction::create(const FreePath& length,
                         const int thetaSize,
                         const int phiSize,
                         const int zSize,
                         const Float maxZ,
                         const int thetaIterations /*= 1000*/,
                         const int phiIterations /*= 1000*/)
{
    m_array = allocate3dArray<Float>(zSize, phiSize, thetaSize);

    m_maxZ = maxZ;
    m_thetaStep = M_PI / thetaSize;
    m_phiStep = 2*M_PI / phiSize;
    m_zStep = maxZ / zSize;

    const Float kThetaIterStep = 0.5*M_PI / thetaIterations;
    const Float kPhiIterStep = M_PI / phiIterations;

    for (int i = 0; i < thetaSize; ++i) {

        Float t_i = i*m_thetaStep;
        Float cost_i = cos(t_i);
        Float sint_i = sin(t_i);

        for (int j = 0; j < phiSize; ++j) {

            Float p_i = j*m_phiStep;

            Vector3 s_i = Vector3(sint_i*cos(p_i), sint_i*sin(p_i), cost_i);
            Indicatrix ind = Indicatrix(s_i, Optics::n);

            for (int k = 0; k < zSize; ++k) {

                Float z = k*m_zStep;

                //double integration
                Float res = 0;

                for (int l = 1; l < thetaIterations; ++l) { //skip t_s == 0

                    Float t_s = 0.5*M_PI + l*kThetaIterStep;
                    Float cost_s = cos(t_s);
                    Float sint_s = sin(t_s);

                    for (int m = 0; m < phiIterations; ++m) {

                        Float p_s = m*kPhiIterStep;
                        Vector3 s_s = Vector3(sint_s*cos(p_s), sint_s*sin(p_s), cost_s);

                        Float dist = z / s_s.z();
                        res += ind(s_s)*sint_s*exp(dist/length(Angle(s_s, Optics::n)));
                    }
                }
            }
        }
    }

    return true;
}

Float EscFunction::operator()(const Float theta, const Float phi, const Float z) const
{
    if (z > m_maxZ)
        return 0;

    return m_array[(int)(z/m_zStep)][(int)(phi/m_phiStep)][(int)(theta/m_thetaStep)];
}

bool EscFunction::load(const std::string& name)
{
    FILE *file = fopen(name.c_str(), "r");
    if (!file)
        return false;

    fclose(file);

    return true;
}

bool EscFunction::save(const std::string& name)
{
    FILE *file = fopen(name.c_str(), "w");
    if (!file)
        return false;

    fclose(file);

    return true;
}




