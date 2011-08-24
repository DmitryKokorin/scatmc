#include "escfunction.h"

EscFunction::EscFunction() :
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

EscFunction::~EscFunction()
{
//    free3dArray(m_array);
}

Float EscFunction::operator()(const Float theta, const Float phi, const Float z) const
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

bool EscFunction::load(const std::string& name)
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

bool EscFunction::save(const std::string& name)
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

void EscFunction::recalcSteps()
{
    m_thetaStep = M_PI / m_thetaSize;
    m_phiStep   = M_PI / m_phiSize;
    m_zStep     = m_maxZ / m_zSize;
}
