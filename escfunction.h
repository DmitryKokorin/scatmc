#ifndef _ESCFUNCTION_H_
#define _ESCFUNCTION_H_

#include <string>

class FreePath;

class EscFunction
{
public:

    EscFunction();
    ~EscFunction();

    bool create(const FreePath& lengths,
                const int thetaSize,
                const int phiSize,
                const int zSize,
                const Float maxZ,
                const int thetaIterations = 360,
                const int phiIterations = 360);

    Float operator()(const Float theta, const Float phi, const Float z) const;

    bool load(const std::string& name);
    bool save(const std::string& name);


private:

    void recalcSteps();

    Float m_maxZ;

    size_t m_thetaSize;
    size_t m_phiSize;
    size_t m_zSize;

    Float m_thetaStep;
    Float m_phiStep;
    Float m_zStep;

    Float ***m_array;

    EscFunction(const EscFunction&);
    EscFunction& operator=(const EscFunction&);
};


#endif /* end of include guard: _ESCFUNCTION_H_ */
