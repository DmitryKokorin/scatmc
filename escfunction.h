#ifndef _ESCFUNCTION_H_
#define _ESCFUNCTION_H_

#include <string>

//class FreePathEE;

class EscFunction
{
public:

    EscFunction();
    ~EscFunction();

    bool create(const LinearInterpolation& lengths,
                const ULong thetaSize,
                const ULong phiSize,
                const ULong zSize,
                const Float maxZ,
                const ULong thetaIterations = 1000,
                const ULong phiIterations = 1000);

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


#endif /* end of include guard: _ESCFUNCTION_H_ */
