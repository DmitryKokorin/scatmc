#ifndef _LINTERPOL_H_

#define _LINTERPOL_H_

#include <vector>
#include <string>

#include "common.h"

class LinearInterpolation
{
public:

    LinearInterpolation(const Float min = 0., const Float max = 0., const size_t size = 0);
    
    bool load(const std::string& name);
    bool save(const std::string& name);

    Float operator()(const Float& x) const;
    Float& operator[](const size_t idx) {return m_data[idx];}

    Float min() const {return m_min;}
    Float max() const {return m_max;}
    Float resolution() const {return m_resolution;}

    void setRange(const Float min_, const Float max_) { m_min = min_; m_max = max_; m_resolution = (m_max - m_min) / (m_data.size() - 1);}
    void resize(const size_t size) { m_data.resize(size); }

protected:

    Float m_min;
    Float m_max;
    Float m_resolution;

    std::vector<Float> m_data;
};



#endif /* end of include guard: _LINTERPOL_H_ */
