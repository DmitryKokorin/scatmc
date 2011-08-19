#ifndef _MATHCOMPAT_H_
#define _MATHCOMPAT_H_

#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#ifdef WIN32 


inline double expm1(double x)
{
    return std::exp(x) - 1.0;
}

inline double log1p(double x)
{
    return std::log(1.0 + x);
}

inline bool isnan(double x)
{
	return x != x;
}
#endif

#endif
