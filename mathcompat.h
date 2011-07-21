#ifndef _MATHCOMPAT_H_
#define _MATHCOMPAT_H_

#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#ifdef WIN32 
#define M_LN2 0.69314718055994530941723212146 // ln(2)

#include <limits>
// Implementation taken from GSL (http://www.gnu.org/software/gsl/)
inline double expm1(double x)
{
	/* FIXME: this should be improved */
	if (std::fabs(x) < M_LN2) {
		/* Compute the taylor series S = x + (1/2!) x^2 + (1/3!) x^3 + ... */
		double i = 1.0;
		double sum = x;
		double term = x / 1.0;
		do {
			i++;
			term *= x / i;
			sum += term;
		}
		while (std::fabs(term) > (std::fabs(sum)
			* std::numeric_limits<double>::epsilon()));
		return sum;
	}
	else {
		return std::exp(x) - 1.0;
	}
}

inline double log1p(double x) { return std::log(1.0 + x); }

inline bool isnan(double x) {
	return x != x;
}
#endif

#endif
