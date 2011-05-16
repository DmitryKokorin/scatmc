#include "common.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

#if defined DOUBLE_PRECISION
const Float kMachineEpsilon = 1e-17;
#else
const Float kMachineEpsilon = 1e-6;
#endif



int solveQuadric(const Float a, const Float b, const Float c, Float& x1, Float& x2)
{
	Float max = std::max(fabs(a), std::max(fabs(b), fabs(c)));

	if (max == 0)
		return 0;

	Float A = a/max;
	Float B = b/max;
	Float C = c/max;

	if (fabs(A) < kMachineEpsilon) {

		if (fabs(B) < kMachineEpsilon) {
				return -1;
		}
		else {

			x1 = -C / B;
			return 1;
		}
	}


	Float discr = B*B - 4.*A*C;

	if (discr < 0.)
		return 0;

	if (discr < kMachineEpsilon) {

		x1 = -0.5*B/A;
		return 1;
	}

	Float discrSqrt = sqrt(discr);

	x1 = 0.5*(-B - discrSqrt)/A;
	x2 = 0.5*(-B + discrSqrt)/A;

	return 2;
}
