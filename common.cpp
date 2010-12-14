#include "common.h"
#include "cmath"
#include <cstdio>

#if defined DOUBLE_PRECISION
const Float kMachineEpsilon = 1e-17;
#else
const Float kMachineEpsilon = 1e-6;
#endif



int solveQuadric(const Float a, const Float b, const Float c, Float& x1, Float& x2)
{
//	fprintf(stderr, "a=%f\tb=%f\tc=%f\n", a, b, c);

	Float max = MAX(fabs(a), MAX(fabs(b), fabs(c)));

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

//	fprintf(stderr, "x1=%f\tx2=%f\n", x1, x2);

	return 2;
}
