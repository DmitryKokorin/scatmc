#include "freepath.h"

Float symmetrizeTheta(const Float theta)
{
    Float res = theta;
    res = fabs(res);
    res = (res >= M_PI) ? res - M_PI : res;
	res = (res < M_PI * 0.5) ? res : (M_PI - res);

	return res;
}

