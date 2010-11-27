#include "common.h"

#if defined DOUBLE_PRECISION
const Float kMachineEpsilon = 1e-17;
#else
const Float kMachineEpsilon = 1e-6;
#endif

