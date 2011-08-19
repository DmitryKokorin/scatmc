#include "integrate.h"


Adapt::Adapt(const Float tol) :
    TOL(tol),
    toler(0.),
    terminate(true),
    out_of_tolerance(false)
{
	TOL = (TOL < 10.*kMachineEpsilon) ? 10.*kMachineEpsilon : TOL;
}


const Float Adapt::alpha = sqrt(2.0/3.0);
const Float Adapt::beta  = 1.0/sqrt(5.0);
const Float Adapt::x1    = 0.942882415695480;
const Float Adapt::x2    = 0.641853342345781;
const Float Adapt::x3    = 0.236383199662150;
const Float Adapt::x[12] = {0, -x1, -alpha, -x2, -beta, -x3, 0.0, x3, beta, x2, alpha, x1};

