#include "optics.h"


namespace Optics {

//physical constants

const Vector3 n = Vector3(1., 0., 0.);	

const Float eps_par = 3.0;
const Float eps_perp = 2.2;
const Float eps_a = eps_par - eps_perp;

const Float K3 = 6.1e-7;
const Float t1 = 0.79;
const Float t2 = 0.43;
const Float K1 = t1*K3;
const Float K2 = t2*K3;

const Float lambda = 4.88e-5; //cm
const Float k0 = 2.*M_PI / lambda;
const Float xi = 4.2e-4;

const Float T = 301.;
const Float Kb = 1.38e-16;

const Float H = 0.5e+4;
const Float Hi_alpha = 5e-6;

//precalculated constants

const Float s0 = (0.25*Kb*T*eps_a*eps_a)/(lambda*lambda*K3);
const Float add = 0.25*lambda*lambda/(M_PI*M_PI*xi*xi);



#ifdef EXPERIMENTAL

const Float g = 0.80188954242914923;

const Float l 		= 2.5e-2;				//free path in cm
const Float ls     = l / (1-g);

#endif


}  //namespace Optics
