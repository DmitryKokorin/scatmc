#ifndef _INDICATRIX_APP_H_
#define _INDICATRIX_APP_H_

#include "photon.h"

class ScatMCApp
{
public:
    ScatMCApp();

	void run();
	void processScattering(const Photon& ph);
	void output();

private:

	static const int   maxPhotons     = 100;
	static const int   maxScatterings = 1000;

	static const int   phiSize    = 16;
	static const int   thetaSize  = 50;

	static const Float thetaMax;

	Float det1[phiSize][thetaSize];
	Float det2[phiSize][thetaSize];
	Float det5[phiSize][thetaSize];
	Float det100[phiSize][thetaSize];
	Float det10000[phiSize][thetaSize];
};

#endif /* _INDICATRIX_APP_H_ */
