#ifndef _INDICATRIX_APP_H_
#define _INDICATRIX_APP_H_

#include "photon.h"

class ScatMCApp
{
public:
    ScatMCApp();

	void run();
	bool getOpts(int argc, char ** argv);
	int  getSeed() {return seed;}

private:

	void processScattering(const Photon& ph);
	void output();
	bool checkResultsReady();

	int seed;

	static const int   maxPhotons     = 101;
	static const int   maxScatterings = 1000;

	static const int   phiSize    = 16;
	static const int   thetaSize  = 50;

	static const Float thetaMax;

	//TODO: list of arrays for detected intensity

	Float det1[phiSize][thetaSize];
	Float det2[phiSize][thetaSize];
	Float det5[phiSize][thetaSize];
	Float det100[phiSize][thetaSize];
	Float detall[phiSize][thetaSize];

	Float lastdet[phiSize][thetaSize];
};

#endif /* _INDICATRIX_APP_H_ */
