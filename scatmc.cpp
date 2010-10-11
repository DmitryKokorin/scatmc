#include <cstdio>
#include <cmath>
#include <memory.h>

#include "common.h"

#include "direction.h"
#include "extlength.h"
#include "photon.h"
#include "indicatrix.h"
#include "partition.h"

#include "scatmc.h"



int main(int /*argc*/, char ** /*argv*/)
{
    ScatMCApp app;
    app.run();
    
	return 0;
}

const Float ScatMCApp::thetaMax = 1e-5;

ScatMCApp::ScatMCApp() 
{
	memset(&det1,     0, sizeof(det1));
	memset(&det2,     0, sizeof(det2));
	memset(&det5,     0, sizeof(det5));
	memset(&det100,   0, sizeof(det100));
	memset(&det10000, 0, sizeof(det10000));
}




void ScatMCApp::run()
{
	fprintf(stderr, "calculating extinction lengths...\n");
	ExtLength length;

	fprintf(stderr, "preparing partition...\n");
	Partition p;

	fprintf(stderr, "scattering...\n");
	Photon::init(&length, &p); 

	int cnt = 0;

	#pragma omp parallel for
	for (int i = 0; i < maxPhotons; ++i) {
		
		Photon ph;
		
		while (ph.pos.z() >= 0 && ph.scatterings < maxScatterings) {

			ph.move();
			ph.scatter();

			processScattering(ph);
		}

		#pragma omp critical
		{
			++cnt;
			fprintf(stderr, "%d\t%d\n", cnt, ph.scatterings);
		}
	}

	output();
}


void ScatMCApp::processScattering(const Photon& ph)
{
    #pragma omp critical
	{
	
	for (int i = 0; i < thetaSize; ++i)
		for (int j = 0; j < phiSize; ++j) {

			Float res = 0;

			if (1 == ph.scatterings)
				det1[j][i] += res;

			if (2 <= ph.scatterings)
				det2[j][i] += res;

			if (5 <= ph.scatterings)
				det5[j][i] += res;

			if (100 <= ph.scatterings)
				det100[j][i] += res;

			det10000[j][i] += res;
		}
	}
}


void ScatMCApp::output()
{
	FILE *det1file     = 0,
		 *det2file     = 0,
		 *det5file     = 0,
		 *det100file   = 0,
		 *det10000file = 0;

	det1file     = fopen("output/peak1.txt", "w");
	det2file     = fopen("output/peak2.txt", "w");
	det5file     = fopen("output/peak5.txt", "w");
	det100file   = fopen("output/peak100.txt", "w");
	det10000file = fopen("output/peak1000.txt", "w");

	for (int i = 0; i < thetaSize; ++i)
		for (int j = 0; j < phiSize; ++j) {

			Float phi   = 2.*M_PI  / phiSize;
			Float theta = thetaMax / thetaSize;


			fprintf(det1file,     "%f\t%f\t%f\n", theta, phi, det1[j][i]);
			fprintf(det2file,     "%f\t%f\t%f\n", theta, phi, det2[j][i]);
			fprintf(det5file,     "%f\t%f\t%f\n", theta, phi, det5[j][i]);
			fprintf(det100file,   "%f\t%f\t%f\n", theta, phi, det100[j][i]);
			fprintf(det10000file, "%f\t%f\t%f\n", theta, phi, det10000[j][i]);
		}

	fclose(det1file);
	fclose(det2file);
	fclose(det5file);
	fclose(det100file);
	fclose(det10000file);
}