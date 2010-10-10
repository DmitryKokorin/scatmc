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


void preparePartitionTree(Partition& p);
void createRegionsList(Partition& p);



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
	preparePartitionTree(p); //TODO: move this code to Partition c'tor

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

void preparePartitionTree(Partition& p)
{
	Float ** data = allocate2dArray<Float>(Partition::size, Partition::size);

	const Float phiStep   = 2*M_PI/Partition::size;
	const Float thetaStep =   M_PI/Partition::size;

	const Float idxIterations = 100;
	const Float idxThetaStep  = M_PI / idxIterations;

	for (int k = 0; k < idxIterations; ++k) {

		Direction d_i = Direction(k*idxThetaStep, 0.);
		Indicatrix ind(d_i);

		Vector3 oX = Vector3(1, 0, 0);
		Vector3 iVector = d_i.toVector3();
				
		Float   angle = acos(oX*iVector);

		Float m[3][3];  //rotation matrix

		if (angle) {

			Vector3 axis  = crossProduct(iVector, oX).normalize(); 

			{
				Float rcos = cos(angle);
				Float rsin = sin(angle);

				Float u = axis.x(); 
				Float v = axis.y(); 
				Float w = axis.z(); 

				m[0][0] =      rcos + u*u*(1-rcos);
				m[1][0] =  w * rsin + v*u*(1-rcos);
				m[2][0] = -v * rsin + w*u*(1-rcos);
				m[0][1] = -w * rsin + u*v*(1-rcos);
				m[1][1] =      rcos + v*v*(1-rcos);
				m[2][1] =  u * rsin + w*v*(1-rcos);
				m[0][2] =  v * rsin + u*w*(1-rcos);
				m[1][2] = -u * rsin + v*w*(1-rcos);
				m[2][2] =      rcos + w*w*(1-rcos);
			}

		}
		else {

			memset(&m, 0, sizeof(m));
			m[0][0] = m[1][1] = m[2][2] = 1.;
		}

		for (int i = 0; i < Partition::size; ++i)
			for (int j = 0; j < Partition::size; ++j) {

				/*
				Direction d_s = Direction(i*thetaStep, j*phiStep);
				data[j][i] = ind(d_s)*d_s.sintheta;
				*/

				Float t = i*thetaStep;
				Float p = j*phiStep;

				Vector3 rel = Direction(t,p).toVector3();
				Vector3 abs(m[0][0]*rel.x() + m[0][1]*rel.y() + m[0][2]*rel.z(),
							m[1][0]*rel.x() + m[1][1]*rel.y() + m[1][2]*rel.z(),
							m[2][0]*rel.x() + m[2][1]*rel.y() + m[2][2]*rel.z());

				Direction d_s = Direction(abs);
				data[j][i] = ind(d_s)*d_s.sintheta;
			}

		p.setData(data, (M_PI/Partition::size) * (2*M_PI/Partition::size));
		p.refine();

		//printf("%d nodes...\n", p.rectCount);
	}

	free2dArray(data);

}

void createRegionsList(Partition& /*p*/)
{
}
