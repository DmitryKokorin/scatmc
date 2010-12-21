#include <cstdio>
#include <cmath>
#include <memory.h>

#include "common.h"

#include "extlength.h"
#include "photon.h"
#include "indicatrix.h"
#include "partition.h"
#include "optics.h"

#include "scatmc.h"



int main(int argc, char ** argv)
{
    ScatMCApp app;

    int res = 0;

	if (app.getOpts(argc, argv)) {

    	res = app.run();
	}
	else {

		app.printHelp();
		res = -1;
	}
    
	return res;
}

const Float ScatMCApp::kThetaMax = 1e-4;


ScatMCApp::ScatMCApp() :
	m_length(),
	m_executableFileName(),
	m_extLengtsFileName(),
	m_partitionFileName(),
	m_loadExtLengths(false),
	m_saveExtLengths(false),
	m_loadPartition(false),
	m_savePartition(false),
	m_seed(1000),
	m_maxPhotons(100),
	m_maxScatterings(10)
{
	memset(&det1,     0, sizeof(det1));
	memset(&det2,     0, sizeof(det2));
	memset(&det5,     0, sizeof(det5));
	memset(&det100,   0, sizeof(det100));
	memset(&detall,   0, sizeof(detall));
	memset(&lastdet,  0, sizeof(lastdet));
}

bool ScatMCApp::getOpts(int argc, char ** argv)
{
	m_executableFileName = argv[0];

	for (int i = 1; i < argc; ++i) {

		if (!strcmp(argv[i], "--seed")) {

			if (++i == argc)
				return false;

			m_seed = atoi(argv[i]);
		}
		else if(!strcmp(argv[i], "--loadlengths")) {

			if (m_saveExtLengths)
				return false;

			if (++i == argc)
				return false;

			m_loadExtLengths    = true;
			m_extLengtsFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--loadpartition")) {

			if (m_savePartition)
				return false;

			if (++i == argc)
				return false;

			m_loadPartition     = true;
			m_partitionFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--savelengths")) {

			if (m_loadExtLengths)
				return false;

			if (++i == argc)
				return false;

			m_saveExtLengths    = true;
			m_extLengtsFileName = argv[i];

		}
		else if (!strcmp(argv[i], "--savepartition")) {

			if (m_loadPartition)
				return false;

			if (++i == argc)
				return false;

			m_savePartition     = true;
			m_partitionFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--photons")) {

			if (++i == argc)
				return false;

			m_maxPhotons = atoi(argv[i]);
		}
		else if (!strcmp(argv[i], "--scatterings")) {

			if (++i == argc)
				return false;

			m_maxScatterings = atoi(argv[i]);
		}
		else 
			return false;
	}

	return true;
}


int ScatMCApp::run()
{
	int res = 0;

	res = prepareExtinctionLengths(m_length);

	if (0 != res)
		return res;

	Partition p;
	res = preparePartition(p);

	if (0 != res)
		return res;

	fprintf(stderr, "scattering...\n");
	Photon::init(&m_length, &p, getSeed()); 

	int cnt = 0;
	bool ready = false;

	#pragma omp parallel for
	for (int i = 0; i < m_maxPhotons; ++i) 
		if (!ready) {

			Photon ph;
		
			while (ph.pos.z() >= 0 && ph.scatterings < m_maxScatterings) {

				ph.move();
				ph.scatter();

				//fprintf(stderr, "x=%f\ty=%f\tz=%f\tsx=%f\tsy=%f\tsz=%f\n", ph.pos.x(), ph.pos.y(), ph.pos.z(), ph.s_i.x(), ph.s_i.y(), ph.s_i.z());
				//fprintf(stderr, "%.17f\t%.17f\t%.17f\n", ph.pos.x(), ph.pos.y(), ph.pos.z());


				processScattering(ph);
			}

			#pragma omp critical
			{
				++cnt;
				fprintf(stderr, "%d\t%d\n", cnt, ph.scatterings);
				if (0 == cnt % 100)
					ready = checkResultsReady();
			}
		}

	output();

	return 0;
}


void ScatMCApp::processScattering(const Photon& ph)
{
	Indicatrix ind(ph.s_i, Optics::n);

	Float thetaStep = kThetaMax / kThetaSize;
	Float phiStep   = 2*M_PI / kThetaSize;

    #pragma omp critical
	{
	
	for (int i = 0; i < kThetaSize; ++i)
		for (int j = 0; j < kPhiSize; ++j) {

			Float theta_s = i*thetaStep;
			Float phi_s   = j*phiStep;

			Float sintheta_s = sin(theta_s);

			Vector3 s_s = Vector3(sintheta_s*cos(phi_s),  //wrong parity, but who cares?
								  sintheta_s*sin(phi_s),
								  -cos(theta_s));

			Angle a_s = Angle(s_s, Optics::n);

			Float scale = ph.pos.z() / s_s.z();
			Float x     = ph.pos.x() + scale*s_s.x();
			Float y     = ph.pos.y() + scale*s_s.y();

			Float length = fabs(scale);
			Float lengthFactor = exp(-length/m_length(a_s));
			
			Vector3 R = Vector3(x, y, 0);
			Vector3 q = Optics::ke(s_s, Optics::n) - Optics::ke(ph.s_i, Optics::n);

			//TODO: think more about it
			Float res = lengthFactor*cos(q*R)/**ind(s_s)/ph.fullIntegral*sintheta_s*/;

			if (1 == ph.scatterings)
				det1[j][i] += res;

			if (2 <= ph.scatterings)
				det2[j][i] += res;

			if (5 <= ph.scatterings)
				det5[j][i] += res;

			if (100 <= ph.scatterings)
				det100[j][i] += res;

			detall[j][i] += res;
		}
	}
}


void ScatMCApp::output()
{
	FILE *det1file     = 0,
		 *det2file     = 0,
		 *det5file     = 0,
		 *det100file   = 0,
		 *detallfile   = 0;

	det1file     = fopen("output/peak1.txt", "w");
	det2file     = fopen("output/peak2.txt", "w");
	det5file     = fopen("output/peak5.txt", "w");
	det100file   = fopen("output/peak100.txt", "w");
	detallfile   = fopen("output/peakall.txt", "w");

	

	Float thetaStep = kThetaMax / kThetaSize;
	Float phiStep   = 2*M_PI / kPhiSize;

	for (int i = 0; i < kThetaSize; ++i)
		for (int j = 0; j < kPhiSize; ++j) {

			Float phi   = j*phiStep;
			Float theta = i*thetaStep;

			fprintf(det1file,     "%e\t%e\t%.17e\n", theta, phi, det1[j][i]);
			fprintf(det2file,     "%e\t%e\t%.17e\n", theta, phi, det2[j][i]);
			fprintf(det5file,     "%e\t%e\t%.17e\n", theta, phi, det5[j][i]);
			fprintf(det100file,   "%e\t%e\t%.17e\n", theta, phi, det100[j][i]);
			fprintf(detallfile,   "%e\t%e\t%.17e\n", theta, phi, detall[j][i]);
			
		}

	fclose(det1file);
	fclose(det2file);
	fclose(det5file);
	fclose(det100file);
	fclose(detallfile);
}

bool ScatMCApp::checkResultsReady()
{
	int size = kPhiSize*kThetaSize;
	for (int i = 0; i < size; ++i) {

		if (fabs((*lastdet)[i] < kMachineEpsilon))
			return false;

		Float err =  fabs(((*detall)[i] - (*lastdet)[i]) / (*lastdet)[i]);
		if (err > 0.01) {

			fprintf(stderr, "relative error: %f\n", err);
			memcpy(lastdet, detall, sizeof(lastdet));
			
			return false;
		}
	}
		
	return true;
}

void ScatMCApp::printHelp()
{
	fprintf(stderr, "Usage: %s [options]", m_executableFileName.c_str());
	fprintf(stderr, "\n\nAvailable options:");
	fprintf(stderr, "\n--seed [seed]\t\t\t\tseed for random numbers generator");
	fprintf(stderr, "\n--loadlengths [filename]\t\tload extinction lengths from file");
	fprintf(stderr, "\n--loadpartition [filename]\t\tload partition from file");
	fprintf(stderr, "\n--savelengths [filename]\t\tsave extinction lengths to file");
	fprintf(stderr, "\n--savepartition [filename]\t\tsave partition to file");
	fprintf(stderr, "\n--photons [photons]\t\t\tnumber of photons to scatter");
	fprintf(stderr, "\n--scatterings [scatterings]\t\tmax scatterings for each photon");
	fprintf(stderr, "\n");
}

int ScatMCApp::prepareExtinctionLengths(ExtLength& l)
{
	if (isLoadExtLengths()) {

		fprintf(stderr, "loading extinction lengths...");

		if (!l.load(getExtLengthsFileName())) {

			fprintf(stderr, "can't load extinction lengths\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}
	else {
		
		fprintf(stderr, "calculating extinction lengths...\n");
		l.create();
	}

	if (isSaveExtLengths()) {

		fprintf(stderr, "saving extinction lengths...");

		if (!l.save(getExtLengthsFileName())) {

			fprintf(stderr, "can't save extinction lengths\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}

	return 0;
}

int ScatMCApp::preparePartition(Partition& p)
{
	if (isLoadPartition()) {

		fprintf(stderr, "loading partition...");

		if (!p.load(getPartitionFileName())) {

			fprintf(stderr, "can't load partition\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}
	else {
		
		fprintf(stderr, "creating partition...\n");
		p.create();
	}

	if (isSavePartition()) {

		fprintf(stderr, "saving partition...");

		if (!p.save(getPartitionFileName())) {

			fprintf(stderr, "can't save partition\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}

	return 0;
}
