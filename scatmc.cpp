#include <cstdio>
#include <memory.h>
#include <string.h>

#include "mathcompat.h"
#include "common.h"

#include "freepath.h"
#include "partition.h"
#include "escfunction.h"
#include "photon.h"
#include "indicatrix.h"
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
const Float ScatMCApp::kThetaStep = ScatMCApp::kThetaMax / ScatMCApp::kThetaSize;
const Float ScatMCApp::kPhiStep   = 2*M_PI / ScatMCApp::kPhiSize;



ScatMCApp::ScatMCApp() :
    m_workDir(),
	m_executableFileName(),
	m_freePathFileName(),
	m_partitionFileName(),
	m_escFunctionFileName(),
	m_loadFreePath(false),
	m_saveFreePath(false),
	m_loadPartition(false),
	m_savePartition(false),
	m_loadEscFunction(false),
	m_saveEscFunction(false),
	m_length(),
	m_seed(1000),
	m_maxPhotons(1000),
	m_maxScatterings(100000),
	m_minPhotonWeight(1e-8),
    m_photonCnt(0),
	m_chunkParams()
{
	memset(&det1,     0, sizeof(det1));
	memset(&det2,     0, sizeof(det2));
	memset(&det5,     0, sizeof(det5));
	memset(&det100,   0, sizeof(det100));
	memset(&det5000,   0, sizeof(det100));
    memset(&det100000,   0, sizeof(det100));

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
		else if(!strcmp(argv[i], "--workdir")) {

			if (++i == argc)
				return false;

			m_workDir = argv[i];
		}
		else if(!strcmp(argv[i], "--loadfreepath")) {

			if (m_saveFreePath)
				return false;

			if (++i == argc)
				return false;

			m_loadFreePath     = true;
			m_freePathFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--loadpartition")) {

			if (m_savePartition)
				return false;

			if (++i == argc)
				return false;

			m_loadPartition     = true;
			m_partitionFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--loadescfunction")) {

			if (m_saveEscFunction)
				return false;

			if (++i == argc)
				return false;

			m_loadEscFunction     = true;
			m_escFunctionFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--savefreepath")) {

			if (m_loadFreePath)
				return false;

			if (++i == argc)
				return false;

			m_saveFreePath     = true;
			m_freePathFileName = argv[i];

		}
		else if (!strcmp(argv[i], "--savepartition")) {

			if (m_loadPartition)
				return false;

			if (++i == argc)
				return false;

			m_savePartition     = true;
			m_partitionFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--saveescfunction")) {

			if (m_loadEscFunction)
				return false;

			if (++i == argc)
				return false;

			m_saveEscFunction     = true;
			m_escFunctionFileName = argv[i];
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
    fprintf(stderr, "# seed = %d\n", getSeed());

	int res = 0;

    //free path
	res = prepareFreePath(m_length);

	if (0 != res)
		return res;

    //partition
#ifdef EXPERIMENTAL
    int numChunks = 1;
#else
    int numChunks = 100;
#endif
    Float chunkStep = 0.5*M_PI / numChunks;

    for (int i = 1; i <= numChunks; ++i)
        m_chunkParams.push_back(ChunkParam(i*chunkStep, 100));
    
	Partition p;
	res = preparePartition(p);

	if (0 != res)
		return res;

	//escape function
	
    EscFunction escFunction;
    res = prepareEscFunction(escFunction);

    if (0 != res)
        return res;



	fprintf(stderr, "scattering...\n");
	Photon::init(&m_length, &p, &escFunction, getSeed());


    //main loop

	#pragma omp parallel for
	for (int i = 0; i < m_maxPhotons; ++i) {

		Photon ph;

		while (ph.pos.z() >= 0.
		        && ph.scatterings < m_maxScatterings
		        && ph.weight > m_minPhotonWeight) {

			ph.move();

			processScattering(ph);

			ph.scatter();

		}


		#pragma omp critical
		{
			++m_photonCnt;
			fprintf(stderr, "Photon: %d\tScatterings: %d\n", m_photonCnt, ph.scatterings);

			if (0 == m_photonCnt % kSaveRate)
                output();
		}
	}

	output();

	return 0;
}


void ScatMCApp::processScattering(const Photon& ph)
{
	if (0 == ph.scatterings)
		return;

	IndicatrixEE ind(ph.s_i, Optics::director);

    #pragma omp critical
	{

	for (int i = 0; i < kThetaSize; ++i)
		for (int j = 0; j < kPhiSize; ++j) {

			Float theta_s = i*kThetaStep;
			Float phi_s   = j*kPhiStep;

			Float sintheta_s = sin(theta_s);
			Float costheta_s = cos(theta_s);

			Vector3 s_s = Vector3(sintheta_s*cos(phi_s),  //wrong parity, but who cares?
								  sintheta_s*sin(phi_s),
								  -costheta_s);

			Angle a_s = Angle(s_s, Optics::director);

			Float dist   = fabs(ph.pos.z() / s_s.z());
			Float x      = ph.pos.x() + dist*s_s.x();
			Float y      = ph.pos.y() + dist*s_s.y();

			Float lengthFactor = exp(-dist/m_length(a_s));

			Vector3 R = Vector3(x, y, 0);
#ifndef EXPERIMENTAL
			Vector3 q = Optics::EBeam::k(s_s, a_s)*Optics::k0;
#else
            Vector3 q = Optics::k0*s_s;
#endif

			Float probFactor = ind(s_s)/ph.fullIntegral;

			Float res = ph.weight*lengthFactor*probFactor*cos(q*R);

			if (!isnan(res)) {
				if (0 == ph.scatterings)
					det1[j][i] += res;
				else {

					if (1 == ph.scatterings)
						det2[j][i] += res;

					if (4 == ph.scatterings)
						det5[j][i] += res;

					if (100 == ph.scatterings)
						det100[j][i] += res;

					if (5000 == ph.scatterings)
					    det5000[j][i] += res;

					if (100000 == ph.scatterings)
					    det100000[j][i] += res;

					detall[j][i] += res;
				}
			}
			else {

				fprintf(stderr, "%f\n", ph.fullIntegral);
			}
		}
	}
}


void ScatMCApp::output()
{
	FILE *det1file     = 0,
		 *det2file     = 0,
		 *det5file     = 0,
		 *det100file   = 0,
		 *det5000file  = 0,
		 *det100000file= 0,
		 *detallfile   = 0;

	// I'll fix this one day, I promise

    std::string fileName;

    fileName = m_workDir + "/peak1.txt";
	det1file     = fopen(fileName.c_str(), "w");

    fileName = m_workDir + "/peak2.txt";
   	det2file     = fopen(fileName.c_str(), "w");

    fileName = m_workDir + "/peak5.txt";
   	det5file     = fopen(fileName.c_str(), "w");

    fileName = m_workDir + "/peak100.txt";
   	det100file     = fopen(fileName.c_str(), "w");

    fileName = m_workDir + "/peak5000.txt";
   	det5000file     = fopen(fileName.c_str(), "w");

    fileName = m_workDir + "/peak100000.txt";
   	det100000file     = fopen(fileName.c_str(), "w");


    fileName = m_workDir + "/peakall.txt";
	detallfile   = fopen(fileName.c_str(), "w");


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
			fprintf(det5000file,  "%e\t%e\t%.17e\n", theta, phi, det5000[j][i]);
			fprintf(det100000file,  "%e\t%e\t%.17e\n", theta, phi, det100000[j][i]);
			fprintf(detallfile,   "%e\t%e\t%.17e\n", theta, phi, detall[j][i]);
		}

	fclose(det1file);
	fclose(det2file);
	fclose(det5file);
	fclose(det100file);
	fclose(det5000file);
	fclose(det100000file);
	fclose(detallfile);
}


void ScatMCApp::printHelp()
{
	fprintf(stderr, "Usage: %s [options]", m_executableFileName.c_str());
	fprintf(stderr, "\n\nAvailable options:");
	fprintf(stderr, "\n--seed [seed]\t\t\t\tseed for random numbers generator");
	fprintf(stderr, "\n--workdir [path]\t\t\toutput path");
	fprintf(stderr, "\n--loadfreepath [filename]\t\tload extinction lengths from file");
	fprintf(stderr, "\n--loadpartition [filename]\t\tload partition from file");
	fprintf(stderr, "\n--loadescfunction [filename]\t\tload escape function from file");
	fprintf(stderr, "\n--savefreepath [filename]\t\tsave extinction lengths to file");
	fprintf(stderr, "\n--savepartition [filename]\t\tsave partition to file");
	fprintf(stderr, "\n--saveescfunction [filename]\t\tsave escape function to file");
	fprintf(stderr, "\n--photons [photons]\t\t\tnumber of photons to scatter");
	fprintf(stderr, "\n--scatterings [scatterings]\t\tmax scatterings for each photon");
	fprintf(stderr, "\n");
}

int ScatMCApp::prepareFreePath(FreePathEE& l)
{
	if (isLoadFreePath()) {

		fprintf(stderr, "loading free path file...");

		if (!l.load(getFreePathFileName())) {

			fprintf(stderr, "can't load free path data\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}
	else {

		fprintf(stderr, "calculating free path data...\n");
		l.create();
	}

	if (isSaveFreePath()) {

		fprintf(stderr, "saving free path data to file...");

		if (!l.save(getFreePathFileName())) {

			fprintf(stderr, "can't save free path data\n");
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

        Float leftBorder = 0.;
        ChunkParamsList::iterator i;
        for (i = m_chunkParams.begin(); i != m_chunkParams.end(); ++i) {

            p.addChunk(leftBorder, i->first, i->second);
            leftBorder = i->first;
        }

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


int ScatMCApp::prepareEscFunction(EscFunction& esc)
{
	if (isLoadEscFunction()) {

		fprintf(stderr, "loading escape function...");

		if (!esc.load(getEscFunctionFileName())) {

			fprintf(stderr, "can't load escape function\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}
	else {

		fprintf(stderr, "creating escape function...\n");
#ifdef EXPERIMENTAL
        esc.create(m_length, 180, 1, 301, 15.*Optics::l, 3600, 100);
#else
        esc.create(m_length, 90, 90, 200, 1., 800, 800);
#endif
	}

	if (isSaveEscFunction()) {

		fprintf(stderr, "saving escape function...");

		if (!esc.save(getEscFunctionFileName())) {

			fprintf(stderr, "can't save escape function\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}

	return 0;
}
