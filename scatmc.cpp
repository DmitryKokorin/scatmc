#include <cstdio>
#include <memory.h>
#include <string.h>
#include <sstream>

#include "mathcompat.h"
#include "common.h"

#include "channel.h"
#include "partition.h"
#include "escfunction.h"
#include "idxnorm.h"
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
	m_oFreePathFileName(),
	m_eFreePathFileName(),
	m_eChannelProbFileName(),
	m_oePartitionFileName(),
	m_eoPartitionFileName(),
	m_eePartitionFileName(),
	m_oEscFunctionFileName(),
	m_eEscFunctionFileName(),
	m_oNormFileName(),
	m_eNormFileName(),
	m_loadOFreePath(false),
	m_saveOFreePath(false),
	m_loadEFreePath(false),
	m_saveEFreePath(false),
	m_loadEChannelProb(false),
	m_saveEChannelProb(false),
	m_loadOEPartition(false),
	m_saveOEPartition(false),
	m_loadEOPartition(false),
	m_saveEOPartition(false),
	m_loadEEPartition(false),
	m_saveEEPartition(false),
	m_loadOEscFunction(false),
	m_saveOEscFunction(false),
	m_loadEEscFunction(false),
	m_saveEEscFunction(false),
	m_loadONorm(false),
	m_saveONorm(false),
	m_loadENorm(false),
	m_saveENorm(false),
	m_eLength(),
	m_oLength(),
	m_eChannelProb(),
	m_oNorm(),
	m_eNorm(),
	m_seed(1000),
	m_maxPhotons(1000),
	m_maxScatterings(100000),
	m_minPhotonWeight(1e-8),
    m_photonCnt(0),
	m_chunkParams(),
	m_ladderFiles(),
	m_cyclicFiles(),
	m_dataLadder(),
	m_dataCyclic()
{
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

			if (!m_workDir.empty())
			    m_workDir = m_workDir + '/';
		}
		else if(!strcmp(argv[i], "--loadofreepath")) {

			if (m_saveOFreePath)
				return false;

			if (++i == argc)
				return false;

			m_loadOFreePath     = true;
			m_oFreePathFileName = argv[i];
		}
		else if(!strcmp(argv[i], "--loadefreepath")) {

			if (m_saveEFreePath)
				return false;

			if (++i == argc)
				return false;

			m_loadEFreePath     = true;
			m_eFreePathFileName = argv[i];
		}
		else if(!strcmp(argv[i], "--loadechannelprob")) {

			if (m_saveEChannelProb)
				return false;

			if (++i == argc)
				return false;

			m_loadEChannelProb     = true;
			m_eChannelProbFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--loadoepartition")) {

			if (m_saveOEPartition)
				return false;

			if (++i == argc)
				return false;

			m_loadOEPartition     = true;
			m_oePartitionFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--loadeopartition")) {

			if (m_saveEOPartition)
				return false;

			if (++i == argc)
				return false;

			m_loadEOPartition     = true;
			m_eoPartitionFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--loadeepartition")) {

			if (m_saveEEPartition)
				return false;

			if (++i == argc)
				return false;

			m_loadEEPartition     = true;
			m_eePartitionFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--loadoescfunction")) {

			if (m_saveOEscFunction)
				return false;

			if (++i == argc)
				return false;

			m_loadOEscFunction     = true;
			m_oEscFunctionFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--loadeescfunction")) {

			if (m_saveEEscFunction)
				return false;

			if (++i == argc)
				return false;

			m_loadEEscFunction     = true;
			m_eEscFunctionFileName = argv[i];
		}
		else if(!strcmp(argv[i], "--loadonorm")) {

			if (m_saveONorm)
				return false;

			if (++i == argc)
				return false;

			m_loadONorm     = true;
			m_oNormFileName = argv[i];
		}
		else if(!strcmp(argv[i], "--loadenorm")) {

			if (m_saveENorm)
				return false;

			if (++i == argc)
				return false;

			m_loadENorm     = true;
			m_eNormFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--saveofreepath")) {

			if (m_loadOFreePath)
				return false;

			if (++i == argc)
				return false;

			m_saveOFreePath     = true;
			m_oFreePathFileName = argv[i];

		}
		else if (!strcmp(argv[i], "--saveefreepath")) {

			if (m_loadEFreePath)
				return false;

			if (++i == argc)
				return false;

			m_saveEFreePath     = true;
			m_eFreePathFileName = argv[i];

		}
		else if (!strcmp(argv[i], "--saveechannelprob")) {

			if (m_loadEChannelProb)
				return false;

			if (++i == argc)
				return false;

			m_saveEChannelProb     = true;
			m_eChannelProbFileName = argv[i];

		}
		else if (!strcmp(argv[i], "--saveoepartition")) {

			if (m_loadOEPartition)
				return false;

			if (++i == argc)
				return false;

			m_saveOEPartition     = true;
			m_oePartitionFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--saveeopartition")) {

			if (m_loadEOPartition)
				return false;

			if (++i == argc)
				return false;

			m_saveEOPartition     = true;
			m_eoPartitionFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--saveeepartition")) {

			if (m_loadEEPartition)
				return false;

			if (++i == argc)
				return false;

			m_saveEEPartition     = true;
			m_eePartitionFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--saveoescfunction")) {

			if (m_loadOEscFunction)
				return false;

			if (++i == argc)
				return false;

			m_saveOEscFunction     = true;
			m_oEscFunctionFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--saveeescfunction")) {

			if (m_loadEEscFunction)
				return false;

			if (++i == argc)
				return false;

			m_saveEEscFunction     = true;
			m_eEscFunctionFileName = argv[i];
		}
		else if (!strcmp(argv[i], "--saveonorm")) {

			if (m_loadONorm)
				return false;

			if (++i == argc)
				return false;

			m_saveONorm     = true;
			m_oNormFileName = argv[i];

		}
		else if (!strcmp(argv[i], "--saveenorm")) {

			if (m_loadENorm)
				return false;

			if (++i == argc)
				return false;

			m_saveENorm     = true;
			m_eNormFileName = argv[i];
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


    //allocate arrays of full data
    m_dataCyclic = allocate2dArray<Float>(kPhiSize, kThetaSize);
    m_dataLadder = allocate2dArray<Float>(kPhiSize, kThetaSize);

    //allocate arrays for individual scattering orders data
    int orders[] = {1, 2, 3, 4, 5, 10, 50, 100, 300, 500, 1000, 3000, 5000, 10000, 30000, 50000, 100000};
    int ordersLength = sizeof(orders)/sizeof(orders[0]);

    for (int i = 0; i < ordersLength; ++i) {

        std::string       s1;
        std::stringstream ss1(s1);

        ss1 << m_workDir << "ladder" << orders[i] << ".txt";
        m_ladderFiles.push_back(ScatteringOrder(orders[i], ss1.str(), allocate2dArray<Float>(kPhiSize, kThetaSize))); 

        std::string       s2;
        std::stringstream ss2(s2);


        ss2 << m_workDir << "cyclic" << orders[i] << ".txt";
        m_cyclicFiles.push_back(ScatteringOrder(orders[i], ss2.str(), allocate2dArray<Float>(kPhiSize, kThetaSize))); 
    }


	int res = 0;

    //free path
   	res = prepareOFreePath(m_oLength);

	if (0 != res)
		return res;

	res = prepareEFreePath(m_eLength);

	if (0 != res)
		return res;

	//channel probabilities
	res = prepareEChannelProb(m_eChannelProb);

	if (0 != res)
		return res;


    res = prepareONorm(m_oNorm);

    if (0 != res)
        return res;

    res = prepareENorm(m_eNorm);

    if (0 != res)
        return res;
	

    //partition
    int numChunks = 100;
    Float chunkStep = 0.5*M_PI / numChunks;

    for (int i = 1; i <= numChunks; ++i)
        m_chunkParams.push_back(ChunkParam(i*chunkStep, 20));
    
	Partition pOE, pEO, pEE;

	res = prepareOEPartition(pOE);

	if (0 != res)
		return res;

	res = prepareEOPartition(pEO);

	if (0 != res)
		return res;

	res = prepareEEPartition(pEE);

	if (0 != res)
		return res;



	//escape function
	
    EscFunction eEscFunction;
    res = prepareEEscFunction(eEscFunction);
    if (0 != res)
        return res;

    EscFunction oEscFunction;
    res = prepareOEscFunction(oEscFunction);
    if (0 != res)
        return res;

	fprintf(stderr, "scattering...\n");
	Photon::init(&m_oLength, &m_eLength, &pOE, &pEO, &pEE, &oEscFunction, &eEscFunction, &m_eChannelProb, getSeed());


    //main loop

	#pragma omp parallel for
	for (int i = 0; i < m_maxPhotons; ++i) {

		Photon ph;

		while (ph.pos.z() >= 0.
		        && ph.scatterings < m_maxScatterings
		        && ph.weight > m_minPhotonWeight) {

			ph.move();

            if (Optics::OCHANNEL == ph.channel)
    			processScattering<Optics::OBeam>(ph);
            else
               	processScattering<Optics::EBeam>(ph);


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



template <class T>
void ScatMCApp::processScattering(const Photon& ph)
{
    Indicatrix<T, Optics::OBeam> indO(ph.s_i, Optics::director);
	Indicatrix<T, Optics::EBeam> indE(ph.s_i, Optics::director);
	
	Float norm;
	Float otheta = symmetrizeTheta(Angle(ph.s_i, Optics::director).theta);

	if (Optics::OCHANNEL == ph.channel)
	    norm = m_oNorm(otheta);
    else
        norm = m_eNorm(otheta);


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

			Float symmetrizedTheta = symmetrizeTheta(a_s.theta);

            Float oLengthFactor = exp(-dist/m_oLength(symmetrizedTheta));
			Float eLengthFactor = exp(-dist/m_eLength(symmetrizedTheta));

			Vector3 R  = Vector3(x, y, 0);
			Vector3 qe = Optics::EBeam::k(s_s, a_s)*Optics::k0;
			Vector3 qo = Optics::OBeam::k(s_s, a_s)*Optics::k0;

			Float oProbFactor = indO(s_s)/norm;
			Float eProbFactor = indE(s_s)/norm;

			Float ladderRes = ph.weight*(oLengthFactor*oProbFactor + eLengthFactor*eProbFactor);

			Float cyclicRes = ph.weight*(oLengthFactor*oProbFactor*cos(qo*R)
			                           + eLengthFactor*eProbFactor*cos(qe*R));

            ScatteringOrderFiles::iterator iter = m_ladderFiles.begin();
            ScatteringOrderFiles::iterator end  = m_ladderFiles.end();

            while (iter != end) {

                if (ph.scatterings + 1 == iter->order)
                    iter->data[j][i] += ladderRes;
             
                ++iter;
            }

            iter = m_cyclicFiles.begin();
            end  = m_cyclicFiles.end();

            while (iter != end) {

                if (ph.scatterings + 1 == iter->order)
                    iter->data[j][i] += cyclicRes;

                ++iter;
            }

            m_dataLadder[j][i] += ladderRes;

            if (0 != ph.scatterings)
                m_dataCyclic[j][i] += cyclicRes;
		}
	}
}


void ScatMCApp::output()
{
    typedef ScatteringOrderFiles::iterator Iter;


    for (Iter iter = m_ladderFiles.begin(); iter != m_ladderFiles.end(); ++iter) {

        iter->file = fopen(iter->fileName.c_str(), "w");
    }

    for (Iter iter = m_cyclicFiles.begin(); iter != m_cyclicFiles.end(); ++iter) {

        iter->file = fopen(iter->fileName.c_str(), "w");
    }


    std::string fileName = m_workDir + "ladder.txt";
	FILE* ladderFile     = fopen(fileName.c_str(), "w");

    fileName = m_workDir + "cyclic.txt";
	FILE* cyclicFile     = fopen(fileName.c_str(), "w");



	Float thetaStep = kThetaMax / kThetaSize;
	Float phiStep   = 2*M_PI / kPhiSize;

	for (int i = 0; i < kThetaSize; ++i)
		for (int j = 0; j < kPhiSize; ++j) {

			Float phi   = j*phiStep;
			Float theta = i*thetaStep;

            for (Iter iter = m_ladderFiles.begin(); iter != m_ladderFiles.end(); ++iter) {

                fprintf(iter->file,     "%e\t%e\t%.17e\n", theta, phi, iter->data[j][i]);
            }

            for (Iter iter = m_cyclicFiles.begin(); iter != m_cyclicFiles.end(); ++iter) {

                fprintf(iter->file,     "%e\t%e\t%.17e\n", theta, phi, iter->data[j][i]);
            }

			fprintf(ladderFile,   "%e\t%e\t%.17e\n", theta, phi, m_dataLadder[j][i]);
			fprintf(cyclicFile,   "%e\t%e\t%.17e\n", theta, phi, m_dataCyclic[j][i]);
		}


    for (Iter iter = m_ladderFiles.begin(); iter != m_ladderFiles.end(); ++iter) {

        fclose(iter->file);
    }

    for (Iter iter = m_cyclicFiles.begin(); iter != m_cyclicFiles.end(); ++iter) {

        fclose(iter->file);
    }


    fclose(ladderFile);
    fclose(cyclicFile);
}


void ScatMCApp::printHelp()
{
	fprintf(stderr, "Usage: %s [options]", m_executableFileName.c_str());
	fprintf(stderr, "\n\nAvailable options:");
	fprintf(stderr, "\n--seed [seed]\t\t\t\tseed for random numbers generator");
	fprintf(stderr, "\n--workdir [path]\t\t\toutput path");
	fprintf(stderr, "\n--loadofreepath [filename]\t\tload o-beam free path from file");
	fprintf(stderr, "\n--loadefreepath [filename]\t\tload e-beam free path from file");
	fprintf(stderr, "\n--loadechannelprob [filename]\t\tload e-e probability from file");
	fprintf(stderr, "\n--loadoepartition [filename]\t\tload o-e partition from file");
	fprintf(stderr, "\n--loadeopartition [filename]\t\tload e-o partition from file");
	fprintf(stderr, "\n--loadeepartition [filename]\t\tload e-e partition from file");
	fprintf(stderr, "\n--loadoescfunction [filename]\t\tload o-escape function from file");
	fprintf(stderr, "\n--loadeescfunction [filename]\t\tload e-escape function from file");
	fprintf(stderr, "\n--loadonorm [filename]\t\tload o-beam indicatrix norm from file");
	fprintf(stderr, "\n--loadenorm [filename]\t\tload e-beam indicatrix from file");
	fprintf(stderr, "\n--saveofreepath [filename]\t\tsave o-beam free path to file");
	fprintf(stderr, "\n--saveefreepath [filename]\t\tsave e-beam free path to file");
	fprintf(stderr, "\n--saveechannelprob [filename]\t\tsave e-e probability to file");
	fprintf(stderr, "\n--saveoepartition [filename]\t\tsave o-e partition to file");
	fprintf(stderr, "\n--saveeopartition [filename]\t\tsave e-o partition to file");
    fprintf(stderr, "\n--saveeepartition [filename]\t\tsave e-e partition to file");
	fprintf(stderr, "\n--saveoescfunction [filename]\t\tsave o-escape function to file");
	fprintf(stderr, "\n--saveeescfunction [filename]\t\tsave e-escape function to file");
	fprintf(stderr, "\n--saveonorm [filename]\t\tsave o-beam indicatrix norm to file");
	fprintf(stderr, "\n--saveenorm [filename]\t\tsave e-beam indicatrix norm to file");
	fprintf(stderr, "\n--photons [photons]\t\t\tnumber of photons to scatter");
	fprintf(stderr, "\n--scatterings [scatterings]\t\tmax scatterings for each photon");
	fprintf(stderr, "\n");
}

int ScatMCApp::prepareOFreePath(LinearInterpolation& l)
{
	if (isLoadOFreePath()) {

		fprintf(stderr, "loading o-beam free path file...");

		if (!l.load(getOFreePathFileName())) {

			fprintf(stderr, "can't load o-beam free path data\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}
	else {

		fprintf(stderr, "calculating o-beam free path data...\n");
		createFreePath<Optics::OBeam>(l);
	}

	if (isSaveOFreePath()) {

		fprintf(stderr, "saving o-beam free path data to file...");

		if (!l.save(getOFreePathFileName())) {

			fprintf(stderr, "can't save o-beam free path data\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}

	return 0;
}


int ScatMCApp::prepareEFreePath(LinearInterpolation& l)
{
	if (isLoadEFreePath()) {

		fprintf(stderr, "loading e-beam free path file...");

		if (!l.load(getEFreePathFileName())) {

			fprintf(stderr, "can't load e-beam free path data\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}
	else {

		fprintf(stderr, "calculating e-beam free path data...\n");
		createFreePath<Optics::EBeam>(l);
	}

	if (isSaveEFreePath()) {

		fprintf(stderr, "saving e-beam free path data to file...");

		if (!l.save(getEFreePathFileName())) {

			fprintf(stderr, "can't save e-beam free path data\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}

	return 0;
}



int ScatMCApp::prepareEChannelProb(LinearInterpolation& l)
{
	if (isLoadEChannelProb()) {

		fprintf(stderr, "loading e-e probability file...");

		if (!l.load(getEChannelProbFileName())) {

			fprintf(stderr, "can't load e-e probability data\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}
	else {

		fprintf(stderr, "calculating e-e probability data...\n");
		createEChannelProb<Optics::EBeam>(l);
	}

	if (isSaveEChannelProb()) {

		fprintf(stderr, "saving e-e probability data to file...");

		if (!l.save(getEChannelProbFileName())) {

			fprintf(stderr, "can't save e-e probability data\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}

	return 0;
}


int ScatMCApp::prepareEEPartition(Partition& p)
{
	if (isLoadEEPartition()) {

		fprintf(stderr, "loading e-e partition...");

		if (!p.load(getEEPartitionFileName())) {

			fprintf(stderr, "can't load e-e partition\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}
	else {

		fprintf(stderr, "creating e-e partition...\n");

        Float leftBorder = 0.;
        ChunkParamsList::iterator i;
        for (i = m_chunkParams.begin(); i != m_chunkParams.end(); ++i) {

            p.addChunk<IndicatrixEE>(leftBorder, i->first, i->second);
            leftBorder = i->first;
        }

	}

	if (isSaveEEPartition()) {

		fprintf(stderr, "saving e-e partition...");

		if (!p.save(getEEPartitionFileName())) {

			fprintf(stderr, "can't save e-e partition\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}

	return 0;
}

int ScatMCApp::prepareOEPartition(Partition& p)
{
	if (isLoadOEPartition()) {

		fprintf(stderr, "loading o-e partition...");

		if (!p.load(getOEPartitionFileName())) {

			fprintf(stderr, "can't load o-e partition\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}
	else {

		fprintf(stderr, "creating o-e partition...\n");

        Float leftBorder = 0.;
        ChunkParamsList::iterator i;
        for (i = m_chunkParams.begin(); i != m_chunkParams.end(); ++i) {

            p.addChunk<IndicatrixOE>(leftBorder, i->first, i->second);
            leftBorder = i->first;
        }

	}

	if (isSaveOEPartition()) {

		fprintf(stderr, "saving o-e partition...");

		if (!p.save(getOEPartitionFileName())) {

			fprintf(stderr, "can't save o-e partition\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}

	return 0;
}


int ScatMCApp::prepareEOPartition(Partition& p)
{
	if (isLoadEOPartition()) {

		fprintf(stderr, "loading e-o partition...");

		if (!p.load(getEOPartitionFileName())) {

			fprintf(stderr, "can't load e-o partition\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}
	else {

		fprintf(stderr, "creating e-o partition...\n");

        Float leftBorder = 0.;
        ChunkParamsList::iterator i;
        for (i = m_chunkParams.begin(); i != m_chunkParams.end(); ++i) {

            p.addChunk<IndicatrixEO>(leftBorder, i->first, i->second);
            leftBorder = i->first;
        }

	}

	if (isSaveEOPartition()) {

		fprintf(stderr, "saving e-o partition...");

		if (!p.save(getEOPartitionFileName())) {

			fprintf(stderr, "can't save e-o partition\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}

	return 0;
}



int ScatMCApp::prepareEEscFunction(EscFunction& esc)
{
	if (isLoadEEscFunction()) {

		fprintf(stderr, "loading e-escape function...");

		if (!esc.load(getEEscFunctionFileName())) {

			fprintf(stderr, "can't load e-escape function\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}
	else {

		fprintf(stderr, "creating e-escape function...\n");
#ifdef EXPERIMENTAL
        esc.create<Optics::EBeam>(m_oLength, m_eLength, 180, 1, 301, 15.*Optics::l);
#else
        esc.create<Optics::EBeam>(m_oLength, m_eLength, 90, 90, 200, 1.);
#endif
	}

	if (isSaveEEscFunction()) {

		fprintf(stderr, "saving e-escape function...");

		if (!esc.save(getEEscFunctionFileName())) {

			fprintf(stderr, "can't save e-escape function\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}

	return 0;
}

int ScatMCApp::prepareOEscFunction(EscFunction& esc)
{
	if (isLoadOEscFunction()) {

		fprintf(stderr, "loading o-escape function...");

		if (!esc.load(getOEscFunctionFileName())) {

			fprintf(stderr, "can't load o-escape function\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}
	else {

		fprintf(stderr, "creating o-escape function...\n");
#ifdef EXPERIMENTAL
        esc.create<Optics::OBeam>(m_oLength, m_eLength, 180, 1, 301, 15.*Optics::l);
#else
        esc.create<Optics::OBeam>(m_oLength, m_eLength, 90, 90, 200, 1.);
#endif
	}

	if (isSaveOEscFunction()) {

		fprintf(stderr, "saving o-escape function...");

		if (!esc.save(getOEscFunctionFileName())) {

			fprintf(stderr, "can't save o-escape function\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}

	return 0;
}

int ScatMCApp::prepareONorm(LinearInterpolation& l)
{
	if (isLoadONorm()) {

		fprintf(stderr, "loading o-beam norm file...");

		if (!l.load(getONormFileName())) {

			fprintf(stderr, "can't load o-beam norm data\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}
	else {

		fprintf(stderr, "calculating o-beam norm data...\n");
		createIndicatrixNorm<Optics::OBeam>(l);
	}

	if (isSaveONorm()) {

		fprintf(stderr, "saving o-beam norm data to file...");

		if (!l.save(getONormFileName())) {

			fprintf(stderr, "can't save o-beam norm data\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}

	return 0;
}


int ScatMCApp::prepareENorm(LinearInterpolation& l)
{
	if (isLoadENorm()) {

		fprintf(stderr, "loading e-beam norm file...");

		if (!l.load(getENormFileName())) {

			fprintf(stderr, "can't load e-beam norm data\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}
	else {

		fprintf(stderr, "calculating e-beam norm data...\n");
		createIndicatrixNorm<Optics::EBeam>(l);
	}

	if (isSaveENorm()) {

		fprintf(stderr, "saving e-beam norm data to file...");

		if (!l.save(getENormFileName())) {

			fprintf(stderr, "can't save e-beam norm data\n");
			return -1;
		}
		else {

			fprintf(stderr, "\tdone\n");
		}
	}

	return 0;
}

