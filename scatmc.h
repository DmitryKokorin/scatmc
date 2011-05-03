#ifndef _INDICATRIX_APP_H_
#define _INDICATRIX_APP_H_

#include <string>

#include "photon.h"

class ScatMCApp
{
public:

    ScatMCApp();

	int  run();
	bool getOpts(int argc, char ** argv);
	void printHelp();


private:

	int  getSeed() { return m_seed; }

	bool isLoadFreePath() const {return m_loadFreePath;}
	bool isSaveFreePath() const {return m_saveFreePath;}
	bool isLoadPartition() const {return m_loadPartition;}
	bool isSavePartition() const {return m_savePartition;}

	const std::string& getFreePathFileName() const {return m_freePathFileName;}
	const std::string& getPartitionFileName() const {return m_partitionFileName;}

	int  prepareFreePath(FreePath& length);
	int  preparePartition(Partition& partition);


	void processScattering(const Photon& ph);
	void output();
	bool checkResultsReady();

	FreePath m_length;

	std::string m_executableFileName;
	std::string m_freePathFileName;
	std::string m_partitionFileName;

	bool m_loadFreePath;
	bool m_saveFreePath;
	bool m_loadPartition;
	bool m_savePartition;

	int m_seed;
	int m_maxPhotons;
	int m_maxScatterings;

	Float m_minPhotonWeight;

	static const int   kPhiSize    = 32;
	static const int   kThetaSize  = 50;

	static const Float kThetaMax;

	//TODO: list of arrays for detected intensity

	Float det1[kPhiSize][kThetaSize];
	Float det2[kPhiSize][kThetaSize];
	Float det5[kPhiSize][kThetaSize];
	Float det100[kPhiSize][kThetaSize];
	Float det5000[kPhiSize][kThetaSize];
    Float det100000[kPhiSize][kThetaSize];
	Float detall[kPhiSize][kThetaSize];

	Float lastdet[kPhiSize][kThetaSize];
};

#endif /* _INDICATRIX_APP_H_ */
