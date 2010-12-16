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

	bool isLoadExtLengths() const {return m_loadExtLengths;}
	bool isSaveExtLengths() const {return m_saveExtLengths;}
	bool isLoadPartition() const {return m_loadPartition;}
	bool isSavePartition() const {return m_savePartition;}

	const std::string& getExtLengthsFileName() const {return m_extLengtsFileName;}
	const std::string& getPartitionFileName() const {return m_partitionFileName;}

	int  prepareExtinctionLengths(ExtLength& length);
	int  preparePartition(Partition& partition);


	void processScattering(const Photon& ph);
	void output();
	bool checkResultsReady();

	ExtLength m_length;

	std::string m_executableFileName;
	std::string m_extLengtsFileName;
	std::string m_partitionFileName;

	bool m_loadExtLengths;
	bool m_saveExtLengths;
	bool m_loadPartition;
	bool m_savePartition;

	int m_seed;
	int m_maxPhotons;
	int m_maxScatterings;

	static const int   kPhiSize    = 32;
	static const int   kThetaSize  = 50;

	static const Float kThetaMax;

	//TODO: list of arrays for detected intensity

	Float det1[kPhiSize][kThetaSize];
	Float det2[kPhiSize][kThetaSize];
	Float det5[kPhiSize][kThetaSize];
	Float det100[kPhiSize][kThetaSize];
	Float detall[kPhiSize][kThetaSize];

	Float lastdet[kPhiSize][kThetaSize];
};

#endif /* _INDICATRIX_APP_H_ */
