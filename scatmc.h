#ifndef _INDICATRIX_APP_H_
#define _INDICATRIX_APP_H_

#include <string>
#include <utility>
#include <list>


//right region border and number of iterations for some partition chunk
typedef std::pair<Float, size_t> ChunkParam;
typedef std::list<ChunkParam> ChunkParamsList;



class ScatMCApp
{
public:

    ScatMCApp();

	int  run();
	bool getOpts(int argc, char ** argv);
	void printHelp();

	static const int   kPhiSize    = 32;
	static const int   kThetaSize  = 50;
	static const Float kThetaMax;  //max interested value, peak is supposed to be in [0, kThetaMAx]


private:

	int  getSeed() { return m_seed; }

	bool isLoadFreePath() const {return m_loadFreePath;}
	bool isSaveFreePath() const {return m_saveFreePath;}
	bool isLoadPartition() const {return m_loadPartition;}
	bool isSavePartition() const {return m_savePartition;}
	bool isLoadEscFunction() const {return m_loadEscFunction;}
	bool isSaveEscFunction() const {return m_saveEscFunction;}


	const std::string& getFreePathFileName() const {return m_freePathFileName;}
	const std::string& getPartitionFileName() const {return m_partitionFileName;}
	const std::string& getEscFunctionFileName() const {return m_escFunctionFileName;}



	int  prepareFreePath(FreePath& length);
	int  preparePartition(Partition& partition);
    int  prepareEscFunction(EscFunction& escFunction);

	void processScattering(const Photon& ph);
	void output();
	bool checkResultsReady();

	std::string m_executableFileName;
	std::string m_freePathFileName;
	std::string m_partitionFileName;
    std::string m_escFunctionFileName;

	bool m_loadFreePath;
	bool m_saveFreePath;
	bool m_loadPartition;
	bool m_savePartition;
	bool m_loadEscFunction;
	bool m_saveEscFunction;

	FreePath m_length;

	int m_seed;
	int m_maxPhotons;
	int m_maxScatterings;
	Float m_minPhotonWeight;

    ChunkParamsList    m_chunkParams;
	

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
