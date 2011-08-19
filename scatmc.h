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

	static const Float kThetaStep;
	static const Float kPhiStep;

	static const int   kSaveRate = 2; //save results each on each kSaveRate photons simulated



private:

	int  getSeed() const { return m_seed; }

	bool isLoadOFreePath() const {return m_loadOFreePath;}
	bool isSaveOFreePath() const {return m_saveOFreePath;}
	bool isLoadEFreePath() const {return m_loadEFreePath;}
	bool isSaveEFreePath() const {return m_saveEFreePath;}

	bool isLoadOChannelProb() const {return m_loadOChannelProb;}
	bool isSaveOChannelProb() const {return m_saveOChannelProb;}
	bool isLoadEChannelProb() const {return m_loadEChannelProb;}
	bool isSaveEChannelProb() const {return m_saveEChannelProb;}


	bool isLoadPartition() const {return m_loadPartition;}
	bool isSavePartition() const {return m_savePartition;}
	bool isLoadEscFunction() const {return m_loadEscFunction;}
	bool isSaveEscFunction() const {return m_saveEscFunction;}


	const std::string& getOFreePathFileName() const {return m_oFreePathFileName;}
	const std::string& getEFreePathFileName() const {return m_eFreePathFileName;}

	const std::string& getOChannelProbFileName() const {return m_oChannelProbFileName;}
	const std::string& getEChannelProbFileName() const {return m_eChannelProbFileName;}

	const std::string& getPartitionFileName() const {return m_partitionFileName;}
	const std::string& getEscFunctionFileName() const {return m_escFunctionFileName;}
	const std::string& getWorkDir() const {return m_workDir;}


    int  prepareOFreePath(LinearInterpolation& length);
    int  prepareEFreePath(LinearInterpolation& length);

    int  prepareOChannelProb(LinearInterpolation& prob);
    int  prepareEChannelProb(LinearInterpolation& prob);

	int  preparePartition(Partition& partition);
    int  prepareEscFunction(EscFunctionEE& escFunction);

	void processScattering(const Photon& ph);
	void output();

    std::string m_workDir;
	std::string m_executableFileName;

	std::string m_oFreePathFileName;
	std::string m_eFreePathFileName;

    std::string m_oChannelProbFileName;
    std::string m_eChannelProbFileName;

	std::string m_partitionFileName;
    std::string m_escFunctionFileName;

   	bool m_loadOFreePath;
	bool m_saveOFreePath;
	bool m_loadEFreePath;
	bool m_saveEFreePath;

   	bool m_loadOChannelProb;
	bool m_saveOChannelProb;
	bool m_loadEChannelProb;
	bool m_saveEChannelProb;


	bool m_loadPartition;
	bool m_savePartition;

	bool m_loadEscFunction;
	bool m_saveEscFunction;

	LinearInterpolation m_eLength;
	LinearInterpolation m_oLength;

	LinearInterpolation m_oChannelProb;
	LinearInterpolation m_eChannelProb;

	int m_seed;
	int m_maxPhotons;
	int m_maxScatterings;
	Float m_minPhotonWeight;

    int m_photonCnt;

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
