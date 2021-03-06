#ifndef _INDICATRIX_APP_H_
#define _INDICATRIX_APP_H_

#include <string>
#include <utility>
#include <list>


//right region border and number of iterations for some partition chunk
typedef std::pair<Float, size_t>    ChunkParam;
typedef std::list<ChunkParam>       ChunkParamsList;


class MeasuredData
{
public:

    MeasuredData(const int order_ = 0);
    ~MeasuredData();

    MeasuredData& operator+=(const MeasuredData& other);
    void clear();

    Float** data;
    int     order;  //scattering order. unused for total intensity data

//private:

    MeasuredData(const MeasuredData& other);
    MeasuredData& operator=(const MeasuredData& other);
};


class DataFile
{
public:

    DataFile(const int order_, const std::string& fileName_);

    void output();

    MeasuredData data;
    std::string fileName;
};

typedef std::list<MeasuredData> DataList;
typedef std::list<DataFile>     DataFilesList;


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

	static const int   kSeedIncrement = 1000; //random generator seeds in threads differ by this number



private:

	int  getSeed() const { return m_seed; }

	bool isLoadOFreePath() const {return m_loadOFreePath;}
	bool isSaveOFreePath() const {return m_saveOFreePath;}
	bool isLoadEFreePath() const {return m_loadEFreePath;}
	bool isSaveEFreePath() const {return m_saveEFreePath;}

	bool isLoadEChannelProb() const {return m_loadEChannelProb;}
	bool isSaveEChannelProb() const {return m_saveEChannelProb;}


	bool isLoadEEPartition() const {return m_loadEEPartition;}
	bool isSaveEEPartition() const {return m_saveEEPartition;}
	bool isLoadOEPartition() const {return m_loadOEPartition;}
	bool isSaveOEPartition() const {return m_saveOEPartition;}
	bool isLoadEOPartition() const {return m_loadEOPartition;}
	bool isSaveEOPartition() const {return m_saveEOPartition;}

	bool isLoadEEscFunction() const {return m_loadEEscFunction;}
	bool isSaveEEscFunction() const {return m_saveEEscFunction;}
	bool isLoadOEscFunction() const {return m_loadOEscFunction;}
	bool isSaveOEscFunction() const {return m_saveOEscFunction;}

	bool isLoadONorm() const {return m_loadONorm;}
	bool isSaveONorm() const {return m_saveONorm;}
	bool isLoadENorm() const {return m_loadENorm;}
	bool isSaveENorm() const {return m_saveENorm;}


	const std::string& getOFreePathFileName() const {return m_oFreePathFileName;}
	const std::string& getEFreePathFileName() const {return m_eFreePathFileName;}

	const std::string& getEChannelProbFileName() const {return m_eChannelProbFileName;}

	const std::string& getOEPartitionFileName() const {return m_oePartitionFileName;}
	const std::string& getEOPartitionFileName() const {return m_eoPartitionFileName;}
    const std::string& getEEPartitionFileName() const {return m_eePartitionFileName;}

	const std::string& getOEscFunctionFileName() const {return m_oEscFunctionFileName;}
	const std::string& getEEscFunctionFileName() const {return m_eEscFunctionFileName;}

    const std::string& getONormFileName() const {return m_oNormFileName;}
    const std::string& getENormFileName() const {return m_eNormFileName;}

	const std::string& getWorkDir() const {return m_workDir;}


    int  prepareOFreePath(LinearInterpolation& length);
    int  prepareEFreePath(LinearInterpolation& length);

    int  prepareEChannelProb(LinearInterpolation& prob);

	int  prepareOEPartition(Partition& partition);
    int  prepareEOPartition(Partition& partition);
    int  prepareEEPartition(Partition& partition);

    int  prepareOEscFunction(EscFunction& escFunction);
    int  prepareEEscFunction(EscFunction& escFunction);

    int  prepareONorm(LinearInterpolation& norm);
    int  prepareENorm(LinearInterpolation& norm);

    template <class T>
	void processScattering(const Photon& ph,
	        MeasuredData& oLadder, MeasuredData& oCyclic, DataList& oLadderList, DataList& oCyclicList,
	        MeasuredData& eLadder, MeasuredData& eCyclic, DataList& eLadderList, DataList& eCyclicList);

	void flushBuffers(   const int scatteredCount,
                         const MeasuredData& oLadderData, 
                         const MeasuredData& oCyclicData,
                         const DataList& oLadderDataList,
                         const DataList& oCyclicDataList,
                         const MeasuredData& eLadderData, 
                         const MeasuredData& eCyclicData,
                         const DataList& eLadderDataList,
                         const DataList& eCyclicDataList);

	void output();

    std::string m_workDir;
	std::string m_executableFileName;

	std::string m_oFreePathFileName;
	std::string m_eFreePathFileName;

    std::string m_eChannelProbFileName;

	std::string m_oePartitionFileName;
	std::string m_eoPartitionFileName;
	std::string m_eePartitionFileName;

    std::string m_oEscFunctionFileName;
    std::string m_eEscFunctionFileName;

    std::string m_oNormFileName;
    std::string m_eNormFileName;


   	bool m_loadOFreePath;
	bool m_saveOFreePath;
	bool m_loadEFreePath;
	bool m_saveEFreePath;

 	bool m_loadEChannelProb;
	bool m_saveEChannelProb;

	bool m_loadOEPartition;
	bool m_saveOEPartition;
	bool m_loadEOPartition;
	bool m_saveEOPartition;
	bool m_loadEEPartition;
	bool m_saveEEPartition;

	bool m_loadOEscFunction;
	bool m_saveOEscFunction;
	bool m_loadEEscFunction;
	bool m_saveEEscFunction;

	bool m_loadONorm;
	bool m_saveONorm;
	bool m_loadENorm;
	bool m_saveENorm;


	LinearInterpolation m_eLength;
	LinearInterpolation m_oLength;

	LinearInterpolation m_eChannelProb;

	LinearInterpolation m_oNorm;
	LinearInterpolation m_eNorm;

	int m_seed;
	int m_maxPhotons;
	int m_maxScatterings;
	Float m_minPhotonWeight;

    int m_photonCnt;
    int m_saveRate;

    ChunkParamsList    m_chunkParams;

    DataFilesList m_oLadderFiles;
    DataFilesList m_oCyclicFiles;
    DataFilesList m_eLadderFiles;
    DataFilesList m_eCyclicFiles;

    DataFile m_oLadder;
    DataFile m_oCyclic;
    DataFile m_eLadder;
    DataFile m_eCyclic;

private:

    ScatMCApp& operator=(const ScatMCApp& other);
    ScatMCApp(const ScatMCApp& other);
};

#endif /* _INDICATRIX_APP_H_ */
