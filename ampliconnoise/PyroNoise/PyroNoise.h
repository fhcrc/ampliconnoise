/*This software and documentation is copyright © 2009 by Christopher Quince.*/

/*Permission is granted for anyone to copy, use, or modify these programs and documents for purposes of research or education, provided this copyright notice is retained, and note is made of any changes that have been made.*/ 

/* These programs and documents are distributed without any warranty, express or implied. As the programs were written for research purposes only, they have not been tested to the degree that would be advisable in any important application. All use of these programs is entirely at the user's own risk.*/

#ifndef CLASSIFY_H
#define CLASSIFY_H

/**Constants**/
#define MAX_SEQUENCE_LENGTH  1024
#define ID_LENGTH       1024
#define BIG_LINE_LENGTH 1048576
#define MAX_LINE_LENGTH 65536
#define DELIM           " \n"

#define FALSE 0
#define TRUE  1

/* Input Definitions */
#define OPTION  0      /* optional */
#define ALWAYS  1      /* required */

#define UNSEEN_PENALTY 10.0
#define BINS  1000
#define MAX_S 10
#define MAX_7 7
#define MAX_DIST 100.0
#define MAX_DBL  1.0e10
#define DELTA 1.0e-7
#define PRECISION 0.01
#define MIN_DELTA 1.0e-6
#define MIN_ITER 10
#define MAX_ITER 1000
#define MIN_COUNT 0.1
#define MIN_TAU   1.0e-4

#define OVERCALL_THRESH  6.0 /*?*/

#define LOOKUP_FILE "../Data/LookUp_E123.dat"

#define OUT_FILE_STUB            "-out"
#define INIT_FILE                "-lin"
#define DATA_FILE                "-din"
#define NO_INDEX                 "-ni"
#define SIGMA                    "-s"
#define INIT_CUT                 "-c"
#define VERBOSE                  "-v"
#define LOOKUP_FILE_FLAG         "-rin"

#define Z_SUFFIX                 ".z"
#define ALIGN_SUFFIX             ".align"
#define LIST_SUFFIX              ".list"
#define OTU_SUFFIX               ".otu"
#define FASTA_SUFFIX             ".fa"
#define ALIGNED_FASTA_SUFFIX     ".afa"
#define QUAL_SUFFIX              ".qual"
#define TAU_SUFFIX               ".tau"
#define CEN_SUFFIX               ".cen"
#define MASTER_SUFFIX            ".master"
#define STATE_SUFFIX             ".state"
#define MAP_SUFFIX               ".mapping"

#define DEF_CLUSTERS 10
#define DEF_SIGMA    60.0
#define INIT_SIZE    1024
#define DEF_CUT      0.01

#define DIAG  0
#define LEFT  1
#define UP    2

#define NOT_SET -1
#define GAP     -2
#define UPPER_S  2

#define MIN_WEIGHT 0.1  
/* i sequence k cluster j pos*/
#define aindex(i,k,j)   ((i)*nM*nK + nM*(k) + (j)) 

#ifndef min
#define min(x, y)	((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define	max(x, y)	((x) > (y) ? (x) : (y))
#endif

/* User defined structures */

#define INIT_MASTER_SIZE 1024

typedef struct s_Letter
{
  int nSize;

  int nM;

  int nTotal;

  int    *anK;

  int    *anI;

  double *adT;

  double *adD;

} t_Letter;

typedef struct s_Master
{
  double   *adTau;  /*probs*/
  double   *adDist; /*distances*/
  
  int      **aanP;   /*pointers*/
  int      **aanI;   /*members*/

  int      nTotal;

  int      *anN;  
  int      *anCN;
  int      *anP;
  int      *anI;
  int      nSize;
  
  int      nK;
  int      nM;
  int      nN;
} t_Master;

typedef struct s_Node
{
  int left; 

  int right; 

  double distance;
} 
t_Node;

typedef struct s_Flows
{
  /*number of flows*/
  int nN;
  /*max flow length*/
  int nM;
  /*flows*/
  double* adData;

  short * asData;

  char  **aszID;

  /*clean lengths*/
  int*    anLengths;
} t_Flows;

typedef struct s_Unique
{
  int   nN;

  int   nM;

  int   nU;

  int   nSize;

  int   nMSize;
  
  int   *anF;

  int   *anMap;

  short *asU;

  int   *anLenU;

  int   *anWeights;
} t_Unique;


typedef struct s_Params
{
  char *szDataFile;

  char *szOutFileStub;

  char *szInitFile;

  char *szLookUpFile;

  int bNoIndex;

  int nK;

  double dSigma;

  double dInitCut;
} t_Params;

/*User defined functions*/

void getCommandLineParams(t_Params *ptParams, int argc, char *argv[]);

void readData(char* szDataFile, t_Flows *ptFlows, t_Params *ptParams);

int imin(int nA, int nB);

int getMove(double *pdMin, double dA, double dB, double dC);

double dist(double dSignal, double dSignal2);

double find_closest_pair(int n, double** distmatrix, int* ip, int* jp);

void outputCluster(t_Params *ptParams, t_Node* tree, char **aszID, int nN);

void initLookUp(t_Params *ptParams);

void writeCentroids(t_Flows *ptFlows, int nK,  int *anZ, double* adWeights, short *asCentroids, t_Master *ptMaster, t_Params *ptParams);

void writeZ(t_Params *ptParams, t_Flows *ptFlows, int *anZ, double* adNorm);

char* flowToSeq(int *pnSeqLength, double *adFlows, int nLength);

void writeSequenceUnaligned(FILE* ofp, char *szLabel, char *szSequence, int nLen);

void readInitFile(t_Params *ptParams, int* anZ);

void writeClustersD(t_Params *ptParams, int nK, int nM, int *anT, int* anCentroids, t_Unique *ptUnique, char** aszIDs);

void writeClustersF(t_Params *ptParams, int nN, int nM, int nK, double *adWeight, short* asCentroids);

void outputClusters(t_Params *ptParams, int nK, int *anZ, int nN, int nM, int *anT, t_Flows *ptFlows, char **aszLabels, int* anCentroids, t_Unique *ptUnique, char** aszIDs);

void allocLookUp();

void calcCentroidsMaster(int *anChange, int nKStart, int nKFinish, t_Flows *ptFlows, t_Master *ptMaster, int *anCentroids, t_Unique *ptUnique);

void allocateMaster(t_Master *ptMaster, int nM, int nN, int nK, int nSize);

void reallocateMaster(t_Master *ptMaster);

void updateMasterI(int nI, t_Master *ptMaster, double *adT, double *adD);

void initAlignment(t_Master *ptMaster, t_Flows *ptFlows, int nK, int *anZ, int *anChange, int nI0);

void fillMaster(t_Master *ptMaster);

void destroyMaster(t_Master *ptMaster);

void broadcastMaster(t_Master *ptMaster);

void receiveMaster(t_Master *ptMaster);

void receiveFlows(t_Flows *ptFlows);

void broadcastUnique(t_Unique *ptUnique);

void receiveUnique(t_Unique *ptUnique);

void broadcastFlows(t_Flows *ptFlows);

void writeTau(double *adTau, int nN, int nK, t_Params* ptParams);

void allocateLetter(t_Letter *ptLetter, int nM, int nSize);

void reallocateLetter(t_Letter *ptLetter);

void destroyLetter(t_Letter *ptLetter);

void updateMasterLetter(t_Master *ptMaster, t_Letter *ptLetter);

void minimiseMaster(t_Master *ptMaster);

void minimiseLetter(t_Letter *ptLetter);

double alignX(short* asA, short* asB, int nLenS, int nLenF);

void setMasterZ(t_Master *ptMaster, t_Flows *ptFlows, int nN, int nK, int* anZ);

void writeState(int nN, int* anZ, int nK, int nM, double* adWeight, short* asCentroids, t_Params *ptParams);

void calcQualitiesMaster(int nK, t_Flows *ptFlows, t_Master *ptMaster, char*** paacQualities, int *anCentroids, t_Unique *ptUnique);

void writeQualitiesD(t_Params *ptParams, int nK, int *anT, char** aacQualities, int* anCentroids, t_Unique *ptUnique, char** aszIDs);

void writeMaster(t_Master *ptMaster, t_Unique *ptUnique, int* anCentroids, t_Flows *ptFlows, t_Params *ptParams);

void checkCentroidUniqueness(double *adWeight, int *anCentroids, int nK);

double calcNewWeights(int nK, double *adWeight, t_Master *ptMaster);

void calcDistancesMaster(double* adDist, int* anChange, double *adNorm, t_Master *ptMaster, double dSigma,
			 int* anCentroids, t_Flows *ptFlows, double* adWeight, int nI0, t_Unique *ptUnique);

void calcDistancesSlave(double *adDist, int nM, int nK, int* anChange, double *adLocalNorm, t_Letter *ptLetter, double dSigma,
			int* anCentroids, t_Flows *ptFlows, double* adWeight, int nIStart, int nIFinish, t_Unique *ptUnique);

void receiveLetters(int nTag, double *adNorm, int numtasks, t_Master *ptMaster, t_Letter *ptLetter, int nI0, int nI);

void setZ(int *anZ, double *adP, t_Master* ptMaster, t_Params *ptParams);

double calcNLLikelihood(t_Master *ptMaster, int* pnKEff, double* adWeights, double dSigma);

void outputMap(t_Params *ptParams, int nK, int *anZ, int nN, int nM, int* anT, int *anCentroids, t_Unique *ptUnique, t_Flows *ptFlows);

void allocateUnique(int nN, int nM, t_Unique *ptUnique);

void destroyUnique(t_Unique *ptUnique);

void reallocateUnique(t_Unique *ptUnique, int nNewSize);

void calcUnique(t_Unique *ptUnique, t_Flows *ptFlows);

#endif
