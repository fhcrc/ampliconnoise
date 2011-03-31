
/*This software and documentation is copyright © 2009 by Christopher Quince.*/

/*Permission is granted for anyone to copy, use, or modify these programs and documents for purposes of research or education, provided this copyright notice is retained, and note is made of any changes that have been made.*/ 

/* These programs and documents are distributed without any warranty, express or implied. As the programs were written for research purposes only, they have not been tested to the degree that would be advisable in any important application. All use of these programs is entirely at the user's own risk.*/

#ifndef SCLUSTER_H
#define SCLUSTER_H

/*User defined structures*/

typedef struct s_Data
{
  char*    acSequences;

  char**   aszID;

  int      nSeq;

  int      nMaxLen;

  int*     anLen;

  double*  adW;
} t_Data;

typedef struct s_Align
{
  char* acA;

  char* acB;

  int   nLen;

  int   nComp;

  double dDist;
} t_Align;

typedef struct s_Params
{
  char *szDataFile;

  char *szDistFile;

  char *szOutFileStub;

  char *szInitFile;

  char *szMappingFile;

  char *szLookUpFile;

  double dSigma;

  double dInitCut;

  int nK;
} t_Params;

typedef struct s_Letter
{
  int nSize;

  int nTotal;

  int    *anK;

  int    *anI;

  double *adT;

} t_Letter;

typedef struct s_Master
{
  double   *adTau; /*probs*/
  
  int      **aanP;   /*pointers*/
  int      **aanI;   /*members*/

  int      nTotal;

  /*flat arrays*/
  int      *anN;  
  int      *anCN;
  int      *anP;
  int      *anI;
  int      nSize;
  
  int      nK;
  int      nN;
} t_Master;

/* Input Definitions */
#define OPTION  0      /* optional */
#define ALWAYS  1      /* required */

/*****Input Constants**************/
#define OUT_FILE_STUB            "-out"
#define INIT_FILE                "-lin"
#define DATA_FILE                "-in"
#define MAPPING_FILE             "-min"
#define LOOKUP_FILE_FLAG         "-rin"
#define DISTANCE_FILE            "-din"
#define SIGMA                    "-s"
#define INIT_CUT                 "-c"
#define VERBOSE                  "-v"
/**********************************/


/********DEFAULTS***********/
#define DEF_SIGMA    30.0
#define DEF_CUT      0.08
#define LOOKUP_FILE  "../Data/Tran.dat"
/***************************/

/******Programming Constants************/
#define MAX_LINE_LENGTH 65536
#define BIG_LINE_LENGTH 1048576
#define DELIM           " \n"
#define DELIMT          "\t\n"
#define COMMA           ","

#define FALSE 0
#define TRUE  1

#define INIT_MASTER_SIZE 1024

#define MAX_DBL  1.0e10

#define NOT_SET  -1
/*******Alignment constants*************/
#define GAP_PENALTY         15.0
#define HOMOPOLYMER_PENALTY 4.0
#define GAP             '-'
#define T_GAP           '.'

#define DIAG  0
#define LEFT  1
#define UP    2

#define NOTSET_CHAR     '\0'
/***************************************/

/*Clustering constants*/
#define MIN_DELTA 1.0e-6
#define MIN_ITER 20
#define MAX_ITER 1000

#define MIN_COUNT 0.1
#define MIN_COUNT_LATE 0.1

#define MIN_TAU   1.0e-4
#define MIN_TAU_LATE 1.0e-4

#define MIN_WEIGHT 0.1
#define MIN_WEIGHT_LATE 0.1
/*********************/

/****File Suffices****/
#define FASTA_SUFFIX             ".fa"
#define TAU_SUFFIX               ".tau"
#define Z_SUFFIX                 ".z"
#define MASTER_SUFFIX            ".master"
#define MAP_SUFFIX               ".mapping"
/***** macros**********************************/
#define aindex(i,k,j)   ((i)*nM*nK + nM*(k) + (j))
#ifndef min
#define min(x, y)	((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define	max(x, y)	((x) > (y) ? (x) : (y))
#endif
/**********************************************/

#define TYPE_A   0
#define TYPE_C   1
#define TYPE_T   2
#define TYPE_G   3
#define TYPE_GAP 4
#define N_BASES  4 /*4 bases order ACTG*/
#define N_TYPES  5 /*4 bases order ACTG then gap*/

void getCommandLineParams(t_Params *ptParams, int argc, char *argv[]);

void readData(t_Data *ptData, t_Params *ptParams);

void destroyData(t_Data *ptData);

double needlemanWunsch(int nM, const char* acA, const char* acB, int nLenA, int nLenB);

void allocateMaster(t_Master *ptMaster, int nN, int nK, int nSize);

void reallocateMaster(t_Master *ptMaster);

void minimiseLetter(t_Letter *ptLetter);

void minimiseMaster(t_Master *ptMaster);

void updateMasterI(int nI, t_Master *ptMaster, double *adT);

void updateMasterLetter(t_Master *ptMaster, t_Letter *ptLetter);

void initAlignment(t_Master *ptMaster, t_Data *ptData, int nK, int *anZ, int *anChange, int nI0);

void fillMaster(t_Master *ptMaster);

void destroyMaster(t_Master *ptMaster);

void broadcastData(t_Data *ptData);

void broadcastMaster(t_Master *ptMaster);

void receiveMaster(t_Master *ptMaster);

void writeTau(double *adTau, int nN, int nK, t_Params* ptParams);

void allocateLetter(t_Letter *ptLetter, int nSize);

void reallocateLetter(t_Letter *ptLetter);

void destroyLetter(t_Letter *ptLetter);

void calcCentroidsMaster(int *anChange, int nKStart, int nKFinish, t_Master *ptMaster, int *anCentroids, double* adW, float *afDist);

void receiveData(t_Data *ptData);

void readInitFile(t_Params *ptParams, int* anZ);

void writeZ(t_Params *ptParams, t_Data *ptData, int *anZ);

void writeClustersF(t_Params *ptParams, int nK, double *adWeight, t_Data *ptData, int* anCentroids);

void writeClustersD(t_Params *ptParams, int nN, int nK, int *anT, int *anZ, t_Data *ptData, int* anCentroids, double* adW);

void outputClusters(t_Params *ptParams, int nK, int* anT, int *anZ, int nN, int nM, t_Data *ptData, int* anCentroids);

double calcNewWeights(int nK, double *adWeight, t_Master *ptMaster, double* adW);

void setZ(int *anZ, double *adP, t_Master* ptMaster, t_Params *ptParams);

void writeSequenceUnaligned(FILE* ofp, char *szLabel, char *szSequence, int nLen);

void checkCentroidUniqueness(double *adWeight, int *anCentroids, int nK);

void initLookUp(t_Params *ptParams);

void setMasterZ(t_Master *ptMaster, int nN, int nK, int* anZ);

void writeMaster(int nIter, t_Master *ptMaster, int* anCentroids, t_Data *ptData, t_Params *ptParams, double *adW);

void outputMap(t_Params *ptParams, int nK, int* anT, int *anZ, int nN, int nM, char **aszLabels, double *adW, int* anCentroids);

void outputMapping(t_Params *ptParams, int nK, int *anT, int *anZ, int nN, double *adW, int *anCentroids, char** aszID);

void readDistanceMatrix(char *szDistFile, float **pafDist, int nN);

#endif
