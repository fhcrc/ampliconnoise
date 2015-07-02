/*This software and documentation is copyright © 2009 by Christopher Quince.*/

/*Permission is granted for anyone to copy, use, or modify these programs and documents for purposes of research or education, provided this copyright notice is retained, and note is made of any changes that have been made.*/ 

/* These programs and documents are distributed without any warranty, express or implied. As the programs were written for research purposes only, they have not been tested to the degree that would be advisable in any important application. All use of these programs is entirely at the user's own risk.*/

#ifndef NDIST_MPI_H
#define NDIST_MPI_H

#define FASTA_SUFFIX ".fa"
#define ALIGN_SUFFIX ".afa"
#define TEMP_ERROR_FILE "Temp.err"

#define MAX_DIFF 1000

#define GOOD    0
#define CHIMERA 1
#define TRIMERA 2
#define QUAMERA 3

typedef struct s_Params
{
  char *szSeqInputFile;
  
  char *szRefInputFile;

  int  bImbalance;

  int bOutputAlignments;

  char *szLookUpFile;

  int nSkew;
} t_Params;

typedef struct s_Data
{
  char*  acSequences;
  
  char** aszID;

  double* adFreq;

  int    nSeq;

  int    nMaxLen;

  int*   anLen;
} t_Data;




/**Constants**/
#define MAX_LINE_LENGTH 65536
#define DELIM           " \n"

#define FALSE 0
#define TRUE  1

/* Input Definitions */
#define OPTION  0      /* optional */
#define ALWAYS  1      /* required */

#define LOOKUP_FILE_FLAG     "-rin"
#define REF_INPUT_FILE       "-tin"
#define SEQ_INPUT_FILE       "-sin"
#define USE_IMBALANCE        "-d"
#define OUTPUT_ALIGNMENTS    "-a"
#define SKEW                 "-s"

#define GAP_PENALTY     1.5
#define GAP             '-'
#define T_GAP           '.'
#define COMMA           ","

#define DEFAULT_SKEW          1
#define GAP_PENALTY_N         15.0
#define HOMOPOLYMER_PENALTY   4.0
#define TERMINAL_PENALTY      1.39
#define MIS_MATCH       '#'

#define DIAG  0
#define LEFT  1
#define UP    2

#define MAX_PACKET_SIZE 1048576
#define BIG_DBL         1.0e12
#define BIG_INT         100000
#define WEIGHTDELIM     '_'

#define N_BASES 4
#define LOOKUP_FILE "../Data/Tran.dat"

#ifndef min
#define min(x, y)	((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define	max(x, y)	((x) > (y) ? (x) : (y))
#endif

/* User defined structures */


typedef struct s_Align
{
  char* acA;

  char* acB;

  int   nLen;

  int   nComp;

  int   nDiff;

  double dDiff;

  double dDist;

  int *anMapD;

  double *adD;

  int *anMapR;

  double* adR;

  int *anD;

  int *anR;
} t_Align;


/*User defined functions*/

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[]);

void readData(char* szInputFile, t_Data *ptData);

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[]);

void needlemanWunsch(t_Align *ptAlign, const char* acA, const char* acB, int nLenA, int nLenB);

void receiveData(t_Data *ptData);

void broadcastData(t_Data *ptData);

double needlemanWunschN(const char* acA, const char* acB, int nLenA, int nLenB, int nM);

void initLookUp();

char* getChimera(int *pnCLength, t_Align* ptA, t_Align* ptB, int nSplit, int nLenI);

double calcCIndex(int nI, int nP1, int nP2, char *acChimera, int nCLength, t_Data* ptRefData, t_Data* ptSeqData);

double calcEIndex(int nX1, int nY1, int nX2, int nY2, int nSplit, int nLenI);

int getDifferenceRight(t_Align* ptA, t_Align* ptB, int nSplit, int nLenI);

void allocateMatrices(int nLenI, int ***paanT, int ***paanBestT, int ***paanT2, int ***paanBestT2);

void allocateMatricesD(int nLenI, double ***paadT, int ***paanBestT, double ***paadT2, int ***paanBestT2);

char* getTrimera(int *pnCLength, t_Align* ptA, t_Align* ptB, t_Align* ptC, int nSplit1, int nSplit2, int nLenI);

char* getQuamera(int *pnCLength, t_Align* ptA, t_Align* ptB, t_Align* ptC, t_Align* ptD, int nSplit1, int nSplit2, int nSplit3, int nLenI);

double calcLoonIndex(t_Data *ptSeqData, t_Data *ptRefData, int nI, int nP1, int nP2, int* pnSplit, t_Params *ptParams);

void destroyData(t_Data *ptData);

double distN(char cA, char cB);

int alignAll(int nI, int nLenI, int *pnBest, int *pnBestJ, int *anRestrict, int nK, t_Data *ptSeqData, t_Data *ptRefData, t_Align* atAlign, t_Params *ptParams);

int getBestChimera(int nK, t_Data *ptRefData, int* pnP1, int* pnP2, int *pnSplit, int* anRestrict, int nLenI, t_Align* atAlign, int* anD, int* anR, int* anBestD, int* anBestR);

int getBestTrimera(int nK, t_Data *ptRefData, int* pnT1, int* pnT2, int* pnT3, int *pnSplit1, int *pnSplit2, int* anRestrict, int nLenI, t_Align* atAlign, int* anD, int* anR, int* anBestD, int *anBestR);

double getBestTrimeraD(int nK, t_Data *ptRefData, int* pnT1, int* pnT2, int* pnT3, int *pnSplit1, int *pnSplit2, int* anRestrict, int nLenI, t_Align* atAlign, double* adD, double* adR, int* anBestD, int *anBestR);

double getBestChimeraD(int nK, t_Data *ptRefData, int* pnP1, int* pnP2, int *pnSplit, int* anRestrict, int nLenI, t_Align* atAlign, double* adD, double* adR, int* anBestD, int* anBestR);

int setDR(int nK, t_Data *ptRefData, int* anRestrict, int nLenI, t_Align* atAlign, int* anD, int* anR, int* anBestD, int* anBestR);

#endif
