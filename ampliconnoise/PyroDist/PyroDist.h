/*This software and documentation is copyright © 2009 by Christopher Quince.*/

/*Permission is granted for anyone to copy, use, or modify these programs and documents for purposes of research or education, provided this copyright notice is retained, and note is made of any changes that have been made.*/ 

/* These programs and documents are distributed without any warranty, express or implied. As the programs were written for research purposes only, they have not been tested to the degree that would be advisable in any important application. All use of these programs is entirely at the user's own risk.*/

#ifndef CLASSIFY_H
#define CLASSIFY_H

/**Constants**/
#define MAX_LINE_LENGTH 65536
#define DELIM           " \n"

#define FALSE 0
#define TRUE  1

/* Input Definitions */
#define OPTION  0      /* optional */
#define ALWAYS  1      /* required */

#define DATA_FILE       "-in"
#define NO_INDEX        "-ni"
#define OUT_FILE_STUB   "-out"
#define LOOKUP_FILE_FLAG "-rin"

#define GAP_PENALTY     15.0
#define MAX_DIST        100.0
#define MAX_DBL         1.0e10
#define GAP             -2.0

#define SUFFIX          ".fdist"

#define DIAG  0
#define LEFT  1
#define UP    2

#define LOOKUP_FILE "../Data/LookUp_E123.dat"
#define MAX_PACKET_SIZE 1048576

#define BINS  1000
#define MAX_S 10
#define DELTA 1.0e-7
#define PRECISION 0.01

#ifndef min
#define min(x, y)	((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define	max(x, y)	((x) > (y) ? (x) : (y))
#endif

/* User defined structures */

typedef struct s_Flows
{
  /*number of flows*/
  int nN;
  /*max flow length*/
  int nM;
  /*flows*/
  int   * anF;
  double* adF;
  double* adData;
  /*clean lengths*/
  int*    anLengths;
} t_Flows;

typedef struct s_Params
{
  char *szOutFileStub;

  char *szDataFile;

  char *szLookUpFile;

  int bNoIndex;
} t_Params;


typedef struct s_Align
{
  double *adA;
  double *adB;

  int nLen;

  int nComp;

  double dDist;
} t_Align;


/*User defined functions*/

void getCommandLineParams(t_Params *ptParams, int argc, char *argv[]);

void readData(double *adLookUp, char* szDataFile, t_Flows *ptFlows, t_Params *ptParams);

int imin(int nA, int nB);

double dmin3(double dA, double dB, double dC);

int getMove(double dA, double dB, double dC);

void initLookUp(t_Params *ptParams, double* adLookUp);

void initFookUp(double* adFookUp, double *adLookUp);

void needlemanWunsch(double *adFookUp, double* adF1, double* adF2, t_Align *ptAlign, int* anA, int* anB, int nLenA, int nLenB);

void broadcastFlows(t_Flows* ptFlows);

void receiveFlows(t_Flows* ptFlows);

double distI(double *adLookUp, int nFlow1, int nFlow2);

double distM(double *adLookUp, int nFlow1, int nFlow2);

double distS(double* adLookUp, int nS, double dFlow);

double distF1(double *adLookUp, double dF1);

double distM1(double *adLookUp, double dF1);

double FDist(double *adFookUp, double* adF1, double* adF2, int* anA, int* anB, int nLenA, int nLenB);

#endif
