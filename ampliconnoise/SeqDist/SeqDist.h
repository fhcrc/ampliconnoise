/*This software and documentation is copyright © 2009 by Christopher Quince.*/

/*Permission is granted for anyone to copy, use, or modify these programs and documents for purposes of research or education, provided this copyright notice is retained, and note is made of any changes that have been made.*/ 

/* These programs and documents are distributed without any warranty, express or implied. As the programs were written for research purposes only, they have not been tested to the degree that would be advisable in any important application. All use of these programs is entirely at the user's own risk.*/

#ifndef NDIST_MPI_H
#define NDIST_MPI_H

typedef struct s_Params
{
  char *szInputFile;

  char *szLookUpFile;

  int  bIdent;

  int  bPhylip;
} t_Params;

typedef struct s_Data
{
  char*  acSequences;
  
  char** aszID;

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

#define LOOKUP_FILE_FLAG "-rin"
#define INPUT_FILE       "-in"
#define IDENT            "-i"
#define PHYLIP           "-p"

#define GAP_PENALTY         15.0
#define HOMOPOLYMER_PENALTY 4.0
#define GAP             '-'
#define T_GAP           '.'

#define COMMA           ","

#define DIAG  0
#define LEFT  1
#define UP    2

#define MAX_PACKET_SIZE 1048576
#define N_BASES 4

#ifndef min
#define min(x, y)	((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define	max(x, y)	((x) > (y) ? (x) : (y))
#endif

/* User defined structures */

#define LOOKUP_FILE "../Data/Tran.dat"


typedef struct s_Align
{
  //  char* acA;

  //char* acB;

  int   nLen;

  int   nComp;

  double dDist;
} t_Align;


/*User defined functions*/

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[]);

void readData(t_Data *ptData, t_Params *ptParams);

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[]);

double needlemanWunsch(const char* acA, const char* acB, int nLenA, int nLenB, int nM);

void initLookUp(t_Params *ptParams);

#endif
