/*
August 2012 Code has been optimised for speed by Tessella, funded by Unilever.
See comments in code for details of optimisation changes.
Output has been verified to be unaffected by these changes.
*/

#ifndef NDIST_MPI_H
#define NDIST_MPI_H

typedef struct s_Params
{
  char *szInputFile;

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

#define INPUT_FILE       "-in"
#define IDENT            "-i"
#define PHYLIP           "-p"

#define GAP_PENALTY     1.5
#define GAP             '-'

#define DIAG  0
#define LEFT  1
#define UP    2

#define MAX_PACKET_SIZE 1048576


#ifndef min
#define min(x, y)	((x) < (y) ? (x) : (y))
#endif

#ifndef max
#define	max(x, y)	((x) > (y) ? (x) : (y))
#endif

/* August 2012
   Speed optimisation 5B: inline loop functions */
#ifndef dist
#define	dist(a, b)	((a) == (b) ? 0 : 1)
#endif

#ifndef getMove
#define	getMove(a, b, c)	((a) < (b) ? ((a) >= (c) ? UP : DIAG) : ((b) < (c) ? LEFT : UP))
#endif

#ifndef dmin3
#define dmin3(a, b, c)	((a) < (b) ? ((c) < (a) ? (c) : (a)) : ((c) < (b) ? (c) : (b)))
#endif

/* User defined structures */


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

void needlemanWunsch(t_Align *ptAlign, const char* acA, const char* acB, int nLenA, int nLenB, double** aadFMatrix, int** aanMoves);

#endif
