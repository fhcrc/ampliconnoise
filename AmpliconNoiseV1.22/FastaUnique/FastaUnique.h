typedef struct s_Params
{
  char *szMapFile;

  char *szInputFile;

  int  nMin;

  int  bTrim;
} t_Params;

typedef struct s_Data
{
  char** aacSequences;
  
  char** aszID;

  int    nSeq;

  int    nMaxLen;

  int*   anLen;
} t_Data;

typedef struct s_Map
{
  int nSeq;

  int **aanMap;

  int *anN;

  int *anSize;
} t_Map;

#define MAX_LINE_LENGTH 1024
#define FALSE 0
#define TRUE  1
#define INIT_MAP_SIZE 128

#define FASTA_SUFFIX    ".fa"
#define MAP_SUFFIX      ".map"


/* Input Definitions */
#define OPTION  0      /* optional */
#define ALWAYS  1      /* required */

#define INPUT_FILE       "-in"
#define TRIM             "-t"
#define MIN              "-m"

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[]);

void allocMap(t_Map *ptMap, int nSeq);

void addEntry(t_Map *ptMap, int nI, int nJ);

void writeMap(FILE *ofp, t_Data *ptData, t_Map *ptMap);

void destroyMap(t_Map *ptMap);

char* szCreateMapFileName(char *szFastaFile);
