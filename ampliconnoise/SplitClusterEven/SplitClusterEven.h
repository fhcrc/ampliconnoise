typedef struct s_Node
{
  struct s_Node *ptLeft;
  double dLeft;


  struct s_Node *ptRight;
  double dRight;

  int nId;

  double dDepth;

  //struct s_Node *ptUp;

  int nR;

  int nN;

  int nLeaves;

  int *anLeaves;
  
} t_Node;

typedef struct s_Params
{
  char* szDatFile;

  char* szTreeFile;

  char* szMapFile;

  int nSplit;

  int nMinSize;
} t_Params;


typedef struct s_Map
{
  int nU;

  int *anNM;

  int **aanMap;

} t_Map;

typedef struct s_Data
{
  char **aszLines;
  
  int nN;

  int nM;
  
} t_Data;

/* Input Definitions */
#define OPTION  0      /* optional */
#define ALWAYS  1      /* required */

#define MAP_FILE                "-min"
#define TREE_FILE               "-tin"
#define SPLIT                   "-s"
#define DAT_FILE                "-din"
#define MIN_SIZE                "-m"

#define LIST_SUFFIX             ".list"
#define TREE_SUFFIX             ".tree"
#define DAT_SUFFIX              ".dat"

#define INTERNAL        -1
#define MAX_WORD_LENGTH 1024
#define MAX_LINE_LENGTH 1000000
#define MAX_SPLIT       10000
#define DELIM           " "
#define BIG_LINE_LENGTH 400000000
#define DELIM2          ",:"


void getCommandLineParams(t_Params *ptParams,int argc,char *argv[]);

char addElement(t_Node** pptTree, FILE* ifp);

void writeTree(t_Node *ptTree, FILE* ofp);

void treeSplit(t_Node* ptTree, double* pdDepth, double dSplit, t_Node **aptSplit,int* pnSplit);

double maxDepth(t_Node* ptTree, double* pdMaxDepth);

void writeIndices(t_Node* ptTree, FILE* ofp);

void countLeaves(t_Node* ptTree, int *pnLeaves);

void getLeaves(t_Node* ptTree, int* anIDs, int *pnLeaves);

void getLeavesR(t_Node* ptTree, int* anIDs, int *pnLeaves);

int compNode(const void *pvA, const void* pvB);

void readData(t_Data *ptData, char* szDatFile);

void writeData(FILE* ofp, t_Data *ptData, int nEntries, int *anEntries, t_Map *ptMap);

void renumberLeaves(t_Node* ptTree, int* pnCount);

void writeList(FILE* ofp, t_Node *ptTree, double dMaxSplit);

void readMapFile(t_Map *ptMap, char *szMapFile);

//void treeSplit(t_Node* ptTree, double dDepth, double dSplit, t_Node **aptSplit, int Split);

void setLeaves(t_Node* ptTree);

void treeSplitEven(t_Node* ptTree, int nSplit, t_Node **aptSplit,int* pnSplit);

double setDepth(t_Node* ptTree, double dDepth);
