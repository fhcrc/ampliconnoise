/*This software and documentation is copyright © 2009 by Christopher Quince.*/

/*Permission is granted for anyone to copy, use, or modify these programs and documents for purposes of research or education, provided this copyright notice is retained, and note is made of any changes that have been made.*/ 

/* These programs and documents are distributed without any warranty, express or implied. As the programs were written for research purposes only, they have not been tested to the degree that would be advisable in any important application. All use of these programs is entirely at the user's own risk.*/

#define MAX_LINE_LENGTH 256
#define TRUE  1
#define FALSE 0

#define INIT_SIZE  1024
#define RESOLUTION 0.01

/* Input Definitions */
#define OPTION  0      /* optional */
#define ALWAYS  1      /* required */

#define OUT_FILE_STUB            "-out"
#define INPUT_FILE               "-in"
#define SCALE                    "-s"
#define IDENT                    "-i"
#define WEIGHTS                  "-w"
#define AVERAGE_LINK             "-a"
#define RES                      "-r"

#define LIST_SUFFIX              ".list"
#define OTU_SUFFIX               ".otu"
#define TREE_SUFFIX              ".tree"

typedef struct s_Params
{
  char *szInputFile;

  char *szOutFileStub;

  double dRes;

  int bScale;

  int bIdent;

  int bAverageLink;

  int bWeight;
} t_Params;

typedef struct s_Node
{
  int left; 

  int right; 

  double distance;
} 
t_Node;

#ifndef min
#define min(x, y)	((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define	max(x, y)	((x) > (y) ? (x) : (y))
#endif

t_Node* AvCluster (int nE, float** aafDistMatrix);

t_Node* MaxCluster (int nE, float** aafDistMatrix);

t_Node* AvClusterW (int nE, float** aafDistMatrix, float *afW);

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[]);

void readDistanceMatrix(char ***paszID, char *szDistFile, float *** paafDist, int *pnN, t_Params *ptParams, float** pafW);

int outputCluster(t_Params *ptParams, t_Node* tree, char **aszID, int nN);

void scaleDistances(float **aafDist, int nN);

double writeNodesPhylip(FILE* ofp, t_Node* tree, char **aszID, int nNode);
