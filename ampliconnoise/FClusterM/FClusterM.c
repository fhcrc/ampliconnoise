/*This software and documentation is copyright © 2009 by Christopher Quince.*/

/*Permission is granted for anyone to copy, use, or modify these programs and documents for purposes of research or education, provided this copyright notice is retained, and note is made of any changes that have been made.*/ 

/* These programs and documents are distributed without any warranty, express or implied. As the programs were written for research purposes only, they have not been tested to the degree that would be advisable in any important application. All use of these programs is entirely at the user's own risk.*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "FClusterM.h"

static char *usage[] = {"FCluster\n",
			"-in    string            distance input file name\n",
			"-out   string            output file stub\n",
			"Options:\n",
			"-r                       resolution\n",
			"-a                       average linkage\n",
			"-w                       use weights\n",
			"-i                       read identifiers\n",
                        "-s                       scale dist.\n"};

static int nLines = 9;

int main(int argc, char* argv[])
{
  float **aafDist = NULL;
  int    i = 0, j = 0, nN        = 0;
  t_Node *atTree   = NULL;
  char   **aszID   = NULL;
  float  *afW      = NULL;
  t_Params tParams;
  char   *szTreeFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE   *ofp = NULL;

  getCommandLineParams(&tParams, argc, argv);

  readDistanceMatrix(&aszID, tParams.szInputFile, &aafDist, &nN, &tParams, &afW);

  for(i = 0; i < nN; i++){
    for(j = 0; j < i; j++){
      aafDist[i][j] = 0.5*(aafDist[i][j] + aafDist[j][i]);
    }
  }

  if(tParams.bScale){
    scaleDistances(aafDist, nN);
  }
  else{
    printf("Not scaling dist\n");
  }
  printf("Read data n = %d\n",nN); fflush(stdout);
  
  if(tParams.bAverageLink == FALSE){
    atTree = MaxCluster(nN, aafDist);
  }
  else{
    if(tParams.bWeight){
      atTree = AvClusterW(nN, aafDist, afW);
    }
    else{
      atTree = AvCluster(nN, aafDist);
    }
  }
  sprintf(szTreeFile, "%s%s", tParams.szOutFileStub, TREE_SUFFIX);

  ofp = fopen(szTreeFile, "w");

  if(ofp){
    writeNodesPhylip(ofp, atTree, aszID, nN - 2);
    fprintf(ofp,";\n");
    fclose(ofp);
  }
  else{
    fprintf(stderr, "Failed to open %s for writing\n", szTreeFile);
    fflush(stderr);
  }
  outputCluster(&tParams, atTree, aszID, nN);

  /*free up memory*/
  free(atTree);
  for(i = 0; i < nN; i++){
    if(i > 0){
	free(aafDist[i]);
    }
    free(aszID[i]);
  }
  free(aafDist);
  free(afW);
  free(aszID);
  free(szTreeFile);

  exit(EXIT_SUCCESS);
}

void writeUsage(FILE* ofp)
{
  int i = 0;
  char *line;

  for(i = 0; i < nLines; i++){
    line = usage[i];
    fputs(line,ofp);
  }
}

char *extractParameter(int argc, char **argv, char *param,int when)
{
  int i = 0;

  while((i < argc) && (strcmp(param,argv[i]))){
    i++;
  }

  if(i < argc - 1){
    return(argv[i + 1]);
  }

  if((i == argc - 1) && (when == OPTION)){
    return "";
  }

  if(when == ALWAYS){
    fprintf(stdout,"Can't find asked option %s\n",param);
  }

  return (char *) NULL;
}

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[])
{
  char *szTemp = NULL;
  char *pcError = NULL;

  /*get parameter file name*/
  ptParams->szInputFile  = extractParameter(argc,argv,INPUT_FILE,ALWAYS);  
  if(ptParams->szInputFile == NULL)
    goto error;

  /*get output filestub*/
  ptParams->szOutFileStub  = extractParameter(argc,argv, OUT_FILE_STUB, ALWAYS);  
  if(ptParams->szOutFileStub == NULL)
    goto error;

  if(extractParameter(argc,argv, SCALE, OPTION)){
    ptParams->bScale = TRUE;
  }
  else{
    ptParams->bScale = FALSE;
  }

  if(extractParameter(argc,argv, IDENT, OPTION)){
    ptParams->bIdent = TRUE;
  }
  else{
    ptParams->bIdent = FALSE;
  }

  if(extractParameter(argc,argv, WEIGHTS, OPTION)){
    ptParams->bWeight = TRUE;
  }
  else{
    ptParams->bWeight = FALSE;
  }

  if(extractParameter(argc,argv, AVERAGE_LINK, OPTION)){
    ptParams->bAverageLink = TRUE;
  }
  else{
    ptParams->bAverageLink = FALSE;
  }

  szTemp = extractParameter(argc,argv, RES, OPTION);
  if(szTemp != NULL){
    ptParams->dRes = strtod(szTemp,&pcError);
    if(*pcError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->dRes = RESOLUTION;
  }

  return;

 error:
  writeUsage(stdout);
  exit(EXIT_FAILURE);
}

float getClosest(int nN, float** aafDistMatrix, int* nI, int* nJ)
{ 
  int i, j;
  float  fTemp;
  float  fDist = aafDistMatrix[1][0];
  *nI = 1;
  *nJ = 0;
  for (i = 1; i < nN; i++){ 
    for (j = 0; j < i; j++){ 
      fTemp = aafDistMatrix[i][j];
      
      if(fTemp < fDist){ 
	fDist = fTemp;
        *nI = i;
        *nJ = j;
      }
    }
  }

  return fDist;
}

char* revstr(char *szID, char cD)
{
  int nLen = strlen(szID);
  int i = nLen - 1;

  while(i >= 0 && szID[i] != cD){
    i--;
  }

  if(i < 0){
    return NULL;
  }
  else{
    return &szID[i];
  }
}

double getWeight(char *szID)
{
  char   *szBreak  = revstr(szID, '_'); 
  double dWeight   = 0.0;
  char   *pcError  = NULL;
  
 
  if(szBreak != NULL){
    dWeight = strtod(szBreak + 1, &pcError);
    if(*pcError != '\0'){
      fprintf(stderr, "Sequence weight format error %s %s\n",szID,szBreak + 1);
      fflush(stderr);
    }
  }
  else{
    dWeight = 1.0;
  }

  return dWeight;
}

void readDistanceMatrix(char*** paszID, char *szDistFile, float *** paafDist, int *pnN, t_Params *ptParams, float** pafW)
{
  FILE *ifp = NULL;
  char szLine[MAX_LINE_LENGTH];
  char *szTok = NULL, *pcError = NULL, *szBrk = NULL;
  float **aafDist = NULL;
  int    i = 0, j = 0, nN = 0;
  char **aszID;
  float* afW = NULL;

  ifp = fopen(szDistFile, "r");

  if(ifp){
    fgets(szLine, MAX_LINE_LENGTH, ifp);
    
    szBrk = strpbrk(szLine, "\n");

    (*szBrk) = '\0';

    nN = strtol(szLine, &pcError, 10);
    if(*pcError != '\0')
      goto fileFormatError;
    
    afW = (float *) malloc(nN*sizeof(float));
    if(!afW)
      goto memoryError;

    for(i = 0; i < nN; i++){
      afW[i] = 1.0;
    }

    aafDist = (float **) malloc(nN*sizeof(float *));
    if(!aafDist)
      goto memoryError;

    aszID = (char **) malloc(nN*sizeof(char*));
    if(!aszID)
      goto memoryError;

    for(i = 0; i < nN; i++){
      aszID[i] = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
      if(!aszID[i])
	goto memoryError;
    }

    if(ptParams->bIdent == TRUE || ptParams->bWeight == TRUE){
      for(i = 0; i < nN; i++){
	fgets(szLine, MAX_LINE_LENGTH, ifp);
	szBrk = strpbrk(szLine, "\n");
	(*szBrk) = '\0';

	strcpy(aszID[i],szLine);
      }
    }
    else{
      for(i = 0; i < nN; i++){
	sprintf(aszID[i],"%d",i);
      }
    }

    if(ptParams->bWeight == TRUE){
      for(i = 0; i < nN; i++){
	afW[i] = (float) getWeight(aszID[i]);
	
	if(ptParams->bIdent != TRUE){
	  sprintf(aszID[i],"%d",i);
	}
      }
    }

    for(i = 0; i < nN; i++){
      
      aafDist[i] = (float *) malloc(nN*sizeof(float));
      if(!aafDist[i])
	goto memoryError;

      for(j = 0; j < nN; j++){
	fgets(szLine, MAX_LINE_LENGTH, ifp);
    
	szBrk = strpbrk(szLine, "\n");
	(*szBrk) = '\0';

	aafDist[i][j] = (float) strtod(szLine, &pcError);
	
	if(*pcError != '\0')
	  goto fileFormatError;
      }
    }

    fclose(ifp);
  }
  else{
    fprintf(stderr, "Failed to open file %s in readDistanceMatrix",szDistFile);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  (*pafW) = afW;
  (*pnN) = nN;
  (*paafDist) = aafDist;
  (*paszID) = aszID;
  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in readDistanceMatrix");
  fflush(stderr);
  exit(EXIT_FAILURE);

 fileFormatError:
  fprintf(stderr, "Incorrectly formatted input file in readDistanceMatrix");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

t_Node* MaxCluster(int nE, float** aafDistMatrix)
{ 
  int j;
  int n;
  int* anCI;
  t_Node* atResult;

  anCI = malloc(nE*sizeof(int));
  if(!anCI) return NULL;
  atResult = malloc((nE-1)*sizeof(t_Node));
  if (!atResult)
  { free(anCI);
    return NULL;
  }

  for (j = 0; j < nE; j++) anCI[j] = j;

  for (n = nE; n > 1; n--)
  { 
    int is = 1;
    int js = 0;
    atResult[nE-n].distance = getClosest(n, aafDistMatrix, &is, &js);

    for (j = 0; j < js; j++)
      aafDistMatrix[js][j] = max(aafDistMatrix[is][j],aafDistMatrix[js][j]);
    for (j = js+1; j < is; j++)
      aafDistMatrix[j][js] = max(aafDistMatrix[is][j],aafDistMatrix[j][js]);
    for (j = is+1; j < n; j++)
      aafDistMatrix[j][js] = max(aafDistMatrix[j][is],aafDistMatrix[j][js]);

    for (j = 0; j < is; j++) aafDistMatrix[is][j] = aafDistMatrix[n-1][j];
    for (j = is+1; j < n-1; j++) aafDistMatrix[j][is] = aafDistMatrix[n-1][j];

    atResult[nE-n].left = anCI[is];
    atResult[nE-n].right = anCI[js];
    anCI[js] = n-nE-1;
    anCI[is] = anCI[n-1];
    if(n % 100 == 0){
      printf("%d ",n);
      fflush(stdout);
    }
  }
  free(anCI);
  printf("\n"); fflush(stdout);
  return atResult;
}

t_Node* AvCluster(int nE, float** aafDistMatrix)
{ 
  int j;
  int n;
  int* anCI;
  int* anNumber;
  t_Node* atResult;

  anCI = malloc(nE*sizeof(int));
  if(!anCI) return NULL;
  anNumber = malloc(nE*sizeof(int));
  if(!anNumber)
  { 
    free(anCI);
    return NULL;
  }
  atResult = malloc((nE-1)*sizeof(t_Node));
  if (!atResult)
  { 
    free(anCI);
    free(anNumber);
    return NULL;
  }

  for (j = 0; j < nE; j++)
  { 
    anNumber[j] = 1;
    anCI[j] = j;
  }

  for (n = nE; n > 1; n--)
  { 
    int nSum;
    int is = 1;
    int js = 0;
    
    atResult[nE-n].distance = getClosest(n, aafDistMatrix, &is, &js);

    atResult[nE-n].left = anCI[is];
    atResult[nE-n].right = anCI[js];

    nSum = anNumber[is] + anNumber[js];
    
    for (j = 0; j < js; j++)
    { 
      aafDistMatrix[js][j] = aafDistMatrix[is][j]*anNumber[is]
                        + aafDistMatrix[js][j]*anNumber[js];
      aafDistMatrix[js][j] /= nSum;
    }
    
    for (j = js+1; j < is; j++)
    { 
      aafDistMatrix[j][js] = aafDistMatrix[is][j]*anNumber[is]
                        + aafDistMatrix[j][js]*anNumber[js];
      aafDistMatrix[j][js] /= nSum;
    }
    
    for (j = is+1; j < n; j++)
    { 
      aafDistMatrix[j][js] = aafDistMatrix[j][is]*anNumber[is]
                        + aafDistMatrix[j][js]*anNumber[js];
      aafDistMatrix[j][js] /= nSum;
    }

    for (j = 0; j < is; j++) aafDistMatrix[is][j] = aafDistMatrix[n-1][j];
    for (j = is+1; j < n-1; j++) aafDistMatrix[j][is] = aafDistMatrix[n-1][j];

    anNumber[js] = nSum;
    anNumber[is] = anNumber[n-1];

    anCI[js] = n-nE-1;
    anCI[is] = anCI[n-1];
  }
  free(anCI);
  free(anNumber);

  return atResult;
}

t_Node* AvClusterW(int nE, float** aafDistMatrix, float* afW)
{ 
  int j;
  int n;
  int* anCI;
  float* afNumber;
  t_Node* atResult;

  anCI = malloc(nE*sizeof(int));
  
  if(!anCI) return NULL;
  
  afNumber = malloc(nE*sizeof(float));
  
  if(!afNumber)
  { 
    free(anCI);
    return NULL;
  }
  
  atResult = malloc((nE-1)*sizeof(t_Node));
  
  if (!atResult)
  { 
    free(anCI);
    free(afNumber);
    return NULL;
  }

  for (j = 0; j < nE; j++)
  { afNumber[j] = afW[j];
    anCI[j] = j;
  }

  for (n = nE; n > 1; n--)
  { 
    float fSum;
    int is = 1;
    int js = 0;
    
    atResult[nE-n].distance = getClosest(n, aafDistMatrix, &is, &js);
    atResult[nE-n].left = anCI[is];
    atResult[nE-n].right = anCI[js];

    fSum = afNumber[is] + afNumber[js];
    
    for (j = 0; j < js; j++)
    { 
      aafDistMatrix[js][j] = aafDistMatrix[is][j]*afNumber[is]
                        + aafDistMatrix[js][j]*afNumber[js];
      aafDistMatrix[js][j] /= fSum;
    }
    
    for (j = js+1; j < is; j++)
    { 
      aafDistMatrix[j][js] = aafDistMatrix[is][j]*afNumber[is]
                        + aafDistMatrix[j][js]*afNumber[js];
      aafDistMatrix[j][js] /= fSum;
    }
    
    for (j = is+1; j < n; j++)
    { 
      aafDistMatrix[j][js] = aafDistMatrix[j][is]*afNumber[is]
                        + aafDistMatrix[j][js]*afNumber[js];
      aafDistMatrix[j][js] /= fSum;
    }

    for (j = 0; j < is; j++) aafDistMatrix[is][j] = aafDistMatrix[n-1][j];
    for (j = is+1; j < n-1; j++) aafDistMatrix[j][is] = aafDistMatrix[n-1][j];

    afNumber[js] = fSum;
    afNumber[is] = afNumber[n-1];

    anCI[js] = n-nE-1;
    anCI[is] = anCI[n-1];
  }
  free(anCI);
  free(afNumber);

  return atResult;
}

int outputCluster(t_Params *ptParams, t_Node* tree, char **aszID, int nN)
{
  int i = 0, j = 0, k = 0, l = 0;
  int nNodes = nN - 1;
  double dRes = ptParams->dRes;
  int    nIter = 0; 
  char ***aaszClusters = NULL;
  int  *anCounts = (int *) malloc(sizeof(int)*nN);
  int  *anSizes  = (int *) malloc(sizeof(int)*nN);
  int  *anIndex  = (int *) malloc(sizeof(int)*nNodes);
  char *szOtuFile  = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  char *szListFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL, *lfp = NULL;

  sprintf(szOtuFile, "%s%s", ptParams->szOutFileStub, OTU_SUFFIX);
  sprintf(szListFile, "%s%s", ptParams->szOutFileStub, LIST_SUFFIX);

  aaszClusters = (char ***) malloc(nN*sizeof(char**));

  for(i = 0; i < nN; i++){
    anCounts[i] = 1;
    anSizes[i]  = INIT_SIZE; 
    aaszClusters[i] = (char **) malloc(anSizes[i]*sizeof(char *));
    for(j = 0; j < anSizes[i]; j++){
      aaszClusters[i][j] = NULL;
    }
    aaszClusters[i][0] = aszID[i];
  }

  ofp = fopen(szOtuFile,"w");
  lfp = fopen(szListFile, "w");

  nIter = (int) floor(tree[nNodes - 1].distance/ptParams->dRes + 1.0e-7);

  j = 0;
  for(i = 0; i < nIter; i++){
    double dPrint = i*dRes;

    while(j < nNodes && tree[j].distance < dPrint){
      int left = tree[j].left, right = tree[j].right;
      int put  = -1, destroy = -1;

      if(left < 0){
	left = anIndex[-left - 1];
      }

      if(right < 0){
	right = anIndex[-right - 1];
      }

      put = right; destroy = left;
      if(left < right){
	put     = left;
	destroy = right;
      }
   
      anIndex[j] = put;
    

      for(k = 0; k < anCounts[destroy]; k++){
	if(anCounts[put] == anSizes[put]){
	  anSizes[put] *= 2;
	  aaszClusters[put] = (char **) realloc(aaszClusters[put],anSizes[put]*sizeof(char *));
	}
	aaszClusters[put][anCounts[put]++] = aaszClusters[destroy][k];
      }
      anCounts[destroy] = 0;
      free(aaszClusters[destroy]);

      j++;
    }

    fprintf(lfp, "%f %d ",dPrint, nN - j);
    fprintf(ofp, "%f %d ",dPrint, nN - j);
    dPrint += 0.01;

    for(k = 0; k < nN; k++){
      if(anCounts[k] > 0){
	for(l = 0; l < anCounts[k] - 1; l++){
	  fprintf(lfp, "%s,",aaszClusters[k][l]);
	}

	fprintf(lfp,"%s ",aaszClusters[k][anCounts[k] - 1]);

	fprintf(ofp, "%d ", anCounts[k]);
      }
    }
     
    fprintf(lfp, "\n"); fprintf(ofp, "\n");
    fflush(lfp); fflush(ofp);
  }
  
  fclose(lfp);
  fclose(ofp);

  for(i = 0; i < nN; i++){
    if(anCounts[i] > 0){
      free(aaszClusters[i]);
    }
  }
  free(aaszClusters);
  /*free up memory*/
  free(anCounts); 
  free(anSizes);
  free(anIndex); 
  free(szOtuFile);
  free(szListFile);
  return;
  
}

double writeNodesPhylip(FILE* ofp, t_Node* tree, char **aszID, int nNode)
{
  int nLeft = tree[nNode].left, nRight = tree[nNode].right;
  double dist = tree[nNode].distance;
  double lastdist = 0.0;
  fprintf(ofp, "(");

  if(nLeft >= 0){
     fprintf(ofp, "%s", aszID[nLeft]);
  }
  else{
    lastdist = writeNodesPhylip(ofp, tree, aszID, -nLeft - 1);
  }

  fprintf(ofp, ":%.3f,",dist - lastdist);

  lastdist = 0.0;
  if(nRight >=0){
    fprintf(ofp, "%s", aszID[nRight]);
  }
  else{
    lastdist = writeNodesPhylip(ofp, tree, aszID, -nRight - 1);
  }
 
  fprintf(ofp, ":%.3f)",dist - lastdist);
  
  return dist;
}   

void scaleDistances(float **aafDist, int nN)
{
  int i = 0, j = 0;
  float fMin = aafDist[1][0], fMax = aafDist[1][0];

  for(i = 2; i < nN; i++){
    for(j = 0; j < i; j++){
      if(aafDist[i][j] < fMin){
	fMin = aafDist[i][j];
      }

      if(aafDist[i][j] > fMax){
	fMax = aafDist[i][j];
      }
    }
  }

  for(i = 1; i < nN; i++){
    for(j = 0; j < i; j++){
      aafDist[i][j] = (aafDist[i][j] - fMin)/(fMax - fMin);
    }
  }
}
