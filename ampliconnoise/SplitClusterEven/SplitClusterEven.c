#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>

#include "SplitClusterEven.h"

static char *usage[] = {"SplitClusterEven\n",
			"-din    string  dat filename\n",
			"-min    string  map filename\n",
			"-tin    string  tree filename\n",
                        "-s      split size\n",
			"-m      min   size\n"};

static int nLines = 6;

int main(int argc, char* argv[])
{
  t_Node    *ptTree = NULL;
  FILE      *ifp    = NULL;
  t_Params  tParams;
  t_Data    tData;
  t_Map     tMap;
  t_Node     **aptSplit = (t_Node **) malloc(MAX_SPLIT*sizeof(t_Node *));
  int        nSplit = 0;
  double     dMaxDepth = 0.0, dSplitDepth = 0.0, dDepth = 0.0;
  int        iL = 0, nLast = 0, nCount = 0, i = 0, j = 0;
  int*       anLast = NULL;
  char szDir[MAX_WORD_LENGTH];
  char szTreeFile[MAX_WORD_LENGTH];
  char szDatFile[MAX_WORD_LENGTH];
  char szListFile[MAX_WORD_LENGTH];
  FILE* tfp = NULL, *dfp = NULL, *lfp = NULL;

  getCommandLineParams(&tParams, argc, argv);

  readData(&tData, tParams.szDatFile);

  readMapFile(&tMap, tParams.szMapFile);

  ifp = fopen(tParams.szTreeFile, "r");

  if(ifp){
    addElement(&ptTree, ifp);

    fclose(ifp);
  }
  else{
    printf("Failed to open tree file\n");
  }
  
  setLeaves(ptTree);

  treeSplitEven(ptTree, tParams.nSplit, aptSplit, &nSplit);

  for(i = 0; i < nSplit; i++){

    countLeaves(aptSplit[i],&(aptSplit[i]->nN));

    if(aptSplit[i]->nN < tParams.nMinSize){
      nLast += aptSplit[i]->nN;
    }

    aptSplit[i]->anLeaves = (int *) malloc(sizeof(int)*aptSplit[i]->nN);
    
    nCount = 0;
    
    getLeaves(aptSplit[i],aptSplit[i]->anLeaves, &nCount);
  }

  maxDepth(ptTree, &dMaxDepth);

  setDepth(ptTree, 0.0);

  /*sort on number of leaves*/
  //void qsort(void* field, size_t nElements, size_t sizeOfAnElement,
  //                 int(_USERENTRY *cmpFunc)(const void*, const void*));


  qsort(aptSplit,nSplit,sizeof(t_Node*),compNode);

  i = 0;
  while(i < nSplit && aptSplit[i]->nN >= tParams.nMinSize){

    sprintf(szDir, "C%03d",i);
    
    mkdir(szDir, S_IRWXU);

    sprintf(szTreeFile,"%s/%s%s",szDir,szDir,TREE_SUFFIX);
    sprintf(szDatFile,"%s/%s%s",szDir,szDir,DAT_SUFFIX);
    sprintf(szListFile,"%s/%s%s",szDir,szDir,LIST_SUFFIX);

    printf("%d %d %f\n",i,aptSplit[i]->nN,dMaxDepth - aptSplit[i]->dDepth);


    tfp = fopen(szTreeFile, "w");

    if(tfp){
      writeTree(aptSplit[i], tfp);
      fprintf(tfp, ";\n");
      fclose(tfp);
    }

    dfp = fopen(szDatFile, "w");

    if(dfp){
      writeData(dfp, &tData, aptSplit[i]->nN, aptSplit[i]->anLeaves, &tMap);
    
      fclose(dfp);
    }

    nCount=0;
    renumberLeaves(aptSplit[i], &nCount);

    lfp = fopen(szListFile, "w");

    if(lfp){
      writeList(lfp, aptSplit[i], dMaxDepth - aptSplit[i]->dDepth);
    
      fclose(lfp);
    }
    i++;
  }
  
  if(nLast > 0){
    anLast = (int *) malloc(sizeof(int)*nLast);
    nCount = 0;
    printf("%d %d\n",i,nLast);
    iL = i;
    for(; i < nSplit; i++){
      for(j = 0; j < aptSplit[i]->nN; j++){
	anLast[nCount + j] = aptSplit[i]->anLeaves[j];
      }
      nCount += aptSplit[i]->nN;
    }

    if(nCount > 0){
      sprintf(szDir, "C%03d+",iL);
  
      mkdir(szDir, S_IRWXU);

      sprintf(szTreeFile,"%s/%s%s",szDir,szDir,TREE_SUFFIX);
      sprintf(szDatFile,"%s/%s%s",szDir,szDir,DAT_SUFFIX);
      sprintf(szListFile,"%s/%s%s",szDir,szDir,LIST_SUFFIX);

      dfp = fopen(szDatFile, "w");

      if(dfp){
	writeData(dfp, &tData, nLast, anLast, &tMap);
    
	fclose(dfp);
      }
    }

    free(anLast);
  }
  exit(EXIT_SUCCESS);
}

char addElement(t_Node** pptTree, FILE* ifp)
{
  char    cNext;
  t_Node* ptTree = (t_Node *) malloc(sizeof(t_Node));

  *pptTree = ptTree;
  ptTree->nId = INTERNAL;
  ptTree->nLeaves = -1;
  ptTree->dDepth = -1.0;

  cNext = fgetc(ifp);

  if(cNext != '('){
    char szIndex[MAX_WORD_LENGTH];
    char *pcError = NULL;
    int  nCount = 0;

    while(cNext != ':'){
      szIndex[nCount] = cNext;
      cNext = fgetc(ifp);
      nCount++;
    }

    szIndex[nCount] = '\0';
    ptTree->ptLeft  = NULL; ptTree->dLeft = 0.0;
    ptTree->nId  = strtol(szIndex, &pcError,10);
    ptTree->ptRight = NULL; ptTree->dRight = 0.0;
  }
  else{

    cNext = addElement(&ptTree->ptLeft, ifp);

    //cNext = fgetc(ifp);
    if(cNext == ':'){
      char szDist[MAX_WORD_LENGTH];
      char *pcError = NULL;
      int  nCount = 0;
      cNext = fgetc(ifp);
      while(cNext != ','){
	szDist[nCount] = cNext;
	cNext = fgetc(ifp);
	nCount++;
      }
      
      szDist[nCount] = '\0';
      ptTree->dLeft = strtod(szDist,&pcError);
    }
    cNext = addElement(&ptTree->ptRight, ifp);
    if(cNext == ':'){
      char szDist[MAX_WORD_LENGTH];
      char *pcError = NULL;
      int  nCount = 0;

      cNext = fgetc(ifp);
      while(cNext != ')'){
	szDist[nCount] = cNext;
	cNext = fgetc(ifp);
	nCount++;
      }
      
      szDist[nCount] = '\0';
      ptTree->dRight = strtod(szDist,&pcError);
    }

    cNext = fgetc(ifp);
  }

  return cNext;
}

void writeTree(t_Node *ptTree, FILE* ofp)
{
  if(ptTree->nId != INTERNAL){
    fprintf(ofp,"%d", ptTree->nId);
  }
  else{
    fprintf(ofp,"(");

    writeTree(ptTree->ptLeft, ofp);

    fprintf(ofp,":%.3f,",ptTree->dLeft);

    writeTree(ptTree->ptRight, ofp);

    fprintf(ofp,":%.3f)",ptTree->dRight);
  }
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
  ptParams->szDatFile  = extractParameter(argc,argv,DAT_FILE,ALWAYS);  
  if(ptParams->szDatFile == NULL)
    goto error;

  /*get output filestub*/
  ptParams->szTreeFile  = extractParameter(argc,argv,TREE_FILE, ALWAYS);  
  if(ptParams->szTreeFile == NULL)
    goto error;

  ptParams->szMapFile  = extractParameter(argc,argv,MAP_FILE, ALWAYS);  
  if(ptParams->szMapFile == NULL)
    goto error;

  szTemp = extractParameter(argc,argv,SPLIT,ALWAYS);
  if(szTemp != NULL){
    ptParams->nSplit = strtol(szTemp,&pcError,10);
    if(*pcError != '\0'){
      goto error;
    }
  }
  else{
    goto error;
  }

  szTemp = extractParameter(argc,argv,MIN_SIZE,ALWAYS);
  if(szTemp != NULL){
    ptParams->nMinSize = strtol(szTemp,&pcError,10);
    if(*pcError != '\0'){
      goto error;
    }
  }
  else{
    goto error;
  }

  return;

 error:
  writeUsage(stdout);
  exit(EXIT_FAILURE);

}


void treeSplit(t_Node* ptTree, double* pdDepth, double dSplit, t_Node **aptSplit,int* pnSplit)
{
  double dOldDepth = *pdDepth;

  if(*pdDepth > dSplit || ptTree->nId != INTERNAL){
    aptSplit[*pnSplit] = ptTree;

    (*pnSplit)++;
    return;
  }
  else{

    (*pdDepth) = dOldDepth + ptTree->dLeft;
    treeSplit(ptTree->ptLeft, pdDepth, dSplit, aptSplit,pnSplit);

    (*pdDepth) = dOldDepth + ptTree->dRight;
    treeSplit(ptTree->ptRight, pdDepth, dSplit, aptSplit,pnSplit);
  }
}

void treeSplitEven(t_Node* ptTree, int nSplit, t_Node **aptSplit,int* pnSplit)
{

  if(ptTree->nLeaves <= nSplit){
    aptSplit[*pnSplit] = ptTree;

    (*pnSplit)++;
    return;
  }
  else{
    treeSplitEven(ptTree->ptLeft, nSplit, aptSplit,pnSplit);
    
    treeSplitEven(ptTree->ptRight, nSplit, aptSplit,pnSplit);
  }
  
}

double maxDepth(t_Node* ptTree, double* pdDepth)
{
  double dOldDepth = *pdDepth;

  if(ptTree->nId == INTERNAL){
    (*pdDepth) += ptTree->dLeft;
    maxDepth(ptTree->ptLeft, pdDepth);
  }

  return;
}

double setDepth(t_Node* ptTree, double dDepth)
{
  double dOldDepth = dDepth;

  ptTree->dDepth = dDepth;
  
  if(ptTree->nId == INTERNAL){
    setDepth(ptTree->ptLeft, dDepth + ptTree->dLeft);

    setDepth(ptTree->ptRight, dDepth + ptTree->dRight);
  }

  return;
}

void writeIndices(t_Node* ptTree, FILE* ofp)
{
  if(ptTree->nId == INTERNAL){
    
    writeIndices(ptTree->ptLeft, ofp);

    writeIndices(ptTree->ptRight, ofp);
  }
  else{
    fprintf(ofp,"%d,",ptTree->nId);
  }
}

void countLeaves(t_Node* ptTree, int *pnLeaves)
{
  if(ptTree->nId != INTERNAL){
    (*pnLeaves)++;
  } 
  else{
    countLeaves(ptTree->ptLeft, pnLeaves);

    countLeaves(ptTree->ptRight, pnLeaves);
  }
}

void setLeaves(t_Node* ptTree)
{
  if(ptTree->nId != INTERNAL){
    ptTree->nLeaves = 1;
  } 
  else{

    if(ptTree->ptLeft->nLeaves < 0){
      setLeaves(ptTree->ptLeft);
    }

    if(ptTree->ptRight->nLeaves < 0){
      setLeaves(ptTree->ptRight);
    }
    
    ptTree->nLeaves = ptTree->ptLeft->nLeaves + ptTree->ptRight->nLeaves;
    
  }
}


void getLeaves(t_Node* ptTree, int* anIDs, int *pnLeaves)
{
  if(ptTree->nId != INTERNAL){
    anIDs[*pnLeaves] = ptTree->nId;

    (*pnLeaves)++;
  } 
  else{
    getLeaves(ptTree->ptLeft, anIDs, pnLeaves);

    getLeaves(ptTree->ptRight, anIDs, pnLeaves);
  }
}

void getLeavesR(t_Node* ptTree, int* anIDs, int *pnLeaves)
{
  if(ptTree->nId != INTERNAL){
    anIDs[*pnLeaves] = ptTree->nR;

    (*pnLeaves)++;
  } 
  else{
    getLeavesR(ptTree->ptLeft, anIDs, pnLeaves);

    getLeavesR(ptTree->ptRight, anIDs, pnLeaves);
  }
}

void renumberLeaves(t_Node* ptTree, int* pnCount)
{
  if(ptTree->nId != INTERNAL){
    ptTree->nR = *pnCount;
    (*pnCount)++;
  } 
  else{
    renumberLeaves(ptTree->ptLeft, pnCount);

    renumberLeaves(ptTree->ptRight, pnCount);
  }
}

int compNode(const void *pvA, const void* pvB)
{
  t_Node* ptA = *((t_Node **) pvA);
  t_Node* ptB = *((t_Node **) pvB);

  return ptA->nN > ptB->nN ? -1 : +1;
}

void readData(t_Data *ptData, char* szDatFile)
{
  char szLine[MAX_LINE_LENGTH];
  char *szTok = NULL, *pcError = NULL, *szBrk = NULL;;
  FILE *ifp = NULL;
  int i = 0, nM = 0, nN = 0;

  ifp = fopen(szDatFile, "r");
  
  if(ifp){
    fgets(szLine, MAX_LINE_LENGTH, ifp);

    //szBrk = strpbrk(szLine, "\n");

    szTok = strtok(szLine, DELIM);

    nN = strtol(szTok, &pcError, 10);
    if(*pcError != '\0')
      goto fileFormatError;

    szTok = strtok(NULL, DELIM);

    nM = strtol(szTok, &pcError, 10);
    if(*pcError != '\n')
      goto fileFormatError;

    ptData->nN = nN; ptData->nM = nM;

    ptData->aszLines = (char **) malloc(sizeof(char *)*3*nN);

    for(i = 0; i < 3*nN; i++){
      
      fgets(szLine, MAX_LINE_LENGTH, ifp);

      ptData->aszLines[i] = strdup(szLine);
    }

    fclose(ifp);
  }
  else{
    printf("Failed to open %s\n",szDatFile);
  }
  
  return;

 fileFormatError:
  fprintf(stderr, "Incorrectly formatted input file in readMeans");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void writeData(FILE* ofp, t_Data *ptData, int nEntries, int *anEntries, t_Map *ptMap)
{
  int i = 0, j = 0, index = 0;
  int nTotal = 0;

  for(i = 0; i < nEntries; i++){
    index = anEntries[i];
    nTotal += ptMap->anNM[index];
  }

  fprintf(ofp, "%d %d\n",nTotal,ptData->nM);

  for(i = 0; i < nEntries; i++){
    int index = anEntries[i];

    for(j = 0; j < ptMap->anNM[index]; j++){
      int nI = ptMap->aanMap[index][j];

      //fprintf(ofp, "%d ", nI);
      fputs (ptData->aszLines[nI], ofp);
    }
  }
}




void writeList(FILE* ofp, t_Node *ptTree, double dMaxSplit)
{
  double dSplitDepth = 0.0, dSplit = 0.0, dDepth = 0.0, dMaxDepth = 0.0;
  t_Node* aptSplit[ptTree->nN];
  int    i = 0, j = 0, nCount = 0, nSplit = 0, nLeaves = 0; 
  int *anLeaves = NULL;

  maxDepth(ptTree, &dMaxDepth);

  while(dSplit < dMaxSplit){
 
    dSplitDepth = dMaxDepth - dSplit;

    nSplit = 0;
    dDepth = 0.0;

    treeSplit(ptTree, &dDepth, dSplitDepth, aptSplit, &nSplit);

    fprintf(ofp,"%f %d ",dSplit, nSplit);

    for(i = 0; i < nSplit; i++){
   
      nLeaves = 0;

      countLeaves(aptSplit[i],&nLeaves);
      
      if(anLeaves != NULL){
	free(anLeaves);
	anLeaves = NULL;
      }

      anLeaves = (int *) malloc(sizeof(int)*nLeaves);
    
      nCount = 0;
    
      getLeavesR(aptSplit[i], anLeaves, &nCount);
      
      for(j = 0; j < nLeaves - 1; j++){
	fprintf(ofp, "%d,",anLeaves[j]);
      }
      
      fprintf(ofp,"%d",anLeaves[nLeaves - 1]);
      if(i < nSplit - 1){
	fprintf(ofp, " ");
      }
    }
    fprintf(ofp, "\n");
    dSplit += 0.01;
  }
}

void readMapFile(t_Map *ptMap, char *szMapFile)
{
  int i = 0, j = 0;
  char *szLine  = (char *) malloc(BIG_LINE_LENGTH*sizeof(char));
  FILE *ifp     = NULL;
  char *pcError = NULL, *szTemp  = NULL, *szTok = NULL;


  ifp = fopen(szMapFile, "r");
  
  if(ifp){
    
    ptMap->nU = 0;
    while(fgets(szLine, BIG_LINE_LENGTH, ifp) != NULL){
      ptMap->nU++;
    }
    fclose(ifp);

    ifp = fopen(szMapFile, "r");

    ptMap->anNM  = (int *) malloc(ptMap->nU*sizeof(int));
    ptMap->aanMap = (int **) malloc(ptMap->nU*sizeof(int *));
    
    if(!ptMap->anNM || !ptMap->aanMap)
      goto memoryError;

    for(i = 0; i < ptMap->nU; i++){
      fgets(szLine, BIG_LINE_LENGTH, ifp);

      szTok = strtok(szLine, DELIM2);

      for(j = 0; j < 3; j++){
	szTok = strtok(NULL, DELIM2);
      }

      ptMap->anNM[i] = strtol(szTok,&pcError, 10);
     
      ptMap->aanMap[i] = (int *) malloc(ptMap->anNM[i]*sizeof(int));
      if(!ptMap->aanMap[i])
	goto memoryError;

      for(j = 0; j < ptMap->anNM[i]; j++){
	szTok = strtok(NULL, DELIM2);
	
	ptMap->aanMap[i][j] = strtol(szTok,&pcError, 10);
      }
      
    }


    fclose(ifp);
  }
  else{
    goto fileError;
  }
  
  free(szLine);

  return;

 fileError:
  fprintf(stderr, "Failed to open map file\n", szMapFile);
  exit(EXIT_FAILURE);
 memoryError:
  fprintf(stderr, "Failed to allocate memory\n");
  exit(EXIT_FAILURE);
}


//void readTree(t_Node **pptTree, FILE* ifp)
//{
//char cNext;

//cNext = fgetc(ifp);

  /*parse of until we start new element*/
  //while(cNext != '(' && cNext != EOF){
//cNext = fgetc(ifp);
//}
  
//if(cNext == EOF){
//  return;
//}

//addElement(pptTree, ifp);

//}
