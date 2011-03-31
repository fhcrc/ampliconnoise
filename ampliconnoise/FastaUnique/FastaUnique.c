#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FastaUnique.h"

static char *usage[] = {"FastaUnique - dereplicates fasta file\n",
			"-in     string            input file name\n",
			"Options:\n"};

static int nLines = 3;

static char szSequence[] = "acgturyndbwsmkhvACGTURYNDBWSMKHV-.";

static char szNoisy[]    = "Nn";

void readData(t_Data *ptData, t_Params *ptParams)
{
  FILE *ifp = NULL;
  char szLine[MAX_LINE_LENGTH];
  int  nPos = 0, i = 0, j = 0, nSequences = 0;
  char *szBrk;  
  char *szRet;

  /*first count sequences and get length*/
  ptData->nSeq    = 0;
  ptData->nMaxLen = 0;
  
  ifp = fopen(ptParams->szInputFile, "r");

  if(ifp){
    while(fgets(szLine, MAX_LINE_LENGTH, ifp)){
      if(szLine[0] == '>'){
	if(nPos > ptData->nMaxLen){
	  ptData->nMaxLen = nPos;
	}
	
	ptData->nSeq++;
	nPos = 0;
      }
      else{
	i = 0;
	while(strrchr(szSequence,szLine[i]) != NULL){
	  i++;
	  nPos++;
	}
      }
    }

    fclose(ifp);
  }
  else{
    fprintf(stderr, "Can't open input file %s\n", ptParams->szInputFile);
    exit(EXIT_FAILURE);
  }

  ptData->aszID        = (char **) malloc(ptData->nSeq*sizeof(char *));
  ptData->aacSequences = (char **) malloc(ptData->nSeq*sizeof(char *));
  ptData->anLen        = (int *)   malloc(ptData->nSeq*sizeof(int));

  ifp = fopen(ptParams->szInputFile, "r");

  if(ifp){
    while(szRet = fgets(szLine, MAX_LINE_LENGTH, ifp)){
      if(szLine[0] == '>'){
	if(nSequences > 0){
	  ptData->anLen[nSequences - 1] = nPos;
	}

	ptData->aacSequences[nSequences] = (char *) malloc(ptData->nMaxLen*sizeof(char));
	szBrk = strpbrk(szLine, " \n");
	(*szBrk) = '\0';
	ptData->aszID[nSequences] = strdup(szLine + 1);
	nPos = 0;
	nSequences++;
      }
    
      i = 0;
      while(szLine[i] != '\0' && strrchr(szSequence,szLine[i]) != NULL){
	ptData->aacSequences[nSequences - 1][nPos] = szLine[i];
	nPos++; i++;
      }
    }

    ptData->anLen[nSequences - 1] = nPos;

    fclose(ifp);
  }
  else{
    fprintf(stderr, "Can't open input file %s\n", ptParams->szInputFile);
    exit(EXIT_FAILURE);
  }
}

void writeSequence(FILE *ofp, t_Data *ptData, int nCopies, int nI)
{
  int nPos = 0;
  
  fprintf(ofp,">%s_%d\n",ptData->aszID[nI],nCopies);

  while(nPos < ptData->anLen[nI]){
    if(nPos == 80){
      fputc('\n',ofp);
    }

    fputc(ptData->aacSequences[nI][nPos],ofp);

    nPos++;
  }

  fputc('\n',ofp);
}


int main(int argc, char** argv){
  t_Params tParams;
  t_Data   tData;
  int      i = 0, j = 0, k = 0;
  double   **aadDMatrix = NULL;
  int      *anUnique    = NULL;
  t_Map    tMap;
  FILE     *ofp         = NULL;

  getCommandLineParams(&tParams, argc, argv);

  readData(&tData, &tParams);

  allocMap(&tMap, tData.nSeq);

  tParams.szMapFile = szCreateMapFileName(tParams.szInputFile);

  anUnique = (int *) malloc(sizeof(int)*tData.nSeq);

  for(i = 0; i < tData.nSeq; i++){
    anUnique[i] = TRUE;
  }
  
  
  for(i = 0; i < tData.nSeq; i++){
  
    if(anUnique[i] == TRUE){
      int nCopies = 1;
      
      addEntry(&tMap,i,i);

      for(j = i + 1; j < tData.nSeq; j++){
	if(anUnique[j] == TRUE){
	  int nComp =  tData.anLen[i] > tData.anLen[j] ? tData.anLen[j]: tData.anLen[i]; 

	  k = 0;
	  while(k < nComp){
	    if(tData.aacSequences[i][k] != tData.aacSequences[j][k]){
	      break;
	    }
	    k++;
	  }

	  if(k == nComp){
	    if(tData.anLen[i] > tData.anLen[j]){
	      addEntry(&tMap,i,j);
	      anUnique[j] = FALSE;
	    }
	    else{
	      anUnique[i] = FALSE;
	      for(k = 0; k < tMap.anN[i]; k++){
		addEntry(&tMap,j,tMap.aanMap[i][k]);
	      }
	      tMap.anN[i] = 0;
	      break;
	    }
	  }
	}
      } /*loop j*/
    }
  } /*loop i*/

  for(i = 0; i < tData.nSeq; i++){
    if(anUnique[i] == TRUE){
      writeSequence(stdout, &tData, tMap.anN[i], i);
    }
  }


  ofp = fopen(tParams.szMapFile, "w");

  if(ofp){
    writeMap(ofp,&tData,&tMap);
    
    fclose(ofp);
  }
  else{
    fprintf(stderr, "Failed to open %s for writing\n", tParams.szMapFile);
    fflush(stderr);
  }
  
  /*free up allocated memory*/
  for(i = 0; i < tData.nSeq; i++){
    free(tData.aacSequences[i]);
    free(tData.aszID[i]);
  }

  destroyMap(&tMap);
  free(anUnique);
  free(aadDMatrix);
  free(tData.anLen);
  free(tData.aacSequences);
  free(tData.aszID);
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
  char *cError = NULL;

  /*get parameter file name*/
  ptParams->szInputFile  = extractParameter(argc,argv, INPUT_FILE,ALWAYS);  
  if(ptParams->szInputFile == NULL)
    goto error;
  
  szTemp = extractParameter(argc,argv,MIN,OPTION);  
  if(szTemp != NULL){
    ptParams->nMin = strtol(szTemp, &cError, 10);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->nMin = -1;
  }

  ptParams->bTrim = FALSE;
  if(extractParameter(argc,argv, TRIM, OPTION)){  
    ptParams->bTrim = TRUE;
  }

  return;

 error:
  writeUsage(stdout);
  exit(EXIT_FAILURE);
}

void allocMap(t_Map *ptMap, int nSeq)
{
  int i = 0;

  ptMap->nSeq = nSeq;

  ptMap->anN = (int *) malloc(nSeq*sizeof(int));
  if(!ptMap->anN)
    goto memoryError;

  ptMap->anSize = (int *) malloc(nSeq*sizeof(int));
  if(!ptMap->anSize)
    goto memoryError;

  ptMap->aanMap = (int **) malloc(nSeq*sizeof(int*));
  if(!ptMap->aanMap)
    goto memoryError;

  for(i = 0; i < nSeq; i++){
    ptMap->anN[i]    = 0;
    ptMap->anSize[i] = INIT_MAP_SIZE;
    ptMap->aanMap[i] = (int *) malloc(INIT_MAP_SIZE*sizeof(int));
    if(!ptMap->aanMap[i])
      goto memoryError;
  }

  return;
 memoryError:
  fprintf(stderr, "Failed allocating memory in allocMap\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void addEntry(t_Map *ptMap, int nI, int nJ)
{
  ptMap->anN[nI]++;
  if(ptMap->anN[nI] > ptMap->anSize[nI]){
    ptMap->anSize[nI] *= 2;

    ptMap->aanMap[nI] = (int *) realloc(ptMap->aanMap[nI],ptMap->anSize[nI]*sizeof(int));
    if(!ptMap->aanMap[nI])
      goto memoryError;
  }

  ptMap->aanMap[nI][ptMap->anN[nI] - 1] = nJ;

  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in allocMap\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void writeMap(FILE *ofp, t_Data *ptData, t_Map *ptMap)
{
  int i = 0, j = 0;
  int nUnique = 0;

  for(i = 0; i < ptMap->nSeq; i++){
    if(ptMap->anN[i] > 0){
      fprintf(ofp, "%d,%d,%s,%d:", i, nUnique, ptData->aszID[i], ptMap->anN[i]);
      for(j = 0; j < ptMap->anN[i]; j++){
	fprintf(ofp, "%d",ptMap->aanMap[i][j]);
	if(j < ptMap->anN[i] - 1){
	  fprintf(ofp, ",");
	}
	else{
	  fprintf(ofp, ":");
	}
      }

      for(j = 0; j < ptMap->anN[i]; j++){
	fprintf(ofp, "%s",ptData->aszID[ptMap->aanMap[i][j]]);
	if(j < ptMap->anN[i] - 1){
	  fprintf(ofp, ",");
	}
      }
    
      nUnique++;
      fprintf(ofp, "\n");
    }
  }
}

void destroyMap(t_Map *ptMap)
{
  int i = 0;

  for(i = 0; i < ptMap->nSeq; i++){
    free(ptMap->aanMap[i]);
  }

  free(ptMap->aanMap);

  free(ptMap->anN);

  free(ptMap->anSize);

}

char* szCreateMapFileName(char *szFastaFile)
{
  char *szMapFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  char *szStub    = strdup(szFastaFile);
  char *szTemp    = NULL;

  szTemp = strstr(szStub, FASTA_SUFFIX);

  if(szTemp != NULL){
    *szTemp = '\0';
  }
  
  sprintf(szMapFile,"%s%s",szStub,MAP_SUFFIX);

  free(szStub);
  return szMapFile;
}
