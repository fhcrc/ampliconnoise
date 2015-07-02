/*This software and documentation is copyright © 2009 by Christopher Quince.*/

/*Permission is granted for anyone to copy, use, or modify these programs and documents for purposes of research or education, provided this copyright notice is retained, and note is made of any changes that have been made.*/ 

/* These programs and documents are distributed without any warranty, express or implied. As the programs were written for research purposes only, they have not been tested to the degree that would be advisable in any important application. All use of these programs is entirely at the user's own risk.*/

/**System includes****/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>
#include "Perseus.h"

/*global constants*/
static char *usage[] = {"Perseus2 - slays monsters\n",
			"-sin     string            seq file name\n",
			"Options:\n",
			"-s       integer\n",
			"-tin     string            reference sequence file\n",
			"-a                         output alignments\n",
			"-d                         use imbalance\n",
			"-rin     string            lookup file name\n"};

static int  nLines = 8;

static char szSequence[] = "ACGTNacgtn-";

static double* adLookUp = NULL;



int main(int argc, char* argv[]){
  int a = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, nN = 0, nK = 0;
  int nNLen = 0, nKLen = 0;
  t_Params tParams;
  t_Data   tSeqData, tRefData;
  char *acTest = NULL;
  int    nTestMatch = -1, nTestLength = -1;
  double dTest = 0.0;

  /*get command line params*/
  getCommandLineParams(&tParams, argc, argv);
  
  /*read sequences to chimera check*/
  readData(tParams.szSeqInputFile, &tSeqData);

  /*read sequences to compare too*/
  readData(tParams.szRefInputFile, &tRefData);
 
  /*set parameters for sequence distances*/
  initLookUp(&tParams);

  /*number of sequences*/
  nN = tSeqData.nSeq; 

  /*number of reference sequences*/
  nK = tRefData.nSeq;
  
  /*max length of sequences*/
  nNLen = tSeqData.nMaxLen; 

  /*max length of references*/
  nKLen = tRefData.nMaxLen;

  printf("%d %d\n",nN, nK);

  for(i = 0; i < nN; i++){
    t_Align atAlign[nK];
    int nLenI = tSeqData.anLen[i];
    double adD[nLenI], adR[nLenI];
    int anD[nLenI], anR[nLenI];
    int anBestD[nLenI], anBestR[nLenI];
    double dBestD = BIG_DBL, dBestTriD = BIG_DBL, dBest  = BIG_DBL, dBestChi = 0.0, dBestTri = 0.0;
    int    nBestJ = -1;
    int    nBest = BIG_INT, nBestChi = BIG_INT, nBestTri = BIG_INT;
    int    nSplit = -1, nSplit1 = -1, nSplit2 = -1, nP1 = -1, nP2 = -1, nT1 = -1, nT2 = -1, nT3 = -1;
    int    anRestrict[nK];
    int    nCompare = 0;
    double dLoon = 0.0, dCIndex = 0.0;

    printf("%d %s ",i, tSeqData.aszID[i]);

    /*do pairwise alignments and get best hit for each sequence i*/
    nCompare = alignAll(i, nLenI, &nBest, &nBestJ, anRestrict, nK, &tSeqData, &tRefData, atAlign, &tParams);

    if(nCompare >= 2){

      dBestD = getBestChimeraD(nK, &tRefData, &nP1, &nP2, &nSplit, anRestrict, nLenI, atAlign, adD, adR, anBestD, anBestR);
      nBestChi = atAlign[nP1].anD[nSplit] + atAlign[nP2].anR[nLenI - nSplit - 2];
      if(nBestChi >= 3 && nCompare >= 3){
	//setDR(nK, &tRefData, anRestrict, nLenI, atAlign, anD, anR, anBestD, anBestR);
	dBestTriD = getBestTrimeraD(nK, &tRefData, &nT1, &nT2, &nT3, &nSplit1, &nSplit2, anRestrict, nLenI, atAlign, adD, adR, anBestD, anBestR);
	nBestTri = atAlign[nT1].anD[nSplit1] - atAlign[nT2].anD[nSplit1] + atAlign[nT2].anD[nSplit2] + atAlign[nT3].anR[nLenI - nSplit2 - 2];
	
      }
      
      dBestChi = ((double) nBestChi)/((double) nLenI);
      dBestTri = ((double) nBestTri)/((double) nLenI);
	
      dBest = needlemanWunschN(&tSeqData.acSequences[i*nNLen],&tRefData.acSequences[nBestJ*nKLen] , nLenI, tRefData.anLen[nBestJ], nKLen);
	
      printf("%d %d %s ",nBest, nBestJ, tRefData.aszID[nBestJ]);
      printf("%d %d %d %s %s ", nBestChi, nP1, nP2, tRefData.aszID[nP1], tRefData.aszID[nP2]);
	
      if(nBestChi - nBestTri >= 3){
	
	nTestMatch = TRIMERA;

	acTest =  getTrimera(&nTestLength, &atAlign[nT1], &atAlign[nT2], &atAlign[nT3],nSplit1, nSplit2,nLenI);
	  /*Trimera*/
      }
      else{
	/*Chimera*/
	nTestMatch = CHIMERA;

	acTest =  getChimera(&nTestLength, &atAlign[nP1], &atAlign[nP2], nSplit, nLenI);
      }
    	
      dTest = needlemanWunschN(&tSeqData.acSequences[i*nNLen], acTest, nLenI, nTestLength, nKLen);

      dCIndex = calcCIndex(i, nP1, nP2, acTest, nTestLength, &tRefData, &tSeqData);

      dLoon = calcLoonIndex(&tSeqData, &tRefData, i, nP1, nP2, &nSplit, &tParams);
      
      printf("%f %f %f %f ",dBest, dCIndex, dCIndex - dBest, dLoon);
      
      printf("%d %d %d ", nBestChi, nBestTri, nSplit);

      switch(nTestMatch){
      case GOOD:
	printf("Good\n");
	break;
      case CHIMERA:
	printf("Chimera\n");
	break;
      case TRIMERA:
	printf("Trimera\n");
	break;
      case QUAMERA:
	printf("Quamera\n");
	break;
      }

      free(acTest);
    }
    else{
      printf("0 0 Null 0 0 0 Null Null 0.0 0.0 0.0 0.0 0 0 0 Null\n");
    }

    for(j = 0; j < nK; j++){
      if(anRestrict[j] == FALSE){
	free(atAlign[j].acA);
	free(atAlign[j].acB);
	free(atAlign[j].anD);
	free(atAlign[j].anR);
	free(atAlign[j].adD);
	free(atAlign[j].adR);
	free(atAlign[j].anMapD);
	free(atAlign[j].anMapR);
      }
    }
    
  }

  /*free allocated memory*/
  free(adLookUp);

  destroyData(&tSeqData);
  destroyData(&tRefData);

  exit(EXIT_SUCCESS);

 memoryError:
  fprintf(stderr, "Failed allocating memory in main\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
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
  ptParams->szSeqInputFile  = extractParameter(argc,argv, SEQ_INPUT_FILE,ALWAYS);  
  if(ptParams->szSeqInputFile == NULL)
    goto error;
  
    /*get parameter file name*/
  ptParams->szRefInputFile  = extractParameter(argc,argv, REF_INPUT_FILE,OPTION);  
  if(ptParams->szRefInputFile == NULL){
    ptParams->szRefInputFile = ptParams->szSeqInputFile;
  }

  szTemp  = extractParameter(argc,argv,SKEW,OPTION);  
  if(szTemp != NULL){
    ptParams->nSkew = strtol(szTemp,&cError,10);
    if(*cError!='\0'){
      goto error;
    }
  }
  else{
    ptParams->nSkew = DEFAULT_SKEW;
  }

  if(extractParameter(argc,argv,USE_IMBALANCE,OPTION)){  
    ptParams->bImbalance = TRUE;
  }
  else{
    ptParams->bImbalance = FALSE;
  }

  if(extractParameter(argc,argv,OUTPUT_ALIGNMENTS,OPTION)){  
    ptParams->bOutputAlignments = TRUE;
  }
  else{
    ptParams->bOutputAlignments = FALSE;
  }

  if(szTemp = extractParameter(argc,argv, LOOKUP_FILE_FLAG, OPTION)){
    ptParams->szLookUpFile = szTemp;
  }
  else{
    ptParams->szLookUpFile = getenv("SEQ_LOOKUP_FILE");
    if(ptParams->szLookUpFile == NULL){
        ptParams->szLookUpFile = LOOKUP_FILE;
    }  
  }

  return;

 error:
  writeUsage(stdout);
  exit(EXIT_FAILURE);
}

double dist(char cA, char cB)
{
  if(cA == cB){
    return 0;
  }
  else{
    return 1;
  }
}

double dmin3(double dA, double dB, double dC)
{
  double dAB = dA < dB ? dA : dB;
  
  return dC < dAB? dC : dAB;
}

int getMove(double dA, double dB, double dC)
{

  if(dA < dB){
    if(dA < dC){
      return DIAG;
    }
    else{
      return UP;
    }
  }
  else{
    if(dB < dC){
      return LEFT;
    }
    else{
      return UP;
    }
  }
}

double calcDistance(char* acA, char* acB, int nLen)
{
  int gap_alert1 = 0, gap_alert2  = 0;
  int start_seq1 = 0;
  int start_seq2 = 0;
  int end_seq1 = nLen;
  int end_seq2 = nLen;
  double residuecount = 0.0, distance = 0.0;
  int k = 0;
  int start_seq = 0, end_seq = 0;
 
  /* get left side gaps */
  k = 0;
  while(acA[k] == GAP && k < nLen) k++;
  start_seq1 = k;
	
  k = 0;
  while(acB[k] == GAP && k < nLen) k++;
  start_seq2 = k;
	
  /* get right side gaps */
  k = nLen - 1;
  while(acA[k] == GAP && k >= 0) k--;
  end_seq1 = k + 1; /* add one on for < not <= further on down */
	
  k = nLen - 1;
  while(acB[k] == GAP && k >= 0) k--;
  end_seq2 = k + 1;

  start_seq = start_seq1;
  if(start_seq < start_seq2){
    start_seq = start_seq2;
  }

  end_seq = end_seq1;
  if(end_seq > end_seq2){
    end_seq = end_seq2;
  }

  gap_alert1 = 0;
  gap_alert2 = 0;

  if(end_seq - start_seq <= 0){
    return 1.0;
  }

  for(k = start_seq; k < end_seq; k++){

    if (acA[k] == GAP && acB[k] == GAP){
      continue;
    }
	  
    if(acA[k] == GAP){
      if(gap_alert1 == 0){
	residuecount+=1.0;
	distance+=1.0;
      }
      gap_alert1 = 1;
      gap_alert2 = 0;
      
      continue;
    }
	
    if(acB[k] == GAP){
      if(gap_alert2 == 0){
	residuecount+=1.0;
	distance+=1.0;
      }
      gap_alert1 = 0;
      gap_alert2 = 1;
	    
      continue;
    }

    /* ok no gaps , check for mismatch */
    if(acA[k] != 'N' && acB[k] != 'N'){

      if (acA[k] != acB[k]){
	distance += 1.0;
	residuecount+=1.0;
	gap_alert1 = 0;
	gap_alert2 = 0;
		      
	continue;
      }

      if(acA[k] == acB[k]){
	residuecount+= 1.0;
	gap_alert1 = 0;
	gap_alert2 = 0;
	
	continue;
      }
    }
    else{
      residuecount+= 1.0;
	    
      gap_alert1 = 0;
      gap_alert2 = 0;

      continue;
    }
	
    printf("FELL THROUGH \n");
  } 

  if (residuecount > 0) {
    distance = distance / residuecount;
  }
  else {
    distance = 1.0;
  }

  return distance;
}

void needlemanWunsch(t_Align *ptAlign, const char* acA, const char* acB, int nLenA, int nLenB)
{
  double **aadFMatrix = NULL;
  int    **aanMoves   = NULL;
  char    *acAlignA   = NULL, *acAlignB = NULL;
  int    nCount = 0;
  int    i = 0, j = 0, nD = 0, nR = 0; 
  double dR = 0.0, dD = 0.0;
  
  aadFMatrix = (double **) malloc((nLenA + 1)*sizeof(double *));
  aanMoves   = (int **)    malloc((nLenA + 1)*sizeof(int *));
  if(!aadFMatrix || !aanMoves)
    goto memoryError;

  for(i = 0; i < nLenA + 1; i++){
    aadFMatrix[i] = (double *)  malloc((nLenB + 1)*sizeof(double));
    aanMoves[i]   = (int    *)  malloc((nLenB + 1)*sizeof(int));
    if(!aadFMatrix[i] || !aanMoves[i])
      goto memoryError;

    for(j = 0; j < nLenB + 1; j++){
      aadFMatrix[i][j] = 0.0;
      aanMoves[i][j]   = -1;
    }
  }

  for(i = 0; i <= nLenA; i++){
    aadFMatrix[i][0] = GAP_PENALTY*i;
    aanMoves[i][0] = UP;
  }
  for(j = 0; j <= nLenB; j++){ 
    aadFMatrix[0][j] = GAP_PENALTY*j;
    aanMoves[0][j] = LEFT;
  }

  for(i = 1; i <= nLenA; i++){
    for(j = 1; j <= nLenB; j++){
      double dChoice1, dChoice2, dChoice3;
      
      dChoice1 = aadFMatrix[i-1][j-1] + dist(acA[i - 1], acB[j - 1]);
      
      if(i == nLenA){
	dChoice2 = aadFMatrix[i][j-1];
      } 
      else{
	dChoice2 = aadFMatrix[i][j-1] + GAP_PENALTY;
      }

      if(j == nLenB){
	dChoice3 = aadFMatrix[i-1][j];
      }
      else{
	dChoice3 = aadFMatrix[i-1][j] + GAP_PENALTY;
      }
      
      aanMoves[i][j]   = getMove(dChoice1, dChoice2, dChoice3);
      aadFMatrix[i][j] = dmin3(dChoice1, dChoice2, dChoice3);
    }
  }

  ptAlign->dDist = aadFMatrix[nLenA][nLenB];
  
  acAlignA = (char *) malloc((nLenA + nLenB)*sizeof(char));
  if(!acAlignA)
    goto memoryError;

  acAlignB = (char *) malloc((nLenA + nLenB)*sizeof(char));
  if(!acAlignB)
    goto memoryError;

  for(i = 0; i < nLenA + nLenB; i++){
    acAlignA[i] = GAP; acAlignB[i] = GAP;
  }

  nCount = 0;
  i = nLenA; j = nLenB;
  while(i > 0 && j > 0){
    
    switch(aanMoves[i][j]){
    case DIAG:
     
      acAlignA[nCount] = acA[i - 1];
      acAlignB[nCount] = acB[j - 1];
      
      i--;
      j--;
      break;
    case UP:
      acAlignA[nCount] = acA[i - 1];
   
      if(j == nLenB){
	acAlignB[nCount] = T_GAP;
      }
      else{
	acAlignB[nCount] = GAP;
      }

      i--;
      break;
    case LEFT:
      if(i == nLenA){
	acAlignA[nCount] = T_GAP;
      }
      else{
	acAlignA[nCount] = GAP;
      }

      acAlignB[nCount] = acB[j - 1];

      j--;
      break;
    }
    
    nCount++;
  }
  
  while(i > 0){
    acAlignA[nCount] = acA[i - 1];
    acAlignB[nCount] = GAP;
    i--;
    nCount++;
  }
    
  while(j > 0){
    acAlignA[nCount] = GAP;
    acAlignB[nCount] = acB[j - 1];
    j--;
    nCount++;
  }

  ptAlign->nLen  = nCount;
  
  ptAlign->acA = (char *) malloc(ptAlign->nLen*sizeof(char));
  if(!ptAlign->acA)
    goto memoryError;

  ptAlign->acB = (char *) malloc(ptAlign->nLen*sizeof(char));
  if(!ptAlign->acB)
    goto memoryError;

  for(i = 0; i < ptAlign->nLen; i++){
    ptAlign->acA[ptAlign->nLen - 1 -i] = acAlignA[i];
    ptAlign->acB[ptAlign->nLen - 1 -i] = acAlignB[i];
  }

  ptAlign->dDist = calcDistance(acAlignA, acAlignB, ptAlign->nLen);

  ptAlign->nComp = ptAlign->nLen;

  while(ptAlign->nComp > 0 && ptAlign->acA[ptAlign->nComp - 1] == T_GAP || ptAlign->acB[ptAlign->nComp - 1] == T_GAP){
    ptAlign->nComp--;
  }

  ptAlign->anD = (int *) malloc(nLenA*sizeof(int));
  if(!ptAlign->anD)
    goto memoryError;

  ptAlign->anR = (int *) malloc(nLenA*sizeof(int));
  if(!ptAlign->anR)
    goto memoryError;

  ptAlign->adD = (double *) malloc(nLenA*sizeof(double));
  if(!ptAlign->adD)
    goto memoryError;

  ptAlign->adR = (double *) malloc(nLenA*sizeof(double));
  if(!ptAlign->adR)
    goto memoryError;


  ptAlign->anMapD = (int *) malloc(nLenA*sizeof(int));
  if(!ptAlign->anMapD)
    goto memoryError;

  ptAlign->anMapR = (int *) malloc(nLenA*sizeof(int));
  if(!ptAlign->anMapR)
    goto memoryError;

  nCount = 0; dD = 0; 
  for(i = 0; i < ptAlign->nLen; i++){
    if(ptAlign->acA[i] == GAP){
      dD+=GAP_PENALTY_N;
    }
    else if(ptAlign->acA[i] == T_GAP){

    }
    else{
      if(ptAlign->acB[i] == GAP){
	dD+=GAP_PENALTY_N;
      }
      else if(ptAlign->acB[i] == T_GAP){
	dD+=TERMINAL_PENALTY;
      }
      else{
	dD+=distN(ptAlign->acA[i],ptAlign->acB[i]);
      }
     
      ptAlign->adD[nCount] = dD;
      nCount++; 
    }
    
  }

  nCount = 0; nD = 0; 
  for(i = 0; i < ptAlign->nLen; i++){
    if(ptAlign->acA[i] == GAP){
      nD++;
    }
    else if(ptAlign->acA[i] == T_GAP){

    }
    else{
      if(ptAlign->acB[i] == GAP){
	nD++;
      }
      else if(ptAlign->acB[i] == T_GAP){

      }
      else if(ptAlign->acA[i] != ptAlign->acB[i]){
	nD++;
      }
     
      ptAlign->anD[nCount] = nD;
      ptAlign->anMapD[nCount] = i;
      nCount++; 
    }
    
  }

  ptAlign->nDiff = nD;

  ptAlign->dDiff = dD;

  nCount = 0; dR = 0.0; 
  for(i = ptAlign->nLen - 1; i >= 0; i--){
    
    if(ptAlign->acA[i] == GAP){
      dR+=GAP_PENALTY_N;
    }
    else if(ptAlign->acA[i] == T_GAP){

    }
    else{
      if(ptAlign->acB[i] == GAP){
	dR+=GAP_PENALTY_N;
      }
      else if(ptAlign->acB[i] == T_GAP){
	dR+=TERMINAL_PENALTY;
      }
      else{
	dR+=distN(ptAlign->acA[i],ptAlign->acB[i]);
      }
      
      ptAlign->adR[nCount] = dR;
      
      nCount++; 
    }
  }

  nCount = 0; nR = 0; 
  for(i = ptAlign->nLen - 1; i >= 0; i--){
    
    if(ptAlign->acA[i] == GAP){
      nR++;
    }
    else if(ptAlign->acA[i] == T_GAP){

    }
    else{
      if(ptAlign->acB[i] == GAP){
	nR++;
      }
      else if(ptAlign->acB[i] == T_GAP){
	//nR++;
      }
      else if(ptAlign->acB[i] != ptAlign->acA[i]){
	nR++;
      }
      
      ptAlign->anR[nCount] = nR;
      ptAlign->anMapR[nCount] = i;
      nCount++; 
    }
  }


  for(i = 0; i <= nLenA; i++){
    free(aadFMatrix[i]);
    free(aanMoves[i]);
  }
  free(aadFMatrix);
  free(aanMoves);
  free(acAlignA);
  free(acAlignB);
  return;

 memoryError:
  fprintf(stderr, "Failed to allocate memory in needlemanWunsch\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
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
  char   *szBreak  = revstr(szID, WEIGHTDELIM); 
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

static char szNoisy[]    = "N";

void readData(char* szInputFile, t_Data *ptData)
{
  FILE *ifp = NULL;
  char szLine[MAX_LINE_LENGTH];
  int  nPos = 0, i = 0, j = 0, nM = 0, nSequences = 0;
  char *szBrk;  
  char *szRet;

  /*first count sequences and get length*/
  ptData->nSeq    = 0;
  ptData->nMaxLen = 0;
  
  ifp = fopen(szInputFile, "r");

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
    fprintf(stderr, "Can't open input file %s\n", szInputFile);
    exit(EXIT_FAILURE);
  }

  if(nPos > ptData->nMaxLen){
    ptData->nMaxLen = nPos;
  }
  ptData->aszID        = (char **) malloc(ptData->nSeq*sizeof(char *));
  
  nM = ptData->nMaxLen;
  ptData->acSequences = (char *) malloc(ptData->nSeq*nM*sizeof(char));
  ptData->anLen       = (int *)  malloc(ptData->nSeq*sizeof(int));
  ptData->adFreq      = (double *) malloc(ptData->nSeq*sizeof(double));

  ifp = fopen(szInputFile, "r");

  if(ifp){
    while(szRet = fgets(szLine, MAX_LINE_LENGTH, ifp)){
      if(szLine[0] == '>'){
	if(nSequences > 0){
	  ptData->anLen[nSequences - 1] = nPos;
	}

	szBrk = strpbrk(szLine, " \n");
	(*szBrk) = '\0';
	ptData->aszID[nSequences] = strdup(szLine + 1);
	ptData->adFreq[nSequences] = getWeight(ptData->aszID[nSequences]);
	nPos = 0;
	nSequences++;
      }
    
      i = 0;
      while(szLine[i] != '\0' && strrchr(szSequence,szLine[i]) != NULL){
	ptData->acSequences[(nSequences - 1)*nM + nPos] = toupper(szLine[i]);
	
	nPos++; i++;
      }
    }

    ptData->anLen[nSequences - 1] = nPos;

    fclose(ifp);
  }
  else{
    fprintf(stderr, "Can't open input file %s\n", szInputFile);
    exit(EXIT_FAILURE);
  }
}

double distN(char cA, char cB)
{
  int nA = 0, nB = 0;

  switch(cA){
  case 'A':
    nA = 0;
    break;
  case 'C':
    nA = 1;
    break;
  case 'T':
    nA = 2;
    break;
  case 'G':
    nA = 3;
    break;
    //default:
    //fprintf(stderr, "Non standard base\n");
  }

  switch(cB){
  case 'A':
    nB = 0;
    break;
  case 'C':
    nB = 1;
    break;
  case 'T':
    nB = 2;
    break;
  case 'G':
    nB = 3;
    break;
    //default:
    //fprintf(stderr, "Non standard base\n");
  }

  return adLookUp[nA*N_BASES + nB];
}

char getLastMatch(int nMove, int nI, int nJ, const char *acA, const char *acB)
{
  switch(nMove){
  case DIAG:
    if((nI > 0 && nJ > 0) && acA[nI - 1] == acB[nJ - 1]){
      return acA[nI - 1];
    }
    else{
      return '\0';
    }
    break;
  case LEFT:
    if(nJ > 0){
      return acB[nJ - 1];
    }
    
    break;
  case UP:
    if(nI > 0){
      return acA[nI - 1];
    }
    break;
  }

  return '\0';
}

void updateHomopolymers(int** aanHA, int** aanHB, int** aanMoves, int nI, int nJ, const char* acA, const char *acB)
{
  int nMove = aanMoves[nI][nJ];

  switch(nMove){
  case DIAG:
    if(acA[nI - 1] == acB[nJ - 1]){
      if(getLastMatch(aanMoves[nI - 1][nJ - 1], nI - 1, nJ - 1, acA, acB) == acA[nI - 1]){
	aanHA[nI][nJ] = aanHA[nI - 1][nJ - 1] + 1;
	aanHB[nI][nJ] = aanHB[nI - 1][nJ - 1] + 1;
      }
      else{
	aanHA[nI][nJ] = 1;
	aanHB[nI][nJ] = 1;
      }
    }
    else{
      aanHA[nI][nJ] = 0;
      aanHB[nI][nJ] = 0;
    }

    break;
  case LEFT:
    if(acB[nJ - 1] == getLastMatch(aanMoves[nI][nJ - 1], nI, nJ - 1, acA, acB)){
      aanHB[nI][nJ] = aanHB[nI][nJ - 1] + 1;
      aanHA[nI][nJ] = aanHA[nI][nJ - 1] + 1;
    }
    else{
      aanHA[nI][nJ] = 0;
      aanHB[nI][nJ] = 0;
    }
    break;

  case UP:
    /*GAP in B*/
    if(acA[nI - 1] == getLastMatch(aanMoves[nI - 1][nJ], nI - 1, nJ, acA, acB)){
      aanHA[nI][nJ] = aanHA[nI - 1][nJ] + 1;
      aanHB[nI][nJ] = aanHB[nI - 1][nJ] + 1;
    }
    else{
      aanHA[nI][nJ] = 0;
      aanHB[nI][nJ] = 0;
    }
    break;
  }
}

int returnHomopolymerA(int nMove, int** aanHA, int** aanMoves, int nI, int nJ, const char* acA, const char *acB)
{
  int retA = 0;

  switch(nMove){
  case DIAG:
    if(acA[nI - 1] == acB[nJ - 1]){
      if(getLastMatch(aanMoves[nI - 1][nJ - 1], nI - 1, nJ - 1, acA, acB) == acA[nI - 1]){
	retA = aanHA[nI - 1][nJ - 1] + 1;
      }
      else{
	retA = 1;
      }
    }
    else{
      retA = 0;
    }

    break;
  case LEFT:
    if(acB[nJ - 1] == getLastMatch(aanMoves[nI][nJ - 1], nI, nJ - 1, acA, acB)){
      retA = aanHA[nI][nJ - 1] + 1;
    }
    else{
      retA = 0;
    }
    break;

  case UP:
    /*GAP in B*/
    if(acA[nI - 1] == getLastMatch(aanMoves[nI - 1][nJ], nI - 1, nJ, acA, acB)){
      retA = aanHA[nI - 1][nJ] + 1;
    }
    else{
      retA = 0;
    }
    break;
  }
  
  return retA;
}


double needlemanWunschN(const char* acA, const char* acB, int nLenA, int nLenB, int nM)
{
  double **aadFMatrix = NULL;
  int    **aanMoves   = NULL;
  int    **aanHA      = NULL;
  int    **aanHB      = NULL;
  int    anHA[nLenA + nLenB];
  int    anHB[nLenA + nLenB];
  char    *acAlignA   = NULL, *acAlignB = NULL;
  int    nCount = 0, nLen = 0, nComp = 0;
  int    i = 0, j = 0; 
  double dDist = 0.0;

  aadFMatrix = (double **) malloc((nLenA + 1)*sizeof(double *));
  aanMoves   = (int **)    malloc((nLenA + 1)*sizeof(int *));
  aanHA      = (int **)    malloc((nLenA + 1)*sizeof(int *));
  aanHB      = (int **)    malloc((nLenA + 1)*sizeof(int *));
  if(!aadFMatrix || !aanMoves)
    goto memoryError;

  for(i = 0; i < nLenA + 1; i++){
    aadFMatrix[i] = (double *)  malloc((nLenB + 1)*sizeof(double));
    aanMoves[i]   = (int    *)  malloc((nLenB + 1)*sizeof(int));
    aanHA[i]      = (int    *)  malloc((nLenB + 1)*sizeof(int));
    aanHB[i]      = (int    *)  malloc((nLenB + 1)*sizeof(int));
    if(!aadFMatrix[i] || !aanMoves[i])
      goto memoryError;

    for(j = 0; j < nLenB + 1; j++){
      aadFMatrix[i][j] = 0.0;
      aanMoves[i][j]   = -1;
      aanHA[i][j] = -1;
      aanHB[i][j] = -1;
    }
  }

  for(i = 0; i <= nLenA; i++){
    aadFMatrix[i][0] = GAP_PENALTY_N*i;
    aanMoves[i][0] = UP;
    aanHA[i][0] = 0;
    aanHB[i][0] = 0;
  }
  for(j = 0; j <= nLenB; j++){ 
    aadFMatrix[0][j] = GAP_PENALTY_N*j;
    aanMoves[0][j] = LEFT;
    aanHA[0][j] = 0;
    aanHB[0][j] = 0;
  }

  for(i = 1; i <= nLenA; i++){
    for(j = 1; j <= nLenB; j++){
      double dChoice1, dChoice2, dChoice3;
     

      dChoice1 = aadFMatrix[i-1][j-1] + distN(acA[i - 1], acB[j - 1]);
      
      if(i == nLenA){
	dChoice2 = aadFMatrix[i][j-1];
      } 
      else{
	double dGap = 0.0;
	int    nCurrH = aanHA[i][j - 1];
	int    nNewH  = returnHomopolymerA(LEFT, aanHA, aanMoves, i, j, acA, acB);
	/*Left gap in A*/
	if(nNewH == 0){
	  dGap = GAP_PENALTY_N;
	}
	else{
	  dGap = HOMOPOLYMER_PENALTY;
	}
	
	dChoice2 = aadFMatrix[i][j-1] + dGap;
      }

      if(j == nLenB){
	dChoice3 = aadFMatrix[i-1][j];
      }
      else{
	double dGap = 0.0;
	int    nCurrH = aanHA[i][j - 1];
	int    nNewH  = returnHomopolymerA(LEFT, aanHA, aanMoves, i, j, acA, acB);
	/*Left gap in A*/
	if(nNewH == 0){
	  dGap = GAP_PENALTY_N;
	}
	else{
	  dGap = HOMOPOLYMER_PENALTY;
	}
	/*Up gap in B*/
	dChoice3 = aadFMatrix[i-1][j] + dGap;
      }
      
      aanMoves[i][j]   = getMove(dChoice1, dChoice2, dChoice3);
      aadFMatrix[i][j] = dmin3(dChoice1, dChoice2, dChoice3);
      updateHomopolymers(aanHA, aanHB, aanMoves, i, j, acA, acB);
    }
  }

  dDist = aadFMatrix[nLenA][nLenB];
  
  acAlignA = (char *) malloc((nLenA + nLenB)*sizeof(char));
  if(!acAlignA)
    goto memoryError;

  acAlignB = (char *) malloc((nLenA + nLenB)*sizeof(char));
  if(!acAlignB)
    goto memoryError;

  for(i = 0; i < nLenA + nLenB; i++){
    acAlignA[i] = GAP; acAlignB[i] = GAP;
  }

  nCount = 0;
  i = nLenA; j = nLenB;
  while(i > 0 && j > 0){
    anHA[nCount] = aanHA[i][j];
    anHB[nCount] = aanHB[i][j];
    switch(aanMoves[i][j]){
    case DIAG:
     
      acAlignA[nCount] = acA[i - 1];
      acAlignB[nCount] = acB[j - 1];
      
      i--;
      j--;
      break;
    case UP:
      acAlignA[nCount] = acA[i - 1];
   
      if(j == nLenB){
	acAlignB[nCount] = T_GAP;
      }
      else{
	acAlignB[nCount] = GAP;
      }
      i--;
      break;
    case LEFT:
      if(i == nLenA){
	acAlignA[nCount] = T_GAP;
      }
      else{
	acAlignA[nCount] = GAP;
      }
      acAlignB[nCount] = acB[j - 1];

      j--;
      break;
    }
    
    nCount++;
  }
  
  while(i > 0){
    acAlignA[nCount] = acA[i - 1];
    acAlignB[nCount] = GAP;
    anHA[nCount] = aanHA[i][j];
    anHB[nCount] = aanHB[i][j];
    i--;
    nCount++;
  }
    
  while(j > 0){
    acAlignA[nCount] = GAP;
    acAlignB[nCount] = acB[j - 1];
    anHA[nCount] = aanHA[i][j];
    anHB[nCount] = aanHB[i][j];
    j--;
    nCount++;
  }

  nLen  = nCount;

  i = 0;
  while(acAlignA[i] == T_GAP || acAlignB[i] == T_GAP){
    i++;
  }
  
  nComp = nLen - i;
  
  dDist = aadFMatrix[nLenA][nLenB]/((double) nComp); //normalise by M * true length not with terminal gaps

  free(acAlignA); free(acAlignB);
  for(i = 0; i <= nLenA; i++){
    free(aadFMatrix[i]);
    free(aanMoves[i]);
    free(aanHA[i]);
    free(aanHB[i]);
  }
  free(aadFMatrix);
  free(aanMoves);
  free(aanHA);
  free(aanHB);
  return dDist;

 memoryError:
  fprintf(stderr, "Failed to allocate memory in needlemanWunsch\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void initLookUp(t_Params *ptParams)
{
  int i = 0, j = 0;
  FILE *ifp = NULL;

  adLookUp = (double *) malloc(N_BASES*N_BASES*sizeof(double));
  if(!adLookUp)
    goto memoryError;

  ifp = fopen(ptParams->szLookUpFile,"r");
  if(ifp){
    char szLine[MAX_LINE_LENGTH];
    char *pcError = NULL;
    char *szRet = NULL;
    char *szTok = NULL;

    for(i = 0; i < N_BASES; i++){
      fgets(szLine, MAX_LINE_LENGTH, ifp);
	
      szRet = strpbrk(szLine, "\n");

      (*szRet) = '\0';
      
      for(j = 0; j < N_BASES; j++){
	
	if(j == 0){
	  szTok = strtok(szLine, COMMA);
	}
	else{
	  szTok = strtok(NULL, COMMA);
	}

	adLookUp[i*N_BASES + j] = strtod(szTok,&pcError);
	if(*pcError != '\0')
	  goto formatError;
      }
    }
   
    fclose(ifp);
  }
  else{
    fprintf(stderr,"Failed to open %s\n",ptParams->szLookUpFile);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  return;
 memoryError:
  fprintf(stderr, "Failed allocating memory in initLookUp\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
 formatError:
  fprintf(stderr, "Format error LookUp.dat\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
} 





char* getChimera(int *pnCLength, t_Align* ptA, t_Align* ptB, int nSplit, int nLenI)
{
  int   nCLength = 0, i = 0;
  int   nCMaxLength = ptA->nLen + ptB->nLen; 
  char* acChimera   = (char *) malloc(sizeof(char)*nCMaxLength);

  for(i = 0; i < nCMaxLength; i++){
    acChimera[i] = '\0';
  }

  nCLength = 0;
  for(i = 0; i <= ptA->anMapD[nSplit]; i++){
    if(ptA->acB[i] != GAP && ptA->acB[i] != T_GAP){
      acChimera[nCLength] = ptA->acB[i];
      nCLength++;
    }
  }

  //printf("split = %d Left %d %d Right %d %d\n",nSplit,ptA->anMapD[nSplit],ptA->nLen,ptB->anMapD[nSplit],ptB->nLen); 

  for(i = ptB->anMapR[nLenI - nSplit - 2]; i < ptB->nLen; i++){
    if(ptB->acB[i] != GAP && ptB->acB[i] != T_GAP){
      acChimera[nCLength] = ptB->acB[i];
      nCLength++;
    }
  }


  //printf("%d %d %d %d\n", nCLength, nCMaxLength, ptA->nLen,ptB->nLen);

  (*pnCLength) = nCLength;
  return acChimera;
}

double calcCIndex(int nI, int nP1, int nP2, char *acChimera, int nCLength, t_Data* ptRefData, t_Data* ptSeqData)
{
  double dCIndex = 0.0; //log(ptRefData->adFreq[nP1]) + log(ptRefData->adFreq[nP2]);
  int    nM = ptSeqData->nMaxLen;

  dCIndex = needlemanWunschN(&(ptSeqData->acSequences[nI*nM]), acChimera, ptSeqData->anLen[nI], nCLength, nM);

  return dCIndex;
}

double calcEIndex(int nX1, int nY1, int nX2, int nY2, int nSplit, int nLenI)
{
  int i = 0, nM = 0, nN = 0, nL = 0, nH = 0, nO = 0;
  double dLN1 = 0.0, dLN2 = 0.0, dLD = 0.0;

  if(nX1 <= nX2){
    nM = nX1 - nY1;

    nH = nY1;

    nN = nX1;

    nL = nLenI - nSplit - 1;

    nO = nX1 + nY1;
  }
  else{
    nM = nX2 - nY2;

    nH = nY2;

    nN = nX2;

    nL = nSplit + 1;

    nO = nX2 + nY2;
  }

  if(nL - nO >= 0 && nH <= nL - nO){
    dLN1 = gsl_sf_lnchoose(nL - nO, nH); 
  }
  else{
    dLN1 = 0.0;
  }
  dLN2 = gsl_sf_lnchoose(nO, nM);

  if(nL >= nN){
    dLD = gsl_sf_lnchoose(nL, nN);
  }
  else{
    dLD = 0.0;
  }
  return dLD - dLN1 - dLN2;
}
      


int getDifferenceRight(t_Align* ptA, t_Align* ptB, int nSplit, int nLenI)
{
  int   nDiff = 0, i = 0, a = 0, b = 0;


  for(i = nSplit + 1; i < nLenI ;i++){
    a = ptA->anMapD[i];
    b = ptB->anMapD[i];

    if(ptA->acB[a] != T_GAP && ptB->acB[b] != T_GAP){

      if(ptA->acB[a] != ptB->acB[b]){
	nDiff++;
      }
    }
  }

  return nDiff;
}

void allocateMatricesD(int nLenI, double ***paadT, int ***paanBestT, double ***paadT2, int ***paanBestT2)
{
  double **aadT     = NULL;
  int **aanBestT = NULL;
  double **aadT2     = NULL;
  int **aanBestT2 = NULL;
  int k = 0, l = 0, m = 0;

  aadT     = (double **) malloc(nLenI*sizeof(double *));
  aanBestT = (int **) malloc(nLenI*sizeof(int *));

  for(k = 0; k < nLenI; k++){
    aadT[k]     = (double *) malloc(nLenI*sizeof(double));
    aanBestT[k] = (int *) malloc(nLenI*sizeof(int));
    
    for(l = 0; l < nLenI; l++){
      aadT[k][l] = 0.0;
      aanBestT[k][l] = 0;
    }
  }

  aadT2     = (double **) malloc(nLenI*sizeof(double *));
  aanBestT2 = (int **) malloc(nLenI*sizeof(int *));

  for(k = 0; k < nLenI; k++){
    aadT2[k]     = (double *) malloc(nLenI*sizeof(double));
    aanBestT2[k] = (int *) malloc(nLenI*sizeof(int));
    
    for(l = 0; l < nLenI; l++){
      aadT2[k][l] = 0.0;
      aanBestT2[k][l] = 0;
    }
  }

  *paadT = aadT;
  *paanBestT = aanBestT; 
  *paadT2 = aadT2;
  *paanBestT2 = aanBestT2;

  return;
}

void allocateMatrices(int nLenI, int ***paanT, int ***paanBestT, int ***paanT2, int ***paanBestT2)
{
  int **aanT     = NULL;
  int **aanBestT = NULL;
  int **aanT2     = NULL;
  int **aanBestT2 = NULL;
  int k = 0, l = 0, m = 0;

  aanT     = (int **) malloc(nLenI*sizeof(int *));
  aanBestT = (int **) malloc(nLenI*sizeof(int *));

  for(k = 0; k < nLenI; k++){
    aanT[k]     = (int *) malloc(nLenI*sizeof(int));
    aanBestT[k] = (int *) malloc(nLenI*sizeof(int));
    
    for(l = 0; l < nLenI; l++){
      aanT[k][l] = 0;
      aanBestT[k][l] = 0;
    }
  }

  aanT2     = (int **) malloc(nLenI*sizeof(int *));
  aanBestT2 = (int **) malloc(nLenI*sizeof(int *));

  for(k = 0; k < nLenI; k++){
    aanT2[k]     = (int *) malloc(nLenI*sizeof(int));
    aanBestT2[k] = (int *) malloc(nLenI*sizeof(int));
    
    for(l = 0; l < nLenI; l++){
      aanT2[k][l] = 0;
      aanBestT2[k][l] = 0;
    }
  }

  *paanT = aanT;
  *paanBestT = aanBestT; 
  *paanT2 = aanT2;
  *paanBestT2 = aanBestT2;

  return;
}

char* getTrimera(int *pnCLength, t_Align* ptA, t_Align* ptB, t_Align* ptC, int nSplit1, int nSplit2, int nLenI)
{
  int   nCLength = 0, i = 0; 
  int   nCMaxLength = ptA->nLen + ptB->nLen + ptC->nLen;
  char* acChimera   = (char *) malloc(sizeof(char)*nCMaxLength);

  for(i = 0; i < nCMaxLength; i++){
    acChimera[i] = '\0';
  }

  nCLength = 0;
  for(i = 0; i <= ptA->anMapD[nSplit1]; i++){
    if(ptA->acB[i] != GAP && ptA->acB[i] != T_GAP){
      acChimera[nCLength] = ptA->acB[i];
      nCLength++;
    }
  }
  
  for(i = ptB->anMapD[nSplit1] + 1; i <= ptB->anMapD[nSplit2]; i++){
    if(ptB->acB[i] != GAP && ptB->acB[i] != T_GAP){
      acChimera[nCLength] = ptB->acB[i];
      nCLength++;
    }
  }

  for(i = ptC->anMapR[nLenI - nSplit2 - 2]; i < ptC->nLen; i++){
    if(ptC->acB[i] != GAP && ptC->acB[i] != T_GAP){
      acChimera[nCLength] = ptC->acB[i];
      nCLength++;
    }
  }


  (*pnCLength) = nCLength;
  return acChimera;
}

char* getQuamera(int *pnCLength, t_Align* ptA, t_Align* ptB, t_Align* ptC, t_Align* ptD, int nSplit1, int nSplit2, int nSplit3, int nLenI)
{
  int   nCLength = 0, i = 0;
  int   nCMaxLength =  ptA->nLen + ptB->nLen + ptC->nLen + ptD->nLen;
  char* acChimera   = (char *) malloc(sizeof(char)*nCMaxLength);

  for(i = 0; i < nCMaxLength; i++){
    acChimera[i] = '\0';
  }

  nCLength = 0;
  for(i = 0; i <= ptA->anMapD[nSplit1]; i++){
    if(ptA->acB[i] != GAP && ptA->acB[i] != T_GAP){
      acChimera[nCLength] = ptA->acB[i];
      nCLength++;
    }
  }
  
  for(i = ptB->anMapD[nSplit1] + 1;i <= ptB->anMapD[nSplit2]; i++){
    if(ptB->acB[i] != GAP && ptB->acB[i] != T_GAP){
      acChimera[nCLength] = ptB->acB[i];
      nCLength++;
    }
  }

  for(i = ptC->anMapD[nSplit2] + 1;i <= ptC->anMapD[nSplit3]; i++){
    if(ptC->acB[i] != GAP && ptC->acB[i] != T_GAP){
      acChimera[nCLength] = ptC->acB[i];
      nCLength++;
    }
  }

  for(i = ptD->anMapR[nLenI - nSplit3 - 2]; i < ptD->nLen; i++){
    if(ptD->acB[i] != GAP && ptD->acB[i] != T_GAP){
      acChimera[nCLength] = ptD->acB[i];
      nCLength++;
    }
  }


  (*pnCLength) = nCLength;
  return acChimera;
}

void writeSequenceI(FILE* ofp, t_Data *ptData, int nI)
{
  int i = 0, nLenI = ptData->anLen[nI], nMaxLen = ptData->nMaxLen;

  fprintf(ofp, ">%s\n",ptData->aszID[nI]);
  for(i = 0; i < nLenI; i++){
    fprintf(ofp, "%c",ptData->acSequences[nI*nMaxLen + i]);
  }
  fprintf(ofp,"\n");
}

char getC(char cA, char cB, char cC){

  if(cB == cC){
    return cB;
  }

  return cA;
}

int calcOffset(char *acAlign, int nLen){
  int i = 0;

  while(acAlign[nLen - i - 1] == GAP){
    i++;
  }

  return i;
}

double calcLoonIndex(t_Data *ptSeqData, t_Data *ptRefData, int nI, int nP1, int nP2, int* pnSplit, t_Params *ptParams)
{
  int i = 0, nLenI = ptSeqData->anLen[nI], nLenP1 = ptRefData->anLen[nP1], nLenP2 = ptRefData->anLen[nP2];
  FILE *ofp = NULL;
  char szCommand[MAX_LINE_LENGTH];
  t_Data tAlign;
  int nDiff1 = 0, anDiff1[MAX_DIFF + 1];
  int nDiff2 = 0, anDiff2[MAX_DIFF + 1];
  int anSplits[2*MAX_DIFF], s = 0, s1 = 0, s2 = 0;
  int nSplit = -1, sMin = 0, nMinD = BIG_INT;
  int anD[2*MAX_DIFF];
  int nMaxLen = 0, nTLen = 0;
  char* acSequences = NULL;
  double dRet = 0;
  char szTempFasta[MAX_LINE_LENGTH], szTempAlign[MAX_LINE_LENGTH];
  int nTGap1 = 0, nTGap2 = 0, nTGap3 = 0,nMaxTGap = -1;

  /*create sequence filenames*/
  if(ptParams->bOutputAlignments == FALSE){
  	sprintf(szTempFasta, "Temp%s",FASTA_SUFFIX);
  	sprintf(szTempAlign, "Temp%s",ALIGN_SUFFIX);
  }
  else{
  	sprintf(szTempFasta, "Temp%d%s",nI,FASTA_SUFFIX);
  	sprintf(szTempAlign, "Temp%d%s",nI,ALIGN_SUFFIX);
  }

  /*write sequences*/
  ofp = fopen(szTempFasta, "w");
  if(ofp){
    writeSequenceI(ofp, ptSeqData, nI);
    writeSequenceI(ofp, ptRefData, nP1);
    writeSequenceI(ofp, ptRefData, nP2);
    fclose(ofp);
  }
  else{
    fprintf(stderr, "Failed to open %s for writing ... abort\n",szTempFasta);
    exit(EXIT_FAILURE);
  }
  /*run mafft remotely*/
  sprintf(szCommand,"mafft-linsi %s > %s 2> %s",szTempFasta,szTempAlign, TEMP_ERROR_FILE);

  system(szCommand);

  /*read in mafft output - three sequence alignment*/
  readData(szTempAlign, &tAlign);
 
  acSequences = tAlign.acSequences;  nMaxLen = tAlign.nMaxLen;

  /*alignment length*/
  nTLen = tAlign.nMaxLen;

  /*find largest terminal gap*/
  while(acSequences[nTLen - 1 - nTGap1] == GAP && nTLen - nTGap1> 1){
    nTGap1++;
  }

  while(acSequences[nMaxLen + nTLen - 1 - nTGap2] == GAP && nTLen - nTGap2> 1){
    nTGap2++;
  }

  while(acSequences[2*nMaxLen + nTLen - 1 - nTGap3] == GAP && nTLen - nTGap3> 1){
    nTGap3++;
  }

  nMaxTGap = nTGap1 > nTGap2 ? nTGap1 : nTGap2;		
  nMaxTGap = nTGap3 > nMaxTGap ? nTGap3 : nMaxTGap;
  /*remove from alignment*/
  nTLen -= nMaxTGap;	

  /*find all differences between chimera and parents and positions*/
  for(i = 0; i < nTLen; i++){
    if(acSequences[i] != acSequences[nMaxLen + i]){
      if(nDiff1 < MAX_DIFF){
	anDiff1[nDiff1] = i;
	nDiff1++;
      }
      else{
	fprintf(stderr,"Max diff reached in calcLoon\n");
      }
    }

    if(acSequences[i] != acSequences[2*nMaxLen + i]){
      if(nDiff2 < MAX_DIFF){
	anDiff2[nDiff2] = i;
	nDiff2++;
      }
      else{
	fprintf(stderr,"Max diff reached in calcLoon\n");
      }
    }
  }

  anDiff1[nDiff1] = nTLen;
  anDiff2[nDiff2] = nTLen;

  s = 0; s1 = 0; s2 = 0;
  
  anSplits[s] = -1;
  anD[s] = nDiff2;
  s++;
  /*loop differences to find optimal split point*/
  while(s1 < nDiff1 || s2 < nDiff2){
    if(anDiff1[s1] <= anDiff2[s2]){
      anD[s] = anD[s - 1] + 1;
      anSplits[s] = anDiff1[s1];

      s++;
      s1++;
    }
    else if(anDiff1[s1] > anDiff2[s2]){
      anD[s] = anD[s - 1] - 1;
      anSplits[s] = anDiff2[s2];

      s++;
      s2++;
    }
  
  }
  

  for(i = 0; i < s; i++){
    if(anD[i] < nMinD){
      nMinD = anD[i];
      sMin = i;
    }
  }
  
  if(sMin < s - 1){
    nSplit = (anSplits[sMin] + anSplits[sMin + 1])/2;
  }
  else{
    nSplit = nTLen - 1;
  }

  /*dummy if, as we always do this*/
  if(TRUE){
    int nA = -1, nB = -1;
    int nDLP1 = 0, nDRP1 = 0, nDLP2 = 0, nDRP2 = 0;
    int nX = 0, nY = 0, nZ = 0, nXZ = 0;
    int nXA = 0, nXB = 0, nYA = 0, nYB = 0;
    char cC = '\0';
    double pA = 0.0, pB = 0.0;
    double dP = 0.0;

    for(i = 0; i < nTLen; i++){
      char cI = acSequences[i], cP1 = acSequences[nMaxLen + i], cP2 = acSequences[2*nMaxLen + i];
      /*get consenus base*/
      cC = getC(cI, cP1, cP2);

      /*count number of differences between chimera and consensus*/
      if(cI != cC){
	nZ++;
      }
      /*count number of differences between chimera and parent 1*/
      if(cP1 != cC){
	/*count number to left*/
	if(i <= nSplit){
	  nDLP1++;
	}
	/*count number to right*/
	else{
	  nDRP1++;
	}
      }
      /*count number of differences between chimera and parent 2*/
      if(cP2 != cC){
	if(i <= nSplit){
	  nDLP2++;
	}
	else{
	  nDRP2++;
	}
      }
    }
    /*if Parent1 is left*/
    if(nDiff1 <= nDiff2){
      nA = nP1;
      nB = nP2;

      nX = nDLP1 + nDRP1;
      nXA = nDLP1; nXB = nDRP1;

      nY = nDLP2 + nDRP2;
      nYA = nDLP2; nYB = nDRP2;

      pA = ((double) nSplit + 1)/((double) nTLen); /*probability of change to left*/
      pB = ((double) nTLen - nSplit - 1)/((double) nTLen);/*prob. to right*/
    }
    else{
      /*if parent 1 is right*/
      nA = nP2;
      nB = nP1;

      nX = nDLP2 + nDRP2;
      nXA = nDRP2; nXB = nDLP2;

      nY = nDLP1 + nDRP1;
      nYA = nDRP1; nYB = nDLP1;

      pB = ((double) nSplit + 1)/((double) nTLen);
      pA = ((double) nTLen - nSplit - 1)/((double) nTLen);
    }

    nXZ = nX + nZ;
    
    dRet = 0.0;
    dP = 0.0;
    
    /*adds extra factor for tree imbalance - not generally used*/
    if(ptParams->bImbalance){
      if(nXZ > 0){
	for(i = nX; i <= nXZ; i++){
	  dP += gsl_ran_binomial_pdf (i, 0.5, nXZ);
	}

	dRet += -log(dP);
      }
    }

    /*contribution from right parent*/
    if(nY > 0){
      dP = 0.0;
      for(i = nYA; i <= nY; i++){
	dP += gsl_ran_binomial_pdf(i, pA, nY);
      }
      dRet += -log(dP);
    }

    /*contribution from left*/
    if(nX > 0){
      dP = 0.0;
      for(i = nXB; i <= nX; i++){
	dP += gsl_ran_binomial_pdf(i, pB, nX);
      }
      dRet += -log(dP);
    }
  }

  (*pnSplit) = nSplit;

  /*free up memory*/
  destroyData(&tAlign);
  return dRet;
} 

void destroyData(t_Data *ptData)
{
  int i = 0;

  if(ptData->nSeq > 0){
    free(ptData->acSequences);
    free(ptData->anLen);
    for(i = 0; i < ptData->nSeq; i++){
      free(ptData->aszID[i]);
    }
    free(ptData->aszID);
    free(ptData->adFreq);
  }
}

int alignAll(int nI, int nLenI, int *pnBest, int *pnBestJ, int *anRestrict, int nK, t_Data *ptSeqData, t_Data *ptRefData, t_Align* atAlign, t_Params *ptParams)
{
  int j = 0, nCompare = 0, nBest = BIG_INT, nBestJ = -1;
  int nNLen = ptSeqData->nMaxLen, nKLen = ptRefData->nMaxLen;

  for(j = 0; j < nK; j++){
    int     nDist = 0;

    if(strcmp(ptRefData->aszID[j],ptSeqData->aszID[nI]) != 0 && ptRefData->adFreq[j] >= ptParams->nSkew*ptSeqData->adFreq[nI]){

      needlemanWunsch(&atAlign[j], &ptSeqData->acSequences[nI*nNLen], &ptRefData->acSequences[j*nKLen], nLenI, ptRefData->anLen[j]);
      
      nDist = atAlign[j].nDiff;
	
      if(nDist < nBest){
	nBest  = nDist;
	nBestJ = j;
      }
      
      nCompare++;
      anRestrict[j] = FALSE;
    }
    else{
      anRestrict[j] = TRUE;
    }
  }

  (*pnBest) = nBest;
  (*pnBestJ) = nBestJ;

  return nCompare;
}

int getBestChimera(int nK, t_Data *ptRefData, int* pnP1, int* pnP2, int *pnSplit, int* anRestrict, int nLenI, t_Align* atAlign, int* anD, int* anR, int* anBestD, int* anBestR)
{
  int j = 0, k = 0;
  int nP1 = -1, nP2 = -1;
  int nSplit;
  int nBestChi = BIG_INT;

  for(k = 0; k < nLenI; k++){
    anD[k] = BIG_INT; anBestD[k] = -1;

    for(j = 0; j < nK; j++){
      if(anRestrict[j] == FALSE){

	if(atAlign[j].anD[k]  < anD[k] || atAlign[j].anD[k] == anD[k] && ptRefData->adFreq[j] > ptRefData->adFreq[anBestD[k]]){
	  anD[k]     = atAlign[j].anD[k];
	  anBestD[k] = j;
	}
	    
      }
    }
  }
      
  for(k = 0; k < nLenI; k++){
    anR[k] = BIG_INT; anBestR[k] = -1;

    for(j = 0; j < nK; j++){
      if(anRestrict[j] == FALSE){
	if(atAlign[j].anR[k] < anR[k] || atAlign[j].anR[k] == anR[k] && ptRefData->adFreq[j] > ptRefData->adFreq[anBestR[k]]){
	  anR[k]     = atAlign[j].anR[k];
	  anBestR[k] = j;
	}
      }
    }
  }
    
  for(k = 0; k < nLenI - 1; k++){
    int nChi = anD[k] + anR[nLenI - k - 2];
    if(nChi < nBestChi){
      nBestChi = nChi;
      nSplit = k;
      nP1 = anBestD[k];
      nP2 = anBestR[nLenI - k - 2];
    }
  }

  *pnSplit = nSplit;
  *pnP1 = nP1;
  *pnP2 = nP2;

  return nBestChi;
}

int setDR(int nK, t_Data *ptRefData, int* anRestrict, int nLenI, t_Align* atAlign, int* anD, int* anR, int* anBestD, int* anBestR)
{
  int j = 0, k = 0;
  int nP1 = -1, nP2 = -1;
  int nSplit;
  int nBestChi = BIG_INT;

  for(k = 0; k < nLenI; k++){
    anD[k] = BIG_INT; anBestD[k] = -1;

    for(j = 0; j < nK; j++){
      if(anRestrict[j] == FALSE){

	if(atAlign[j].anD[k]  < anD[k] || atAlign[j].anD[k] == anD[k] && ptRefData->adFreq[j] > ptRefData->adFreq[anBestD[k]]){
	  anD[k]     = atAlign[j].anD[k];
	  anBestD[k] = j;
	}
	    
      }
    }
  }
      
  for(k = 0; k < nLenI; k++){
    anR[k] = BIG_INT; anBestR[k] = -1;

    for(j = 0; j < nK; j++){
      if(anRestrict[j] == FALSE){
	if(atAlign[j].anR[k] < anR[k] || atAlign[j].anR[k] == anR[k] && ptRefData->adFreq[j] > ptRefData->adFreq[anBestR[k]]){
	  anR[k]     = atAlign[j].anR[k];
	  anBestR[k] = j;
	}
      }
    }
  }
    
  for(k = 0; k < nLenI - 1; k++){
    int nChi = anD[k] + anR[nLenI - k - 2];
    if(nChi < nBestChi){
      nBestChi = nChi;
      nSplit = k;
      nP1 = anBestD[k];
      nP2 = anBestR[nLenI - k - 2];
    }
  }

  return nBestChi;
}

double getBestChimeraD(int nK, t_Data *ptRefData, int* pnP1, int* pnP2, int *pnSplit, int* anRestrict, int nLenI, t_Align* atAlign, double* adD, double* adR, int* anBestD, int* anBestR)
{
  int j = 0, k = 0;
  int nP1 = -1, nP2 = -1;
  int nSplit;
  double dBestChi = BIG_DBL;
  int nBestChi = -1;

  for(k = 0; k < nLenI; k++){
    adD[k] = BIG_DBL; anBestD[k] = -1;

    for(j = 0; j < nK; j++){
      if(anRestrict[j] == FALSE){

	if(atAlign[j].adD[k]  < adD[k] || atAlign[j].adD[k] == adD[k] && ptRefData->adFreq[j] > ptRefData->adFreq[anBestD[k]]){
	  adD[k]     = atAlign[j].adD[k];
	  anBestD[k] = j;
	}
	    
      }
    }
  }
      
  for(k = 0; k < nLenI; k++){
    adR[k] = BIG_DBL; anBestR[k] = -1;

    for(j = 0; j < nK; j++){
      if(anRestrict[j] == FALSE){
	if(atAlign[j].adR[k] < adR[k] || atAlign[j].adR[k] == adR[k] && ptRefData->adFreq[j] > ptRefData->adFreq[anBestR[k]]){
	  adR[k]     = atAlign[j].adR[k];
	  anBestR[k] = j;
	}
      }
    }
  }
    
  for(k = 0; k < nLenI - 1; k++){
    double dChi = adD[k] + adR[nLenI - k - 2];
    if(dChi < dBestChi){
      dBestChi = dChi;
      nSplit = k;
      nP1 = anBestD[k];
      nP2 = anBestR[nLenI - k - 2];
    }
  }

  *pnSplit = nSplit;
  *pnP1 = nP1;
  *pnP2 = nP2;
  
  return dBestChi;
}

int getBestTrimera(int nK, t_Data *ptRefData, int* pnT1, int* pnT2, int* pnT3, int *pnSplit1, int *pnSplit2, int* anRestrict, int nLenI, t_Align* atAlign, int* anD, int* anR, int* anBestD, int *anBestR)
{
  int j = 0, k = 0, l = 0;
  int ** aanT = NULL, **aanBestT = NULL;
  int ** aanT2 = NULL, **aanBestT2 = NULL;
  int nBestTri = BIG_INT;
  int nT1 = -1, nT2 = -1, nT3 = -1, nSplit1 = -1, nSplit2 = -1;

  allocateMatrices(nLenI, &aanT,&aanBestT,&aanT2,&aanBestT2);

  for(k = 0; k < nLenI; k++){
    for(l = k; l < nLenI - 1; l++){
      aanT[k][l] = BIG_INT; aanBestT[k][l] = -1;
	
      for(j = 0; j < nK; j++){
	if(anRestrict[j] == FALSE){

	  int nX = atAlign[j].anD[l] - atAlign[j].anD[k];
		
	  if(nX < aanT[k][l] || nX == aanT[k][l] && ptRefData->adFreq[j] > ptRefData->adFreq[aanBestT[k][l]]){
	    aanT[k][l] = nX;
	    aanBestT[k][l] = j;
	  }
	  
	}
      }

      aanT[k][l] += anD[k] + anR[nLenI - l - 2];

      if(aanT[k][l] < nBestTri){
	nBestTri = aanT[k][l];
	nSplit1 = k;
	nSplit2 = l;
	nT1 = anBestD[k];
	nT2 = aanBestT[k][l];
	nT3 = anBestR[nLenI -l -2];
      }
    }
  }

  (*pnT1) = nT1;
  (*pnT2) = nT2;
  (*pnT3) = nT3;

  (*pnSplit1) = nSplit1;
  (*pnSplit2) = nSplit2;

  for(j = 0; j < nLenI; j++){
    free(aanT[j]);
    free(aanBestT[j]);
    free(aanT2[j]);
    free(aanBestT2[j]);
  }

  free(aanT);
  free(aanBestT);
  free(aanT2);
  free(aanBestT2);

  return nBestTri;
}

double getBestTrimeraD(int nK, t_Data *ptRefData, int* pnT1, int* pnT2, int* pnT3, int *pnSplit1, int *pnSplit2, int* anRestrict, int nLenI, t_Align* atAlign, double* adD, double* adR, int* anBestD, int *anBestR)
{
  int j = 0, k = 0, l = 0;
  double ** aadT = NULL;
  int **aanBestT = NULL;
  double **aadT2 = NULL;
  int  **aanBestT2 = NULL;
  double dBestTri = BIG_DBL;
  int nT1 = -1, nT2 = -1, nT3 = -1, nSplit1 = -1, nSplit2 = -1;

  allocateMatricesD(nLenI, &aadT,&aanBestT,&aadT2,&aanBestT2);

  for(k = 0; k < nLenI; k++){
    for(l = k; l < nLenI - 1; l++){
      aadT[k][l] = BIG_DBL; aanBestT[k][l] = -1;
	
      for(j = 0; j < nK; j++){
	if(anRestrict[j] == FALSE){
	  double dX = atAlign[j].adD[l] - atAlign[j].adD[k];
		
	  if(dX < aadT[k][l] || dX == aadT[k][l] && ptRefData->adFreq[j] > ptRefData->adFreq[aanBestT[k][l]]){
	    aadT[k][l] = dX;
	    aanBestT[k][l] = j;
	  }
	  
	}
      }

      aadT[k][l] += adD[k] + adR[nLenI - l - 2];

      if(aadT[k][l] < dBestTri){
	dBestTri = aadT[k][l];
	nSplit1 = k;
	nSplit2 = l;
	nT1 = anBestD[k];
	nT2 = aanBestT[k][l];
	nT3 = anBestR[nLenI -l -2];
      }
    }
  }

  (*pnT1) = nT1;
  (*pnT2) = nT2;
  (*pnT3) = nT3;

  (*pnSplit1) = nSplit1;
  (*pnSplit2) = nSplit2;

  for(j = 0; j < nLenI; j++){
    free(aadT[j]);
    free(aanBestT[j]);
    free(aadT2[j]);
    free(aanBestT2[j]);
  }

  free(aadT);
  free(aanBestT);
  free(aadT2);
  free(aanBestT2);

  return dBestTri;
}

/*void writeAlignment(int nLenI)
{
  int a = 0, a1 = 0, a2 = 0;

  for(a = 0; a < nLenI; a++){
    printf("%d %c ", a, tSeqData.acSequences[nNLen*i + a]);

    while((atAlign[nP1].acA[a1] == GAP || atAlign[nP1].acA[a1] == T_GAP) && a1 < atAlign[nP1].nLen){
      a1++;
    }

    printf("%c %d %d ", atAlign[nP1].acB[a1], atAlign[nP1].anD[a], atAlign[nP1].anR[nLenI - a - 2]);
    a1++;

    while((atAlign[nP2].acA[a2] == GAP || atAlign[nP2].acA[a2] == T_GAP) && a2 < atAlign[nP2].nLen ){
	    a2++;
    }

    printf("%c %d %d ", atAlign[nP2].acB[a2], atAlign[nP2].anD[a], atAlign[nP2].anR[nLenI - a - 2]);
    a2++;

    printf("\n");
  }
	
 }
 }*/
