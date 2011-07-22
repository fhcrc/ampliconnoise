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
#include <mpi.h>

#include "SeqDist.h"

/*global constants*/
static char *usage[] = {"SeqDist - pairwise distance matrix from a fasta file\n",
			"-in     string            fasta file name\n",
			"Options:\n",
			"-i output identifiers\n",
			"-rin    string            lookup file name\n"};

static int  nLines = 5;

static char szSequence[] = "ACGT";

static double* adLookUp = NULL;

void broadcastData(t_Data *ptData)
{
  int i = 0;

  MPI_Bcast((void *) &ptData->nSeq, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast((void *) &ptData->nMaxLen, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast((void *) ptData->anLen, ptData->nSeq, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast((void *) ptData->acSequences, ptData->nSeq*ptData->nMaxLen, MPI_CHAR, 0, MPI_COMM_WORLD);

}

void receiveData(t_Data *ptData)
{
 
  MPI_Bcast((void *) &ptData->nSeq, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast((void *) &ptData->nMaxLen, 1, MPI_INT, 0, MPI_COMM_WORLD);

  ptData->anLen = (int *) malloc(ptData->nSeq*sizeof(int));
  if(!ptData->anLen)
    goto memoryError;

  MPI_Bcast((void *) ptData->anLen, ptData->nSeq, MPI_INT, 0, MPI_COMM_WORLD);

  ptData->acSequences = (char *) malloc(ptData->nSeq*ptData->nMaxLen*sizeof(char));
  if(!ptData->acSequences)
    goto memoryError;

  MPI_Bcast((void *) ptData->acSequences, ptData->nSeq*ptData->nMaxLen, MPI_CHAR, 0, MPI_COMM_WORLD);
 
  return;
 memoryError:
  fprintf(stderr, "Failed allocating memory in receiveData\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}


int main(int argc, char* argv[]){
  int a = 0, i = 0, j = 0, nN = 0, nM = 0;
  t_Params tParams;
  t_Data   tData;
  double *adDists = NULL;
  int    numtasks, rank, rc; 
  int    offset, nA, nA0, nSize, nCount, nStart, nFinish, nTag = 1;
  int    nPackets = 0, nPacketSize = 0, nPacketCurr = 0;
  MPI_Status   status;

  fflush(stdout);

  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  /*get command line params*/
  getCommandLineParams(&tParams, argc, argv);
  
  fflush(stdout);

  if(rank == 0){ /*head node reads data*/
    
    //printf("%d read data\n",rank); fflush(stdout);

    initLookUp(&tParams);

    readData(&tData, &tParams);
 
    MPI_Bcast(adLookUp, N_BASES*N_BASES, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //printf("%d broadcast data\n",rank); fflush(stdout);

    /*sends it out to other nodes*/
    broadcastData(&tData);

    /*do my bit*/
    nN = tData.nSeq; nM = tData.nMaxLen;

    /*allocate memory for whole dist matrix*/
    nSize = (nN*(nN - 1))/2;
    
    nA = (int) (floor(nSize / numtasks));

    nA0 = nA + (nSize % numtasks);

    adDists = (double *) malloc(nA*sizeof(double)); 
    if(!adDists)
      goto memoryError;

    printf("%d\n",nN);
    /*output identifiers here*/
    if(tParams.bIdent == TRUE){
      for(i = 0; i < tData.nSeq; i++){
	printf("%s\n", tData.aszID[i]);
      }
    }
    nCount = 0;
    for(i = 0; i < nN; i++){
      for(j = 0; j < i; j++){
	double  dDist = needlemanWunsch(&tData.acSequences[i*nM], &tData.acSequences[j*nM],tData.anLen[i], tData.anLen[j], tData.nMaxLen);

	printf("%.8e\n",dDist);

	nCount++;
	if(nCount == nA0){
	  goto finish;
	}
      }
      fflush(stdout);
    }

  finish:
    nPackets = ((nA - 1)/MAX_PACKET_SIZE) + 1;

    for (i = 1; i < numtasks; i++){
      
      nPacketCurr = 0;
      for(j = 0; j < nPackets; j++){
       	
	if(j < nPackets - 1){
	  nPacketSize = MAX_PACKET_SIZE;
	}
	else{
	  nPacketSize = nA % MAX_PACKET_SIZE;
	}
	
	MPI_Recv((void *) (&adDists[nPacketCurr]), nPacketSize, 
		 MPI_DOUBLE, i, nTag, MPI_COMM_WORLD, &status);
	
	nPacketCurr += nPacketSize;
      }
      

      for(j = 0; j < nA; j++){	
	printf("%.8e\n",adDists[j]); 
      }
    }
    
  }
  else{
    /* receive data*/
    adLookUp = (double *) malloc(N_BASES*N_BASES*sizeof(double));
    if(!adLookUp)
      goto memoryError;

    MPI_Bcast(adLookUp, N_BASES*N_BASES, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    receiveData(&tData);

    nN = tData.nSeq; nM = tData.nMaxLen;
    
    nSize = (nN*(nN - 1))/2;
    
    nA = (int) floor(nSize / numtasks);

    nA0 = nA + (nSize % numtasks);

    adDists = (double *) malloc(nA*sizeof(double)); 
    if(!adDists)
      goto memoryError;

    nCount = 0; nStart = nA0 + (rank - 1)*nA; nFinish = nA0 + rank*nA;
    
    for(i = 0; i < nN; i++){
      for(j = 0; j < i; j++){

	if(nCount >= nStart){
	  double  dDist = needlemanWunsch(&tData.acSequences[i*nM], &tData.acSequences[j*nM],tData.anLen[i], tData.anLen[j], tData.nMaxLen);

	  adDists[nCount - nStart] = dDist;
	}	

  	nCount++;
	if(nCount == nFinish){
	    goto finish2;
	}
	
      }
    }
  
  finish2:
      
    nPackets = ((nA - 1)/MAX_PACKET_SIZE) + 1;
      
    nPacketCurr = 0;
    for(j = 0; j < nPackets; j++){
      if(j < nPackets - 1){
	nPacketSize = MAX_PACKET_SIZE;
      }
      else{
	nPacketSize = nA % MAX_PACKET_SIZE;
      }   

      //printf("Send %d %d %d\n",rank, nPacketCurr, nPacketSize);
      MPI_Send((void *) (&adDists[nPacketCurr]), nPacketSize, 
		 MPI_DOUBLE, 0, nTag, MPI_COMM_WORLD);
      nPacketCurr += nPacketSize;
    }

    if(nPacketCurr != nA){
      printf("Message error\n");
    }
  }
 
  /*free allocated memory*/
  free(adDists);
  free(tData.acSequences);
  free(tData.anLen);
  // Wait for everyone to stop   
  
  MPI_Barrier(MPI_COMM_WORLD);

  // Always use MPI_Finalize as the last instruction of the program

  MPI_Finalize();
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
  ptParams->szInputFile  = extractParameter(argc,argv, INPUT_FILE,ALWAYS);  
  if(ptParams->szInputFile == NULL)
    goto error;
  
  if(szTemp = extractParameter(argc,argv, LOOKUP_FILE_FLAG, OPTION)){
    ptParams->szLookUpFile = szTemp;
  }
  else{
    ptParams->szLookUpFile = getenv("SEQ_LOOKUP_FILE");
    if(ptParams->szLookUpFile == NULL){
        ptParams->szLookUpFile = LOOKUP_FILE;
    }  
  }


  /*identifiers*/
  szTemp = extractParameter(argc, argv, IDENT, OPTION);
  if(szTemp != NULL){
    ptParams->bIdent = TRUE;
  }
  else{
    ptParams->bIdent = FALSE;
  }

  szTemp = extractParameter(argc, argv, PHYLIP, OPTION);
  if(szTemp != NULL){
    ptParams->bPhylip = TRUE;
  }
  else{
    ptParams->bPhylip = FALSE;
  }

  return;

 error:
  writeUsage(stdout);
  exit(EXIT_FAILURE);
}

double dist(char cA, char cB)
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
  default:
    fprintf(stderr, "Non standard base %c\n", cA);
    return 0.0;
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
  default:
    fprintf(stderr, "Non standard base %c\n", cB);
    return 0.0;
  }

  return adLookUp[nA*N_BASES + nB];
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

/*int getHomopolymerLengthA(char cAStart, const char* acA, const char* acB, int **aanMoves, int nI, int nJ)
{
  int i = 0, j = 0, k = 0, nCount = 0;
  char cACurr = '\0', cBEnd = '\0';

  nCount = 0;

  i = nI; j = nJ;
  while(i > 0 && j > 0){
    
    if(aanMoves[i][j] == 0){

      cACurr = acA[i - 1];
      cBEnd = acB[j - 1];
      if(cAStart == cACurr){
	nCount++;
      }
      else{
	goto noMatch;
      }

      break;
    }
    else if(aanMoves[i][j] < 0){
      int nUp = -aanMoves[i][j];

      for(k = 0; k < nUp; k++){
	cACurr = acA[i - k];
	if(cACurr == cAStart){
	  nCount++;
	}
	else{
	  goto noMatch;
	}
	i--;
      }
    }
    else{
      goto noMatch;
    }
     
  }

  if(cBEnd != cAStart){
    nCount = 0;
  }

  return nCount;
  
 noMatch:
  return 0;
}*/

char getLastMatch(int nMove, int nI, int nJ, const char *acA, const char *acB)
{
  switch(nMove){
  case DIAG:
    if(acA[nI - 1] == acB[nJ - 1]){
      return acA[nI - 1];
    }
    else{
      return '\0';
    }
    break;
  case LEFT:
    return acB[nJ - 1];
    break;
  case UP:
    return acA[nI - 1];
    break;
  }
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


double needlemanWunsch(const char* acA, const char* acB, int nLenA, int nLenB, int nM)
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
    aadFMatrix[i][0] = GAP_PENALTY*i;
    aanMoves[i][0] = UP;
    aanHA[i][0] = 0;
    aanHB[i][0] = 0;
  }
  for(j = 0; j <= nLenB; j++){ 
    aadFMatrix[0][j] = GAP_PENALTY*j;
    aanMoves[0][j] = LEFT;
    aanHA[0][j] = 0;
    aanHB[0][j] = 0;
  }

  for(i = 1; i <= nLenA; i++){
    for(j = 1; j <= nLenB; j++){
      double dChoice1, dChoice2, dChoice3;
     

      dChoice1 = aadFMatrix[i-1][j-1] + dist(acA[i - 1], acB[j - 1]);
      
      if(i == nLenA){
	dChoice2 = aadFMatrix[i][j-1];
      } 
      else{
	double dGap = 0.0;
	int    nCurrH = aanHA[i][j - 1];
	int    nNewH  = returnHomopolymerA(LEFT, aanHA, aanMoves, i, j, acA, acB);
	/*Left gap in A*/
	if(nNewH == 0){
	  dGap = GAP_PENALTY;
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
	  dGap = GAP_PENALTY;
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



void readData(t_Data *ptData, t_Params *ptParams)
{
  FILE *ifp = NULL;
  char szLine[MAX_LINE_LENGTH];
  int  nPos = 0, i = 0, j = 0, nM = 0, nSequences = 0;
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
  if(nPos > ptData->nMaxLen){
    ptData->nMaxLen = nPos;
  }
  nM = ptData->nMaxLen;
  ptData->acSequences = (char *) malloc(ptData->nSeq*nM*sizeof(char));
  ptData->anLen       = (int *)  malloc(ptData->nSeq*sizeof(int));

  ifp = fopen(ptParams->szInputFile, "r");

  if(ifp){
    while(szRet = fgets(szLine, MAX_LINE_LENGTH, ifp)){
      if(szLine[0] == '>'){
	if(nSequences > 0){
	  ptData->anLen[nSequences - 1] = nPos;
	}

	szBrk = strpbrk(szLine, " \n");
	(*szBrk) = '\0';
	ptData->aszID[nSequences] = strdup(szLine + 1);
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
    fprintf(stderr, "Can't open input file %s\n", ptParams->szInputFile);
    exit(EXIT_FAILURE);
  }
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
