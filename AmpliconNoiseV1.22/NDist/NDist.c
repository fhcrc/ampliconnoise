/**System includes****/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#include "NDist.h"

/*global constants*/
static char *usage[] = {"NDist - pairwise Needleman-Wunsch sequence distance matrix from a fasta file\n",
			"-in     string            fata file name\n",
			"Options:\n",
			"-i output identifiers\n"};

static int  nLines = 4;

static char szSequence[] = "acgturynACGTURYN-";

void broadcastData(t_Data *ptData)
{
  int i = 0;

  MPI_Bcast((void *) &ptData->nSeq, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast((void *) &ptData->nMaxLen, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast((void *) ptData->anLen, ptData->nSeq, MPI_INT, 0, MPI_COMM_WORLD);

  for(i = 0; i < ptData->nMaxLen;i++){
  	MPI_Bcast((void *) &ptData->acSequences[i*ptData->nSeq], ptData->nSeq, MPI_CHAR, 0, MPI_COMM_WORLD);
  }
}

void receiveData(t_Data *ptData)
{
  int i = 0; 
  MPI_Bcast((void *) &ptData->nSeq, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast((void *) &ptData->nMaxLen, 1, MPI_INT, 0, MPI_COMM_WORLD);

  ptData->anLen = (int *) malloc(ptData->nSeq*sizeof(int));
  if(!ptData->anLen)
    goto memoryError;

  MPI_Bcast((void *) ptData->anLen, ptData->nSeq, MPI_INT, 0, MPI_COMM_WORLD);

  ptData->acSequences = (char *) malloc(ptData->nSeq*ptData->nMaxLen*sizeof(char));
  if(!ptData->acSequences)
    goto memoryError;

  for(i = 0; i < ptData->nMaxLen; i++){
  	MPI_Bcast((void *) &ptData->acSequences[i*ptData->nSeq], ptData->nSeq, MPI_CHAR, 0, MPI_COMM_WORLD);
  }
  return;
 memoryError:
  fprintf(stderr, "Failed allocating memory in receiveData\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}


int main(int argc, char* argv[]){
  int a = 0, nN = 0, nM = 0;
  unsigned long i = 0, j = 0;
  t_Params tParams;
  t_Data   tData;
  double *adDists = NULL;
  int    numtasks, rank, rc; 
  unsigned long   ulNumTasks, ulOffset, ulA, ulA0, ulSize, ulCount, ulStart, ulFinish; 
  int    nTag = 1;
  unsigned long   ulPackets = 0, ulPacketSize = 0, ulPacketCurr = 0;
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
    
    //fprintf(stdout, "%d read data\n",rank); fflush(stdout);

    readData(&tData, &tParams);
 
    //fprintf(stdout, "%d broadcast data\n",rank); fflush(stdout);

    /*sends it out to other nodes*/
    broadcastData(&tData);
    
    //fprintf(stdout, "%d start calculations\n",rank); fflush(stdout);

    /*do my bit*/
    nN = tData.nSeq; nM = tData.nMaxLen;

    /*allocate memory for whole dist matrix*/
    ulSize = ((unsigned long) nN*( (unsigned long) nN - 1))/2;
     	
    ulNumTasks = (unsigned long) numtasks;
    ulA = (unsigned long) (floor(ulSize / ulNumTasks));

    ulA0 = ulA + (ulSize % ulNumTasks);

    //fprintf(stderr, "%d allocate memory %lu %lu\n",rank, ulSize, ulA); fflush(stderr);

    adDists = (double *) malloc(ulA*sizeof(double)); 
    if(!adDists)
      goto memoryError;

    printf("%d\n",nN);
    /*output identifiers here*/
    if(tParams.bIdent == TRUE){
      for(i = 0; i < tData.nSeq; i++){
	printf("%s\n", tData.aszID[i]);
      }
    }
    ulCount = 0;
    for(i = 0; i < nN; i++){
      for(j = 0; j < i; j++){
	t_Align tAlign;
	double  dDist = 0.0;

	needlemanWunsch(&tAlign, &tData.acSequences[i*nM], &tData.acSequences[j*nM], 
			tData.anLen[i], tData.anLen[j]);


	dDist = tAlign.dDist;
	printf("%.8e\n",dDist);

	ulCount++;
	if(ulCount == ulA0){
	  goto finish;
	}
      }
      fflush(stdout);
    }

  finish:
    ulPackets = ((ulA - 1)/MAX_PACKET_SIZE) + 1;

    for (i = 1; i < numtasks; i++){
      
      ulPacketCurr = 0;
      for(j = 0; j < ulPackets; j++){
       	
	if(j < ulPackets - 1){
	  ulPacketSize = MAX_PACKET_SIZE;
	}
	else{
	  ulPacketSize = ulA % MAX_PACKET_SIZE;
	}
	
	MPI_Recv((void *) (&adDists[ulPacketCurr]), ulPacketSize, 
		 MPI_DOUBLE, i, nTag, MPI_COMM_WORLD, &status);
	
	ulPacketCurr += ulPacketSize;
      }
      

      for(j = 0; j < ulA; j++){	
	printf("%.8e\n",adDists[j]); 
      }
    }
    
  }
  else{
    /* receive data*/

    //fprintf(stdout,"%d receive data\n",rank); fflush(stderr);
    
    receiveData(&tData);

    nN = tData.nSeq; nM = tData.nMaxLen;
    
    ulSize = ((unsigned long) nN*( (unsigned long) nN - 1))/2;
    ulNumTasks = (unsigned long) numtasks;
    ulA = (unsigned long) floor(ulSize / ulNumTasks);

    ulA0 = ulA + (ulSize % ulNumTasks);

    //fprintf(stdout,"%d allocate memory %lu %lu\n",rank, ulSize, ulA); fflush(stdout);
    adDists = (double *) malloc(ulA*sizeof(double)); 
    if(!adDists)
      goto memoryError;

    ulCount = 0; ulStart = ulA0 + (rank - 1)*ulA; ulFinish = ulA0 + rank*ulA;
    
    //fprintf(stdout,"%d begin calculations\n",rank); fflush(stdout);
    for(i = 0; i < nN; i++){
      for(j = 0; j < i; j++){
	t_Align tAlign;
	double  dDist = 0.0;

	if(ulCount >= ulStart){

	  needlemanWunsch(&tAlign, &tData.acSequences[i*nM], &tData.acSequences[j*nM], 
			  tData.anLen[i], tData.anLen[j]);

	  dDist = tAlign.dDist;

	  adDists[ulCount - ulStart] = dDist;
	}	

  	ulCount++;
	if(ulCount == ulFinish){
	    goto finish2;
	}
	
      }
    }
  
  finish2:
      
    ulPackets = ((ulA - 1)/MAX_PACKET_SIZE) + 1;
      
    ulPacketCurr = 0;
    for(j = 0; j < ulPackets; j++){
      if(j < ulPackets - 1){
	ulPacketSize = MAX_PACKET_SIZE;
      }
      else{
	ulPacketSize = ulA % MAX_PACKET_SIZE;
      }   

      //printf("Send %d %d %d\n",rank, nPacketCurr, nPacketSize);
      MPI_Send((void *) (&adDists[ulPacketCurr]), ulPacketSize, 
		 MPI_DOUBLE, 0, nTag, MPI_COMM_WORLD);
      ulPacketCurr += ulPacketSize;
    }

    if(ulPacketCurr != ulA){
      fprintf(stderr, "Message error\n");
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
  int    i = 0, j = 0; 
  
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
      acAlignB[nCount] = GAP;

      i--;
      break;
    case LEFT:
      acAlignA[nCount] = GAP;
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

  //for(i = 0; i < ptAlign->nLen; i++){
  //fprintf(stderr,"%c %c\n",acAlignA[i],acAlignB[i]);
  //}
  //fprintf(stderr,"\n\n");
  ptAlign->dDist = calcDistance(acAlignA, acAlignB, ptAlign->nLen);

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

static char szNoisy[]    = "N";

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
