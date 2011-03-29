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

#include "PyroDist.h"

/*global constants*/
static char *usage[] = {"PyroDist - pairwise distance matrix from flowgrams\n",
			"-in     string            flow file name\n",
			"-out    stub              out file stub\n",
			"Options:\n",
                        "-ni                       no index in dat file\n",
                        "-rin    string            lookup file name\n"};


static int  nLines = 6;

int main(int argc, char* argv[]){
  int a = 0, i = 0, j = 0, nN = 0, nM = 0;
  t_Params tParams;
  t_Flows tFlows;
  double *adDists = NULL;
  int    numtasks, rank, rc; 
  int    offset, nA, nA0, nSize, nCount, nStart, nFinish, nTag = 1;
  int    nPackets = 0, nPacketSize = 0, nPacketCurr = 0;
  MPI_Status   status;
  double *adLookUp = NULL, *adFookUp = NULL;

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
  adLookUp = (double *) malloc(sizeof(double)*MAX_S*BINS);
  adFookUp = (double *) malloc(sizeof(double)*BINS*BINS);

  if(!adLookUp || !adFookUp)
    goto memoryError;

  if(rank == 0){ /*head node reads data*/
    char   szOutFile[MAX_LINE_LENGTH];
    FILE*  ofp = NULL;

    sprintf(szOutFile,"%s%s",tParams.szOutFileStub,SUFFIX);    

    initLookUp(&tParams,adLookUp);

    initFookUp(adFookUp, adLookUp);

    printf("%d: Read data\n", rank);
    readData(adLookUp, tParams.szDataFile, &tFlows, &tParams);
    
    /*sends it out to other nodes*/
    printf("%d: Broadcast data\n", rank);
    MPI_Bcast((void *) adLookUp, MAX_S*BINS, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast((void *) adFookUp, BINS*BINS, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    printf("%d: Broadcast flows\n", rank);

    broadcastFlows(&tFlows);
    
    /*do my bit*/
    nN = tFlows.nN; nM = tFlows.nM;

    /*allocate memory for whole dist matrix*/
    nSize = (nN*(nN - 1))/2;
    
    nA = (int) (floor(nSize / numtasks));

    nA0 = nA + (nSize % numtasks);

    adDists = (double *) malloc(nA*sizeof(double)); 
    if(!adDists)
      goto memoryError;

    ofp = fopen(szOutFile, "w");
    if(!ofp)
      goto fileError;

    fprintf(ofp, "%d\n",nN);
    nCount = 0;
    for(i = 0; i < nN; i++){
      for(j = 0; j < i; j++){
       
	double  dDist = FDist(adFookUp, &(tFlows.adF[i*nM]),&(tFlows.adF[j*nM]), &(tFlows.anF[i*nM]), 
			 &(tFlows.anF[j*nM]), tFlows.anLengths[i], tFlows.anLengths[j]);

	fprintf(ofp,"%.8e\n",dDist); 
	nCount++;
	if(nCount == nA0){
	  goto finish;
	}
      }
      fflush(ofp);
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
	fprintf(ofp, "%.8e\n",adDists[j]); 
      }
    }
    
  }
  else{
    /* receive data*/

    MPI_Bcast((void *) adLookUp, MAX_S*BINS, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast((void *) adFookUp, BINS*BINS, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    printf("%d: Receive flows\n",rank); fflush(stdout);
    
    receiveFlows(&tFlows);

    printf("%d: Perform calculations\n",rank); fflush(stdout);

    /*allocate memory for my part of distance matrix*/
  /*allocate memory for whole dist matrix*/
    nN = tFlows.nN; nM = tFlows.nM;
    
    nSize = (nN*(nN - 1))/2;
    
    nA = (int) floor(nSize / numtasks);

    nA0 = nA + (nSize % numtasks);

    adDists = (double *) malloc(nA*sizeof(double)); 
    if(!adDists)
      goto memoryError;


    nCount = 0; nStart = nA0 + (rank - 1)*nA; nFinish = nA0 + rank*nA;
    
    for(i = 0; i < nN; i++){
      for(j = 0; j < i; j++){
	t_Align tAlign;
	double  dDist = 0.0;

	if(nCount >= nStart){

	  dDist = FDist(adFookUp, &(tFlows.adF[i*nM]),&(tFlows.adF[j*nM]), 
			   &(tFlows.anF[i*nM]), &(tFlows.anF[j*nM]), 
			   tFlows.anLengths[i], tFlows.anLengths[j]);
	  
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
  if(rank == 0){
    free(tFlows.adData);
  }
  free(tFlows.adF);
  free(tFlows.anF);
  free(tFlows.anLengths);
  free(adDists);
  free(adLookUp);
  free(adFookUp);
  // Wait for everyone to stop   
  
  MPI_Barrier(MPI_COMM_WORLD);

  // Always use MPI_Finalize as the last instruction of the program

  MPI_Finalize();
  exit(EXIT_SUCCESS);

 fileError:
  fprintf(stderr, "Failed opening output file in main\n");
  fflush(stderr);
  exit(EXIT_FAILURE);

 memoryError:
  fprintf(stderr, "Failed allocating memory in main\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void broadcastFlows(t_Flows* ptFlows)
{
  double* adF = ptFlows->adF;
  int   * anF = ptFlows->anF;
  int   * anLengths = ptFlows->anLengths;
  int i = 0;
  int nPackets = 0, nPacketSize = 0, nPacketCurr = 0;
  int nM = ptFlows->nM, nN = ptFlows->nN;
  int nSize = nM*nN;

  MPI_Bcast((void *) &ptFlows->nN,1,MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast((void *) &ptFlows->nM,1,MPI_INT, 0, MPI_COMM_WORLD);

 
  nPackets = ((nSize - 1)/MAX_PACKET_SIZE) + 1;
      
  printf("nN=%d nM=%d nSize=%d\n",nN,nM,nSize);fflush(stdout);
  nPacketCurr = 0;
  for(i = 0; i < nPackets; i++){
       	
    if(i < nPackets - 1){
      nPacketSize = MAX_PACKET_SIZE;
    }
    else{
      nPacketSize = nSize % MAX_PACKET_SIZE;
    }
	
    //   printf("Broadcast Total=%d Packet=%d Size=%d\n",nPackets,i,nPacketSize);fflush(stdout);

    MPI_Bcast((void *) (&adF[nPacketCurr]),nPacketSize,MPI_DOUBLE,0, MPI_COMM_WORLD);

    MPI_Bcast((void *) (&anF[nPacketCurr]),nPacketSize,MPI_INT,0, MPI_COMM_WORLD);
	
    nPacketCurr += nPacketSize;
  }

  nSize = ptFlows->nN;

  nPackets = ((nSize - 1)/MAX_PACKET_SIZE) + 1;
      
  nPacketCurr = 0;
  for(i = 0; i < nPackets; i++){
       	
    if(i < nPackets - 1){
      nPacketSize = MAX_PACKET_SIZE;
    }
    else{
      nPacketSize = nSize % MAX_PACKET_SIZE;
    }
       
    //printf("Broadcast Total=%d Packet=%d Size=%d\n",nPackets,i,nPacketSize);fflush(stdout);

    MPI_Bcast((void *) (&anLengths[nPacketCurr]),nPacketSize,MPI_INT,0, MPI_COMM_WORLD);
	
    nPacketCurr += nPacketSize;
  }
      
}

void receiveFlows(t_Flows* ptFlows)
{
  int i = 0, nSize = 0, nPackets = 0, nPacketSize = 0, nPacketCurr = 0;

  MPI_Bcast((void *) &ptFlows->nN,1,MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast((void *) &ptFlows->nM,1,MPI_INT, 0, MPI_COMM_WORLD);

  ptFlows->adF = (double *) malloc(ptFlows->nM*ptFlows->nN*sizeof(double));
  if(!ptFlows->adF)
    goto memoryError;

  ptFlows->anF = (int *) malloc(ptFlows->nM*ptFlows->nN*sizeof(int));
  if(!ptFlows->anF)
    goto memoryError;

  ptFlows->anLengths = (int *) malloc(ptFlows->nN*sizeof(int));
  if(!ptFlows->anLengths)
    goto memoryError;

  nSize = ptFlows->nN*ptFlows->nM;

  nPackets = ((nSize - 1)/MAX_PACKET_SIZE) + 1;
      
  printf("nN=%d nM=%d nSize=%d\n",ptFlows->nN,ptFlows->nM,nSize);fflush(stdout);
  nPacketCurr = 0;
  for(i = 0; i < nPackets; i++){
       	
    if(i < nPackets - 1){
      nPacketSize = MAX_PACKET_SIZE;
    }
    else{
      nPacketSize = nSize % MAX_PACKET_SIZE;
    }
	
    //printf("Receive Total=%d Packet=%d Size=%d\n",nPackets,i,nPacketSize);fflush(stdout);

    MPI_Bcast((void *) (&ptFlows->adF[nPacketCurr]),nPacketSize,MPI_DOUBLE,0, MPI_COMM_WORLD);

    MPI_Bcast((void *) (&ptFlows->anF[nPacketCurr]),nPacketSize,MPI_INT,0, MPI_COMM_WORLD);
	
    nPacketCurr += nPacketSize;
  }

  nSize = ptFlows->nN;

  nPackets = ((nSize - 1)/MAX_PACKET_SIZE) + 1;
      
  nPacketCurr = 0;
  for(i = 0; i < nPackets; i++){
       	
    if(i < nPackets - 1){
      nPacketSize = MAX_PACKET_SIZE;
    }
    else{
      nPacketSize = nSize % MAX_PACKET_SIZE;
    }
       
    //printf("Receive Total=%d Packet=%d Size=%d\n",nPackets,i,nPacketSize);fflush(stdout);

    MPI_Bcast((void *) (&ptFlows->anLengths[nPacketCurr]),nPacketSize,MPI_INT,0, MPI_COMM_WORLD);
	
    nPacketCurr += nPacketSize;
  }

  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in receiveFlows\n"); fflush(stderr);
  exit(EXIT_FAILURE);
}


void readData(double *adLookUp, char* szDataFile, t_Flows* ptFlows, t_Params *ptParams)
{
  char szLine[MAX_LINE_LENGTH];
  char *szTok;
  int nN = 0, nM = 0;
  double *adData = NULL;
  double *adF    = NULL;
  int    *anF    = NULL;
  int    *anLengths = NULL;
  char *pcError = NULL;
  int i = 0, j = 0, nE = 0;
  FILE *ifp = NULL;

  ifp = fopen(szDataFile, "r");

  if(ifp){
    /*read header line*/
    fgets(szLine, MAX_LINE_LENGTH, ifp);
    /*get read number*/
    szTok = strtok(szLine, DELIM);
    nN = strtol(szTok,&pcError, 10);
    if(*pcError != '\0')
      goto formatError;

    /*get flow number*/
    szTok = strtok(NULL,DELIM);
    nM = strtol(szTok,&pcError, 10);
    if(*pcError != '\0')
      goto formatError;
    
    adData = (double *) malloc(nN*nM*sizeof(double));
    if(!adData)
      goto memoryError;

    adF = (double *) malloc(nN*nM*sizeof(double));
    if(!adF)
      goto memoryError;

    anF = (int *) malloc(nN*nM*sizeof(int));
    if(!anF)
      goto memoryError;

    anLengths = (int *) malloc(nN*sizeof(int));
    if(!anLengths)
      goto memoryError;

    for(i = 0; i < nN; i++){
      fgets(szLine, MAX_LINE_LENGTH, ifp);
      
      szTok = strtok(szLine, DELIM);

      if(ptParams->bNoIndex == FALSE){ 
	szTok = strtok(NULL, DELIM);
      }

      anLengths[i] = strtol(szTok, &pcError,10);
      if(*pcError != '\0')
	goto formatError;
    
      nE = anLengths[i] % 4;
      anLengths[i] -= nE;
      //      printf("%d %d %d\n",anLengths[i], nE, nM);
      for(j = 0; j < anLengths[i]; j++){
	double dF = 0.0;

	szTok = strtok(NULL, DELIM);
	dF = strtod(szTok, &pcError);
	//	printf("%d %s %f\n",j, szTok,dF);
	if(dF >= 9.99){
	  //fprintf(stderr,"%d %d %f too large max flow 9.99\n",i,j,dF);
	  dF = 9.99;
	  fflush(stderr);
	}

	adData[i*nM + j] = dF;
	anF[i*nM + j] =  (int) floor((dF + DELTA)/PRECISION);
	if(anF[i*nM + j] < 0){
	  anF[i*nM + j] = 0;
	}
	if(anF[i*nM + j] > BINS - 1){
	  anF[i*nM + j] = BINS - 1;
	}

	adF[i*nM + j] = distM1(adLookUp, dF);
	
	if(*pcError != '\0')
	  goto formatError;
      }

      for(; j < nM; j++){
	double dF = 0.0;
	adData[i*nM + j] = dF;
	anF[i*nM + j] =  (int) floor((dF + DELTA)/PRECISION);
	adF[i*nM + j] = distM1(adLookUp, dF);
	
      }
    }

    fclose(ifp);
  }
  else{
    printf("Failed to open file %s in readData\n",szDataFile);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  ptFlows->anLengths = anLengths;
  ptFlows->adData    = adData;
  ptFlows->adF       = adF;
  ptFlows->anF       = anF;
  ptFlows->nN = nN; ptFlows->nM = nM;
  return;
 memoryError:
  fprintf(stderr, "Failed allocating memory in readData\n");
  fflush(stderr);
  exit(EXIT_FAILURE);

 formatError:
  fprintf(stderr, "Incorrect input format error\n");
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
  ptParams->szDataFile  = extractParameter(argc,argv, DATA_FILE, ALWAYS);  
  if(ptParams->szDataFile == NULL)
    goto error;

  if(extractParameter(argc,argv, NO_INDEX, OPTION)){
    ptParams->bNoIndex = TRUE;
  }
  else{
    ptParams->bNoIndex = FALSE;
  }

  if(szTemp = extractParameter(argc,argv, LOOKUP_FILE_FLAG, OPTION)){
    ptParams->szLookUpFile = szTemp;
  }
  else{
    ptParams->szLookUpFile = getenv("PYRO_LOOKUP_FILE");
    if(ptParams->szLookUpFile == NULL){
    	ptParams->szLookUpFile = LOOKUP_FILE;
    }	
  }


  ptParams->szOutFileStub  = extractParameter(argc,argv, OUT_FILE_STUB, ALWAYS);  
  if(ptParams->szOutFileStub == NULL)
    goto error;

  return;

 error:
  writeUsage(stdout);
  exit(EXIT_FAILURE);
}



void needlemanWunsch(double *adFookUp, double* adF1, 
		      double* adF2, t_Align *ptAlign, int* anA, 
		      int* anB, int nLenA, int nLenB)
{
  double **aadFMatrix = NULL;
  int    **aanMoves   = NULL;
  int     *anPath     = NULL;
  int    nCount = 0;
  int    i = 0, j = 0, k = 0; 
  
  aadFMatrix = (double **) malloc((nLenA + 1)*sizeof(double *));
  aanMoves   = (int **)    malloc((nLenA + 1)*sizeof(int *));
  if(!aadFMatrix || !aanMoves)
    goto memoryError;

  anPath = (int *) malloc(sizeof(int)*(nLenA + nLenB));
  if(!anPath)
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

  //printf("allocated memory %d %d\n",nLenA,nLenB); fflush(stdout);

  for(i = 1; i <= nLenA; i++){
    for(j = 1; j <= nLenB; j++){
      double dChoice1, dChoice2, dChoice3;
      
      dChoice1 = aadFMatrix[i-1][j-1] + adFookUp[anA[i - 1]*BINS +  anB[j - 1]] - adF1[i - 1] - adF2[j - 1];
      
      if(j >= 4){
	if(i == nLenA){
	  dChoice2 = aadFMatrix[i][j-4];
	} 
	else{
	  dChoice2 = aadFMatrix[i][j-4] + 4.0*GAP_PENALTY;
	}
      }
      else{
         dChoice2 = MAX_DBL;
      }

      if(i >= 4){
	if(j == nLenB){
	  dChoice3 = aadFMatrix[i-4][j];
	}
	else{
	  dChoice3 = aadFMatrix[i-4][j] + 4.0*GAP_PENALTY;
	}
      }
      else{
	dChoice3 = MAX_DBL;
      }
      aanMoves[i][j]   = getMove(dChoice1, dChoice2, dChoice3);
      aadFMatrix[i][j] = dmin3(dChoice1, dChoice2, dChoice3);
    }
  }

  ptAlign->dDist = aadFMatrix[nLenA][nLenB];
  
  i = nLenA; j = nLenB;
  while(i > 0 && j > 0){
    anPath[nCount] = aanMoves[i][j];
    switch(aanMoves[i][j]){
    case DIAG:
      i--;
      j--;
      nCount++;
      break;
    case UP:
      i--;
      nCount++;
      k = 1;
      while(k < 4 && i >= 0){

	anPath[nCount] = UP;
	i--;

	k++;
	nCount++;
      }
      break;
    case LEFT:
      j--;
      nCount++;
      while(k < 4 && j >= 0){

	anPath[nCount] = LEFT;
	j--;

	k++;
	nCount++;
      }
      break;
    }
  }
  
  while(i > 0){
    anPath[nCount] = aanMoves[i][j];
    i--;
    nCount++;
  }
    
  while(j > 0){
    anPath[nCount] = aanMoves[i][j];
    j--;
    nCount++;
  }

  ptAlign->nLen  = nCount;
  ptAlign->nComp = nCount;
  
  if(anPath[0] != DIAG){
    int nStart = anPath[0];
    i = 0;
    while(anPath[i] == nStart && i < nCount){
      ptAlign->nComp--;
      i++;
    }
  }
  //printf("free memory %d %d\n",nLenA,nLenB); fflush(stdout);
  /*free up memory*/
  free(anPath);

  for(i = 0; i <= nLenA; i++){
    free(aadFMatrix[i]);
    free(aanMoves[i]);
  }
  free(aadFMatrix);
  free(aanMoves);
  return;
 memoryError:
  fprintf(stderr, "Failed to allocate memory in needlemanWunsch\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

double FDist(double *adFookUp, double* adF1, double* adF2, int* anA, int* anB, int nLenA, int nLenB)
{
  int    nCount = 0;
  int    i = 0; 
  int    nComp = min(nLenA, nLenB);
  double dDist = 0.0;

  for(i = 0; i < nComp; i++){
    dDist += adFookUp[anA[i]*BINS +  anB[i]] - adF1[i] - adF2[i];
  }
    
  dDist /= (double) nComp;
    
  return dDist;
}



int imin(int nA, int nB)
{
  return nA > nB ? nB : nA;
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

void initFookUp(double* adFookUp, double *adLookUp)
{
  int i = 0, j = 0;

  for(i = 0; i < BINS; i++){
    for(j = 0; j < BINS; j++){
      adFookUp[i*BINS + j] = distM(adLookUp, i, j);
    }
  }
}

void initLookUp(t_Params *ptParams, double* adLookUp)
{
  int i = 0, j = 0;
  FILE *ifp = NULL;
  double dTemp = 0.0;

  ifp = fopen(ptParams->szLookUpFile,"r");
  if(ifp){
    char szLine[MAX_LINE_LENGTH];
    char *pcError = NULL;
    char *szRet = NULL;

    for(i = 0; i < MAX_S; i++){
      fgets(szLine, MAX_LINE_LENGTH, ifp);
	
      szRet = strpbrk(szLine, "\n");

      (*szRet) = '\0';

      dTemp = strtod(szLine,&pcError);

      for(j = 0; j < BINS; j++){
	fgets(szLine, MAX_LINE_LENGTH, ifp);
	
	szRet = strpbrk(szLine, "\n");

	(*szRet) = '\0';

	adLookUp[i*BINS + j] = strtod(szLine,&pcError);
	if(*pcError != '\0')
	  goto formatError;
      }
    }
   
    fclose(ifp);
  }
  else{
    fprintf(stderr,"Failed to open %s", ptParams->szLookUpFile);
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

double distS(double* adLookUp, int nS, double dFlow)
{
  int nFlow = (int) floor((dFlow + DELTA)/PRECISION);

  if (nS < MAX_S && nFlow < BINS){
    return adLookUp[nS*BINS + nFlow];
  }

  return MAX_DIST;
}

double distI(double *adLookUp, int nFlow1, int nFlow2)
{                          
  int    s = 0;                                      
  double distF = 0.0; 
  double dNLL = 0.0;
  int    nSB = 0;
  int    nMax   = max(nFlow1,nFlow2) + 3;
  
  if(nMax > MAX_S){
    nMax = MAX_S;
  }

  for(s = 0; s < nMax; s++){//loop signal strength
    nSB = s*BINS;

    dNLL = adLookUp[nSB + nFlow1] + adLookUp[nSB + nFlow2]; 
    
    distF += exp(-dNLL);
  }
	  
  return -log(distF);
}

double distM(double *adLookUp, int nFlow1, int nFlow2)
{                          
  int    s     = 0;                                      
  double distF = MAX_DBL; 
  double dNLL  = 0.0;
  int    nSB   = 0;
  
  for(s = 0; s < MAX_S; s++){//loop signal strength
    nSB = s*BINS;

    dNLL = adLookUp[nSB + nFlow1] + adLookUp[nSB + nFlow2]; 
    
    if(dNLL < distF){
      distF = dNLL;
    }
  }
	  
  return distF;
}

double distF1(double *adLookUp, double dF1)
{                          
  int    s = 0;                                      
  double distF = 0.0; 
  double dNLL = 0.0;
  
  for(s = 0; s < MAX_S; s++){//loop signal strength
    double dT1 = distS(adLookUp, s, dF1);

    dNLL = dT1; 
    
    distF += exp(-dNLL);
  }
	  
  return -log(distF);
}

double distM1(double *adLookUp, double dF1)
{                          
  int    s = 0;                                      
  double distF = MAX_DBL; 
  
  for(s = 0; s < MAX_S; s++){//loop signal strength
    double dT1 = distS(adLookUp, s, dF1);

    if(dT1 < distF){
      distF = dT1; 
    }
 
  }
	  
  return distF;
}
