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
#include <sys/stat.h>
#include <mpi.h>

#include "PyroNoiseM.h"

/*global constants*/
static char *usage[] = {"PyroNoiseM - clusters flowgrams without alignments\n",
			"-din     string            flow file name\n",
			"-out     string            cluster input file name\n",
			"-lin     string            list file\n",
			"Options:\n",
			"-v       verbose\n",
			"-c       double            initial cut-off\n",
			"-ni                        no index in dat file\n",
			"-i       integer           number of iterations (default 1000)\n",
			"-s       double            precision\n",
			"-rin     file              lookup file name\n"};

static int  nLines = 11;

static double *adLookUp = NULL;

static double *adSPrior = NULL;

static int bVerbose = FALSE;

static char szFlows[] = "TACG";
//static char szFlows2[] = "GATC";
static int  nOffSet   = 8;
//static int  nOffSet2  = 9;

int main(int argc, char* argv[]){
  int i = 0, j = 0, k = 0, n = 0, nN = 0, nM = 0, nSize = 0, nPacket = 0;
  t_Params tParams;
  t_Flows  tFlows;
  short    *asData;
  int      *anLengths;
  /*no. of clusters*/
  int      nK = 0, nKEff = 0;
  /*provisional assignments*/
  int      *anCentroids = NULL;
  double   *adDist      = NULL;
  float    *afDistX     = NULL;
  t_Unique tUnique;
  int*     anChange = NULL;
  int      numtasks, rank, rc;   
  MPI_Status   status;
  int      nA = 0, nA0 = 0, nTag = 1;
  int      nU = 0, nI = 0, nI0 = 0;
  int      bCont = TRUE;
  double*  adWeight = NULL;
  int      *anZ = NULL;
  t_Master tMaster;
  double   dSigma;


  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  /*get command line params*/
  getCommandLineParams(&tParams, argc, argv);

  /*everyone allocates memory*/
  if(bVerbose){printf("%d: Allocate memory for lookup\n", rank);}
  allocLookUp();

  if(rank == 0){//head node reads data etc.
    /*variables localised to head node*/
    double   *adP = NULL;
    double   *adLast      = NULL;
    int      nIter = 0;
    double   dDelta = 0.0, dMaxDelta = 0.0;
    char     **aacQualities = NULL;
    double   **aadQualities = NULL;
    double   *adNorm = NULL;

    t_Letter tLetter;

    if(bVerbose){printf("%d: Init lookup\n", rank);}
    fflush(stdout);

    initLookUp(&tParams);

    if(bVerbose){printf("%d: Read data\n", rank);}
    fflush(stdout);

    readData(tParams.szDataFile, &tFlows, &tParams);

    calcUnique(&tUnique, &tFlows);

    nN = tFlows.nN; nM = tFlows.nM;
    asData = tFlows.asData; anLengths = tFlows.anLengths;

    /*broadcast data*/
    if(bVerbose){printf("%d: Broadcast data N = %d M = %d\n", rank, nN, nM);}
    
    MPI_Bcast(adLookUp, MAX_S*BINS, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(adSPrior, MAX_S, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*broadcast flows*/
    if(bVerbose){printf("%d: Broadcast data flows\n", rank);}
    fflush(stdout);

    broadcastFlows(&tFlows);

    broadcastUnique(&tUnique);
    
    /*---------------*/

    anZ = (int *) malloc(nN*sizeof(int));
    if(!anZ)
      goto memoryError;

    adNorm = (double *) malloc(nN*sizeof(double));
    if(!adNorm)
      goto memoryError;
    
    /*read init file*/
    if(tParams.szInitFile != NULL){
      if(bVerbose){printf("%d: Read init file %s\n", rank, tParams.szInitFile);}
      readInitFile(&tParams, anZ);
    }    

    nK     = tParams.nK;

    if(bVerbose){printf("%d: Broadcast cluster number K = %d\n", rank, nK); fflush(stdout);}
    
    MPI_Bcast((void *) &nK, 1,MPI_INT, 0, MPI_COMM_WORLD);    

    dSigma = tParams.dSigma;

    MPI_Bcast((void *) &dSigma, 1,MPI_DOUBLE, 0, MPI_COMM_WORLD);

    nI = (int) (floor(nN / numtasks));

    nI0 = nI + (nN % numtasks);

    nSize = nI0*nM*nK;
    nPacket = nN*nM;

    if(bVerbose){printf("%d: Allocate alignment size = %d packetsize = %d\n", rank, nSize, nPacket); fflush(stdout);}

    adDist = (double *) malloc(nI0*nK*sizeof(double));
    if(!adDist)
      goto memoryError;

    for(i = 0; i < nI0*nK; i++){
      adDist[i] = 0.0;
    }

    nU = tUnique.nU;
    afDistX = (float *) malloc(nN*nU*sizeof(float));
    if(!afDistX)
      goto memoryError;

    /*initialize alignments and tau*/	
    if(bVerbose){printf("%d: Calculate distances %d %d\n", rank, nU, nN);fflush(stdout);}	
    calcDistX(afDistX, 0, nI0, &tFlows, &tUnique);
	
    for (i = 1; i < numtasks; i++){
	int nIStart = (nI0 + (i - 1)*nI);
	if(bVerbose){printf("%d: Receive distances %d %d\n", rank, nIStart, nI*nU);fflush(stdout);}	
	MPI_Recv(&afDistX[nIStart*nU], nI*nU, MPI_FLOAT, i, nTag, MPI_COMM_WORLD, &status);
    }	
	
    if(bVerbose){printf("%d: Broadcast distances\n", rank);fflush(stdout);}
    MPI_Bcast(afDistX, nU*nI0, MPI_FLOAT, 0, MPI_COMM_WORLD);
    for(i = 1; i < numtasks; i++){
	int nIStart = (nI0 + (i - 1)*nI);
	
	MPI_Bcast(&afDistX[nIStart*nU], nU*nI, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    adWeight = (double *) malloc(nK*sizeof(double));
    if(!adWeight)
      goto memoryError;

    for(i = 0; i < nK; i++){
      adWeight[i] = 0.0;
    }

    anChange = (int *) malloc(nK*sizeof(int));
    if(!anChange)
      goto memoryError;

    adP = (double *) malloc(nN*sizeof(double));
    if(!adP)
      goto memoryError;

    anCentroids = (int *) malloc(nK*sizeof(int));
    if(!anCentroids)
      goto memoryError;

    for(i = 0; i < nK; i++){
      anCentroids[i] = NOT_SET;
    }

    /*initialize alignments and tau*/
    if(bVerbose){printf("%d: Initialize alignments and tau\n",rank);fflush(stdout);}

    /*allocate memory for sparse array*/
    if(bVerbose){printf("%d: Allocate master\n",rank);fflush(stdout);}
    allocateMaster(&tMaster, nM, nN, nK, INIT_MASTER_SIZE);

    allocateLetter(&tLetter, nM, INIT_MASTER_SIZE);

    /*Broadcast Z*/
    if(bVerbose){printf("%d: Broadcast Z\n",rank);fflush(stdout);}
    MPI_Bcast(anZ, nN, MPI_INT, 0, MPI_COMM_WORLD);

    if(bVerbose){printf("%d: Init. alignment nI0 = %d nN = %d\n",rank, nI0, tMaster.nN);fflush(stdout);}
    initAlignment(&tMaster, &tFlows, nK, anZ, anChange, nI0);

    nA = (int) (floor(nK / numtasks));

    nA0 = nA + (nK % numtasks);

    if(bVerbose){
      printf("%d: Partition work %d %d %d %d\n",rank, nA0, nA, nI0, nI);
      fflush(stdout);
    }

    while((nIter < MIN_ITER) || ((dMaxDelta > MIN_DELTA) && (nIter < tParams.nMaxIter))){
      /*Broadcast tau and alignments*/
      if(bVerbose){printf("%d-%d: Broadcast master data\n",rank,nIter);fflush(stdout);}
      MPI_Bcast(&bCont, 1, MPI_INT, 0, MPI_COMM_WORLD);

      /*Broadcast master data*/
      fillMaster(&tMaster);

      broadcastMaster(&tMaster);

      /*calculate nK div numtasks centroids*/
      if(bVerbose){printf("%d-%d: Calc centroids\n",rank,nIter); fflush(stdout);}

      calcCentroidsMaster(anChange, 0, nA0, &tFlows, &tMaster, anCentroids, &tUnique, afDistX);
      
      /*then receive the other data*/
      if(bVerbose){printf("%d-%d: Receive data\n",rank,nIter);fflush(stdout);}
      for (i = 1; i < numtasks; i++){
	int nStart = (nA0 + (i - 1)*nA);
	
	MPI_Recv(&anCentroids[nStart], nA, 
		 MPI_INT, i, nTag, MPI_COMM_WORLD, &status);
	
	MPI_Recv(&anChange[nStart], nA, 
		 MPI_INT, i, nTag, MPI_COMM_WORLD, &status);
      }

      printf("%d ",nIter);
      /*calc weights*/

      dMaxDelta = calcNewWeights(nK, adWeight, &tMaster);
    
      printf(" %f",calcNLLikelihood(&tMaster, &nKEff, adWeight, dSigma));

      printf(" %d %f\n",nKEff, dMaxDelta);
      fflush(stdout);
    
      /*remove degenerate centroids*/
      checkCentroidUniqueness(adWeight, anCentroids, nK);

      /*broadcast weights*/
      MPI_Bcast(adWeight, nK, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      /*broadcast centroids*/
      MPI_Bcast(anCentroids, nK, MPI_INT, 0, MPI_COMM_WORLD);

      /*broadcast change*/
      MPI_Bcast(anChange, nK, MPI_INT, 0, MPI_COMM_WORLD);

      if(bVerbose){printf("%d-%d: Calc distances\n",rank,nIter);fflush(stdout);}

      calcDistancesMaster(adDist, anChange, adNorm, &tMaster, dSigma, anCentroids, &tFlows, adWeight, nI0,&tUnique, afDistX);
      
      /*receiveLetters*/
      receiveLetters(nTag, adNorm, numtasks, &tMaster, &tLetter, nI0, nI);
      
      if(bVerbose){
	printf("%d: Updated master total = %d, size = %d\n", rank, tMaster.nTotal, tMaster.nSize); fflush(stdout);
      }
      //writeMasterN(nIter,&tMaster, &tUnique, anCentroids, &tFlows, &tParams);
      fflush(stdout);
      nIter++;
    }

    /*Fill master*/
    fillMaster(&tMaster);

    /*kill off slaves*/
    bCont = FALSE;
    MPI_Bcast((void *) &bCont, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /*find closest centroid for each read*/
    setZ(anZ, adP, &tMaster, &tParams);

    /*output centroids at end of EM algorithm*/
    //writeState(nN, anZ, nK, nM, adWeight, asCentroids, &tParams);
    writeMaster(&tMaster, &tUnique, anCentroids, &tFlows, &tParams);
    //writeMaster(&tMaster, asCentroids, &tFlows, &tParams);
    /*output clusters*/
    if(bVerbose){printf("%d: Output clusters\n",rank);fflush(stdout);}
    //writeClustersF(&tParams, nN, nM, nK, adWeight, asCentroids);

    setMasterZ(&tMaster, &tFlows, nN, nK, anZ);

    {
      int anT[nK];

      for(i = 0; i < nK; i++){
	anT[i] = 0;
      }
      for(i = 0; i < nN; i++){
	anT[anZ[i]]++;
      }

      calcCentroidsMaster(anChange, 0, nK, &tFlows, &tMaster, anCentroids, &tUnique, afDistX);
   
      calcQualitiesMaster(nK, &tFlows, &tMaster, &aacQualities, anCentroids, &tUnique);

      writeClustersD(&tParams, nK, nM, anT, anCentroids, &tUnique, tFlows.aszID);

      writeQualitiesD(&tParams, nK, anT, aacQualities, anCentroids, &tUnique, tFlows.aszID);

      outputClusters(&tParams, nK, anZ, nN, nM, anT, &tFlows, tFlows.aszID, anCentroids, &tUnique, tFlows.aszID); 

      outputMap(&tParams, nK, anZ, nN, nM, anT, anCentroids, &tUnique, &tFlows);
    }

    /*free up memory*/
    destroyLetter(&tLetter);

    destroyUnique(&tUnique);

    free(adP);
    free(adDist);
    free(afDistX);
    free(adWeight);
    free(anChange);

    for(i = 0; i < nN; i++){
      free(tFlows.aszID[i]);
    }

    for(i = 0; i < nK; i++){
      if(tMaster.anN[i] > 0){
	free(aacQualities[i]);
      }
    }

    free(aacQualities);
    free(adLast);

    free(anCentroids);
    free(adNorm);
    free(tFlows.aszID);
    free(tFlows.adData);
    free(tFlows.asData);
    free(tFlows.anLengths);
  }
  else{
    int nIStart = 0, nIFinish = 0;
    double* adLocalNorm = NULL;

    t_Letter tLetter;

    /*recieve data*/
    if(bVerbose){printf("%d: Receive data\n",rank);}

    MPI_Bcast(adLookUp, MAX_S*BINS, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(adSPrior, MAX_S, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(bVerbose){printf("%d: Receive flows\n",rank);}

    receiveFlows(&tFlows);
    nN = tFlows.nN; nM = tFlows.nM;
    asData = tFlows.asData; anLengths = tFlows.anLengths;

    receiveUnique(&tUnique);

    if(bVerbose){printf("%d: Receive dK and dSigma\n",rank);}

    MPI_Bcast(&nK, 1,MPI_INT, 0, MPI_COMM_WORLD);    

    MPI_Bcast((void *) &dSigma, 1,MPI_DOUBLE, 0, MPI_COMM_WORLD);

    tParams.dSigma = dSigma;
    tParams.nK = nK;
    
    /*memory allocation on slave nodes*/
    nI = (int) (floor(nN / numtasks));

    nI0 = nI + (nN % numtasks);

    nSize   = nM*nI*nK;
    nPacket = nM*nN;

    if(bVerbose){printf("%d: Allocate memory nI = %d nK = %d\n",rank, nI, nK);}

    adDist = (double *) malloc(nI*nK*sizeof(double));
    if(!adDist)
      goto memoryError;

    nU = tUnique.nU;
    afDistX = (float *) malloc(nN*nU*sizeof(float));
    if(!afDistX)
      goto memoryError;

    if(bVerbose){printf("%d: Calc distances\n", rank);fflush(stdout);}	
    
    nIStart = (nI0 + (rank - 1)*nI);
	
    calcDistX(afDistX, nIStart, nIStart + nI, &tFlows, &tUnique);
    
    if(bVerbose){printf("%d: Send %d %d distances\n", rank, nIStart*nU, nI*nU);fflush(stdout);}	
    MPI_Send(&afDistX[nIStart*nU], nI*nU, MPI_FLOAT, 0, nTag, MPI_COMM_WORLD);

    if(bVerbose){printf("%d: Receive broadcast distances\n", rank);fflush(stdout);}
    
    MPI_Bcast(afDistX, nU*nI0, MPI_FLOAT, 0, MPI_COMM_WORLD);
    for(i = 1; i < numtasks; i++){
	int nIStart = (nI0 + (i - 1)*nI);
	
	MPI_Bcast(&afDistX[nIStart*nU], nU*nI, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    anChange = (int *) malloc(nK*sizeof(int));
    if(!anChange)
      goto memoryError;

    adWeight = (double *) malloc(nK*sizeof(double));
    if(!adWeight)
      goto memoryError;

    anCentroids = (int *) malloc(nK*sizeof(int));
    if(!anCentroids)
      goto memoryError;

    for(i = 0; i < nK; i++){
      anCentroids[i] = NOT_SET;
    }

    /*--------------------------*/
    nA = (int) (floor(nK / numtasks));

    nA0 = nA + (nK % numtasks);

    nIStart  = nI0 + (rank - 1)*nI;
    nIFinish = nIStart + nI;

    anZ = (int *) malloc(nN*sizeof(int));
    if(!anZ)
      goto memoryError;

    adLocalNorm = (double *) malloc(nI*sizeof(double));
    if(!adLocalNorm)
      goto memoryError;

    if(bVerbose){printf("%d: Receive anZ\n",rank);}

    MPI_Bcast(anZ, nN ,MPI_INT, 0, MPI_COMM_WORLD);

    if(bVerbose){printf("%d: Init. slave\n",rank);}
    
    allocateMaster(&tMaster, nM, nN, nK, INIT_MASTER_SIZE);

    allocateLetter(&tLetter, nM, INIT_MASTER_SIZE);

    fflush(stdout);
    while(TRUE){
      int nStart = nA0 + (rank - 1)*nA;
      int nFinish = nStart + nA;
      /*Receive tau and alignments*/
      if(bVerbose){printf("%d: Receive data\n",rank);}
      MPI_Bcast((void *) &bCont, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(bCont == FALSE)
	goto finish;

      receiveMaster(&tMaster);

      if(bVerbose){printf("%d: Calc centroids\n",rank);}

      calcCentroidsMaster(anChange, nStart, nFinish, &tFlows, &tMaster, anCentroids, &tUnique, afDistX);

      if(bVerbose){printf("%d: Return centroids\n",rank);}
      MPI_Send(&anCentroids[nStart], nA, MPI_INT, 0, nTag, MPI_COMM_WORLD);

      MPI_Send(&anChange[nStart], nA, MPI_INT, 0, nTag, MPI_COMM_WORLD);
      
      /*receive weights*/
      MPI_Bcast(adWeight, nK, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      /*receive centroids*/
      MPI_Bcast(anCentroids, nK, MPI_INT, 0, MPI_COMM_WORLD);

      /*receive change*/
      MPI_Bcast(anChange, nK, MPI_INT, 0, MPI_COMM_WORLD);

      if(bVerbose){printf("%d: Calc distances\n",rank);}

      calcDistancesSlave(adDist, nM, nK, anChange, adLocalNorm, &tLetter, dSigma,
			      anCentroids, &tFlows, adWeight, nIStart, nIFinish, &tUnique, afDistX);

      if(bVerbose){printf("%d: Sending letter size %d\n",rank,tLetter.nTotal);fflush(stdout);}

      minimiseLetter(&tLetter);

      MPI_Send(&tLetter.nTotal, 1, MPI_INT, 0, nTag, MPI_COMM_WORLD);
	
      MPI_Send(tLetter.anK, tLetter.nTotal, MPI_INT, 0, nTag, MPI_COMM_WORLD);

      MPI_Send(tLetter.anI, tLetter.nTotal, MPI_INT, 0, nTag, MPI_COMM_WORLD);

      MPI_Send(tLetter.adT, tLetter.nTotal, MPI_DOUBLE, 0, nTag, MPI_COMM_WORLD);

      MPI_Send(tLetter.adD, tLetter.nTotal, MPI_DOUBLE, 0, nTag, MPI_COMM_WORLD);

      MPI_Send(adLocalNorm, nI, MPI_DOUBLE, 0, nTag, MPI_COMM_WORLD);

      fflush(stdout);

    }

  finish:

    destroyLetter(&tLetter);

    destroyUnique(&tUnique);

    free(asData); asData = NULL;
    
    free(anLengths); anLengths = NULL;

    free(anChange); anChange = NULL;

    free(adWeight); adWeight = NULL;

    free(adLocalNorm); adLocalNorm = NULL;

    free(adDist); adDist = NULL;

    free(afDistX); afDistX = NULL;
  }

  destroyMaster(&tMaster);

  free(anZ); anZ = NULL;

  /*free lookup table*/
  free(adLookUp); adLookUp = NULL;

  free(adSPrior); adSPrior = NULL;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(EXIT_SUCCESS);

 memoryError:
  fprintf(stderr, "Failed allocating memory in main\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void readData(char* szDataFile, t_Flows* ptFlows, t_Params *ptParams)
{
  char szLine[MAX_LINE_LENGTH];
  char *szTok;
  int nN = 0, nM = 0;
  double *adData = NULL;
  short  *asData = NULL;
  int    *anLengths = NULL;
  char   **aszID = NULL;
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

    aszID = (char **) malloc(nN*sizeof(char *));
    if(!aszID)
      goto memoryError;

    adData = (double *) malloc(nN*nM*sizeof(double));
    if(!adData)
      goto memoryError;

    asData = (short *) malloc(nN*nM*sizeof(short));
    if(!asData)
      goto memoryError;

    anLengths = (int *) malloc(nN*sizeof(int));
    if(!anLengths)
      goto memoryError;

    for(i = 0; i < nN; i++){
      fgets(szLine, MAX_LINE_LENGTH, ifp);
      szTok = strtok(szLine, DELIM);

      if(ptParams->bNoIndex == FALSE){
	aszID[i] = strdup(szTok);
      }
      else{
	aszID[i] = (char *) malloc(sizeof(char)*ID_LENGTH);
	sprintf(aszID[i],"%d",i);
      }

      if(ptParams->bNoIndex == FALSE){
	szTok = strtok(NULL, DELIM);
      }
      
      anLengths[i] = strtol(szTok, &pcError,10);
      if(*pcError != '\0')
	goto formatError;
    
      nE = anLengths[i] % 4;
      anLengths[i] -= nE;

      for(j = 0; j < anLengths[i]; j++){
	szTok = strtok(NULL, DELIM);
	adData[i*nM + j] = strtod(szTok, &pcError);
	if(adData[i*nM + j] >= 9.49){
		adData[i*nM + j] = 9.49;
	}
	asData[i*nM + j] = (short) floor((adData[i*nM + j] + DELTA)/PRECISION);
	if(*pcError != '\0')
	  goto formatError;
      }

      for(; j < nM; j++){
	adData[i*nM + j] = 0.0;
	asData[i*nM + j] = (short) floor((adData[i*nM + j] + DELTA)/PRECISION);
	if(*pcError != '\0')
	  goto formatError;
      }
    }
  }
  else{
    printf("Failed to open file %s in readData\n",szDataFile);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  ptFlows->aszID     = aszID;
  ptFlows->anLengths = anLengths;
  ptFlows->adData    = adData;
  ptFlows->asData    = asData;
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

  /*get output filestub*/
  ptParams->szOutFileStub  = extractParameter(argc,argv, OUT_FILE_STUB, ALWAYS);  
  if(ptParams->szOutFileStub == NULL)
    goto error;

  ptParams->szInitFile  = extractParameter(argc,argv, INIT_FILE, ALWAYS);  
  if(ptParams->szInitFile == NULL)
    goto error;

  if(szTemp = extractParameter(argc,argv, LOOKUP_FILE_FLAG, OPTION)){
    ptParams->szLookUpFile = szTemp;
  }
  else{
    ptParams->szLookUpFile = getenv("PYRO_LOOKUP_FILE");
    if(ptParams->szLookUpFile == NULL){
        ptParams->szLookUpFile = LOOKUP_FILE;
    }  
  }

  szTemp  = extractParameter(argc,argv, SIGMA, OPTION);  
  if(szTemp != NULL){
    ptParams->dSigma = strtod(szTemp, &cError);
    if(*cError != '\0')
      goto error;
  }
  else{
    ptParams->dSigma = DEF_SIGMA;
  }

  szTemp  = extractParameter(argc,argv,N_MAX_ITER, OPTION);  
  if(szTemp != NULL){
    ptParams->nMaxIter = strtol(szTemp, &cError, 10);
    if(*cError != '\0')
      goto error;
  }
  else{
    ptParams->nMaxIter = MAX_ITER;
  }

  if(extractParameter(argc,argv, VERBOSE, OPTION)){
    bVerbose = TRUE;
  }
  else{
    bVerbose = FALSE;
  }

  if(extractParameter(argc,argv, NO_INDEX, OPTION)){
    ptParams->bNoIndex = TRUE;
  }
  else{
    ptParams->bNoIndex = FALSE;
  }

  szTemp  = extractParameter(argc,argv, INIT_CUT, OPTION);  
  if(szTemp != NULL){
    ptParams->dInitCut = strtod(szTemp, &cError);
    if(*cError != '\0')
      goto error;
  }
  else{
    ptParams->dInitCut = DEF_CUT;
  }
  
  return;

 error:
  writeUsage(stdout);
  exit(EXIT_FAILURE);
}


double alignX(short* asA, short* asB, int nLenS, int nLenF)
{
  double dDist = 0.0;
  int    i     = 0;
  int    nComp = nLenS;

  if(nLenF < nComp){
	nComp = nLenF;
  }

  for(i = 0; i < nComp; i++){
    dDist += adLookUp[asA[i]*BINS + asB[i]];
  }

  if(nLenF > nLenS){
    dDist += UNSEEN_PENALTY*(nLenF - nLenS);
    nComp = nLenF;
  }
  //else{
  //dDist += 5.0*(nLenS - nLenF);
  //}

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

int getMove(double *pdMin, double dA, double dB, double dC)
{

  if(dA < dB){
    if(dA < dC){
      (*pdMin) = dA;
      return DIAG;
    }
    else{
      (*pdMin) = dC;
      return UP;
    }
  }
  else{
    if(dB < dC){
      (*pdMin) = dB;
      return LEFT;
    }
    else{
      (*pdMin) = dC;
      return UP;
    }
  }
}

double find_closest_pair(int n, double** distmatrix, int* ip, int* jp)
{ 
  int i, j;
  double temp;
  double distance = distmatrix[1][0];
  *ip = 1;
  *jp = 0;
  for (i = 1; i < n; i++)
  { for (j = 0; j < i; j++)
    { temp = distmatrix[i][j];
      if (temp<distance)
      { distance = temp;
        *ip = i;
        *jp = j;
      }
    }
  }
  return distance;
}

void outputCluster(t_Params *ptParams, t_Node* tree, char **aszID, int nN)
{
  int i = 0, j = 0, k = 0;
  int nNodes = nN - 1;
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

  for(i = 0; i < nNodes; i++){
    int left = tree[i].left, right = tree[i].right;
    int put  = -1, destroy = -1;
    double dPrint = tree[i].distance;

    fprintf(lfp, "%f %d ",dPrint, nN - i);
    fprintf(ofp, "%f %d ",dPrint, nN - i);
   
    for(j = 0; j < nN; j++){
      if(anCounts[j] > 0){
	for(k = 0; k < anCounts[j] - 1; k++){
	  fprintf(lfp, "%s,",aaszClusters[j][k]);
	}

	fprintf(lfp,"%s ",aaszClusters[j][anCounts[j] - 1]);

	fprintf(ofp, "%d ", anCounts[j]);
      }
    }

    fprintf(lfp, "\n"); fprintf(ofp, "\n");
    fflush(lfp); fflush(ofp);
  

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
   
    anIndex[i] = put;
    
    for(j = 0; j < anCounts[destroy]; j++){
      if(anCounts[put] == anSizes[put]){
	anSizes[put] *= 2;
	aaszClusters[put] = (char **) realloc(aaszClusters[put],anSizes[put]*sizeof(char *));
      }
      aaszClusters[put][anCounts[put]++] = aaszClusters[destroy][j];
    }
    anCounts[destroy] = 0;
    free(aaszClusters[destroy]);
  
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
  free(anIndex);
  free(anCounts); 
  free(anSizes); 
  free(szOtuFile);
  free(szListFile);
  return;
}

void allocLookUp()
{
  adLookUp = (double *) malloc(BINS*MAX_S*sizeof(double));
  if(!adLookUp)
    goto memoryError;

  adSPrior = (double *) malloc(MAX_S*sizeof(double));
  if(!adSPrior)
    goto memoryError;

  return;

  memoryError:
  fprintf(stderr, "Failed allocating memory in allocLookUp\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void initLookUp(t_Params *ptParams)
{
  int i = 0, j = 0;
  FILE *ifp = NULL;

  ifp = fopen(ptParams->szLookUpFile,"r");
  if(ifp){
    char szLine[MAX_LINE_LENGTH];
    char *pcError = NULL;
    char *szRet = NULL;

    for(i = 0; i < MAX_S; i++){
      fgets(szLine, MAX_LINE_LENGTH, ifp);
	
      szRet = strpbrk(szLine, "\n");

      (*szRet) = '\0';

      adSPrior[i] = strtod(szLine,&pcError);
      
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
    fprintf(stderr,"Failed to open %s\n",ptParams->szLookUpFile);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  return;
 formatError:
  fprintf(stderr, "Format error LookUp.dat\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
} 

double distS(int nS, double dFlow)
{
  int nFlow = (int) floor((dFlow + DELTA)/PRECISION);

  if (nS < MAX_S && nFlow < BINS){
    return adLookUp[nS*BINS + nFlow];
  }

  return MAX_DIST;
}

void setZ(int *anZ, double *adP, t_Master* ptMaster, t_Params *ptParams)
{
  int i = 0, j = 0, k = 0;
  int nN = ptMaster->nN, nK = ptMaster->nK;
  double* adTauMatrix = (double *) malloc(nN*nK*sizeof(double));
  double* adTau = ptMaster->adTau;
  int *anN = ptMaster->anN, *anCN = ptMaster->anCN, *anI = ptMaster->anI, *anP = ptMaster->anP;

  if(!adTauMatrix)
    goto memoryError;
			
  for(i = 0; i < nN*nK; i++){
    adTauMatrix[i] = 0.0;
  }

  for(k = 0; k < nK; k++){ //loop clusters  

    for(i = 0; i < anN[k]; i++){ //loop data points
      int    nIndex = anCN[k] + i;
      double dTau   = adTau[anP[nIndex]];
      int    nI     = anI[nIndex];
  
      adTauMatrix[nI*nK + k] = dTau;
    }
  }  


  for(i = 0; i < nN; i++){
    double dMaxTau = -1.0;
    int    nMaxK   = -1;

    for(j = 0; j < nK; j++){
      if(adTauMatrix[i*nK + j] > dMaxTau){
	dMaxTau = adTauMatrix[i*nK + j];
	nMaxK   = j;
      }
    }

    anZ[i] = nMaxK;
    adP[i] = 1.0 - dMaxTau;
  }

  //writeTau(adTauMatrix, nN, nK, ptParams);

  free(adTauMatrix);
  return;
  
 memoryError:
  fprintf(stderr,"Failed allocating memory in setZ\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

char* flowToSeq(int *pnSeqLength, double *adFlows, int nLength)
{
  char *szSeq = (char *) malloc(MAX_SEQUENCE_LENGTH*sizeof(char));
  int  i = 0;
  int  nCount = 0;

  if(!szSeq)
    goto memoryError;

  for(i = nOffSet; i < nLength; i++){
    char cBase = szFlows[i % 4];
    int  s = 0, nS    =  (int) floor(adFlows[i] + 0.5);
    
    while(s < nS){
      szSeq[nCount] = cBase;
      s++;
      nCount++;
    }
  }

  (*pnSeqLength) = nCount;
  return szSeq;
 memoryError:
  fprintf(stderr, "Failed allocating memory in flowToSeq\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void writeSequenceUnaligned(FILE* ofp, char *szLabel, char *szSequence, int nLen)
{
  int nPos = 0;
  int nUL = 0, i = 0;

  for(i = 0; i < nLen; i++){
    if(szSequence[i] != '-'){
      nUL++;
    } 
  }

  fprintf(ofp,">%s\n",szLabel);

  i = 0;
  while(i < nLen && nPos < nUL){
    if(szSequence[i] != '-'){

      if(nPos > 0 && nPos % 80 == 0){
	fputc('\n',ofp);
      }

      fputc(szSequence[i],ofp);

      nPos++;
    }
    i++;
  }

  fputc('\n',ofp);
}

void writeClustersF(t_Params *ptParams, int nN, int nM, int nK, double *adWeight, short* asCentroids)
{
  char *szClustFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;
  int i = 0, j = 0, s = 0;

  sprintf(szClustFile, "%s_cf%s", ptParams->szOutFileStub, FASTA_SUFFIX);
  ofp = fopen(szClustFile, "w");

  for(i = 0; i < nK; i++){
    if(adWeight[i] > 0){
      fprintf(ofp, ">%s_%d_%.2f\n",ptParams->szOutFileStub,i, adWeight[i]);

      for(j = nOffSet; j < nM; j++){
	char   cBase = szFlows[j % 4];

	for(s = 0; s < asCentroids[i*nM + j]; s++){
	  fprintf(ofp, "%c", cBase);
	}
      }
      fprintf(ofp,"\n");
    }
  }

  fclose(ofp);
  free(szClustFile);
}

void writeClustersD(t_Params *ptParams, int nK, int nM, int *anT, int* anCentroids, t_Unique *ptUnique, char** aszIDs)
{
  char *szClustFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;
  int i = 0, j = 0, s = 0, index = 0;
  int *anF = ptUnique->anF;

  sprintf(szClustFile, "%s_cd%s", ptParams->szOutFileStub, FASTA_SUFFIX);
  ofp = fopen(szClustFile, "w");

  for(i = 0; i < nK; i++){
    int nUI = anCentroids[i];

    if(anT[i] > 0){
      fprintf(ofp, ">%s_%d_%d\n",ptParams->szOutFileStub,index, anT[i]);

      for(j = nOffSet; j < nM; j++){
	char   cBase = szFlows[j % 4];

	for(s = 0; s < ptUnique->asU[nUI*nM + j]; s++){
	  fprintf(ofp, "%c", cBase);
	}
      }
      
      fprintf(ofp,"\n");

      index++;
    }
  }

  fclose(ofp);
  free(szClustFile);
}

void writeQualitiesD(t_Params *ptParams, int nK, int *anT, char** aacQualities, int* anCentroids, t_Unique *ptUnique, char** aszIDs)
{
  char *szClustFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;
  int i = 0, j = 0, s = 0, index = 0;

  sprintf(szClustFile, "%s_cd%s", ptParams->szOutFileStub, QUAL_SUFFIX);
  ofp = fopen(szClustFile, "w");
  
  for(i = 0; i < nK; i++){
    if(anT[i] > 0){
      fprintf(ofp, ">%s_%d_%d\n",aszIDs[ptUnique->anF[anCentroids[i]]],index, anT[i]);

      j = 4;
      while(aacQualities[i][j] != (char) NOT_SET){
	fprintf(ofp, "%d ", (int) aacQualities[i][j]);
	j++;
      }
      fprintf(ofp,"\n");

      index++;
    }
  }

  fclose(ofp);
  free(szClustFile);
}


void writeClustersD2(t_Params *ptParams, int nN, int nM, int nK, int *anZ, short* asCentroids)
{
  char *szClustFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;
  int i = 0, j = 0, s = 0;
  int*  anT = (int * ) malloc(nK*sizeof(int));

  sprintf(szClustFile, "%s_cd2%s", ptParams->szOutFileStub, FASTA_SUFFIX);
  ofp = fopen(szClustFile, "w");
  
  for(i = 0; i < nK; i++){
    anT[i] = 0;
  }
  for(i = 0; i < nN; i++){
    anT[anZ[i]]++;
  }

  for(i = 0; i < nK; i++){
    if(anT[i] > 0){
      fprintf(ofp, ">%s_%d_%d\n",ptParams->szOutFileStub,i, anT[i]);

      for(j = nOffSet; j < nM; j++){
	char   cBase = szFlows[j % 4];

	for(s = 0; s < asCentroids[i*nM + j]; s++){
	  fprintf(ofp, "%c", cBase);
	}
      }
      fprintf(ofp,"\n");
    }
  }

  free(anT);
  fclose(ofp);
  free(szClustFile);
}

void writeClustersD3(t_Params *ptParams, int nN, int nM, int nK, int *anZ, short* asCentroids)
{
  char *szClustFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;
  int i = 0, j = 0, s = 0;
  int*  anT = (int * ) malloc(nK*sizeof(int));

  sprintf(szClustFile, "%s_cd3%s", ptParams->szOutFileStub, FASTA_SUFFIX);
  ofp = fopen(szClustFile, "w");
  
  for(i = 0; i < nK; i++){
    anT[i] = 0;
  }
  for(i = 0; i < nN; i++){
    anT[anZ[i]]++;
  }

  for(i = 0; i < nK; i++){
    if(anT[i] > 0){
      fprintf(ofp, ">%s_%d_%d\n",ptParams->szOutFileStub,i, anT[i]);

      for(j = nOffSet; j < nM; j++){
	char   cBase = szFlows[j % 4];

	for(s = 0; s < asCentroids[i*nM + j]; s++){
	  fprintf(ofp, "%c", cBase);
	}
      }
      fprintf(ofp,"\n");
    }
  }

  free(anT);
  fclose(ofp);
  free(szClustFile);
}

void writeZ(t_Params *ptParams, t_Flows *ptFlows, int *anZ, double *adNorm)
{
  char *szZFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;
  int i = 0;
  
  sprintf(szZFile, "%s%s", ptParams->szOutFileStub, Z_SUFFIX);
  
  ofp = fopen(szZFile, "w");
  
  for(i = 0; i < ptFlows->nN; i++){
    fprintf(ofp, "%d %s %d %f\n", i, ptFlows->aszID[i], anZ[i], adNorm[i]);
  }

  fclose(ofp);
  free(szZFile);

  return;
}

void readInitFile(t_Params *ptParams, int* anZ)
{
  FILE *ifp = NULL;
  char *szLine = NULL; 
  char *szTok, *pcError;
  int  i = 0;

  szLine = (char *) malloc(sizeof(char)*BIG_LINE_LENGTH);
  if(!szLine)
    goto memoryError;

  ifp = fopen(ptParams->szInitFile, "r");
  if(ifp){
     /*read header line*/
    while(fgets(szLine, BIG_LINE_LENGTH, ifp) != NULL){
      /*get read number*/
      szTok = strtok(szLine, DELIM);
      
      if(strtod(szTok,&pcError) == ptParams->dInitCut){
	szTok = strtok(NULL, DELIM);
	ptParams->nK = strtol(szTok, &pcError, 10);
	if(*pcError != '\0')
	  goto fileFormatError;

	for(i = 0; i < ptParams->nK; i++){
	  char *szDup = NULL, *szTemp = NULL, *szStart = NULL;
	  int  nIndex = 0;
	  szTok = strtok(NULL, DELIM);
	  szDup = strdup(szTok);
	  szTemp = szDup;
	  //printf("%s\n",szDup);
	  while(szStart = strstr(szTemp,",")){
	    (*szStart) = '\0';
	    nIndex = strtol(szTemp,&pcError,10);
	    anZ[nIndex] = i;
	    szTemp = szStart + 1;
	    
	  }
	  nIndex = strtol(szTemp,&pcError,10);
	  anZ[nIndex] = i;
	  free(szDup);
	}

      }
    }
    fclose(ifp);
  }
  else{
    fprintf(stderr, "Failed to open %s\n", ptParams->szInitFile);
    fflush(stderr);
  }
  free(szLine);
  return;

 fileFormatError:
  fprintf(stderr, "Format error in %s\n", ptParams->szInitFile);
  fflush(stderr);
  exit(EXIT_FAILURE);

 memoryError:
  fprintf(stderr, "Failed allocating memory in readInitFile\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}


void outputClusters(t_Params *ptParams, int nK, int *anZ, int nN, int nM, int* anT, t_Flows *ptFlows, char **aszLabels, int* anCentroids, t_Unique *ptUnique, char** aszIDs)
{
  FILE *ofp = NULL;
  char *szOutFile = (char *) malloc(sizeof(char)*MAX_LINE_LENGTH);
  int  s = 0, i = 0, j = 0, index = 0;

  mkdir(ptParams->szOutFileStub, S_IRWXU);
  
  if(!szOutFile)
    goto memoryError;

  for(i = 0; i < nK; i++){
    if(anT[i] > 0){
      sprintf(szOutFile, "%s/i_%d%s",ptParams->szOutFileStub,index, FASTA_SUFFIX);
      ofp = fopen(szOutFile, "w");
      if(ofp){
	int nUI = anCentroids[i];
	fprintf(ofp, ">%s_%d_%d\n",ptParams->szOutFileStub,index, anT[i]);

	for(j = nOffSet; j < nM; j++){
	  char   cBase = szFlows[j % 4];

	  for(s = 0; s < ptUnique->asU[nUI*nM + j]; s++){
	    fprintf(ofp, "%c", cBase);
	  }
	}
      
	fprintf(ofp,"\n");

	for(j = 0; j < nN; j++){
	  if(anZ[j] == i){
	    int  nLen     = 0;
	    char *szSeq   = flowToSeq(&nLen, &(ptFlows->adData[j*nM]), ptFlows->anLengths[j]);

	    writeSequenceUnaligned(ofp, aszLabels[j], szSeq, nLen);

	    free(szSeq);
	  }
	}
	
	fclose(ofp);
      }
      else{
	fprintf(stderr,"Can't open file %s for writing\n",szOutFile);
      }

      index++;
    }
  }
				      
  free(szOutFile);

  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in outputClusters\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void outputMap(t_Params *ptParams, int nK, int *anZ, int nN, int nM, int* anT, int *anCentroids, t_Unique *ptUnique, t_Flows *ptFlows)
{
  FILE *ofp = NULL;
  char *szOutFile = (char *) malloc(sizeof(char)*MAX_LINE_LENGTH);
  int  i = 0, j = 0, index = 0;
  
  if(!szOutFile)
    goto memoryError;

  sprintf(szOutFile, "%s%s",ptParams->szOutFileStub,MAP_SUFFIX);

  ofp = fopen(szOutFile, "w");

  if(ofp){
    for(i = 0; i < nK; i++){
      if(anT[i] > 0){
      	int bFirst = TRUE;
	fprintf(ofp, "%s",ptFlows->aszID[ptUnique->anF[anCentroids[i]]]);

	for(j = 0; j < nN; j++){
	  if(anZ[j] == i){
	    if(bFirst == FALSE){
	    	fprintf(ofp, ",%s",ptFlows->aszID[j]);
	    }
	    else{
	    	fprintf(ofp, " %s",ptFlows->aszID[j]);
		bFirst = FALSE;
	    }
	  }
	}

	fprintf(ofp, "\n");
       
	index++;
      }
    }

    fclose(ofp);
  }
  else{
    fprintf(stderr,"Can't open file %s for writing\n",szOutFile);
  }				      
  
  free(szOutFile);

  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in outputClusters\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void allocateMaster(t_Master *ptMaster, int nM, int nN, int nK, int nSize)
{
  int i      = 0, j = 0;

  ptMaster->nM = nM;
  ptMaster->nN = nN;
  ptMaster->nK = nK;

  ptMaster->nSize = nSize;
  ptMaster->nTotal = 0;

  ptMaster->adTau = (double *) malloc(sizeof(double)*nSize);
  if(!ptMaster->adTau)
    goto memoryError;

  ptMaster->adDist = (double *) malloc(sizeof(double)*nSize);
  if(!ptMaster->adDist)
    goto memoryError;
  
  for(i = 0; i < nSize; i++){
    ptMaster->adTau[i]  = 0.0;
    ptMaster->adDist[i] = 0.0;
  }

  ptMaster->anP = (int *) malloc(sizeof(int)*nSize);
  if(!ptMaster->anP)
    goto memoryError;

  for(i = 0; i < nSize; i++){
    ptMaster->anP[i] = 0;
  }

  ptMaster->anI = (int *) malloc(sizeof(int)*nSize);
  if(!ptMaster->anI)
    goto memoryError;

  for(i = 0; i < nSize; i++){
    ptMaster->anI[i] = 0;
  }

  ptMaster->anN = (int *) malloc(sizeof(int)*nK);
  if(!ptMaster->anN)
    goto memoryError;

  ptMaster->anCN = (int *) malloc(sizeof(int)*nK);
  if(!ptMaster->anCN)
    goto memoryError;

  ptMaster->aanP = (int **) malloc(sizeof(int*)*nK);
  ptMaster->aanI = (int **) malloc(sizeof(int*)*nK);
  if(!ptMaster->aanP || !ptMaster->aanI)
    goto memoryError;

  for(i = 0; i < nK; i++){
    ptMaster->anN[i]  = 0;
    ptMaster->anCN[i] = 0;

    ptMaster->aanP[i]  = (int *) malloc(nN*sizeof(int));
    ptMaster->aanI[i]  = (int *) malloc(nN*sizeof(int));
  
    for(j = 0; j < nN; j++){
      ptMaster->aanP[i][j] = 0;
      ptMaster->aanI[i][j] = 0;
    }

  }

  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in allocateMaster\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void reallocateMaster(t_Master *ptMaster)
{
  int nM = ptMaster->nM, i = 0; 
  int nOldSize = ptMaster->nSize;

  while(ptMaster->nSize < ptMaster->nTotal){
    ptMaster->nSize  *= 2;
  }

  ptMaster->adTau = (double *) realloc(ptMaster->adTau, sizeof(double)*ptMaster->nSize);
  if(!ptMaster->adTau)
    goto memoryError;

  ptMaster->adDist = (double *) realloc(ptMaster->adDist, sizeof(double)*ptMaster->nSize);
  if(!ptMaster->adDist)
    goto memoryError;

  for(i = nOldSize; i < ptMaster->nSize; i++){
    ptMaster->adTau[i]  = 0.0;
    ptMaster->adDist[i] = 0.0;
  }

  ptMaster->anP = (int *) realloc(ptMaster->anP,sizeof(int)*ptMaster->nSize);
  if(!ptMaster->anP)
    goto memoryError;

  for(i = nOldSize; i < ptMaster->nSize; i++){
    ptMaster->anP[i] = 0;
  }

  ptMaster->anI = (int *) realloc(ptMaster->anI,sizeof(int)*ptMaster->nSize);
  if(!ptMaster->anI)
    goto memoryError;

  for(i = nOldSize; i < ptMaster->nSize; i++){
    ptMaster->anI[i] = 0;
  }


  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in reallocMaster %d\n",ptMaster->nSize);
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void minimiseLetter(t_Letter *ptLetter)
{
  int nM = ptLetter->nM;

  ptLetter->nSize  = ptLetter->nTotal;

  ptLetter->adT = (double *) realloc(ptLetter->adT, sizeof(double)*ptLetter->nSize);
  if(!ptLetter->adT)
    goto memoryError;

  ptLetter->adD = (double *) realloc(ptLetter->adD, sizeof(double)*ptLetter->nSize);
  if(!ptLetter->adD)
    goto memoryError;


  ptLetter->anI = (int *) realloc(ptLetter->anI,sizeof(int)*ptLetter->nSize);
  if(!ptLetter->anI)
    goto memoryError;

  ptLetter->anK = (int *) realloc(ptLetter->anK,sizeof(int)*ptLetter->nSize);
  if(!ptLetter->anK)
    goto memoryError;

  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in minimiseLetter\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void minimiseMaster(t_Master *ptMaster)
{
  int nM = ptMaster->nM;

  ptMaster->nSize  = ptMaster->nTotal;

  ptMaster->adTau = (double *) realloc(ptMaster->adTau, sizeof(double)*ptMaster->nSize);
  if(!ptMaster->adTau)
    goto memoryError;

  ptMaster->adDist = (double *) realloc(ptMaster->adDist, sizeof(double)*ptMaster->nSize);
  if(!ptMaster->adDist)
    goto memoryError;

  ptMaster->anP = (int *) realloc(ptMaster->anP,sizeof(int)*ptMaster->nSize);
  if(!ptMaster->anP)
    goto memoryError;

  ptMaster->anI = (int *) realloc(ptMaster->anI,sizeof(int)*ptMaster->nSize);
  if(!ptMaster->anI)
    goto memoryError;

  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in minimiseMaster %d\n", ptMaster->nSize);
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void updateMasterI(int nI, t_Master *ptMaster, double *adT, double *adD)
{
  int i = 0, j = 0, nK = ptMaster->nK;
  

  for(i = 0; i < nK; i++){
    if(adT[i] > MIN_TAU){
      int nOldTotal = ptMaster->nTotal;

      ptMaster->nTotal++;

      reallocateMaster(ptMaster);

      ptMaster->adTau[nOldTotal] = adT[i];
      
      ptMaster->adDist[nOldTotal] = adD[i];

      ptMaster->aanP[i][ptMaster->anN[i]] = nOldTotal;
      
      ptMaster->aanI[i][ptMaster->anN[i]] = nI;

      ptMaster->anN[i]++;
      
    }
  }
}

void updateMasterLetter(t_Master *ptMaster, t_Letter *ptLetter)
{
  int i = 0, j = 0, nAdd = 0;
  int nOldTotal = ptMaster->nTotal;

  ptMaster->nTotal += ptLetter->nTotal;

  reallocateMaster(ptMaster);

  for(i = nOldTotal; i < ptMaster->nTotal; i++){
    int nK = ptLetter->anK[nAdd];
    int nI = ptLetter->anI[nAdd];
    //printf("Add %d to %d cluster %d sample %d\n",nAdd,i,nK,nI); fflush(stdout);
    ptMaster->adTau[i]  = ptLetter->adT[nAdd];
    ptMaster->adDist[i] = ptLetter->adD[nAdd];  
    ptMaster->aanP[nK][ptMaster->anN[nK]] = i;
      
    ptMaster->aanI[nK][ptMaster->anN[nK]] = nI;

    ptMaster->anN[nK]++;
    nAdd++;
  }
}



      

void initAlignment(t_Master *ptMaster, t_Flows *ptFlows, int nK, int *anZ, int *anChange, int nI0)
{
/*initialize alignments and tau*/
  int nN = ptFlows->nN, nM = ptFlows->nM;
  int i = 0, j = 0, k = 0;
  double adTau[nK], adDist[nK];
  
  for(k = 0; k < nK; k++){
    adDist[k] = 0.0;
  }

  for(i = 0; i < nI0; i++){
    int nLen = 0;

    for(k = 0; k < nK; k++){
      adTau[k]  = 0.0;
    }

    k = anZ[i];
    adTau[k] = 1.0;

    updateMasterI(i, ptMaster, adTau, adDist);
  }

  for(i = nI0; i < nN; i++){
    int nLen = 0;
    
    for(k = 0; k < nK; k++){
      adTau[k] = 0.0;
    }

    k = anZ[i];
    adTau[k] = 1.0;

    updateMasterI(i, ptMaster, adTau, adDist);
  }

  for(i = 0; i < nK; i++){anChange[i] = TRUE;}

  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in initAlignment\n");
  fflush(stderr);
  exit(EXIT_FAILURE);

}

void fillMaster(t_Master *ptMaster)
{
  int i = 0, j = 0;
  int nCount = 0;

  /*realloc to current size*/
  minimiseMaster(ptMaster);

  for(i = 0; i < ptMaster->nK; i++){
    ptMaster->anCN[i] = nCount;
    for(j = 0; j < ptMaster->anN[i]; j++){
      ptMaster->anP[nCount] = ptMaster->aanP[i][j];
      ptMaster->anI[nCount] = ptMaster->aanI[i][j];
      nCount++;
    }
  }
}

void destroyMaster(t_Master *ptMaster)
{
  int i = 0;

  free(ptMaster->adTau);
  free(ptMaster->adDist);
  
  for(i = 0; i < ptMaster->nK; i++){
    free(ptMaster->aanP[i]);   /*pointers*/
    free(ptMaster->aanI[i]);   /*members*/
  }

  free(ptMaster->aanP);
  free(ptMaster->aanI);

  free(ptMaster->anN);  
  free(ptMaster->anCN);
  free(ptMaster->anP);
  free(ptMaster->anI);
}

void broadcastFlows(t_Flows *ptFlows)
{
  int nN = ptFlows->nN, nM = ptFlows->nM;

  MPI_Bcast(&ptFlows->nN,1,MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(&ptFlows->nM,1,MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptFlows->asData, nN*nM, MPI_SHORT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptFlows->anLengths, nN, MPI_INT, 0, MPI_COMM_WORLD);
}

void broadcastUnique(t_Unique *ptUnique)
{
  int i = 0, nN = ptUnique->nN, nU = ptUnique->nU, nM = ptUnique->nM;

  MPI_Bcast(&ptUnique->nN,1,MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(&ptUnique->nM,1,MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(&ptUnique->nU,1,MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptUnique->anMap, nN, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptUnique->anF, nU, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptUnique->anLenU, nU, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptUnique->asU, nU*nM, MPI_SHORT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptUnique->anWeights, nU, MPI_INT, 0, MPI_COMM_WORLD);
}

void receiveUnique(t_Unique *ptUnique)
{
  int nN = -1, nM = -1, nU = -1, i = 0;

  MPI_Bcast(&ptUnique->nN,1,MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(&ptUnique->nM,1,MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(&ptUnique->nU,1,MPI_INT, 0, MPI_COMM_WORLD);

  nN = ptUnique->nN; nM = ptUnique->nM; nU = ptUnique->nU;
  
  ptUnique->nSize = nU;
  ptUnique->nMSize = nU*nM;
  
  ptUnique->anMap = (int *) malloc(nN*sizeof(int));
  if(!ptUnique->anMap)
    goto memoryError;

  MPI_Bcast(ptUnique->anMap, nN, MPI_INT, 0, MPI_COMM_WORLD);

  ptUnique->anF = (int *) malloc(nU*sizeof(int));
  if(!ptUnique->anF)
    goto memoryError;

  MPI_Bcast(ptUnique->anF, nU, MPI_INT, 0, MPI_COMM_WORLD);

  ptUnique->anLenU = (int *) malloc(nU*sizeof(int));
  if(!ptUnique->anLenU)
    goto memoryError;

  MPI_Bcast(ptUnique->anLenU, nU, MPI_INT, 0, MPI_COMM_WORLD);

  ptUnique->asU = (short *) malloc(nU*nM*sizeof(short));
  if(!ptUnique->asU)
    goto memoryError;

  MPI_Bcast(ptUnique->asU, nU*nM, MPI_SHORT, 0, MPI_COMM_WORLD);

  ptUnique->anWeights = (int *) malloc(nU*sizeof(int));
  if(!ptUnique->anWeights)
    goto memoryError;

  MPI_Bcast(ptUnique->anWeights, nU, MPI_INT, 0, MPI_COMM_WORLD);

  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in receiveUnique\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void broadcastMaster(t_Master *ptMaster)
{
  MPI_Bcast(&ptMaster->nTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptMaster->adTau, ptMaster->nTotal, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
  MPI_Bcast(ptMaster->adDist, ptMaster->nTotal, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptMaster->anN, ptMaster->nK, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptMaster->anCN, ptMaster->nK, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptMaster->anP, ptMaster->nTotal, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptMaster->anI, ptMaster->nTotal, MPI_INT, 0, MPI_COMM_WORLD);
}

void receiveMaster(t_Master *ptMaster)
{
  MPI_Bcast(&ptMaster->nTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);

  minimiseMaster(ptMaster);

  MPI_Bcast(ptMaster->adTau, ptMaster->nTotal, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptMaster->adDist, ptMaster->nTotal, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
  MPI_Bcast(ptMaster->anN, ptMaster->nK, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptMaster->anCN, ptMaster->nK, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptMaster->anP, ptMaster->nTotal, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptMaster->anI, ptMaster->nTotal, MPI_INT, 0, MPI_COMM_WORLD);
}

void writeTau(double *adTau, int nN, int nK, t_Params* ptParams)
{
  int i = 0, j = 0;
  char *szTauFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;

  if(!szTauFile)
    goto memoryError;

  sprintf(szTauFile, "%s%s", ptParams->szOutFileStub, TAU_SUFFIX);
  ofp = fopen(szTauFile, "w");  
  if(ofp){
    for(i = 0; i < nN; i++){
      for(j = 0; j < nK; j++){
	fprintf(ofp, "%4.3f ",adTau[i*nK + j]);
      }
      fprintf(ofp, "\n");
    }
    fclose(ofp);
  }
  else{
    fprintf(stderr, "Failed to open %s for writing\n", szTauFile);
    fflush(stderr);
  }

  free(szTauFile);

  return;

memoryError:
  fprintf(stderr, "Failed allocating memory in writeTau\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void writeCentroids(t_Flows *ptFlows, int nK,  int *anZ, double* adWeights, short *asCentroids, t_Master *ptMaster, t_Params *ptParams)
{
  int nM = ptFlows->nM, nN = ptFlows->nN;
  int i = 0, j = 0, k = 0;
  short s = 0;
  short* asData = ptFlows->asData;
  double dTotal = 0.0, adS[MAX_S], adP[MAX_S];
  double dCount = 0.0;
  char *szCenFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;
  int   *anN = ptMaster->anN, *anCN = ptMaster->anCN;
  int   *anP = ptMaster->anP; 
  double *adT = ptMaster->adTau;
  int **aanI = ptMaster->aanI;
  int *anLengths = ptFlows->anLengths;

  if(!szCenFile)
    goto memoryError;

  sprintf(szCenFile, "%s%s", ptParams->szOutFileStub, CEN_SUFFIX);
  ofp = fopen(szCenFile, "w");  
  if(ofp){
    for(i = 0; i < nK; i++){ //loop clusters
      fprintf(ofp, "c=%d %f\n",i, adWeights[i]);

      for(j = 0; j < nM; j++){//loop position
	char   cBase = szFlows[j % 4];
	double dSMax = MAX_DBL;
	short  sSMax = -1;
	double dH = 0.0;

	for(s = 0; s < MAX_S; s++){//loop signal strength
	  adS[s] = 0.0;
	}

	dCount = 0.0;
	for(k = 0; k < anN[i]; k++){ //loop data points
	  int    nIndex = anCN[i] + k;
	  double dTau   = adT[anP[nIndex]];
	  int    nI     = aanI[i][k];
	  short  sFlow  = asData[nI*nM + j];
	    
	  if(j < anLengths[nI]){
	    dCount += dTau;
	    
	    for(s = 0; s < MAX_S; s++){//loop signal strength
	      adS[s] += dTau*adLookUp[s*BINS + sFlow];
	    }
	  }
	  
	}//loop data

      
	for(s = 0; s < MAX_S; s++){
	  if(adS[s] < dSMax){
	    dSMax = adS[s];
	    sSMax = s;
	  }
	}//loop signal

	dTotal = 0.0;
	for(s = 0; s < MAX_S; s++){
	  adP[s] = exp(-adS[s]+dSMax);
	  dTotal += adP[s];
	}

	for(s = 0; s < MAX_S; s++){
	  adP[s] /= dTotal;
	  if(adP[s] > 0.0){
	    dH += -adP[s]*log(adP[s]);
	  }
	}

	if(dCount > MIN_COUNT){
	  short sOld = asCentroids[i*nM + j];

	  asCentroids[i*nM + j] = sSMax;

	  fprintf(ofp, "%d %c %d %d %3.2f %3.2f %3.2f ",j, cBase,sSMax,sOld, adP[sSMax], dCount, dH);

	  for(k = 0; k < anN[i]; k++){ //loop data points
	    int    nIndex = anCN[i] + k;
	    int    nI     = aanI[i][k];
	    double dTau   = adT[anP[nIndex]];
	    int    sFlow  = asData[nI*nM + j];
	    
	    if(anZ[nI] == i){
	      fprintf(ofp, "%f ", sFlow*PRECISION);
	    }
	  }

	  fprintf(ofp,"\n");
	}
	else{
	  asCentroids[i*nM + j] = 0;
	}

      }//loop pos j

      fprintf(ofp, "\n");
    }//loop centroid i
  
    fclose(ofp);
  }
  else{
    fprintf(stderr, "Can't open szCenFile for writing\n",szCenFile);
    fflush(stderr);
  }

  free(szCenFile);
  return;

memoryError:
  fprintf(stderr, "Failed allocating memory in writeCentroids\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void allocateLetter(t_Letter *ptLetter, int nM, int nSize)
{
  int i = 0;

  ptLetter->nM     = nM;
  ptLetter->nSize  = nSize;
  ptLetter->nTotal = 0;

  ptLetter->anK =(int *) malloc(sizeof(int)*nSize);
  if(!ptLetter->anK)
    goto memoryError;

  ptLetter->anI =(int *) malloc(sizeof(int)*nSize);
  if(!ptLetter->anI)
    goto memoryError;

  ptLetter->adT =(double *) malloc(sizeof(double)*nSize);
  if(!ptLetter->adT)
    goto memoryError;

  ptLetter->adD =(double *) malloc(sizeof(double)*nSize);
  if(!ptLetter->adD)
    goto memoryError;

  for(i = 0; i < nSize; i++){
    ptLetter->anI[i] = 0;
    ptLetter->anK[i] = 0;
    ptLetter->adT[i] = 0.0;
    ptLetter->adD[i] = 0.0;
  }

  return;

memoryError:
  fprintf(stderr, "Failed allocating memory in allocateLetter\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void reallocateLetter(t_Letter *ptLetter)
{
  int nM = ptLetter->nM, i = 0; 
  int nOldSize = ptLetter->nSize;

  while(ptLetter->nSize < ptLetter->nTotal){
    ptLetter->nSize  *= 2;
  }

  ptLetter->anI = (int *) realloc(ptLetter->anI, sizeof(int)*ptLetter->nSize);
  if(!ptLetter->anI)
    goto memoryError;

  ptLetter->anK = (int *) realloc(ptLetter->anK, sizeof(int)*ptLetter->nSize);
  if(!ptLetter->anK)
    goto memoryError;

  ptLetter->adT = (double *) realloc(ptLetter->adT, sizeof(double)*ptLetter->nSize);
  if(!ptLetter->adT)
    goto memoryError;

  ptLetter->adD = (double *) realloc(ptLetter->adD, sizeof(double)*ptLetter->nSize);
  if(!ptLetter->adD)
    goto memoryError;

  for(i = nOldSize; i < ptLetter->nSize; i++){
    ptLetter->adT[i] = 0.0;
    ptLetter->anK[i] = 0;
    ptLetter->anI[i] = 0;
    ptLetter->adD[i] = 0.0;
  }

  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in reallocateLetter\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void destroyLetter(t_Letter *ptLetter)
{
  free(ptLetter->anI); ptLetter->anI = NULL;

  free(ptLetter->adT); ptLetter->adT = NULL;

  free(ptLetter->adD); ptLetter->adD = NULL;

  free(ptLetter->anK); ptLetter->anK = NULL;
}

int lenCmp(const void *pvA, const void* pvB)
{
  t_DoubleInt* ptA = (t_DoubleInt *) pvA;
  t_DoubleInt* ptB = (t_DoubleInt *) pvB;

  if(ptA->nLen < ptB->nLen){
    return -1;
  }
  else if(ptA->nLen > ptB->nLen){
    return +1;
  }
  else{
    if(ptA->nIdx < ptB->nIdx){
      return -1;
    }
    else{
      return +1;
    }
    return 0;
  }
}

void calcCentroidsMaster(int *anChange, int nKStart, int nKFinish, t_Flows *ptFlows, t_Master *ptMaster, int *anCentroids, t_Unique *ptUnique, float* afDistX)
{
  int i = 0, j = 0, k = 0;
  double *adTau = ptMaster->adTau;
  int nM   = ptMaster->nM;
  
  short *asData = ptFlows->asData;
  int   *anLengths = ptFlows->anLengths;
  int   *anN = ptMaster->anN, *anCN = ptMaster->anCN, *anP = ptMaster->anP, *anI = ptMaster->anI;
  int   **aanI = ptMaster->aanI;
  
  int   nU     = ptUnique->nU;
  int   *anMap = ptUnique->anMap;
  short *asU   = ptUnique->asU;
  int   *anLenU = ptUnique->anLenU;

  for(i = nKStart; i < nKFinish; i++){ //loop clusters
    int nIU = 0, nI = 0, nIndex = 0;  
    double dTau = 0.0, dCount = 0.0;
    double dMinF = MAX_DBL;
    int    nMinF = NOT_SET;

    anChange[i] = FALSE;

    for(j = 0; j < anN[i]; j++){ 
      nIndex   = anCN[i] + j;

      dCount += adTau[anP[nIndex]];
    }

    if(anN[i] > 0 && dCount > MIN_COUNT){
      double *adF = (double *) malloc(anN[i]*sizeof(double));
      int    *anL = (int *) malloc(anN[i]*sizeof(int));
      int    nL = 0;

      if(!anL)
	goto memoryError;

      if(!adF)
	goto memoryError;

      for(j = 0; j < anN[i]; j++){ 
	nIndex   = anCN[i] + j;
	nI       = anI[nIndex];
	nIU      = anMap[nI];

	for(k = 0; k < nL; k++){
	  if(nIU == anL[k]){
	    break;
	  }
	}

	if(k == nL){
	  anL[nL] = nIU;
	  adF[nL] = 0.0;
	  nL++;
	}
      }

      for(j = 0; j < anN[i]; j++){ 
	nIndex   = anCN[i] + j;
	nI       = anI[nIndex];
	nIU      = anMap[nI];
	dTau     = adTau[anP[nIndex]];

	for(k = 0; k < nL; k++){
	  double dDistX =  afDistX[nI*nU + anL[k]];

	  adF[k] += dDistX*dTau;
	}
      }

      for(k = 0; k < nL; k++){
	if(adF[k] < dMinF){
	  nMinF = k;
	  dMinF = adF[k];
	}
      }

      if(anCentroids[i] != anL[nMinF]){
	anChange[i] = TRUE;
	anCentroids[i] = anL[nMinF];
      } 
     
      free(adF);
      free(anL);
    }
    else{
      if(anCentroids[i] != NOT_SET){
	anChange[i] = TRUE;
	anCentroids[i] = NOT_SET;
      }
    }
  }

  return;
 memoryError:
  fprintf(stderr, "Failed allocating memory in calcCentroidsMaster\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void receiveFlows(t_Flows *ptFlows)
{
  int nN = 0, nM = 0;


  MPI_Bcast(&ptFlows->nN,1,MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(&ptFlows->nM,1,MPI_INT, 0, MPI_COMM_WORLD);

  nN = ptFlows->nN; nM = ptFlows->nM;

  ptFlows->asData = (short *) malloc(nM*nN*sizeof(short));
  if(!ptFlows->asData)
    goto memoryError;

  ptFlows->anLengths = (int *) malloc(nN*sizeof(int));
  if(!ptFlows->anLengths)
    goto memoryError;

  MPI_Bcast(ptFlows->asData, nN*nM, MPI_SHORT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptFlows->anLengths, nN, MPI_INT, 0, MPI_COMM_WORLD);

  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in receiveFlows\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}


void writeState(int nN, int* anZ, int nK, int nM, double* adWeight, short* asCentroids, t_Params *ptParams)
{
  int i = 0, j = 0, nTotal = 0;
  int  anMap[nK];
  FILE *ofp = NULL;
  char *szStateFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));

  if(!szStateFile)
    goto memoryError;
  
  sprintf(szStateFile, "%s%s", ptParams->szOutFileStub,STATE_SUFFIX);

  for(i = 0; i < nK; i++){
    if(adWeight[i] > MIN_WEIGHT){
      anMap[i] = nTotal;
      nTotal++;
    }
    else{
      anMap[i] = NOT_SET;
    }
  }

  ofp = fopen(szStateFile, "w");

  if(!ofp)
    goto fileError;

  fprintf(ofp, "%d %d\n",nTotal,nM);

  for(i = 0; i < nK; i++){
    if(adWeight[i] > MIN_WEIGHT){
      
      fprintf(ofp,"%d %f ",i, adWeight[i]);

      for(j = 0; j < nM; j++){
	fprintf(ofp, "%d ",asCentroids[i*nM + j]);
      }

      fprintf(ofp, "\n");
    }
  }

  for(i = 0; i < nN; i++){
    int nMap = anMap[anZ[i]];
    fprintf(ofp,"%d %d\n",i,nMap);
  }

  fclose(ofp);

  free(szStateFile);
  return;

 fileError:
  fprintf(stderr, "Failed opening file in writeState\n");
  fflush(stderr);
  exit(EXIT_FAILURE);

 memoryError:
  fprintf(stderr, "Failed alllocating memory in writeTau\n");
  fflush(stderr);
  exit(EXIT_FAILURE);

}

void setMasterZ(t_Master *ptMaster, t_Flows *ptFlows, int nN, int nK, int* anZ)
{
  int i = 0, j = 0;
  
  ptMaster->nTotal = nN;
  
  for(i = 0; i < nK; i++){ptMaster->anN[i] = 0;}

  minimiseMaster(ptMaster);

  for(i = 0; i < nN; i++){
    int k = anZ[i];

    ptMaster->nTotal++;

    ptMaster->adTau[i] = 1.0;
    ptMaster->adDist[i] = 0.0;

    ptMaster->aanP[k][ptMaster->anN[k]] = i;
    
    ptMaster->aanI[k][ptMaster->anN[k]] = i;

    ptMaster->anN[k]++;
  }
    
  fillMaster(ptMaster);
}

void calcQualitiesMaster(int nK, t_Flows *ptFlows, t_Master *ptMaster, char*** paacQualities, int *anCentroids, t_Unique *ptUnique)
{
  int i = 0, j = 0, k = 0;
  short s = 0;
  double adS[MAX_S], *adTau = ptMaster->adTau;
  short  *asData = ptFlows->asData;
  int nM   = ptMaster->nM; 
  int   *anN = ptMaster->anN, *anCN = ptMaster->anCN, *anP = ptMaster->anP;
  int   **aanI = ptMaster->aanI;
  char **aacQualities = NULL;

  aacQualities = (char **) malloc(nK*sizeof(char *));
  if(!aacQualities)
    goto memoryError;

  for(i = 0; i < nK; i++){ //loop clusters
    int j = 0, base  = 0;

    aacQualities[i] = NULL;

    if(anN[i] > 0){
      aacQualities[i] = (char*) malloc(MAX_SEQUENCE_LENGTH*sizeof(char));
      if(!aacQualities[i])
	goto memoryError;

      for(j = 0; j < MAX_SEQUENCE_LENGTH; j++){
	aacQualities[i][j] = NOT_SET;
      }

      j = 0;
      while(j < nM){//loop position
	double dSMax = MAX_DBL;
	short  sSMax = -1;
	double dCount = 0.0;
	
	/*set priors - not used!*/
	for(s = 0; s < MAX_S; s++){//loop signal strength
	  adS[s] = 0.0;
	}

	for(k = 0; k < anN[i]; k++){ //loop data points
	  int    nIndex  = anCN[i] + k;
	  double dTau    = adTau[anP[nIndex]];
	  int    nI      = aanI[i][k];
	  short  sFlow   = asData[nI*nM + j];
	   
	  dCount += dTau;

	  for(s = 0; s < MAX_S; s++){
	    adS[s] += dTau*adLookUp[s*BINS + sFlow];
	  }
	  
	}

	sSMax = ptUnique->asU[anCentroids[i]*nM + j];
	dSMax = adS[sSMax];

	if(dCount > MIN_COUNT){ 
	  double dU = 0.0, dNorm = 0.0;
	 
	  for(s = 0; s < MAX_S; s++){
	    dNorm += exp(-(adS[s] - dSMax));
	  }

	  for(s = 1; s <= sSMax; s++){
	    int nVal = 0;
	    double dTemp = 0.0;

	    dU += exp(-(adS[s - 1] - dSMax))/dNorm; 
	      
	    if(dU > 0.0){
	      dTemp = log10(dU);
	    }
	    else{
	      dTemp = -10.1;
	    }
	    dTemp = floor(-10.0*dTemp);
	    nVal = (int) dTemp;
	    if(nVal > 100){
	      nVal = 100;
	    }

	    aacQualities[i][base] = (char) nVal;
	    base++;
	  }
	}
	
	j++;
      }
    }
  }

  (*paacQualities) = aacQualities;
  return;
 memoryError:
  fprintf(stderr, "Failed allocating memory in calcCentroidsAndQualitiesMaster aborting...\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}


void calcCentroidsAndQualitiesMasterS(int nK, t_Flows *ptFlows, t_Master *ptMaster, char*** paacQualities, short *asCentroids)
{
  int i = 0, j = 0, k = 0;
  short s = 0;
  double adS[MAX_S], *adTau = ptMaster->adTau;
  short *asData = ptFlows->asData;
  int  nM   = ptMaster->nM; 
  int  *anN = ptMaster->anN, *anCN = ptMaster->anCN, *anP = ptMaster->anP;
  int  **aanI = ptMaster->aanI;
  char **aacQualities = NULL;

  aacQualities = (char **) malloc(nK*sizeof(char *));
  if(!aacQualities)
    goto memoryError;

  for(i = 0; i < nK; i++){ //loop clusters
    int j = 0, base  = 0;

    aacQualities[i] = NULL;

    if(anN[i] > 0){
      aacQualities[i] = (char*) malloc(MAX_SEQUENCE_LENGTH*sizeof(char));
      if(!aacQualities[i])
	goto memoryError;

      for(j = 0; j < MAX_SEQUENCE_LENGTH; j++){
	aacQualities[i][j] = NOT_SET;
      }

      j = 0;
      while(j < nM){//loop position
	double dSMax = MAX_DBL;
	short  sSMax = -1;
	double dCount = 0.0;
	
	/*set priors*/
	for(s = 0; s < MAX_S; s++){//loop signal strength
	  adS[s] = 0.0;
	}

	if(anN[i] == 1){
	    int    nIndex  = anCN[i];
	    double dTau    = adTau[anP[nIndex]];
	    int    nI      = aanI[i][0];
	    short  sFlow   = asData[nI*nM + j];

	    dCount = dTau;
	    sSMax  = (short) floor(((double) sFlow + 50.0)/100.0);

	    for(s = 0; s < MAX_S; s++){
	      adS[s] += dTau*adLookUp[s*BINS + sFlow];
	    }

	    dSMax = adS[sSMax];
	}
	else{
	  for(k = 0; k < anN[i]; k++){ //loop data points
	    int    nIndex  = anCN[i] + k;
	    double dTau    = adTau[anP[nIndex]];
	    int    nI      = aanI[i][k];
	    short  sFlow   = asData[nI*nM + j];
	   
	    dCount += dTau;

	    for(s = 0; s < MAX_S; s++){
	      adS[s] += dTau*adLookUp[s*BINS + sFlow];
	    }	 
	  }//loop data points

	  /*find best s*/

	  for(s = 0; s < MAX_S; s++){//loop signal strength  
	    if(adS[s] < dSMax){
	      dSMax = adS[s];
	      sSMax = s;
	    }
	  }//loop signal
	}
      
	if(dCount > MIN_COUNT){ 
	  double dU = 0.0, dNorm = 0.0;
 
	  asCentroids[i*nM + j] = sSMax;
	 
	  for(s = 0; s < MAX_S; s++){
	    dNorm += exp(-(adS[s] - dSMax));
	  }

	  for(s = 1; s <= sSMax; s++){
	    int nVal = 0;
	    double dTemp = 0.0;

	    dU += exp(-(adS[s - 1] - dSMax))/dNorm; 
	      
	    if(dU > 0.0){
	      dTemp = log10(dU);
	    }
	    else{
	      dTemp = -10.1;
	    }
	    dTemp = floor(-10.0*dTemp);
	    nVal = (int) dTemp;
	    if(nVal > 100){
	      nVal = 100;
	    }

	    aacQualities[i][base] = (char) nVal;
	    base++;
	  }
	}
	else{
	  asCentroids[i*nM + j] = 0;
	}
	j++;
      }//loop pos
    }
    else{
      for(j = 0; j < nM; j++){
	asCentroids[i*nM + j] = 0;
      }
    }

  }//loop centroid

  (*paacQualities) = aacQualities;
  return;
 memoryError:
  fprintf(stderr, "Failed allocating memory in calcCentroidsAndQualitiesMaster aborting...\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}
void writeQualitiesD2(t_Params *ptParams, int nN, int nM, int nK, int *anZ, char** aacQualities)
{
  char *szClustFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;
  int i = 0, j = 0, s = 0;
  int*  anT = (int * ) malloc(nK*sizeof(int));

  sprintf(szClustFile, "%s_cd2%s", ptParams->szOutFileStub, QUAL_SUFFIX);
  ofp = fopen(szClustFile, "w");
  
  for(i = 0; i < nK; i++){
    anT[i] = 0;
  }
  for(i = 0; i < nN; i++){
    anT[anZ[i]] ++;
  }

  for(i = 0; i < nK; i++){
    if(anT[i] > 0){
      fprintf(ofp, ">%s_%d_%d\n",ptParams->szOutFileStub,i, anT[i]);

      j = 4;
      while(aacQualities[i][j] != (char) NOT_SET){
	fprintf(ofp, "%d ", (int) aacQualities[i][j]);
	j++;
      }
      fprintf(ofp,"\n");
    }
  }

  free(anT);
  fclose(ofp);
  free(szClustFile);
}

void writeQualitiesD3(t_Params *ptParams, int nN, int nM, int nK, int *anZ, double** aadQualities)
{
  char *szClustFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;
  int i = 0, j = 0, s = 0;
  int*  anT = (int * ) malloc(nK*sizeof(int));

  sprintf(szClustFile, "%s_cd3%s", ptParams->szOutFileStub, QUAL_SUFFIX);
  ofp = fopen(szClustFile, "w");
  
  for(i = 0; i < nK; i++){
    anT[i] = 0;
  }
  for(i = 0; i < nN; i++){
    anT[anZ[i]] ++;
  }

  for(i = 0; i < nK; i++){
    if(anT[i] > 0){
      fprintf(ofp, ">%s_%d_%d\n",ptParams->szOutFileStub,i, anT[i]);

      j = 0;
      while(aadQualities[i][j] != NOT_SET){
	fprintf(ofp, "%f ", aadQualities[i][j]);
	j++;
      }
      fprintf(ofp,"\n");
    }
  }

  free(anT);
  fclose(ofp);
  free(szClustFile);
}

void writeMaster(t_Master *ptMaster, t_Unique *ptUnique, int* anCentroids, t_Flows *ptFlows, t_Params *ptParams)
{
  char *szMasterFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;
  int i = 0, j = 0, l = 0, nK = ptMaster->nK, nM = ptMaster->nM;
  int *anN = ptMaster->anN, *anCN = ptMaster->anCN, *anI = ptMaster->anI, *anP = ptMaster->anP;
  double *adTau = ptMaster->adTau;
  double *adDist = ptMaster->adDist;
  int nI = 0, nIndex = 0;
  double dTau = 0.0, dDist = 0.0, dDistMean = 0.0;
  short *asU = ptUnique->asU;
  short asCentroid[nM]; 
  if(!szMasterFile)
    goto memoryError;

  sprintf(szMasterFile, "%s%s", ptParams->szOutFileStub, MASTER_SUFFIX);
  ofp = fopen(szMasterFile, "w");
  if(ofp){
    fprintf(ofp, "%d %d\n",nK, nM);
    
    for(i = 0; i < nK; i++){
      for(j = 0; j < nM; j++){
	asCentroid[j] = NOT_SET;
      }

      if(anCentroids[i] >= 0){
	for(j = 0; j < nM; j++){
	  asCentroid[j] = asU[anCentroids[i]*nM + j];
	}
      }
      /*calc mean dist in cluster*/
      dDistMean = 0.0;
      for(j = 0; j < anN[i]; j++){
	nIndex = anCN[i] + j;
	dTau   = adTau[anP[nIndex]];
	dDist  = adDist[anP[nIndex]];
	dDistMean += dDist*dTau;
      }

      
      fprintf(ofp, "%d %d %f\n",i, anN[i], dDistMean);
	
      for(j = 0; j < anN[i]; j++){
	nIndex = anCN[i] + j;
	dTau   = adTau[anP[nIndex]];
	dDist  = adDist[anP[nIndex]];
	nI     = anI[nIndex];
        
	fprintf(ofp, "%d-%.5f-%.3f ",nI, dTau, dDist);
      }
      fprintf(ofp, "\n");
      
      
      for(l = 0; l < nM; l++){
	fprintf(ofp, "%d %d ",l,asCentroid[l]);
	for(j = 0; j < anN[i]; j++){
	  nIndex = anCN[i] + j;
	  nI     = anI[nIndex];

	  if(l < ptFlows->anLengths[nI]){
	    fprintf(ofp, "%d ",ptFlows->asData[nI*nM + l]);
	  }
	  else{
	    fprintf(ofp, "-1 ");
	  }
	}
	fprintf(ofp,"\n");
      }
       
      //fprintf(ofp, "\n");
      
    }
    fclose(ofp);
  }
  else{
    fprintf(stderr, "Failed opening %s for writing\n",szMasterFile);
  }
  
  free(szMasterFile);
  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in writeMaster\n");
}

void writeMasterN(int nIter, t_Master *ptMaster, t_Unique *ptUnique, int* anCentroids, t_Flows *ptFlows, t_Params *ptParams)
{
  char *szMasterFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;
  int i = 0, j = 0, l = 0, nK = ptMaster->nK, nM = ptMaster->nM;
  int *anN = ptMaster->anN, *anCN = ptMaster->anCN, *anI = ptMaster->anI, *anP = ptMaster->anP;
  double *adTau = ptMaster->adTau;
  double *adDist = ptMaster->adDist;
  int nI = 0, nIndex = 0;
  double dTau = 0.0, dDist = 0.0, dDistMean = 0.0;
  short *asU = ptUnique->asU;
  short asCentroid[nM]; 
  if(!szMasterFile)
    goto memoryError;

  sprintf(szMasterFile, "%s_%d%s", ptParams->szOutFileStub, nIter, MASTER_SUFFIX);
  ofp = fopen(szMasterFile, "w");
  if(ofp){
    fprintf(ofp, "%d %d\n",nK, nM);
    
    for(i = 0; i < nK; i++){
      for(j = 0; j < nM; j++){
	asCentroid[j] = NOT_SET;
      }

      if(anCentroids[i] >= 0){
	for(j = 0; j < nM; j++){
	  asCentroid[j] = asU[anCentroids[i]*nM + j];
	}
      }
      /*calc mean dist in cluster*/
      dDistMean = 0.0;
      for(j = 0; j < anN[i]; j++){
	nIndex = anCN[i] + j;
	dTau   = adTau[anP[nIndex]];
	dDist  = adDist[anP[nIndex]];
	dDistMean += dDist*dTau;
      }

      
      fprintf(ofp, "%d %d %f\n",i, anN[i], dDistMean);
	
      for(j = 0; j < anN[i]; j++){
	nIndex = anCN[i] + j;
	dTau   = adTau[anP[nIndex]];
	dDist  = adDist[anP[nIndex]];
	nI     = anI[nIndex];
        
	fprintf(ofp, "%d-%.5f-%.3f ",nI, dTau, dDist);
      }
      fprintf(ofp, "\n");
      
      
      for(l = 0; l < nM; l++){
	fprintf(ofp, "%d %d ",l,asCentroid[l]);
	for(j = 0; j < anN[i]; j++){
	  nIndex = anCN[i] + j;
	  nI     = anI[nIndex];

	  if(l < ptFlows->anLengths[nI]){
	    fprintf(ofp, "%d ",ptFlows->asData[nI*nM + l]);
	  }
	  else{
	    fprintf(ofp, "-1 ");
	  }
	}
	fprintf(ofp,"\n");
      }
       
      //fprintf(ofp, "\n");
      
    }
    fclose(ofp);
  }
  else{
    fprintf(stderr, "Failed opening %s for writing\n",szMasterFile);
  }
  
  free(szMasterFile);
  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in writeMaster\n");
}

void checkCentroidUniqueness(double *adWeight, int *anCentroids, int nK)
{
  int i = 0, j = 0;
  int anUnique[nK];
  int nUnique = 0;

  for(i = 0; i < nK; i++){
    anUnique[i] = TRUE;
    
    if(anCentroids[i] == NOT_SET || adWeight[i] < MIN_WEIGHT){
      anUnique[i] = NOT_SET;
    }
  }

  for(i = 0; i < nK; i++){
    
    if(anUnique[i] == TRUE){

      for(j = i + 1; j < nK; j++){
      
	if(anUnique[j] == TRUE){

	  if(anCentroids[j] == anCentroids[i]){
	    anUnique[j] = FALSE;
	    anCentroids[j] = NOT_SET;
	    
	    adWeight[i] += adWeight[j];
	    adWeight[j] = 0.0;
	  }
	}
      }
    }
  }

  return;
}

double calcNewWeights(int nK, double *adWeight, t_Master *ptMaster)
{
  int k = 0, i = 0, j = 0, nM = ptMaster->nM;
  double dMaxChange = 0.0;
  int* anI = ptMaster->anI, *anCN = ptMaster->anCN, *anP = ptMaster->anP, *anN = ptMaster->anN;
  double *adTau = ptMaster->adTau;

  /*calc weights*/
  for(i = 0; i < nK; i++){
    double dChange = adWeight[i];

    adWeight[i] = 0.0;

    for(k = 0; k < anN[i]; k++){
      int    nIndex = anCN[i] + k, nI = anI[nIndex];
      double dTau   = adTau[anP[nIndex]];

      adWeight[i] += dTau;
    }
    
    printf("%5.2f ",adWeight[i]);

    dChange = fabs(adWeight[i] - dChange);
    
    if(dChange > dMaxChange){
      dMaxChange = dChange;
    }
  }

  
  return dMaxChange;
}

void calcDistancesSlave(double* adDist, int nM, int nK, int* anChange, double *adLocalNorm, t_Letter *ptLetter, double dSigma,
			 int* anCentroids, t_Flows *ptFlows, double* adWeight, int nIStart, int nIFinish, t_Unique *ptUnique, float* afDistX)
{
  int k = 0, i = 0, nU = ptUnique->nU;
  double adTau[nK];
  int   *anLengths = ptFlows->anLengths;
  short *asData    = ptFlows->asData;
  short *asU       = ptUnique->asU;

  ptLetter->nTotal = 0;      
  for(i = nIStart; i < nIFinish; i++){
    double dOffset = MAX_DBL;
    int    nO = (i - nIStart)*nK;
	
    adLocalNorm[i - nIStart] = 0.0;

    for(k = 0; k < nK; k++){
      if(adWeight[k] > MIN_WEIGHT && anChange[k] == TRUE){
	adDist[nO + k] = afDistX[i*nU + anCentroids[k]];
      }
      
      if(adWeight[k] > MIN_WEIGHT && adDist[nO + k] < dOffset){
	dOffset = adDist[nO + k];
      }
    }

    for(k = 0; k < nK; k++){
      if(adWeight[k] > MIN_WEIGHT){
	adTau[k] = exp(dSigma*(-adDist[nO + k] + dOffset))*adWeight[k];
	adLocalNorm[i - nIStart] += adTau[k];
      }
      else{
	adTau[k] = 0.0;
      }
    }

    for(k = 0; k < nK; k++){
      adTau[k] /= adLocalNorm[i - nIStart];
	  
      if(adTau[k] > MIN_TAU){
	int nAdd = ptLetter->nTotal;
	    
	ptLetter->nTotal++;

	reallocateLetter(ptLetter);
	    
	ptLetter->anI[nAdd] = i;
	ptLetter->anK[nAdd] = k;
	ptLetter->adT[nAdd] = adTau[k];
	ptLetter->adD[nAdd] = adDist[nO + k];
      }
    }

  } /*loop i*/
}


void calcDistancesMaster(double *adDist, int* anChange, double *adNorm, t_Master *ptMaster, double dSigma,
			 int* anCentroids, t_Flows *ptFlows, double* adWeight, int nI0, t_Unique *ptUnique, float* afDistX)
{
  int k = 0, i = 0, nK = ptMaster->nK, nM = ptMaster->nM, nU = ptUnique->nU;
  double adTau[nK];
  int *anLengths = ptFlows->anLengths;
  short *asData  = ptFlows->asData;
  short *asU     = ptUnique->asU;

  /*clear master*/
  ptMaster->nTotal = 0;
  for(i = 0; i < nK; i++){ptMaster->anN[i] = 0;}

  for(i = 0; i < nI0; i++){
    double dOffset = MAX_DBL;
	
    adNorm[i] = 0.0;

    for(k = 0; k < nK; k++){
      if(adWeight[k] > MIN_WEIGHT && anChange[k] == TRUE){
	adDist[i*nK + k] = afDistX[i*nU + anCentroids[k]];
      }

      if(adWeight[k] > MIN_WEIGHT && adDist[i*nK + k] < dOffset){
	dOffset = adDist[i*nK + k];
      }
    }

    for(k = 0; k < nK; k++){
      if(adWeight[k] > MIN_WEIGHT){
	adTau[k] = exp(dSigma*(-adDist[i*nK + k] + dOffset))*adWeight[k];
	adNorm[i] += adTau[k];
      }
      else{
	adTau[k] = 0.0;
      }
    }

    for(k = 0; k < nK; k++){
      adTau[k] /= adNorm[i];
    }

    updateMasterI(i, ptMaster, adTau, &adDist[i*nK]);
  }
}

void receiveLetters(int nTag, double *adNorm, int numtasks, t_Master *ptMaster, t_Letter *ptLetter, int nI0, int nI)
{
  int n = 0; 
  MPI_Status status;

  for(n = 1; n < numtasks; n++){
    int nStart = nI0 + (n - 1)*nI;	

    MPI_Recv(&ptLetter->nTotal, 1, MPI_INT, n, nTag, MPI_COMM_WORLD, &status);

    minimiseLetter(ptLetter);

    MPI_Recv(ptLetter->anK, ptLetter->nTotal, MPI_INT, n, nTag, MPI_COMM_WORLD, &status);

    MPI_Recv(ptLetter->anI, ptLetter->nTotal, MPI_INT, n, nTag, MPI_COMM_WORLD, &status);

    MPI_Recv(ptLetter->adT, ptLetter->nTotal, MPI_DOUBLE, n, nTag, MPI_COMM_WORLD, &status);

    MPI_Recv(ptLetter->adD, ptLetter->nTotal, MPI_DOUBLE, n, nTag, MPI_COMM_WORLD, &status);

    updateMasterLetter(ptMaster, ptLetter);

    MPI_Recv(&adNorm[nStart], nI, MPI_DOUBLE, n, nTag, MPI_COMM_WORLD, &status);
    
  }
}

double calcNLLikelihood(t_Master *ptMaster, int *pnKEff, double* adWeights, double dSigma)
{
  int i = 0, k = 0, nK = ptMaster->nK, nN = ptMaster->nN;
  double adP[nN];
  double *adDist = ptMaster->adDist;
  int    *anP = ptMaster->anP, *anN = ptMaster->anN, *anCN = ptMaster->anCN, *anI = ptMaster->anI; 
  double dNLL = 0.0;

  for(i = 0; i < nN; i++){
    adP[i] = 0.0;
  }

  (*pnKEff) = 0;
  for(i = 0; i < nK; i++){
    if(adWeights[i] > MIN_WEIGHT){
      (*pnKEff)++;
    }
  }

  for(i = 0; i < nK; i++){

      for(k = 0; k < anN[i]; k++){ //loop data points
	int    nIndex = anCN[i] + k;
	int    nI     = anI[nIndex];
        double dDist  = adDist[anP[nIndex]];

	adP[nI] += adWeights[i]*exp(-dDist*dSigma);
      }
  }

  for(i = 0; i < nN; i++){
    dNLL += -log(adP[i]);
  }
  return dNLL - (((double) nN)*log(dSigma));
}


void calcCentroidsAndQualitiesMaster2(int nK, t_Flows *ptFlows, t_Master *ptMaster, double*** paadQualities, short *asCentroids)
{
  int i = 0, j = 0, k = 0;
  short s = 0;
  double adS[MAX_S], *adTau = ptMaster->adTau;
  short  *asData = ptFlows->asData;
  int nM   = ptMaster->nM; 
  int   *anN = ptMaster->anN, *anCN = ptMaster->anCN, *anP = ptMaster->anP;
  int   **aanI = ptMaster->aanI;
  double **aadQualities = NULL;
  double dOverCallThresh = 1.0 - exp(-OVERCALL_THRESH);


  aadQualities = (double **) malloc(nK*sizeof(double *));
  if(!aadQualities)
    goto memoryError;

  for(i = 0; i < nK; i++){ //loop clusters
    int j = 0, base  = 0;

    aadQualities[i] = NULL;

    if(anN[i] > 0){
      aadQualities[i] = (double *) malloc(MAX_SEQUENCE_LENGTH*sizeof(double));
      if(!aadQualities[i])
	goto memoryError;

      for(j = 0; j < MAX_SEQUENCE_LENGTH; j++){
	aadQualities[i][j] = NOT_SET;
      }

      j = nOffSet;
      while(j < nM){//loop position
	double dSMax = MAX_DBL;
	short  sSMax = -1;
	double dCount = 0.0;
	
	/*set priors*/
	for(s = 0; s < MAX_S; s++){//loop signal strength
	  adS[s] = 0.0;
	}

	for(k = 0; k < anN[i]; k++){ //loop data points
	  int    nIndex  = anCN[i] + k;
	  double dTau    = adTau[anP[nIndex]];
	  int    nI      = aanI[i][k];
	  short  sFlow   = asData[nI*nM + j];
	   
	  dCount += dTau;

	  for(s = 0; s < MAX_S; s++){
	    adS[s] += dTau*adLookUp[s*BINS + sFlow];
	  }
	  
	 
	}//loop data points

	/*find best s*/

	for(s = 0; s < MAX_S; s++){//loop signal strength  
	  if(adS[s] < dSMax){
	    dSMax = adS[s];
	    sSMax = s;
	  }
	}//loop signal

      
	if(dCount > MIN_COUNT){ 
	  double dU = 0.0, dNorm = 0.0;
 
	  asCentroids[i*nM + j] = sSMax;
	 
	  for(s = 0; s < MAX_S; s++){
	    dNorm += exp(-(adS[s] - dSMax));
	  }

	  for(s = 1; s <= sSMax; s++){
	    double dTemp = 0.0;

	    dU += exp(-(adS[s - 1] - dSMax))/dNorm; 

	    aadQualities[i][base] = 1.0 - dU;
	    base++;
	  }


	  dU += exp(-(adS[s - 1] - dSMax))/dNorm; 
	  while(dU < dOverCallThresh && s < MAX_S){
	    double dTemp = 0.0;

	    aadQualities[i][base] = 1.0 - dU;
	    asCentroids[i*nM + j]++;
	    base++;
	    s++;
	    dU += exp(-(adS[s - 1] - dSMax))/dNorm; 
	  }

	}
	else{
	  asCentroids[i*nM + j] = 0;
	}
	j++;
      }//loop pos
    }
    else{
      for(j = 0; j < nM; j++){
	asCentroids[i*nM + j] = 0;
      }
    }

  }//loop centroid

  (*paadQualities) = aadQualities;
  return;
 memoryError:
  fprintf(stderr, "Failed allocating memory in calcCentroidsAndQualitiesMaster aborting...\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void calcCurr(short* asCurr, short* asFlows, int nM)
{
  int i = 0;

  for(i = 0; i < nM; i++){
    asCurr[i] =  floor(((double) asFlows[i] + 50.0)/100.0);
  }
}

void calcUnique(t_Unique *ptUnique, t_Flows *ptFlows)
{
  int nM = ptFlows->nM, nN = ptFlows->nN;
  int i = 0, j = 0, l = 0;
  short  asCurr[nM];
  short* asData = ptFlows->asData;

  allocateUnique(nN, nM, ptUnique);

  for(i = 0; i < nN; i++){
 
    calcCurr(asCurr, &(asData[i*nM]), nM);

    for(j = 0; j < ptUnique->nU; j++){
      int nOffset = j*nM;

      l = 0;
      while(l < nM && asCurr[l] == ptUnique->asU[nOffset + l]){
	l++;
      } 

      if(l == nM){
	ptUnique->anMap[i] = j;
	ptUnique->anWeights[j]++;
	break;
      }
    }

    if(j == ptUnique->nU){
      if(ptUnique->nU == ptUnique->nSize){
	reallocateUnique(ptUnique, ptUnique->nSize*2);
      }
	
      ptUnique->anLenU[ptUnique->nU] = ptFlows->anLengths[i];
 
      ptUnique->anWeights[ptUnique->nU] = 1;

      ptUnique->anMap[i] = ptUnique->nU;
	
      ptUnique->anF[ptUnique->nU] = i;

      for(l = 0; l < nM; l++){
	ptUnique->asU[nM*ptUnique->nU + l] = asCurr[l];
      }

      ptUnique->nU++;
    }
  }

  reallocateUnique(ptUnique, ptUnique->nU);

 
  return;

 memoryError:
  fprintf(stderr,"Failed allocating memory in calcUnique\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void allocateUnique(int nN, int nM, t_Unique *ptUnique)
{
  int i = 0;

  ptUnique->nN = nN;

  ptUnique->nM = nM;

  ptUnique->nU = 0;

  ptUnique->nSize = INIT_SIZE;

  ptUnique->nMSize = nM*ptUnique->nSize;

  ptUnique->anMap = (int *) malloc(ptUnique->nN*sizeof(int));
  if(!ptUnique->anMap)
    goto memoryError;

  for(i = 0; i < ptUnique->nN; i++){
    ptUnique->anMap[i] = NOT_SET;
  }

  ptUnique->asU = (short *) malloc(ptUnique->nMSize*sizeof(short));
  if(!ptUnique->asU)
    goto memoryError;

  for(i = 0; i < ptUnique->nMSize; i++){
    ptUnique->asU[i] = NOT_SET;
  }
  
  ptUnique->anWeights = (int *) malloc(ptUnique->nSize*sizeof(int));
  if(!ptUnique->anWeights)
    goto memoryError;

  for(i = 0; i < ptUnique->nSize; i++){
    ptUnique->anWeights[i] = 0;
  }

  ptUnique->anF = (int *) malloc(ptUnique->nSize*sizeof(int));
  if(!ptUnique->anF)
    goto memoryError;

  for(i = 0; i < ptUnique->nSize; i++){
    ptUnique->anF[i] = NOT_SET;
  }


  ptUnique->anLenU = (int *) malloc(ptUnique->nSize*sizeof(int));
  if(!ptUnique->anLenU)
    goto memoryError;

  for(i = 0; i < ptUnique->nSize; i++){
    ptUnique->anLenU[i] = NOT_SET;
  }

  return;
 memoryError:
  fprintf(stderr, "Failed allocating memory in allocateUnique\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void destroyUnique(t_Unique *ptUnique)
{
  int i = 0;

  ptUnique->nN = 0;

  ptUnique->nM = 0;

  ptUnique->nU = 0;

  ptUnique->nSize = 0;

  free(ptUnique->anMap);

  free(ptUnique->asU);
    
  free(ptUnique->anWeights);

  free(ptUnique->anF);

  free(ptUnique->anLenU);
  return;
 memoryError:
  fprintf(stderr, "Failed allocating memory in allocateUnique\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void reallocateUnique(t_Unique *ptUnique, int nNewSize)
{
  int i = 0;
  int nOldSize  = ptUnique->nSize;
  int nOldMSize = ptUnique->nMSize;

  ptUnique->nSize = nNewSize;

  ptUnique->nMSize = ptUnique->nM*ptUnique->nSize;

  ptUnique->asU = (short *) realloc(ptUnique->asU, ptUnique->nMSize*sizeof(short));
  if(!ptUnique->asU)
    goto memoryError;

  for(i = nOldMSize; i < ptUnique->nMSize; i++){
    ptUnique->asU[i] = NOT_SET;
  }
  
  ptUnique->anWeights = (int *) realloc(ptUnique->anWeights, ptUnique->nSize*sizeof(int));
  if(!ptUnique->anWeights)
    goto memoryError;

  for(i = nOldSize; i < ptUnique->nSize; i++){
    ptUnique->anWeights[i] = 0;
  }

  ptUnique->anF = (int *) realloc(ptUnique->anF, ptUnique->nSize*sizeof(int));
  if(!ptUnique->anF)
    goto memoryError;

  for(i = nOldSize; i < ptUnique->nSize; i++){
    ptUnique->anF[i] = 0;
  }

  ptUnique->anLenU = (int *) realloc(ptUnique->anLenU, ptUnique->nSize*sizeof(int));
  if(!ptUnique->anLenU)
    goto memoryError;

  for(i = nOldSize; i < ptUnique->nSize; i++){
    ptUnique->anLenU[i] = NOT_SET;
  }

  return;
 memoryError:
  fprintf(stderr, "Failed allocating memory in reallocateUnique\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void calcDistX(float* afDistX, int iStart, int iFinish, t_Flows *ptFlows, t_Unique *ptUnique)
{
    int i = 0, j = 0, nN = ptFlows->nN, nU = ptUnique->nU, nM = ptUnique->nM;
    short* asU = ptUnique->asU;
    short* asData = ptFlows->asData;
    int*   anLengths = ptFlows->anLengths;
			
    for(i = iStart; i < iFinish; i++){
	for(j = 0; j < nU; j++){
      		afDistX[i*nU + j] = alignX(&asU[j*nM], &asData[i*nM], ptUnique->anLenU[j], anLengths[i]);;
	}	
    }
}




