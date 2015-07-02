#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <mpi.h>

#include "SeqNoise.h"

/*Global Constants*/

static char szSequence[] = "ACGT";

static int bVerbose      = FALSE;

static int bLate         = FALSE;

static char *usage[] = {"SeqNoise - clusters sequences\n",
			"-in      string            sequence file name\n",
			"-din     string            distance matrix file name\n",
			"-out     string            cluster input file name\n",
			"-lin     string            list file\n",
			"Options:\n",
			"-min       mapping file\n",
			"-v       verbose\n",
			"-c       double            initial cut-off\n",
			"-s       double            precision\n",
			"-rin     string            lookup file name\n"};

static double* adLookUp = NULL;

static int  nLines = 11;

int main(int argc, char* argv[])
{
  int i = 0, j = 0, k = 0, n = 0, nN = 0, nM = 0, nSize = 0, nPacket = 0;
  t_Params tParams;
  t_Data   tData;
  /*number of centroids*/
  int      nK = 0;
  int      *anCentroids = NULL;
  float    *afDist      = NULL;
  int*     anChange     = NULL;
  int      numtasks, rank, rc; 
  MPI_Status   status;
  int      nA = 0, nA0 = 0, nTag = 1;
  int      nI = 0, nI0 = 0;
  int      bCont = TRUE;
  double*  adWeight = NULL;
  int*     anZ = NULL;
  t_Master tMaster;
  double   dSigma;
  char*    acSequences = NULL;
  int*     anLengths   = NULL;
  double*  adW         = NULL;

  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  /*get command line params*/
  getCommandLineParams(&tParams, argc, argv);

  if(rank == 0){//head node reads data etc.
    /*variables localised to head node*/
    double   *adP      = NULL;
    int      nIter = 0;
    double   dDelta = 0.0, dMaxDelta = MAX_DBL;

    t_Letter tLetter;

    initLookUp(&tParams);

    readData(&tData, &tParams);

    nN = tData.nSeq; nM = tData.nMaxLen;
    acSequences = tData.acSequences; anLengths = tData.anLen; adW = tData.adW; 

    /*broadcast data*/
    if(bVerbose){printf("%d: Broadcast data N = %d M = %d\n", rank, nN, nM);}
    
    /*broadcast flows*/
    if(bVerbose){printf("%d: Broadcast data flows\n", rank);}
    fflush(stdout);

    MPI_Bcast(adLookUp, N_BASES*N_BASES, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    broadcastData(&tData);

    anZ = (int *) malloc(nN*sizeof(int));
    if(!anZ)
      goto memoryError;
    
    /*read init file*/
    
    if(bVerbose){printf("%d: Read init file %s\n", rank, tParams.szInitFile);}
    readInitFile(&tParams, anZ);
   
    nK     = tParams.nK;

    if(bVerbose){printf("%d: Broadcast cluster number K = %d\n", rank, nK); fflush(stdout);}
    
    MPI_Bcast((void *) &nK, 1,MPI_INT, 0, MPI_COMM_WORLD);    

    dSigma = tParams.dSigma;

    MPI_Bcast((void *) &dSigma, 1,MPI_DOUBLE, 0, MPI_COMM_WORLD);

    nI = (int) (floor(nN / numtasks));

    nI0 = nI + (nN % numtasks);

    readDistanceMatrix(tParams.szDistFile, &afDist, nN);

    if(bVerbose){printf("%d: Broadcast distance matrix\n", rank, nSize, nPacket); fflush(stdout);}

    for(n = 1; n < numtasks; n++){
      MPI_Send(&afDist[nN*(nI0 + nI*(n-1))], nI*nN, MPI_FLOAT, n, nTag, MPI_COMM_WORLD);
    }

    adWeight = (double *) malloc(nK*sizeof(double));
    if(!adWeight)
      goto memoryError;

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
      anCentroids[i] = -1;
    }

    /*initialize alignments and tau*/
    if(bVerbose){printf("%d: Initialize alignments and tau\n",rank);fflush(stdout);}

    /*allocate memory for sparse array*/
    if(bVerbose){printf("%d: Allocate master\n",rank);fflush(stdout);}
    allocateMaster(&tMaster, nN, nK, INIT_MASTER_SIZE);

    allocateLetter(&tLetter, INIT_MASTER_SIZE);

    if(bVerbose){printf("%d: Init. alignment nI0 = %d nN = %d\n",rank, nI0, tMaster.nN);fflush(stdout);}
    initAlignment(&tMaster, &tData, nK, anZ, anChange, nI0);

    nA = (int) (floor(nK / numtasks));

    nA0 = nA + (nK % numtasks);

    if(bVerbose){
      printf("%d: Partition work %d %d %d %d\n",rank, nA0, nA, nI0, nI);
      fflush(stdout);
    }


    while((nIter < MIN_ITER) || ((dMaxDelta > MIN_DELTA) && (nIter < MAX_ITER))){
      /*Broadcast tau and alignments*/
      if(bVerbose){printf("%d-%d: Broadcast master data\n",rank,nIter);fflush(stdout);}
      MPI_Bcast(&bCont, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(nIter >= MIN_ITER){bLate = TRUE;}

      MPI_Bcast(&bLate, 1, MPI_INT, 0, MPI_COMM_WORLD);
      /*Broadcast master data*/
      fillMaster(&tMaster);

      printf("%d ",nIter);

      dMaxDelta = calcNewWeights(nK, adWeight, &tMaster, adW);
    
      printf("%f\n",dMaxDelta);

      /*calculate nK div numtasks centroids*/
      if(bVerbose){printf("%d-%d: Calc centroids\n",rank,nIter); fflush(stdout);}

      calcCentroidsMaster(anChange, 0, nK, &tMaster, anCentroids, adW, afDist);
    
      /*remove degenerate centroids*/
      checkCentroidUniqueness(adWeight, anCentroids, nK);

      /*broadcast centroids*/
      if(bVerbose){printf("%d-%d: Broadcast centroids\n",rank,nIter); fflush(stdout);}
      MPI_Bcast(anCentroids, nK, MPI_INT, 0, MPI_COMM_WORLD);
      if(bVerbose){printf("%d-%d: Broadcast weights\n",rank,nIter); fflush(stdout);}
      MPI_Bcast(adWeight, nK, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      /*broadcast change*/
      if(nIter == 0){
	for(i = 0; i < nK; i++){
	  anChange[i] = TRUE;
	}
      }

      MPI_Bcast(anChange, nK, MPI_INT, 0, MPI_COMM_WORLD);

      if(bVerbose){printf("%d-%d: Calc distances\n",rank,nIter);fflush(stdout);}

      //writeMaster(nIter, &tMaster, acCentroids, anCLengths, &tData, &tParams, adW);

      tMaster.nTotal = 0;
      for(i = 0; i < nK; i++){tMaster.anN[i] = 0;}

      for(i = 0; i < nI0; i++){
	double dOffset = MAX_DBL;
	double dNorm   = 0.0;
	double dMinWeight = MIN_WEIGHT;
	double adTau[nK];

	if(bLate == TRUE){dMinWeight = MIN_WEIGHT_LATE;} /*still need this?*/

	for(k = 0; k < nK; k++){

	  if(adWeight[k] > dMinWeight && afDist[i*nN + anCentroids[k]] < dOffset){
	    dOffset = afDist[i*nN + anCentroids[k]];
	  }
	}

	for(k = 0; k < nK; k++){
	  if(adWeight[k] > dMinWeight){
	    adTau[k] = exp(dSigma*(-afDist[i*nN + anCentroids[k]] + dOffset))*adWeight[k];
	    dNorm += adTau[k];
	  }
	  else{
	    adTau[k] = 0.0;
	  }
	}

	for(k = 0; k < nK; k++){
	  adTau[k] /= dNorm;
	}

	updateMasterI(i, &tMaster, adTau);
      }
      
      for(n = 1; n < numtasks; n++){
	
	MPI_Recv(&tLetter.nTotal, 1, MPI_INT, n, nTag, MPI_COMM_WORLD, &status);

	minimiseLetter(&tLetter);

	MPI_Recv(tLetter.anK, tLetter.nTotal, MPI_INT, n, nTag, MPI_COMM_WORLD, &status);

	MPI_Recv(tLetter.anI, tLetter.nTotal, MPI_INT, n, nTag, MPI_COMM_WORLD, &status);

	MPI_Recv(tLetter.adT, tLetter.nTotal, MPI_DOUBLE, n, nTag, MPI_COMM_WORLD, &status);

	updateMasterLetter(&tMaster, &tLetter);

      }
      
      if(bVerbose){
	printf("%d: Updated master total = %d, size = %d\n", rank, tMaster.nTotal, tMaster.nSize); fflush(stdout);
      }

      fflush(stdout);
      nIter++;
    }

    /*Fill master*/
    fillMaster(&tMaster);

    /*kill off slaves*/
    bCont = FALSE;
    MPI_Bcast((void *) &bCont, 1, MPI_INT, 0, MPI_COMM_WORLD);

    setZ(anZ, adP, &tMaster, &tParams);

    /*output clusters*/
    if(bVerbose){printf("%d: Output clusters\n",rank);fflush(stdout);}
    
    for(i = 0; i < nK; i++){anChange[i] = TRUE;}

    setMasterZ(&tMaster, nN, nK, anZ);

    calcCentroidsMaster(anChange, 0, nK, &tMaster, anCentroids, adW, afDist);

    if(1){
      int anT[nK];

      for(i = 0; i < nK; i++){
	anT[i] = 0;
      }
      for(i = 0; i < nN; i++){
	anT[anZ[i]] += (int) adW[i];
      }

      writeClustersD(&tParams, nN, nK, anT, anZ, &tData, anCentroids, adW);

      outputClusters(&tParams, nK, anT, anZ, nN, nM, &tData, anCentroids);  

      outputMap(&tParams, nK, anT, anZ, nN, nM, tData.aszID, adW, anCentroids);

      if(tParams.szMappingFile){
	outputMapping(&tParams, nK, anT, anZ, nN, adW, anCentroids, tData.aszID);
      }
    }

    /*free up memory*/
    destroyLetter(&tLetter);
    free(anZ); anZ = NULL;  
    free(adP);
    free(afDist);
    free(adWeight);
    free(anChange);

    destroyData(&tData);
  }
  else{
    int nIStart = 0, nIFinish = 0;
    t_Letter tLetter;

    /*recieve data*/
    if(bVerbose){printf("%d: Get lookup table\n",rank); fflush(stdout);}

    adLookUp = (double *) malloc(N_BASES*N_BASES*sizeof(double));
    if(!adLookUp)
      goto memoryError;

    MPI_Bcast(adLookUp, N_BASES*N_BASES, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(bVerbose){printf("%d: Receive data\n",rank); fflush(stdout);}

    receiveData(&tData);
    nN = tData.nSeq; nM = tData.nMaxLen;
    acSequences = tData.acSequences; anLengths = tData.anLen; adW = tData.adW; 

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

    afDist = (float *) malloc(nI*nN*sizeof(float));
    if(!afDist)
      goto memoryError;

    MPI_Recv(afDist, nI*nN, MPI_FLOAT, 0, nTag, MPI_COMM_WORLD, &status);    

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
      anCentroids[i] = -1;
    }

    /*--------------------------*/
    nA = (int) (floor(nK / numtasks));

    nA0 = nA + (nK % numtasks);

    nIStart  = nI0 + (rank - 1)*nI;
    nIFinish = nIStart + nI;

    if(bVerbose){printf("%d: Init. slave\n",rank);}
    
    allocateMaster(&tMaster, nN, nK, INIT_MASTER_SIZE);

    allocateLetter(&tLetter, INIT_MASTER_SIZE);

    fflush(stdout);
    while(TRUE){
      int nStart = nA0 + (rank - 1)*nA;
      int nFinish = nStart + nA;
      /*Receive tau and alignments*/
      if(bVerbose){printf("%d: Receive data\n",rank);}
      MPI_Bcast((void *) &bCont, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(bCont == FALSE)
	goto finish;

      MPI_Bcast((void *) &bLate, 1, MPI_INT, 0, MPI_COMM_WORLD);

      /*receive centroids*/
      MPI_Bcast(anCentroids, nK, MPI_INT, 0, MPI_COMM_WORLD);

      MPI_Bcast(adWeight, nK, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      /*receive change*/
      MPI_Bcast(anChange, nK, MPI_INT, 0, MPI_COMM_WORLD);

      if(bVerbose){printf("%d: Calc distances\n",rank);}

      tLetter.nTotal = 0;      
      for(i = nIStart; i < nIFinish; i++){
	double dOffset = MAX_DBL;
	double dNorm   = 0.0;
	int    nO      = (i - nIStart)*nN;
	double dMinWeight = MIN_WEIGHT;
	double adTau[nK];

	if(bLate == TRUE){dMinWeight = MIN_WEIGHT_LATE;}	


	for(k = 0; k < nK; k++){
	  
	  if(adWeight[k] > dMinWeight && afDist[nO + anCentroids[k]] < dOffset){
	    dOffset = afDist[nO + anCentroids[k]];
	  }
	}

	for(k = 0; k < nK; k++){
	  if(adWeight[k] > dMinWeight){
	    adTau[k] = exp(dSigma*(-afDist[nO + anCentroids[k]] + dOffset))*adWeight[k];
	    dNorm += adTau[k];
	  }
	  else{
	    adTau[k] = 0.0;
	  }
	}

	for(k = 0; k < nK; k++){
	  double dMinTau = MIN_TAU;

	  if(bLate == TRUE){dMinTau = MIN_TAU_LATE;}

	  adTau[k] /= dNorm;

	  if(adTau[k] > dMinTau){
	    int nAdd = tLetter.nTotal;
	    
	    tLetter.nTotal++;

	    reallocateLetter(&tLetter);
	    
	    tLetter.anI[nAdd] = i;
	    tLetter.anK[nAdd] = k;
	    tLetter.adT[nAdd] = adTau[k];
	  }
	}

      } /*loop i*/

      if(bVerbose){printf("%d: Sending letter size %d\n",rank,tLetter.nTotal);fflush(stdout);}

      minimiseLetter(&tLetter);

      MPI_Send(&tLetter.nTotal, 1, MPI_INT, 0, nTag, MPI_COMM_WORLD);
	
      MPI_Send(tLetter.anK, tLetter.nTotal, MPI_INT, 0, nTag, MPI_COMM_WORLD);

      MPI_Send(tLetter.anI, tLetter.nTotal, MPI_INT, 0, nTag, MPI_COMM_WORLD);

      MPI_Send(tLetter.adT, tLetter.nTotal, MPI_DOUBLE, 0, nTag, MPI_COMM_WORLD);

      fflush(stdout);
    }
  
    finish:
    
    destroyLetter(&tLetter);

    free(adW); adW = NULL;

    free(acSequences); acSequences = NULL;
    
    free(anLengths); anLengths = NULL;

    free(afDist); afDist = NULL;

    free(anChange); anChange = NULL;

    free(adWeight); adWeight = NULL;
  }

  destroyMaster(&tMaster);

  free(adLookUp); adLookUp = NULL;

  free(anCentroids); anCentroids = NULL;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  exit(EXIT_SUCCESS);
 memoryError:
  fprintf(stderr, "Failed allocating memory in main\n");
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
  char   *szBreak  = revstr(szID, '_'); 
  double dWeight   = 1.0;
  char   *pcError  = NULL;
  
 
  if(szBreak != NULL){
    dWeight = strtod(szBreak + 1, &pcError);
    if(*pcError != '\0'){
      fprintf(stderr, "Sequence weight format error %s %s\n",szID,szBreak + 1);
      fflush(stderr);
    }
  }

  return dWeight;
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
  
  ifp = fopen(ptParams->szDataFile, "r");

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
    fprintf(stderr, "Can't open input file %s\n", ptParams->szDataFile);
    exit(EXIT_FAILURE);
  }

  ptData->aszID        = (char **) malloc(ptData->nSeq*sizeof(char *));
  if(nPos > ptData->nMaxLen){
    ptData->nMaxLen = nPos;
  }  
  nM = ptData->nMaxLen;
  ptData->acSequences = (char *) malloc(ptData->nSeq*nM*sizeof(char));
  ptData->anLen       = (int *)  malloc(ptData->nSeq*sizeof(int));
  ptData->adW         = (double *) malloc(ptData->nSeq*sizeof(double));

  if(!ptData->acSequences)
    goto memoryError;

  if(!ptData->anLen)
    goto memoryError;

  if(!ptData->adW)
    goto memoryError;

  for(i = 0; i < ptData->nSeq; i++){
    ptData->anLen[i] = 0;
  }
  
  for(i = 0; i < ptData->nSeq*nM; i++){
    ptData->acSequences[i] = NOTSET_CHAR;
  }


  ifp = fopen(ptParams->szDataFile, "r");

  if(ifp){
    while(szRet = fgets(szLine, MAX_LINE_LENGTH, ifp)){
      if(szLine[0] == '>'){
	if(nSequences > 0){
	  ptData->anLen[nSequences - 1] = nPos;
	}

	szBrk = strpbrk(szLine, " \n");
	(*szBrk) = '\0';
	ptData->aszID[nSequences] = strdup(szLine + 1);
	ptData->adW[nSequences]   = getWeight(ptData->aszID[nSequences]);

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
    fprintf(stderr, "Can't open input file %s\n", ptParams->szDataFile);
    exit(EXIT_FAILURE);
  }

  return;

 memoryError:
  fprintf(stderr,"Failed allocating memory in readData\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}


void destroyData(t_Data *ptData)
{
  int i = 0;
  
  free(ptData->acSequences); ptData->acSequences = NULL;

  for(i = 0; i < ptData->nSeq; i++){
    free(ptData->aszID[i]);
    ptData->aszID[i] = NULL;
  }

  free(ptData->aszID); ptData->aszID = NULL;

  free(ptData->anLen); ptData->anLen = NULL;

  free(ptData->adW); ptData->adW = NULL;
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

void writeUsage(FILE* ofp)
{
  int i = 0;
  char *line;

  for(i = 0; i < nLines; i++){
    line = usage[i];
    fputs(line,ofp);
  }
}

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[])
{
  char *szTemp = NULL;
  char *cError = NULL;

  /*get parameter file name*/
  ptParams->szDataFile  = extractParameter(argc,argv, DATA_FILE, ALWAYS);  
  if(ptParams->szDataFile == NULL)
    goto error;

  ptParams->szDistFile  = extractParameter(argc,argv, DISTANCE_FILE, ALWAYS);  
  if(ptParams->szDistFile == NULL)
    goto error;

  /*get output filestub*/
  ptParams->szOutFileStub  = extractParameter(argc,argv, OUT_FILE_STUB, ALWAYS);  
  if(ptParams->szOutFileStub == NULL)
    goto error;

  ptParams->szInitFile  = extractParameter(argc,argv, INIT_FILE, ALWAYS);  
  if(ptParams->szInitFile == NULL)
    goto error;

  szTemp  = extractParameter(argc,argv, SIGMA, OPTION);  
  if(szTemp != NULL){
    ptParams->dSigma = strtod(szTemp, &cError);
    if(*cError != '\0')
      goto error;
  }
  else{
    ptParams->dSigma = DEF_SIGMA;
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

  szTemp  = extractParameter(argc,argv, MAPPING_FILE, OPTION);  
  if(szTemp != NULL){
    ptParams->szMappingFile = szTemp;
  }
  else{
    ptParams->szMappingFile = NULL;
  }

  if(extractParameter(argc,argv, VERBOSE, OPTION)){
    bVerbose = TRUE;
  }
  else{
    bVerbose = FALSE;
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
    fprintf(stderr, "Non standard base\n");
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
    fprintf(stderr, "Non standard base\n");
  }

  return adLookUp[nA*N_BASES + nB];
}

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


double needlemanWunsch(int nM, const char* acA, const char* acB, int nLenA, int nLenB)
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

  nComp = nLen;

  i = 0;
  while(acAlignA[i] == T_GAP || acAlignB[i] == T_GAP){
    i++;
  }
  
  nComp = nLen - i;
  
  dDist = aadFMatrix[nLenA][nLenB]/((double) nComp); //normalise by true length not with terminal gaps nComp -> nM

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


void allocateMaster(t_Master *ptMaster, int nN, int nK, int nSize)
{
  int i      = 0, j = 0;
  int nSSize = 0;

  ptMaster->nN = nN;
  ptMaster->nK = nK;

  ptMaster->nSize = nSize;
  ptMaster->nTotal = 0;

  ptMaster->adTau = (double *) malloc(sizeof(double)*nSize);
  if(!ptMaster->adTau)
    goto memoryError;
  
  for(i = 0; i < nSize; i++){
    ptMaster->adTau[i] = 0.0;
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
  fprintf(stderr, "Failed allocating memory in main\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void reallocateMaster(t_Master *ptMaster)
{
  int i = 0; 
  int nOldSize = ptMaster->nSize;

  while(ptMaster->nSize < ptMaster->nTotal){
    ptMaster->nSize  *= 2;
  }

  ptMaster->adTau = (double *) realloc(ptMaster->adTau, sizeof(double)*ptMaster->nSize);
  if(!ptMaster->adTau)
    goto memoryError;

  for(i = nOldSize; i < ptMaster->nSize; i++){
    ptMaster->adTau[i] = 0.0;
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
  fprintf(stderr, "Failed allocating memory in main\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void minimiseLetter(t_Letter *ptLetter)
{
  ptLetter->nSize  = ptLetter->nTotal;

  ptLetter->adT = (double *) realloc(ptLetter->adT, sizeof(double)*ptLetter->nSize);
  if(!ptLetter->adT)
    goto memoryError;

  ptLetter->anI = (int *) realloc(ptLetter->anI,sizeof(int)*ptLetter->nSize);
  if(!ptLetter->anI)
    goto memoryError;

  ptLetter->anK = (int *) realloc(ptLetter->anK,sizeof(int)*ptLetter->nSize);
  if(!ptLetter->anK)
    goto memoryError;

  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in main\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void minimiseMaster(t_Master *ptMaster)
{

  ptMaster->nSize  = ptMaster->nTotal;

  ptMaster->adTau = (double *) realloc(ptMaster->adTau, sizeof(double)*ptMaster->nSize);
  if(!ptMaster->adTau)
    goto memoryError;
  
  ptMaster->anP = (int *) realloc(ptMaster->anP,sizeof(int)*ptMaster->nSize);
  if(!ptMaster->anP)
    goto memoryError;

  ptMaster->anI = (int *) realloc(ptMaster->anI,sizeof(int)*ptMaster->nSize);
  if(!ptMaster->anI)
    goto memoryError;

  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in main\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void updateMasterI(int nI, t_Master *ptMaster, double *adT)
{
  int i = 0, j = 0, nK = ptMaster->nK;
  double dMinTau = MIN_TAU;

  if(bLate == TRUE){dMinTau = MIN_TAU_LATE;}

  for(i = 0; i < nK; i++){
    if(adT[i] > dMinTau){
      int nOldTotal = ptMaster->nTotal;

      ptMaster->nTotal++;

      reallocateMaster(ptMaster);

      ptMaster->adTau[nOldTotal] = adT[i];
      
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
    ptMaster->adTau[i] = ptLetter->adT[nAdd];
      
    ptMaster->aanP[nK][ptMaster->anN[nK]] = i;
      
    ptMaster->aanI[nK][ptMaster->anN[nK]] = nI;

    ptMaster->anN[nK]++;
    nAdd++;
  }
}

void initAlignment(t_Master *ptMaster, t_Data *ptData, int nK, int *anZ, int *anChange, int nI0)
{
/*initialize alignments and tau*/
  int nN = ptData->nSeq, nM = ptData->nMaxLen;
  int i = 0, j = 0, k = 0;
  double adTau[nK];

  for(i = 0; i < nN; i++){
    int nLen = 0;
    
    for(j = 0; j < nK; j++){
      adTau[j] = 0.0;
    }

    k = anZ[i];
    adTau[k] = 1.0;
    
    updateMasterI(i, ptMaster, adTau);
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

void broadcastData(t_Data *ptData)
{
  int i = 0;

  MPI_Bcast((void *) &ptData->nSeq, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast((void *) &ptData->nMaxLen, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast((void *) ptData->anLen, ptData->nSeq, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast((void *) ptData->acSequences, ptData->nSeq*ptData->nMaxLen, MPI_CHAR, 0, MPI_COMM_WORLD);

  MPI_Bcast((void *) ptData->adW, ptData->nSeq, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}

void broadcastMaster(t_Master *ptMaster)
{
  MPI_Bcast(&ptMaster->nTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(ptMaster->adTau, ptMaster->nTotal, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
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

void allocateLetter(t_Letter *ptLetter, int nSize)
{
  int i = 0;

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

  for(i = 0; i < nSize; i++){
    ptLetter->anI[i] = 0;
    ptLetter->anK[i] = 0;
    ptLetter->adT[i] = 0.0;
  }

  return;

memoryError:
  fprintf(stderr, "Failed allocating memory in allocateLetter\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void reallocateLetter(t_Letter *ptLetter)
{
  int i = 0; 
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

  for(i = nOldSize; i < ptLetter->nSize; i++){
    ptLetter->adT[i] = 0.0;
    ptLetter->anK[i] = 0;
    ptLetter->anI[i] = 0;
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

  free(ptLetter->anK); ptLetter->anK = NULL;
}

void calcCentroidsMaster(int *anChange, int nKStart, int nKFinish, t_Master *ptMaster, int *anCentroids, double* adW, float *afDist)
{
  int i = 0, j = 0, k = 0;
  double *adTau = ptMaster->adTau;
  int nN = ptMaster->nN; 
  int *anN   = ptMaster->anN, *anCN = ptMaster->anCN, *anP = ptMaster->anP, *anI = ptMaster->anI;
  int **aanI = ptMaster->aanI;
 
  for(i = nKStart; i < nKFinish; i++){ //loop clusters
    int nI = 0, nIndex = 0;  
    double dTau = 0.0, dCount = 0.0, dDistX = 0.0;
    double dMinF = MAX_DBL;
    int    nMinF = NOT_SET;

    anChange[i] = FALSE;

    for(j = 0; j < anN[i]; j++){ 
      nIndex   = anCN[i] + j;

      dCount += adW[anI[nIndex]]*adTau[anP[nIndex]];
    }

    if(anN[i] > 0 && dCount > MIN_COUNT){
      double *adF = (double *) malloc(anN[i]*sizeof(double));
      int    *anL = (int *) malloc(anN[i]*sizeof(int));

      if(!anL)
	goto memoryError;

      if(!adF)
	goto memoryError;

      for(k = 0; k < anN[i]; k++){
	nIndex   = anCN[i] + k;
	nI       = anI[nIndex];
	anL[k]   = nI;
	adF[k]   = 0.0;
      }

      for(j = 0; j < anN[i]; j++){ 
	nIndex   = anCN[i] + j;
	dTau     = adTau[anP[nIndex]];

	for(k = 0; k < anN[i]; k++){
	  dDistX =  afDist[anL[j]*nN + anL[k]]; 

	  adF[k] += dDistX*dTau*adW[anL[j]];
	}
      }

      for(k = 0; k < anN[i]; k++){
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
 
  ptData->adW = (double *) malloc(ptData->nSeq*sizeof(double));
  if(!ptData->adW)
    goto memoryError;

  MPI_Bcast((void *) ptData->adW, ptData->nSeq, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return;
 memoryError:
  fprintf(stderr, "Failed allocating memory in receiveData\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void writeZ(t_Params *ptParams, t_Data *ptData, int *anZ)
{
  char *szZFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;
  int i = 0;
  
  sprintf(szZFile, "%s%s", ptParams->szOutFileStub, Z_SUFFIX);
  
  ofp = fopen(szZFile, "w");
  
  for(i = 0; i < ptData->nSeq; i++){
    fprintf(ofp, "%s %d\n", ptData->aszID[i], anZ[i]);
  }

  fclose(ofp);
  free(szZFile);

  return;
}

void readInitFile(t_Params *ptParams, int* anZ)
{
  FILE *ifp = NULL;
  char* szLine = NULL;
  char *szTok, *pcError;
  int  i = 0;

  szLine = (char *) malloc(BIG_LINE_LENGTH*sizeof(char));
  if(!szLine){
	fprintf(stderr, "Failed allocating memory in readInitFile aborting\n");
	fflush(stderr);
	exit(EXIT_FAILURE);
  }	

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

  return;

 fileFormatError:
  fprintf(stderr, "Format error in %s\n", ptParams->szInitFile);
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void writeClustersF(t_Params *ptParams, int nK, double *adWeight, t_Data *ptData, int* anCentroids)
{
  char *szClustFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;
  int i = 0, j = 0, s = 0;
  int nM = ptData->nMaxLen;

  sprintf(szClustFile, "%s_cf%s", ptParams->szOutFileStub, FASTA_SUFFIX);
  ofp = fopen(szClustFile, "w");

  for(i = 0; i < nK; i++){
    int nI = anCentroids[i];

    if(adWeight[i] > 0){
      fprintf(ofp, ">%s_%d_%.2f\n",ptParams->szOutFileStub,i, adWeight[i]);

      for(j = 0; j < ptData->anLen[nI]; j++){
	char cBase = ptData->acSequences[nI*nM + j];
	if(cBase != NOTSET_CHAR){
	  fprintf(ofp, "%c", cBase);
	}
      }
      fprintf(ofp,"\n");
    }
  }

  fclose(ofp);
  free(szClustFile);
}

void writeClustersD(t_Params *ptParams, int nN, int nK, int* anT, int *anZ, t_Data *ptData, int* anCentroids, double* adW)
{
  char *szClustFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;
  int index =0, i = 0, j = 0, s = 0, nM = ptData->nMaxLen;

  sprintf(szClustFile, "%s_cd%s", ptParams->szOutFileStub, FASTA_SUFFIX);
  ofp = fopen(szClustFile, "w");

  for(i = 0; i < nK; i++){
    if(anT[i] > 0.0){
      int nI = anCentroids[i];

      fprintf(ofp, ">%s_%d_%d\n",ptParams->szOutFileStub,index, anT[i]);
      
      for(j = 0; j < ptData->anLen[nI]; j++){
	char cBase = ptData->acSequences[nI*nM + j];
	if(cBase != NOTSET_CHAR){
	  fprintf(ofp, "%c", cBase);
	}
      }
      index++;
      fprintf(ofp,"\n");
    }
  }

  
  fclose(ofp);
  free(szClustFile);
}

void outputClusters(t_Params *ptParams, int nK, int* anT, int *anZ, int nN, int nM, t_Data *ptData, int* anCentroids)
{
  FILE *ofp = NULL;
  char *szOutFile = (char *) malloc(sizeof(char)*MAX_LINE_LENGTH);
  int  index = 0, i = 0, j = 0;

  mkdir(ptParams->szOutFileStub, S_IRWXU);
  
  if(!szOutFile)
    goto memoryError;

  for(i = 0; i < nK; i++){
    if(anT[i] > 0){
      sprintf(szOutFile, "%s/i_%d%s",ptParams->szOutFileStub,index, FASTA_SUFFIX);
      ofp = fopen(szOutFile, "w");
      if(ofp){
	int nI = anCentroids[i];
	char szTemp[MAX_LINE_LENGTH];

	sprintf(szTemp, ">%s_%d_%d\n",ptParams->szOutFileStub,index, anT[i]);
      
	writeSequenceUnaligned(ofp, szTemp, &ptData->acSequences[nI*nM], ptData->anLen[nI]);	
	for(j = 0; j < nN; j++){
	  if(anZ[j] == i){
	    writeSequenceUnaligned(ofp,ptData->aszID[j], &ptData->acSequences[j*nM],ptData->anLen[j]);
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
  fprintf(stderr, "Failed allocating memory in main\n");
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

double calcNewWeights(int nK, double *adWeight, t_Master *ptMaster, double* adW)
{
  int k = 0, i = 0, j = 0;
  double dMaxChange = 0.0;
  int* anI = ptMaster->anI, *anCN = ptMaster->anCN, *anP = ptMaster->anP, *anN = ptMaster->anN;
  double *adTau = ptMaster->adTau;

  /*calc weights*/
  for(i = 0; i < nK; i++){
    double dChange = adWeight[i];

    adWeight[i] = 0.0;

    for(k = 0; k < anN[i]; k++){
      int    nIndex = anCN[i] + k, nI = anI[nIndex];
      double dTau   = adTau[anP[nIndex]], dW = adW[nI];

      adWeight[i] += dTau*dW;
    }
    
    printf("%5.2f ",adWeight[i]);

    dChange = fabs(adWeight[i] - dChange);
    
    if(dChange > dMaxChange){
      dMaxChange = dChange;
    }
  }

  
  return dMaxChange;
}

void initLookUp(t_Params *ptParams)
{
  int i = 0, j = 0;
  FILE *ifp = NULL;

  adLookUp = (double *) malloc(N_BASES*N_BASES*sizeof(double));
  if(!adLookUp)
    goto memoryError;

  ifp = fopen(ptParams->szLookUpFile,"r");
  
  if(bVerbose){
	printf("Using look up file %s\n",ptParams->szLookUpFile);
  }

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

void checkCentroidUniqueness(double *adWeight, int *anCentroids, int nK)
{
  int i = 0, j = 0, l = 0;
  int anUnique[nK];
  double dMinWeight = MIN_WEIGHT;

  if(bLate == TRUE){dMinWeight = MIN_WEIGHT_LATE;}

  for(i = 0; i < nK; i++){
    anUnique[i] = TRUE;
    
    if(adWeight[i] < dMinWeight){
      anUnique[i] = NOT_SET;
    }
  }

  for(i = 0; i < nK; i++){
    
    if(anUnique[i] == TRUE){

      for(j = i + 1; j < nK; j++){
      
	if(anUnique[j] == TRUE){

	  if(anCentroids[i] == anCentroids[j]){
	    anUnique[j] = FALSE;
	    adWeight[i] += adWeight[j];
	    adWeight[j] = 0.0;
	  }
	}
      }
    }
  }

}

void setMasterZ(t_Master *ptMaster, int nN, int nK, int* anZ)
{
  int i = 0, j = 0;

  ptMaster->nTotal = nN;
  
  for(i = 0; i < nK; i++){ptMaster->anN[i] = 0;}

  minimiseMaster(ptMaster);

  for(i = 0; i < nN; i++){
    int k = anZ[i];

    ptMaster->nTotal++;

    ptMaster->adTau[i] = 1.0;

    ptMaster->aanP[k][ptMaster->anN[k]] = i;
    
    ptMaster->aanI[k][ptMaster->anN[k]] = i;

    ptMaster->anN[k]++;
  }
    
  fillMaster(ptMaster);

  return;
}

void writeMaster(int nIter, t_Master *ptMaster, int* anCentroids, t_Data *ptData, t_Params *ptParams, double *adW)
{
  char *szMasterFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  FILE *ofp = NULL;
  int i = 0, j = 0, l = 0, k = 0, nK = ptMaster->nK;
  char* acSeq = ptData->acSequences;
  int *anLen = ptData->anLen, *anN = ptMaster->anN, *anCN = ptMaster->anCN;
  int *anI = ptMaster->anI, *anP = ptMaster->anP;
  double *adTau = ptMaster->adTau;
  int nM = ptData->nMaxLen, nMap = 0, nI = 0, nIndex = 0;
  double dTau = 0.0;
  if(!szMasterFile)
    goto memoryError;

  sprintf(szMasterFile, "%s_%d%s", ptParams->szOutFileStub, nIter, MASTER_SUFFIX);
  ofp = fopen(szMasterFile, "w");
  if(ofp){
    fprintf(ofp, "%d\n",nK);
    
    for(i = 0; i < nK; i++){
      int nMap = anCentroids[i];
      int nMax = -1;
      double dTotal = 0.0;

      if(nMap >= 0){
	fprintf(ofp, "%d %d %d %s ",i, anN[i], nMap, ptData->aszID[nMap]);
	
	for(j = 0; j < anN[i]; j++){
	  nIndex = anCN[i] + j;
	  dTau   = adTau[anP[nIndex]];
	  nI     = anI[nIndex];
	  if(anLen[nI] > nMax){
	    nMax = anLen[nI];
	  }
	  fprintf(ofp, "%d(%s)-%.5f ",nI, ptData->aszID[nI],dTau);

	  dTotal += adW[nI]*dTau;
	}
	fprintf(ofp, "%f\n", dTotal);

      }
    }

    fclose(ofp);
  }
  
  free(szMasterFile);
  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in writeMaster\n");
}

void outputMap(t_Params *ptParams, int nK, int* anT, int *anZ, int nN, int nM, char **aszLabels, double *adW, int* anCentroids)
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
	//fprintf(ofp, "%s_%d_%d",ptParams->szOutFileStub,index, anT[i]);

	fprintf(ofp, "%s", aszLabels[anCentroids[i]]);

	for(j = 0; j < nN; j++){
	  if(anZ[j] == i){
	    if(bFirst == FALSE){	
	    	fprintf(ofp, ",%s",aszLabels[j]);
	    }
	    else{
	    	fprintf(ofp, " %s",aszLabels[j]);
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

void outputMapping(t_Params *ptParams, int nK, int* anT, int *anZ, int nN, double *adW, int *anCentroids, char** aszID)
{
  FILE *ifp = NULL, *ofp = NULL;
  char *szOutFile = (char *) malloc(sizeof(char)*MAX_LINE_LENGTH);
  int  i = 0, j = 0, index = 0;
  char **aszMapping = (char **) malloc(sizeof(char*)*nN);
  char *szLine = (char *) malloc(sizeof(char)*BIG_LINE_LENGTH);
  char **aszCenters = (char **) malloc(sizeof(char*)*nN);

  if(!szOutFile)
    goto memoryError;

  if(!aszMapping)
    goto memoryError;

  if(!aszCenters)
    goto memoryError;

  ifp = fopen(ptParams->szMappingFile,"r");

  if(ifp){
    for(i = 0; i < nN; i++){
      char *szTemp = NULL, *szStart = NULL, *szEnd = NULL;

      fgets(szLine,BIG_LINE_LENGTH,ifp);

      szStart = strstr(szLine, " ");
      szEnd = strstr(szLine, "\n");

      (*szEnd) = '\0';

      szTemp = strdup(szStart + 1);

      aszMapping[i] = szTemp;
      (*szStart) ='\0';	
      szTemp = strdup(szLine);
      aszCenters[i] = szTemp;		
    }

    fclose(ifp);
  }
  else{
    fprintf(stderr, "Failed to open %s in outputMapping\n", ptParams->szMappingFile);

    fflush(stderr);

    return;
  }

  sprintf(szOutFile, "%s_cd%s",ptParams->szOutFileStub,MAP_SUFFIX);

  ofp = fopen(szOutFile, "w");

  if(ofp){
    for(i = 0; i < nK; i++){
      if(anT[i] > 0){
      	int bFirst = TRUE;

	//fprintf(ofp, "%s_%d_%d",ptParams->szOutFileStub,index, anT[i]);
	fprintf(ofp, "%s",aszCenters[anCentroids[i]]);
	for(j = 0; j < nN; j++){
	  if(anZ[j] == i){
	    if(bFirst == FALSE){	
	    	fprintf(ofp, ",%s",aszMapping[j]);
	    }
	    else{
	    	fprintf(ofp, " %s",aszMapping[j]);
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

  for(i = 0; i < nN; i++){
    free(aszCenters[i]);
    free(aszMapping[i]);
  }
  free(aszCenters);
  free(aszMapping);
  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in outputClusters\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void readDistanceMatrix(char *szDistFile, float **pafDist, int nN)
{
  FILE *ifp = NULL;
  char szLine[MAX_LINE_LENGTH];
  char *szTok = NULL, *pcError = NULL, *szBrk = NULL;
  float *afDist = NULL;
  int    i = 0, j = 0, nTest = 0;
  char **aszID;

  ifp = fopen(szDistFile, "r");

  if(ifp){
    fgets(szLine, MAX_LINE_LENGTH, ifp);
    
    szBrk = strpbrk(szLine, "\n");

    (*szBrk) = '\0';

    nTest = strtol(szLine, &pcError, 10);
    if(*pcError != '\0')
      goto fileFormatError;
    
    if(nN != nTest){
      fprintf(stderr,"Incorrectly sized distance file in SCluster\n");
      fflush(stderr);
      exit(EXIT_FAILURE);
    }

    afDist = (float *) malloc(nN*nN*sizeof(float));
    if(!afDist)
      goto memoryError;
    afDist[0] = 0.0;
    for(i = 1; i < nN; i++){

      afDist[i*nN + i] = 0.0;

      for(j = 0; j < i; j++){
	fgets(szLine, MAX_LINE_LENGTH, ifp);
    
	szBrk = strpbrk(szLine, "\n");
	(*szBrk) = '\0';

	afDist[i*nN + j] = (float) strtod(szLine, &pcError);
	afDist[j*nN + i] = afDist[i*nN + j];
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

  (*pafDist) = afDist;
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

void setZ(int *anZ, double *adP, t_Master* ptMaster, t_Params *ptParams)
{
  int i = 0, j = 0, k = 0;
  int nN = ptMaster->nN, nK = ptMaster->nK;
  double* adBest = (double *) malloc(nN*sizeof(double));
  int*    anBest = (int *)    malloc(nN*sizeof(int));
  double* adTau = ptMaster->adTau;
  int *anN = ptMaster->anN, *anCN = ptMaster->anCN, *anI = ptMaster->anI, *anP = ptMaster->anP;

  if(!adBest)
    goto memoryError;
		
  if(!anBest)
    goto memoryError;
	
  for(i = 0; i < nN; i++){
    adBest[i] = 0.0;
    anBest[i] = NOT_SET;
  }

  for(k = 0; k < nK; k++){ //loop clusters  

    for(i = 0; i < anN[k]; i++){ //loop data points
      int    nIndex = anCN[k] + i;
      double dTau   = adTau[anP[nIndex]];
      int    nI     = anI[nIndex];
  
      if(dTau > adBest[nI]){
	adBest[nI] = dTau;
	anBest[nI] = k;
      }
    }
  }  


  for(i = 0; i < nN; i++){
    double dMaxTau = -1.0;
    int    nMaxK   = -1;

    anZ[i] = anBest[i];
    adP[i] = 1.0 - adBest[i];
  }

  free(adBest); free(anBest);
  return;
  
 memoryError:
  fprintf(stderr,"Failed allocating memory in setZ\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}
