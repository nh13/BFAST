#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include "BLibDefinitions.h"
#include "BError.h"
#include "BLib.h"
#include "RGBinary.h"
#include "RGMatch.h" 
#include "RGMatches.h" 
#include "AlignedRead.h"
#include "AlignedEnd.h" 
#include "AlignedEntry.h" 
#include "ScoringMatrix.h"
#include "Align.h"
#include "RunLocalAlign.h"

/* TODO */
void RunAligner(char *fastaFileName,
		char *matchFileName,
		char *scoringMatrixFileName,
		int32_t ungapped,
		int32_t unconstrained,
		int32_t bestOnly,
		int32_t space,
		int32_t startReadNum,
		int32_t endReadNum,
		int32_t offsetLength,
		int32_t maxNumMatches,
		int32_t avgMismatchQuality,
		int32_t numThreads,
		int32_t queueLength,
		int32_t usePairedEndLength,
		int32_t pairedEndLength,
		int32_t mirroringType,
		int32_t forceMirroring,
		int32_t *totalReferenceGenomeTime,
		int32_t *totalAlignedTime,
		int32_t *totalFileHandlingTime)
{
	char *FnName = "RunAligner";
	gzFile outputFP=NULL;
	gzFile matchFP=NULL;
	int32_t startTime, endTime;
	RGBinary rg;

	startTime = time(NULL);
	RGBinaryReadBinary(&rg,
			NTSpace, // always NT space
			fastaFileName);
	endTime = time(NULL);
	/* Unpack */
	/*
	   RGBinaryUnPack(&rg);
	   */
	(*totalReferenceGenomeTime) = endTime - startTime;

	/* Check rg to make sure it is in NT Space */
	if(rg.space != NTSpace) {
		PrintError(FnName, "rg->space", "The reference genome must be in NT space", Exit, OutOfRange);
	}

	/* Open output file */
	if((outputFP=gzdopen(fileno(stdout), "wb"))==0) {
		PrintError(FnName, "stdout", "Could not open stdout file for writing", Exit, OpenFileError);
	}

	if(0 <= VERBOSE) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Reading match file from %s.\n",
				(NULL == matchFileName) ? "stdin" : matchFileName);
	}

	/* Open current match file */
	if(NULL == matchFileName) {
		if((matchFP=gzdopen(fileno(stdin), "rb"))==0) {
			PrintError(FnName, "stdin", "Could not open stdin for reading", Exit, OpenFileError);
		}
	}
	else {
		if((matchFP=gzopen(matchFileName, "rb"))==0) {
			PrintError(FnName, matchFileName, "Could not open file for reading", Exit, OpenFileError);
		}
	}

	RunDynamicProgramming(matchFP,
			&rg,
			scoringMatrixFileName,
			ungapped,
			unconstrained,
			bestOnly,
			space,
			startReadNum,
			endReadNum,
			offsetLength,
			maxNumMatches,
			avgMismatchQuality,
			numThreads,
			queueLength,
			usePairedEndLength,
			pairedEndLength,
			mirroringType,
			forceMirroring,
			outputFP,
			totalAlignedTime,
			totalFileHandlingTime);

	if(0 <= VERBOSE) {
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Close the match file */
	gzclose(matchFP);

	/* Close output file */
	gzclose(outputFP);

	/* Free the Reference Genome */
	RGBinaryDelete(&rg);
}

/* TODO */
void RunDynamicProgramming(gzFile matchFP,
		RGBinary *rg,
		char *scoringMatrixFileName,
		int32_t ungapped,
		int32_t unconstrained,
		int32_t bestOnly,
		int32_t space,
		int32_t startReadNum,
		int32_t endReadNum,
		int32_t offsetLength,
		int32_t maxNumMatches,
		int32_t avgMismatchQuality,
		int32_t numThreads,
		int32_t queueLength,
		int32_t usePairedEndLength,
		int32_t pairedEndLength,
		int32_t mirroringType,
		int32_t forceMirroring,
		gzFile outputFP,
		int32_t *totalAlignedTime,
		int32_t *totalFileHandlingTime)
{
	char *FnName="RunDynamicProgramming";
	/* local variables */
	ScoringMatrix sm;
	double mismatchScore;
	RGMatches m;
	int32_t i;

	int32_t numAligned=0;
	int32_t numNotAligned=0;
	int32_t startTime, endTime;
	int64_t numLocalAlignments=0;
	/* Thread specific data */
	ThreadData *data;
	pthread_t *threads=NULL;
	int32_t errCode;
	void *status;
	ThreadFileData fdata;

	/* Initialize */
	RGMatchesInitialize(&m);
	ScoringMatrixInitialize(&sm);

	/* Allocate memory for thread arguments */
	data = malloc(sizeof(ThreadData)*numThreads);
	if(NULL==data) {
		PrintError(FnName, "data", "Could not allocate memory", Exit, MallocMemory);
	}
	/* Allocate memory for thread ids */
	threads = malloc(sizeof(pthread_t)*numThreads);
	if(NULL==threads) {
		PrintError(FnName, "threads", "Could not allocate memory", Exit, MallocMemory);
	}

	/* Start file handling timer */
	startTime = time(NULL);

	/* Read in scoring matrix */
	if(NULL != scoringMatrixFileName) {
		ScoringMatrixRead(scoringMatrixFileName, &sm, space); 
	}
	/* End file handling timer */
	endTime = time(NULL);
	(*totalFileHandlingTime) += endTime - startTime;

	/* Calculate mismatch score */
	/* Assumes all match scores are the same and all substitution scores are the same */
	if(space == NTSpace) {
		mismatchScore = sm.ntMatch - sm.ntMismatch;
	}
	else {
		mismatchScore = sm.colorMatch - sm.colorMismatch;
	}

	/* Initialize fdata */
	// matchFP
	pthread_mutex_init(&fdata.matchFP_mutex, NULL);
	pthread_mutex_lock(&fdata.matchFP_mutex);
	fdata.matchFP = matchFP;
	fdata.matchFPctr = 1;
	pthread_mutex_unlock(&fdata.matchFP_mutex);
	// outputFP
	pthread_mutex_init(&fdata.outputFP_mutex, NULL);
	pthread_mutex_lock(&fdata.outputFP_mutex);
	fdata.outputFP = outputFP;
	fdata.outputFPctr = 0;
	pthread_mutex_unlock(&fdata.outputFP_mutex);
	// start/end read numbers
	fdata.startReadNum = startReadNum;
	fdata.endReadNum = endReadNum;

	// Skip matches
	pthread_mutex_lock(&fdata.matchFP_mutex);
	SkipMatches(&fdata);
	pthread_mutex_unlock(&fdata.matchFP_mutex);

	/* Create thread arguments */
	for(i=0;i<numThreads;i++) {
		data[i].fdata = &fdata;
		data[i].rg=rg;
		data[i].space=space;
		data[i].offsetLength=offsetLength;
		data[i].usePairedEndLength = usePairedEndLength;
		data[i].pairedEndLength = pairedEndLength;
		data[i].mirroringType = mirroringType;
		data[i].forceMirroring = forceMirroring;
		data[i].sm = &sm;
		data[i].ungapped = ungapped;
		data[i].unconstrained = unconstrained;
		data[i].bestOnly = bestOnly;
		data[i].numLocalAlignments = 0;
		data[i].avgMismatchQuality = avgMismatchQuality;
		data[i].mismatchScore = mismatchScore;
		data[i].queueLength = queueLength;
		data[i].threadID = i;
		data[i].numAligned = 0;
		data[i].numNotAligned = 0;
		data[i].fileTime = 0;
	}

	if(0 <= VERBOSE) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Performing alignment...\n");
		fprintf(stderr, "Currently on:\n0");
	}

	startTime = time(NULL);

	/* Create threads */
	for(i=0;i<numThreads;i++) {
		/* Start thread */
		errCode = pthread_create(&threads[i], /* thread struct */
				NULL, /* default thread attributes */
				RunDynamicProgrammingThread, /* start routine */
				&data[i]); /* data to routine */
		if(0!=errCode) {
			PrintError(FnName, "pthread_create: errCode", "Could not start thread", Exit, ThreadError);
		}
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Created threadID:%d\n",
					i);
		}
	}

	/* Wait for the threads to finish */
	for(i=0;i<numThreads;i++) {
		/* Wait for the given thread to return */
		errCode = pthread_join(threads[i],
				&status);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Thread returned with errCode:%d\n",
					errCode);
		}
		/* Check the return code of the thread */
		if(0!=errCode) {
			PrintError(FnName, "pthread_join: errCode", "Thread returned an error", Exit, ThreadError);
		}
	}

	if(0 <= VERBOSE) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Alignment complete.\n");
	}

	/* End align timer */
	endTime = time(NULL);
	(*totalAlignedTime) += endTime - startTime;

	/* Sum up statistics */
	for(i=0;i<numThreads;i++) {
		numAligned += data[i].numAligned;
		numNotAligned += data[i].numNotAligned;
		numLocalAlignments += data[i].numLocalAlignments;
		(*totalFileHandlingTime) += data[i].fileTime;
		(*totalAlignedTime) -= data[i].fileTime; // substract time reading and writing
	}

	if(VERBOSE >=0) {
		fprintf(stderr, "Performed %lld local alignments.\n", (long long int)numLocalAlignments);
		fprintf(stderr, "Outputted alignments for %d reads.\n", numAligned);
		fprintf(stderr, "Outputted %d reads for which there were no alignments.\n", numNotAligned); 
		fprintf(stderr, "Outputting complete.\n");
	}

	if(numAligned + numNotAligned != fdata.outputFPctr) {
		PrintError(FnName, "Inconsistent results", "numAligned + numNotAligned != fdata->outputFPctr", Exit, OutOfRange);
	}

	/* Free memory */
	free(data);
	free(threads);
}

/* TODO */
void *RunDynamicProgrammingThread(void *arg)
{
	/* Recover arguments */
	ThreadData *data = (ThreadData *)(arg);
	ThreadFileData *fdata = data->fdata;
	RGBinary *rg=data->rg;
	int32_t space=data->space;
	int32_t offsetLength=data->offsetLength;
	int32_t usePairedEndLength=data->usePairedEndLength;
	int32_t pairedEndLength=data->pairedEndLength;
	int32_t mirroringType=data->mirroringType;
	int32_t forceMirroring=data->forceMirroring;
	ScoringMatrix *sm = data->sm;
	int32_t ungapped=data->ungapped;
	int32_t unconstrained=data->unconstrained;
	int32_t bestOnly=data->bestOnly;
	int32_t threadID=data->threadID;
	int32_t avgMismatchQuality=data->avgMismatchQuality;
	double mismatchScore=data->mismatchScore;
	int32_t queueLength=data->queueLength;
	/* Local variables */
	char *FnName = "RunDynamicProgrammingThread";
	int32_t i, j, wasAligned;
	int32_t numAlignedRead=0;
	int32_t numMatches=0;
	int32_t ctr=0;
	int32_t startTime, endTime, readTime;

	AlignMatrix matrix;
	AlignedRead *alignedQueue=NULL;
	RGMatches *matchQueue=NULL;
	int32_t matchQueueLength=queueLength;
	int32_t numMatchesRead=0;

	AlignMatrixInitialize(&matrix);

	/* Allocate match queue */
	matchQueue = malloc(sizeof(RGMatches)*matchQueueLength);
	if(NULL == matchQueue) {
		PrintError(FnName, "matchQueue", "Could not allocate memory", Exit, MallocMemory);
	}
	alignedQueue = malloc(sizeof(AlignedRead)*matchQueueLength);
	if(NULL == alignedQueue) {
		PrintError(FnName, "alignedQueue", "Could not allocate memory", Exit, MallocMemory);
	}

	/* Initialize */

	/* Go through each read in the match file */
	while(0 != (numMatchesRead = GetMatches(fdata, matchQueue, matchQueueLength, &readTime))) {

		data->fileTime += readTime;

		for(i=0;i<numMatchesRead;i++) {
			AlignedReadInitialize(&alignedQueue[i]);

			if(0 <= VERBOSE && ctr%ALIGN_ROTATE_NUM==0) {
				fprintf(stderr, "\rthread:%d\t[%d]", threadID, ctr);
			}

			wasAligned=0;
			if(1 == IsValidMatch(&matchQueue[i])) {
				numMatches++;
				numAlignedRead = 0;

				/* Update the number of local alignments performed */
				data->numLocalAlignments += AlignRGMatches(&matchQueue[i],
						rg,
						&alignedQueue[i],
						space,
						offsetLength,
						sm,
						ungapped,
						unconstrained,
						bestOnly,
						usePairedEndLength,
						pairedEndLength,
						mirroringType,
						forceMirroring,
						&matrix);

				for(j=wasAligned=0;j<alignedQueue[i].numEnds;j++) {
					if(0 < alignedQueue[i].ends[j].numEntries) {
						wasAligned = 1;
					}
				}
			}

			if(1 == wasAligned) {
				/* Remove duplicates */
				AlignedReadRemoveDuplicates(&alignedQueue[i],
						AlignedEntrySortByAll);
				/* Updating mapping quality */
				AlignedReadUpdateMappingQuality(&alignedQueue[i], 
						mismatchScore, 
						avgMismatchQuality);
			}
			else {
				/* Copy over to alignedQueue[i] */
				AlignedReadAllocate(&alignedQueue[i],
						matchQueue[i].readName,
						matchQueue[i].numEnds,
						space);
				for(j=0;j<matchQueue[i].numEnds;j++) {
					AlignedEndAllocate(&alignedQueue[i].ends[j],
							matchQueue[i].ends[j].read,
							matchQueue[i].ends[j].qual,
							0);
				}
			}

			if(0 == wasAligned) {
				data->numNotAligned++;
			}
			else {
				data->numAligned++;
			}

			/* Free memory */
			RGMatchesFree(&matchQueue[i]);
			ctr++;
		}

		if(0 <= VERBOSE) {
			fprintf(stderr, "\rthread:%d\t[%d]", threadID, ctr);
		}

		/* Print */
		pthread_mutex_lock(&fdata->outputFP_mutex);
		startTime = time(NULL);
		for(i=0;i<numMatchesRead;i++) {
			AlignedReadPrint(&alignedQueue[i],
					fdata->outputFP);

			/* Free memory */
			AlignedReadFree(&alignedQueue[i]);
			fdata->outputFPctr++;
		}
		endTime = time(NULL);
		data->fileTime += endTime - startTime;
		pthread_mutex_unlock(&fdata->outputFP_mutex);
	}
	/* Free the matrix, free your mind */
	AlignMatrixFree(&matrix);

	if(0 <= VERBOSE) {
		fprintf(stderr, "\rthread:%d\t[%d]", threadID, ctr);
	}

	free(matchQueue);

	return arg;
}

int32_t GetMatches(ThreadFileData *fdata, RGMatches *m, int32_t maxToRead, int32_t *readTime) 
{
	int32_t numRead = 0;

	/* Get lock */
	pthread_mutex_lock(&fdata->matchFP_mutex);
	(*readTime) = time(NULL);

	assert(fdata->startReadNum <= fdata->matchFPctr);

	while(numRead < maxToRead && fdata->matchFPctr <= fdata->endReadNum) {
		RGMatchesInitialize(&(m[numRead]));
		if(EOF == RGMatchesRead(fdata->matchFP, &(m[numRead]))) {
			break;
		}
		numRead++;
		fdata->matchFPctr++;
	}
	(*readTime) = time(NULL) - (*readTime);
	pthread_mutex_unlock(&fdata->matchFP_mutex);
	return numRead;
}

void SkipMatches(ThreadFileData *fdata) 
{
	RGMatches m;

	if(fdata->startReadNum <= 1) {
		return;
	}

	if(0 <= VERBOSE) {
		fprintf(stderr, "Skipping matches...\nCurrently on:\n0");
	}

	RGMatchesInitialize(&m);
	while(fdata->matchFPctr < fdata->startReadNum && EOF != RGMatchesRead(fdata->matchFP, &m)) {
		if(0 <= VERBOSE && fdata->matchFPctr%ALIGN_SKIP_ROTATE_NUM==0) {
			fprintf(stderr, "\r%d", fdata->matchFPctr);
		}
		RGMatchesFree(&m);
		fdata->matchFPctr++;
	}
	if(0 <= VERBOSE) {
		fprintf(stderr, "\r%d\n", fdata->matchFPctr);
	}
}
