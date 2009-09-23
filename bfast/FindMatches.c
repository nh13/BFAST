#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>
#include <limits.h>
#include <zlib.h>
#include "BLibDefinitions.h"
#include "BError.h"
#include "BLib.h"
#include "RGBinary.h"
#include "RGIndex.h"
#include "RGReads.h"
#include "RGMatch.h"
#include "RGMatches.h"
#include "MatchesReadInputFiles.h"
#include "FindMatches.h"

/* TODO */
void FindMatches(
		char *fastaFileName,
		char *mainIndexes,
		char *secondaryIndexes,
		char *readFileName, 
		char *offsetsInput,
		int space,
		int startReadNum,
		int endReadNum,
		int keySize,
		int maxKeyMatches,
		int maxNumMatches,
		int whichStrand,
		int numThreads,
		int queueLength,
		char *tmpDir,
		int timing
		)
{
	int numMainIndexes=0;
	char **mainIndexFileNames=NULL;
	int32_t **mainIndexIDs=NULL;

	int numSecondaryIndexes=0;
	char **secondaryIndexFileNames=NULL;
	int32_t **secondaryIndexIDs=NULL;

	int *offsets=NULL;
	int numOffsets=0;

	FILE *seqFP=NULL;
	FILE **tempSeqFPs=NULL;
	char **tempSeqFileNames=NULL;
	gzFile outputFP;
	int i;

	int numMatches;
	int numReads;

	time_t startTime, endTime;
	int seconds, minutes, hours;
	int totalReadRGTime = 0;
	int totalDataStructureTime = 0; /* This will only give the to load and deleted the indexes (excludes searching and other things) */
	int totalSearchTime = 0; /* This will only give the time searching (excludes load times and other things) */
	int totalOutputTime = 0; /* This wll only give the total time to merge and output */

	RGMatches tempMatches;
	RGBinary rg;
	int startChr, startPos, endChr, endPos;

	/* Read in the main RGIndex File Names */
	if(0<=VERBOSE) {
		fprintf(stderr, "Searching for main indexes...\n");
	}
	numMainIndexes=GetIndexFileNames(fastaFileName, space, mainIndexes, &mainIndexFileNames, &mainIndexIDs);
	if(numMainIndexes<=0) {
		PrintError("FindMatches", "numMainIndexes", "Read zero indexes", Exit, OutOfRange);
	}

	/* Read in the secondary RGIndex File Names */
	if(0<=VERBOSE) {
		fprintf(stderr, "%s", BREAK_LINE);
	}
	if(NULL != secondaryIndexes) {
		if(0<=VERBOSE) {
			fprintf(stderr, "Searching for secondary indexes...\n");
		}
		numSecondaryIndexes=GetIndexFileNames(fastaFileName, space, secondaryIndexes, &secondaryIndexFileNames, &secondaryIndexIDs);
	}
	else {
		if(0<=VERBOSE) {
			fprintf(stderr, "Not using secondary indexes.\n");
		}
		numSecondaryIndexes=0;
	}

	/* Check the indexes.
	 * 1. We want the two sets of files to have the same range.
	 * */
	if(numSecondaryIndexes > 0) {
		CheckRGIndexes(mainIndexFileNames, 
				numMainIndexes,
				secondaryIndexFileNames,
				numSecondaryIndexes,
				&startChr,
				&startPos,
				&endChr,
				&endPos,
				space);
	}

	/* Read in the reference genome */
	startTime = time(NULL);
	RGBinaryReadBinary(&rg,
			space,
			fastaFileName);
	assert(rg.space == space);
	endTime = time(NULL);
	totalReadRGTime = endTime - startTime;

	/* Read in the offsets */
	numOffsets = (NULL == offsetsInput) ? 0 : ReadOffsets(offsetsInput, &offsets);

	/* open read file */
	if(NULL == readFileName) {
		if(0 == (seqFP = fdopen(fileno(stdin), "r"))) {
			PrintError("FindMatches", "stdin", "Could not open stdin for reading", Exit, OpenFileError);
		}
	}
	else {
		if(0 == (seqFP = fopen(readFileName, "r"))) {
			PrintError("FindMatches", readFileName, "Could not open readFileName for reading", Exit, OpenFileError);
		}
	}
	/* Allocate memory for the temp file pointers - one for each thread */
	tempSeqFPs=malloc(sizeof(FILE*)*numThreads);
	if(NULL==tempSeqFPs) {
		PrintError("FindMatches", "tempSeqFPs", "Could not allocate memory", Exit, MallocMemory);
	}
	tempSeqFileNames=malloc(sizeof(char*)*numThreads);
	if(NULL==tempSeqFileNames) {
		PrintError("FindMatches", "tempSeqFileNames", "Could not allocate memory", Exit, MallocMemory);
	}
	for(i=0;i<numThreads;i++) {
		tempSeqFPs[i] = NULL;
		tempSeqFileNames[i] = NULL;
	}
	/* Read the reads to the thread temp files */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Reading %s into temp files.\n",
				(readFileName == NULL) ? "stdin" : readFileName);
	}
	/* This will close the reads file */
	WriteReadsToTempFile(seqFP,
			&tempSeqFPs,
			&tempSeqFileNames,
			startReadNum,
			endReadNum,
			numThreads,
			tmpDir,
			&numReads,
			space);
	/* Close the read file */
	fclose(seqFP);
	if(VERBOSE >= 0) {
		fprintf(stderr, "Will process %d reads.\n",
				numReads);
	}
	assert(numReads >= numThreads);

	/* Open output file */
	if(0 == (outputFP=gzdopen(fileno(stdout), "wb"))) {
		PrintError("FindMatches", "stdout", "Could not open stdout for writing", Exit, OpenFileError);
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Processing %d reads using %d main indexes.\n",
				numReads,
				numMainIndexes);
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Do step 1: search the main indexes for all reads */
	numMatches=FindMatchesInIndexes(mainIndexFileNames,
			mainIndexIDs,
			numMainIndexes,
			&rg,
			offsets,
			numOffsets,
			space,
			keySize,
			maxKeyMatches,
			maxNumMatches,
			whichStrand,
			numThreads,
			queueLength,
			&tempSeqFPs,
			&tempSeqFileNames,
			outputFP,
			(0 < numSecondaryIndexes)?CopyForNextSearch:EndSearch,
			MainIndexes,
			tmpDir,
			timing,
			&totalDataStructureTime,
			&totalSearchTime,
			&totalOutputTime
				);

	/* Do secondary index search */

	if(0 < numSecondaryIndexes) { /* Only if there are secondary indexes */
		if(0 < numReads - numMatches) { /* Only if enough reads are left */
			if(VERBOSE >= 0) {
				fprintf(stderr, "%s", BREAK_LINE);
				fprintf(stderr, "%s", BREAK_LINE);
				fprintf(stderr, "Processing remaining %d reads using %d secondary indexes.\n",
						numReads - numMatches,
						numSecondaryIndexes);
				fprintf(stderr, "%s", BREAK_LINE);
			}

			/* Do step 2: search the indexes for all reads */
			numMatches+=FindMatchesInIndexes(secondaryIndexFileNames,
					secondaryIndexIDs,
					numSecondaryIndexes,
					&rg,
					offsets,
					numOffsets,
					space,
					keySize,
					maxKeyMatches,
					maxNumMatches,
					whichStrand,
					numThreads,
					queueLength,
					&tempSeqFPs,
					&tempSeqFileNames,
					outputFP,
					EndSearch,
					SecondaryIndexes,
					tmpDir,
					timing,
					&totalDataStructureTime,
					&totalSearchTime,
					&totalOutputTime
						);
		}
		else {
			/* Output the reads not aligned and close the temporary read files */
			for(i=0;i<numThreads;i++) {
				/* Go to the beginning of the temp file */
				fseek(tempSeqFPs[i], 0, SEEK_SET);

				/* Initialize */
				RGMatchesInitialize(&tempMatches);

				/* Read in the reads */
				while(EOF!=GetRead(tempSeqFPs[i],
							&tempMatches,
							space)) {
					/* Print the match to the output file */
					RGMatchesPrint(outputFP,
							&tempMatches);
					/* Free the matches data structure */
					RGMatchesFree(&tempMatches);
				}

				/* Close the temp file */
				CloseTmpFile(&tempSeqFPs[i],
						&tempSeqFileNames[i]);
			}
			/* Close the output file */
			gzclose(outputFP);
		}
	}

	if(VERBOSE>=0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "In total, found matches for %d out of %d reads.\n", 
				numMatches,
				numReads);
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Free main RGIndex file names */
	for(i=0;i<numMainIndexes;i++) {
		free(mainIndexFileNames[i]);
		free(mainIndexIDs[i]);
	}
	free(mainIndexFileNames);
	free(mainIndexIDs);

	/* Free RGIndex file names */
	for(i=0;i<numSecondaryIndexes;i++) {
		free(secondaryIndexFileNames[i]);
		free(secondaryIndexIDs[i]);
	}
	free(secondaryIndexFileNames);
	free(secondaryIndexIDs);

	/* Free the offsets */
	free(offsets);

	/* Free reference genome */
	RGBinaryDelete(&rg);

	/* Free offsets */
	free(offsets);

	free(tempSeqFPs);
	free(tempSeqFileNames);

	/* Print timing */
	if(timing == 1) {
		/* Read RG time */
		seconds = totalReadRGTime;
		hours = seconds/3600;
		seconds -= hours*3600;
		minutes = seconds/60;
		seconds -= minutes*60;
		fprintf(stderr, "Total time loading the reference genome: %d hour, %d minutes and %d seconds.\n",
				hours,
				minutes,
				seconds);
		/* Data structure time */
		seconds = totalDataStructureTime;
		hours = seconds/3600;
		seconds -= hours*3600;
		minutes = seconds/60;
		seconds -= minutes*60;
		fprintf(stderr, "Total time loading and deleting indexes: %d hour, %d minutes and %d seconds.\n",
				hours,
				minutes,
				seconds);
		/* Search time */
		seconds = totalSearchTime;
		hours = seconds/3600;
		seconds -= hours*3600;
		minutes = seconds/60;
		seconds -= minutes*60;
		fprintf(stderr, "Total time searching indexes: %d hour, %d minutes and %d seconds.\n",
				hours,
				minutes,
				seconds);
		/* Output time */
		seconds = totalOutputTime;
		hours = seconds/3600;
		seconds -= hours*3600;
		minutes = seconds/60;
		seconds -= minutes*60;
		fprintf(stderr, "Total time merging and writing output: %d hour, %d minutes and %d seconds.\n",
				hours,
				minutes,
				seconds);
	}
}

int FindMatchesInIndexes(char **indexFileNames,
		int32_t **indexIDs,
		int numIndexes,
		RGBinary *rg,
		int *offsets,
		int numOffsets,
		int space,
		int keySize,
		int maxKeyMatches,
		int maxNumMatches,
		int whichStrand,
		int numThreads,
		int queueLength,
		FILE ***tempSeqFPs,
		char ***tempSeqFileNames,
		gzFile outputFP,
		int copyForNextSearch,
		int indexesType,
		char *tmpDir,
		int timing,
		int *totalDataStructureTime,
		int *totalSearchTime,
		int *totalOutputTime)
{
	char *FnName = "FindMatchesInIndexes";
	int i;
	gzFile tempOutputFP;
	char *tempOutputFileName=NULL;
	gzFile *tempOutputIndexFPs=NULL;
	char **tempOutputIndexFileNames=NULL;
	int numWritten=0, numReads=0;
	int numMatches = 0;
	time_t startTime, endTime;
	int seconds, minutes, hours;
	FILE *tempSeqFP=NULL;
	char *tempSeqFileName=NULL;

	/* IDEA: for each index, split search into threads generating one output file per thread.
	 * After the threads have finished their searches, merge their output into one output file
	 * specific for each index.  After all indexes have been searched, merge the index specific
	 * output.
	 * */

	/* Optimize if we have only one index, only one thread, or
	 * both one index and one thread */
	/* Allocate memory for the index specific file pointers */
	tempOutputIndexFPs = malloc(sizeof(gzFile)*numIndexes);
	if(NULL == tempOutputIndexFPs) {
		PrintError(FnName, "tempOutputIndexFPs", "Could not allocate memory", Exit, MallocMemory);
	}
	tempOutputIndexFileNames = malloc(sizeof(char*)*numIndexes);
	if(NULL == tempOutputIndexFileNames) {
		PrintError(FnName, "tempOutputIndexFileNames", "Could not allocate memory", Exit, MallocMemory);
	}
	/* If we are ending the search, output to the final output file.  Otherwise,
	 * output to a temporary file.
	 * */
	if(CopyForNextSearch == copyForNextSearch) {
		/* Open temporary file for the entire index search */
		tempOutputFP=OpenTmpGZFile(tmpDir, &tempOutputFileName);
	}
	else {
		assert(EndSearch == copyForNextSearch);
		/* Write directly to the output file */
		tempOutputFP=outputFP;
	}

	/* If we have only one index, output the temp output file */
	if(numIndexes > 1) {
		/* Open tmp files for each index */
		for(i=0;i<numIndexes;i++) {
			tempOutputIndexFPs[i] = OpenTmpGZFile(tmpDir, &tempOutputIndexFileNames[i]); 
		}
	}
	else {
		tempOutputIndexFPs[0] = tempOutputFP;
	}

	/* For each RGIndex, write temporary output */
	for(i=0;i<numIndexes;i++) { /* For each RGIndex */
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Searching index file %d/%d (index #%d, bin #%d)...\n", 
					i+1, numIndexes,
					indexIDs[i][0], indexIDs[i][1]);
		}
		numMatches = FindMatchesInIndex(indexFileNames[i],
				rg,
				offsets,
				numOffsets,
				space,
				keySize,
				maxKeyMatches,
				maxNumMatches,
				whichStrand,
				numThreads,
				queueLength,
				tempSeqFPs,
				tempOutputIndexFPs[i],
				tmpDir,
				timing,
				totalDataStructureTime,
				totalSearchTime,
				totalOutputTime
				);
		if(VERBOSE >= 0) {
			fprintf(stderr, "Searching index file %d/%d (index #%d, bin #%d) complete...\n", 
					i+1, numIndexes,
					indexIDs[i][0], indexIDs[i][1]);
		}
	}

	/* Merge temporary output from each index and output to the output file. */
	if(numIndexes > 1) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "Merging the output from each index...\n");
		}
		for(i=0;i<numIndexes;i++) {
			CloseTmpGZFile(&tempOutputIndexFPs[i], 
					&tempOutputIndexFileNames[i],
					0);
			if(!(tempOutputIndexFPs[i]=gzopen(tempOutputIndexFileNames[i], "rb"))) {
				PrintError(FnName, tempOutputIndexFileNames[i], "Could not re-open file for reading", Exit, OpenFileError);
			}
		}

		startTime=time(NULL);
		/* Merge the temp index files into the all indexes file */
		numWritten=RGMatchesMergeFilesAndOutput(tempOutputIndexFPs,
				numIndexes,
				tempOutputFP,
				maxNumMatches);
		endTime=time(NULL);
		if(VERBOSE >= 0 && timing == 1) {
			seconds = (int)(endTime - startTime);
			hours = seconds/3600;
			seconds -= hours*3600;
			minutes = seconds/60;
			seconds -= minutes*60;
			fprintf(stderr, "Merging matches from the indexes took: %d hours, %d minutes and %d seconds\n",
					hours,
					minutes,
					seconds);
		}
		(*totalOutputTime)+=endTime-startTime;

		/* If we had more than one index, then this is the total merged number of matches */
		numMatches = numWritten;

		/* Close the temporary index files */
		for(i=0;i<numIndexes;i++) {
			CloseTmpGZFile(&tempOutputIndexFPs[i],
					&tempOutputIndexFileNames[i],
					1);
		}
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "Found matches for %d reads.\n", numMatches);
	}


	/* Close the temporary read files */
	for(i=0;i<numThreads;i++) {
		/* Close temporary file */
		CloseTmpFile(&(*tempSeqFPs)[i], &(*tempSeqFileNames)[i]);
	}

	if(CopyForNextSearch == copyForNextSearch) {
		/* Go through the temporary output file and output those reads that have 
		 * at least one match to the final output file.  For those reads that have
		 * zero matches, output them to the temporary read file */

		if(VERBOSE >= 0) {
			fprintf(stderr, "Copying unmatched reads for secondary index search.\n");
		}

		/* Open a new temporary read file */
		tempSeqFP = OpenTmpFile(tmpDir, &tempSeqFileName);

		startTime=time(NULL);
		assert(tempOutputFP != outputFP); // this is very important
		numWritten=ReadTempReadsAndOutput(tempOutputFP,
				tempOutputFileName,
				outputFP,
				tempSeqFP);
		endTime=time(NULL);
		(*totalOutputTime)+=endTime-startTime;

		/* Move to the beginning of the read file */
		fseek(tempSeqFP, 0, SEEK_SET);

		if(VERBOSE >= 0) {
			fprintf(stderr, "Splitting unmatched reads into temp files.\n");
		}
		/* Now apportion the remaining reads into temp files for the threads when 
		 * searching the secondary indexes 
		 * */
		WriteReadsToTempFile(tempSeqFP,
				tempSeqFPs,
				tempSeqFileNames,
				0,
				0,
				numThreads,
				tmpDir,
				&numReads,
				space);
		/* In this case, all the reads should be valid so we should apportion all reads */
		assert(numReads == numWritten);

		/* Close the tempSeqFP */
		CloseTmpFile(&tempSeqFP, &tempSeqFileName);
		/* Close the temporary output file */
		CloseTmpGZFile(&tempOutputFP, &tempOutputFileName, 1);
	}
	else {
		gzclose(tempOutputFP);
	}

	/* Free memory for temporary file pointers */
	free(tempOutputIndexFPs);
	free(tempOutputIndexFileNames);

	if(VERBOSE >= 0) {
		if(MainIndexes == indexesType) {
			fprintf(stderr, "Searching main indexes complete.\n");
		}
		else {
			fprintf(stderr, "Searching secondary indexes complete.\n");
		}
	}

	return numMatches;
}

int FindMatchesInIndex(char *indexFileName,
		RGBinary *rg,
		int *offsets,
		int numOffsets,
		int space,
		int keySize,
		int maxKeyMatches,
		int maxNumMatches,
		int whichStrand,
		int numThreads,
		int queueLength,
		FILE ***tempSeqFPs,
		gzFile indexFP,
		char *tmpDir,
		int timing,
		int *totalDataStructureTime,
		int *totalSearchTime,
		int *totalOutputTime)
{
	char *FnName = "FindMatchesInIndex";
	int i, j;
	gzFile *tempOutputThreadFPs=NULL;
	char **tempOutputThreadFileNames=NULL;
	RGIndex index;
	int numMatches = 0;
	time_t startTime, endTime;
	int seconds, minutes, hours;
	int errCode;
	ThreadIndexData *data=NULL;
	pthread_t *threads=NULL;
	void *status;

	/* Allocate memory for threads */
	threads=malloc(sizeof(pthread_t)*numThreads);
	if(NULL==threads) {
		PrintError(FnName, "threads", "Could not allocate memory", Exit, MallocMemory);
	}
	/* Allocate memory to pass data to threads */
	data=malloc(sizeof(ThreadIndexData)*numThreads);
	if(NULL==data) {
		PrintError(FnName, "data", "Could not allocate memory", Exit, MallocMemory);
	}

	/* Allocate memory for one file pointer per thread */
	tempOutputThreadFPs=malloc(sizeof(gzFile)*numThreads); 
	if(NULL == tempOutputThreadFPs) {
		PrintError(FnName, "tempOutputThreadFPs", "Could not allocate memory", Exit, MallocMemory);
	}
	tempOutputThreadFileNames = malloc(sizeof(char*)*numThreads);
	if(NULL == tempOutputThreadFileNames) {
		PrintError("FindMatchesInThreades", "tempOutputThreadFileNames", "Could not allocate memory", Exit, MallocMemory);
	}
	/* If we have only one thread, output directly to the index fp */
	if(numThreads > 1) {
		/* Open files for thread output */
		for(i=0;i<numThreads;i++) {
			tempOutputThreadFPs[i] = OpenTmpGZFile(tmpDir, &tempOutputThreadFileNames[i]);
			assert(tempOutputThreadFPs[i] != NULL);
		}
	}
	else {
		/* Output directly to the indexe file pointer */
		tempOutputThreadFPs[0] = indexFP; 
	}

	/* Initialize index */
	RGIndexInitialize(&index);

	/* Read in the RG Index */
	startTime = time(NULL);
	ReadRGIndex(indexFileName, &index, space);
	/* Adjust if necessary */
	if(0 < keySize &&
			index.hashWidth <= keySize &&
			keySize < index.keysize) {
		/* Adjust key size and width */
		for(j=i=0;i < index.width && j < keySize;i++) {
			if(1 == index.mask[i]) {
				j++;
			}
		}
		assert(j == keySize);
		index.width = i;
		index.keysize = keySize;
	}
	endTime = time(NULL);
	(*totalDataStructureTime)+=endTime - startTime;	

	/* Set position to read from the beginning of the file */
	assert((*tempSeqFPs)!=NULL);
	for(i=0;i<numThreads;i++) {
		assert((*tempSeqFPs)!=NULL);
		fseek((*tempSeqFPs)[i], 0, SEEK_SET);
	}

	/* Execute */
	startTime = time(NULL);
	/* Initialize arguments to threads */
	for(i=0;i<numThreads;i++) {
		data[i].tempOutputFP = tempOutputThreadFPs[i];
		data[i].index = &index;
		data[i].rg = rg;
		data[i].offsets = offsets;
		data[i].numOffsets = numOffsets;
		data[i].space = space;
		data[i].maxKeyMatches = maxKeyMatches;
		data[i].maxNumMatches = maxNumMatches;
		data[i].queueLength = queueLength;
		data[i].whichStrand = whichStrand;
		data[i].threadID = i;
	}

	/* Open threads */
	for(i=0;i<numThreads;i++) {
		data[i].tempSeqFP = (*tempSeqFPs)[i];
		/* Start thread */
		errCode = pthread_create(&threads[i], /* thread struct */
				NULL, /* default thread attributes */
				FindMatchesInIndexThread, /* start routine */
				&data[i]); /* data to routine */
		if(0!=errCode) {
			PrintError(FnName, "pthread_create: errCode", "Could not start thread", Exit, ThreadError);
		}
	}
	/* Wait for threads to return */
	for(i=0;i<numThreads;i++) {
		/* Wait for the given thread to return */
		errCode = pthread_join(threads[i],
				&status);
		/* Check the return code of the thread */
		if(0!=errCode) {
			PrintError(FnName, "pthread_join: errCode", "Thread returned an error", Exit, ThreadError);
		}
		numMatches += data[i].numMatches;
	}
	endTime = time(NULL);
	(*totalSearchTime)+=endTime - startTime;
	if(VERBOSE >= 0) {
		fprintf(stderr, "\n");
	}

	/* Merge temp thread output into temp index output */
	/* Idea: the reads were apportioned in a given order,
	 * so merge to recover the initial order.
	 *
	 * Only do this if we have do not have one thread.
	 * */
	if(numThreads > 1) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "Merging thread temp files...\n");
		}
		/* Close files and open for reading only */
		for(i=0;i<numThreads;i++) {
			CloseTmpGZFile(&tempOutputThreadFPs[i], 
					&tempOutputThreadFileNames[i],
					0);
			if(!(tempOutputThreadFPs[i]=gzopen(tempOutputThreadFileNames[i], "rb"))) {
				PrintError(FnName, tempOutputThreadFileNames[i], "Could not re-open file for reading", Exit, OpenFileError);
			}
		}
		startTime = time(NULL);
		i=RGMatchesMergeThreadTempFilesIntoOutputTempFile(tempOutputThreadFPs,
				numThreads,
				indexFP);
		endTime = time(NULL);
		if(VERBOSE >= 0 && timing == 1) {
			seconds = (int)(endTime - startTime);
			hours = seconds/3600;
			seconds -= hours*3600;
			minutes = seconds/60;
			seconds -= minutes*60;
			fprintf(stderr, "Merging matches from threads took: %d hours, %d minutes and %d seconds\n",
					hours,
					minutes,
					seconds);
		}
		if(VERBOSE >= 0) {
			fprintf(stderr, "Merged %d reads from threads.\n", i);
		}

		/* Close temp thread output */
		for(i=0;i<numThreads;i++) {
			CloseTmpGZFile(&tempOutputThreadFPs[i],
					&tempOutputThreadFileNames[i],
					1);
		}
	}

	/* Free memory of the RGIndex */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Cleaning up index.\n");
	}
	startTime = time(NULL);
	RGIndexDelete(&index);
	endTime = time(NULL);
	(*totalDataStructureTime)+=endTime - startTime;	

	if(VERBOSE >= 0) {
		fprintf(stderr, "Found %d matches.\n", numMatches);
	}

	/* Free memory for temporary file pointers */
	free(tempOutputThreadFPs);
	free(tempOutputThreadFileNames);
	/* Free thread data */
	free(threads);
	free(data);

	return numMatches;
}

/* TODO */
void *FindMatchesInIndexThread(void *arg)
{
	char *FnName="FindMatchesInIndexThread";
	int32_t i, j;
	int numRead = 0;
	int foundMatch = 0;
	ThreadIndexData *data = (ThreadIndexData*)(arg);
	/* Function arguments */
	FILE *tempSeqFP = data->tempSeqFP;
	gzFile tempOutputFP = data->tempOutputFP;
	RGIndex *index = data->index;
	RGBinary *rg = data->rg;
	int *offsets = data->offsets;
	int numOffsets = data->numOffsets;
	int space = data->space;
	int maxKeyMatches = data->maxKeyMatches;
	int maxNumMatches = data->maxNumMatches;
	int queueLength = data->queueLength;
	int whichStrand = data->whichStrand;
	int threadID = data->threadID;
	data->numMatches = 0;

	RGMatches *matchQueue=NULL;
	int32_t matchQueueLength=queueLength;
	int32_t numMatches=0;

	/* Allocate match queue */
	matchQueue = malloc(sizeof(RGMatches)*matchQueueLength); 
	if(NULL == matchQueue) {
		PrintError(FnName, "matchQueue", "Could not allocate memory", Exit, MallocMemory);
	}

	/* Initialize match structures */
	for(i=0;i<matchQueueLength;i++) {
		RGMatchesInitialize(&matchQueue[i]);
	}

	/* For each read */
	if(VERBOSE >= 0) {
		fprintf(stderr, "\rthreadID:%d\tnumRead:[%d]",
				threadID,
				numRead);
	}
	while(0!=(numMatches = GetReads(tempSeqFP, matchQueue, matchQueueLength, space))) {
		for(i=0;i<numMatches;i++) {
			numRead++;
			if(VERBOSE >= 0 && numRead%FM_ROTATE_NUM==0) {
				fprintf(stderr, "\rthreadID:%d\tnumRead:[%d]",
						threadID,
						numRead);
			}
			/* Read */
			foundMatch = 0;
			for(j=0;j<matchQueue[i].numEnds;j++) {
				RGReadsFindMatches(index,
						rg,
						&matchQueue[i].ends[j],
						offsets,
						numOffsets,
						space,
						0,
						0,
						0,
						0,
						0,
						maxKeyMatches,
						maxNumMatches,
						whichStrand);
				if(0 < matchQueue[i].ends[j].numEntries && 1 != matchQueue[i].ends[j].maxReached) {
					foundMatch = 1;
				}
			}
			if(1 == foundMatch) {
				data->numMatches++;
			}
		}
		if(VERBOSE >= 0) {
			fprintf(stderr, "\rthreadID:%d\tnumRead:[%d]",
					threadID,
					numRead);
		}

		/* Output to file */
		for(i=0;i<numMatches;i++) {
			RGMatchesPrint(tempOutputFP, 
					&matchQueue[i]);

			/* Free matches */
			RGMatchesFree(&matchQueue[i]);
		}
	}
	assert(0 == numMatches);
	if(VERBOSE >= 0) {
		fprintf(stderr, "\rthreadID:%d\tnumRead:[%d]",
				threadID,
				numRead);
	}

	free(matchQueue);

	return arg;
}
