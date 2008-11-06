#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>
#include <limits.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "../blib/RGBinary.h"
#include "../blib/RGIndex.h"
#include "../blib/RGReads.h"
#include "../blib/RGMatch.h"
#include "../blib/RGMatches.h"
#include "Definitions.h"
#include "ReadInputFiles.h"
#include "FindMatches.h"

/* TODO */
void FindMatches(
		int binaryOutput,
		char *rgFileName,
		char *rgIndexMainListFileName,
		char *rgIndexSecondaryListFileName,
		char *readFileName, 
		char *offsetsFileName,
		int space,
		int binaryInput,
		int startReadNum,
		int endReadNum,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions,
		int pairedEnd,
		int maxKeyMatches,
		int maxNumMatches,
		int forwardStrandOnly,
		int numThreads,
		char *outputID,
		char *outputDir,
		char *tmpDir,
		int timing
		)
{
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	char readsFilteredFileName[MAX_FILENAME_LENGTH]="\0";

	int numMainIndexes=0;
	char **mainIndexFileNames=NULL;

	int numSecondaryIndexes=0;
	char **secondaryIndexFileNames=NULL;

	int *offsets=NULL;
	int numOffsets=0;

	FILE *seqFP=NULL;
	FILE *seqFilteredFP=NULL;
	FILE **tempSeqFPs=NULL;
	char **tempSeqFileNames=NULL;
	FILE *outputFP=NULL;
	int i;

	int numMatches;
	int numReads;
	int numReadsFiltered;

	time_t startTime, endTime;
	int seconds, minutes, hours;
	int totalReadRGTime = 0;
	int totalDataStructureTime = 0; /* This will only give the to load and deleted the indexes (excludes searching and other things) */
	int totalSearchTime = 0; /* This will only give the time searching (excludes load times and other things) */
	int totalOutputTime = 0; /* This wll only give the total time to merge and output */

	RGMatches tempMatches;
	RGBinary rg;
	int startChr, startPos, endChr, endPos;

	/* Create output file name */
	sprintf(outputFileName, "%s%s.matches.file.%s.%d.%d.%d.%d.%d.%d.%d.%d.%d.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			space,
			startReadNum,
			endReadNum,
			numMismatches,
			numInsertions,
			numDeletions,
			numGapInsertions,
			numGapDeletions,
			pairedEnd,
			BFAST_MATCHES_FILE_EXTENSION);

	/* Read in the main RGIndex File Names */
	numMainIndexes=ReadFileNames(rgIndexMainListFileName, &mainIndexFileNames);
	if(numMainIndexes<=0) {
		PrintError("FindMatches",
				"numMainIndexes",
				"Read zero indexes",
				Exit,
				OutOfRange);
	}

	/* Read in the secondary RGIndex File Names */
	numSecondaryIndexes=ReadFileNames(rgIndexSecondaryListFileName, &secondaryIndexFileNames);

	/* Check the indexes.
	 * 1. We want the two sets of files to have the same range.
	 * */
	if(numSecondaryIndexes > 0) {
		CheckRGIndexes(mainIndexFileNames, 
				numMainIndexes,
				secondaryIndexFileNames,
				numSecondaryIndexes,
				binaryInput,
				&startChr,
				&startPos,
				&endChr,
				&endPos,
				space);
	}

	/* Read in the reference genome */
	startTime = time(NULL);
	RGBinaryReadBinary(&rg,
			rgFileName);
	assert(rg.space == space);
	endTime = time(NULL);
	totalReadRGTime = endTime - startTime;

	/* Read in the offsets */
	numOffsets=ReadOffsets(offsetsFileName, &offsets);

	/* Since we may be running through many indexes and only look for a small portion
	 * of the reads (see startReadNum and endReadNum), we copy the relevant reads
	 * to a temporary file, thereby elmininating the need to iterate through the 
	 * source read read file for each index. 
	 * */
	/* Create filtered reads file name */
	sprintf(readsFilteredFileName, "%s%s.reads.filtered.file.%s.%d.%d.%d.%d.%d.%d.%d.%d.%d.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			space,
			startReadNum,
			endReadNum,
			numMismatches,
			numInsertions,
			numDeletions,
			numGapInsertions,
			numGapDeletions,
			pairedEnd,
			BFAST_MATCHES_READS_FILTERED_FILE_EXTENSION);
	/* open read file */
	if((seqFP=fopen(readFileName, "r"))==0) {
		PrintError("FindMatches",
				readFileName,
				"Could not open readFileName for reading",
				Exit,
				OpenFileError);
	}
	/* open reads filtered file */
	if((seqFilteredFP=fopen(readsFilteredFileName, "w"))==0) {
		PrintError("FindMatches",
				readsFilteredFileName,
				"Could not open readsFilteredFileName for reading",
				Exit,
				OpenFileError);
	}
	/* Allocate memory for the temp file pointers - one for each thread */
	tempSeqFPs=malloc(sizeof(FILE*)*numThreads);
	if(NULL==tempSeqFPs) {
		PrintError("FindMatches",
				"tempSeqFPs",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	tempSeqFileNames=malloc(sizeof(char*)*numThreads);
	if(NULL==tempSeqFileNames) {
		PrintError("FindMatches",
				"tempSeqFileNames",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	for(i=0;i<numThreads;i++) {
		tempSeqFPs[i] = NULL;
		tempSeqFileNames[i] = NULL;
	}
	/* Read the reads to the thread temp files */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Reading %s into temp files.\n",
				readFileName);
	}
	/* This will close the reads file */
	WriteReadsToTempFile(seqFP,
			seqFilteredFP,
			&tempSeqFPs,
			&tempSeqFileNames,
			startReadNum,
			endReadNum,
			pairedEnd,
			numThreads,
			tmpDir,
			&numReads,
			&numReadsFiltered);
	/* Close the read file */
	fclose(seqFP);
	fclose(seqFilteredFP);
	if(VERBOSE >= 0) {
		fprintf(stderr, "Out of %d reads, will process %d reads and omit %d reads due to filtering.\n",
				numReads+numReadsFiltered,
				numReads,
				numReadsFiltered);
	}
	assert(numReads >= numThreads);


	/* Open output file */
	if((outputFP=fopen(outputFileName, "w"))==0) {
		PrintError("FindMatches",
				outputFileName,
				"Could not open outputFileName for writing",
				Exit,
				OpenFileError);
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Will output to %s.\n",
				outputFileName);
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Processing %d reads using %d main indexes.\n",
				numReads,
				numMainIndexes);
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Do step 1: search the main indexes for all reads */
	numMatches=FindMatchesInIndexes(mainIndexFileNames,
			binaryInput,
			&rg,
			numMainIndexes,
			offsets,
			numOffsets,
			space,
			numMismatches,
			numInsertions,
			numDeletions,
			numGapInsertions,
			numGapDeletions,
			pairedEnd,
			maxKeyMatches,
			maxNumMatches,
			forwardStrandOnly,
			numThreads,
			&tempSeqFPs,
			&tempSeqFileNames,
			outputFP,
			binaryOutput,
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
					binaryInput,
					&rg,
					numSecondaryIndexes,
					offsets,
					numOffsets,
					space,
					numMismatches,
					numInsertions,
					numDeletions,
					numGapInsertions,
					numGapDeletions,
					pairedEnd,
					maxKeyMatches,
					maxNumMatches,
					forwardStrandOnly,
					numThreads,
					&tempSeqFPs,
					&tempSeqFileNames,
					outputFP,
					binaryOutput,
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
				while(EOF!=GetNextRead(tempSeqFPs[i],
							&tempMatches,
							pairedEnd)) {
					/* Print the match to the output file */
					RGMatchesPrint(outputFP,
							&tempMatches,
							binaryOutput);
					/* Free the matches data structure */
					RGMatchesFree(&tempMatches);
				}

				/* Close the temp file */
				CloseTmpFile(&tempSeqFPs[i],
						&tempSeqFileNames[i]);
			}
			/* Close the output file */
			fclose(outputFP);
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
	}
	free(mainIndexFileNames);

	/* Free RGIndex file names */
	for(i=0;i<numSecondaryIndexes;i++) {
		free(secondaryIndexFileNames[i]);
	}
	free(secondaryIndexFileNames);

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
		int binaryInput,
		RGBinary *rg,
		int numIndexes,
		int *offsets,
		int numOffsets,
		int space,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions,
		int pairedEnd,
		int maxKeyMatches,
		int maxNumMatches,
		int forwardStrandOnly,
		int numThreads,
		FILE ***tempSeqFPs,
		char ***tempSeqFileNames,
		FILE *outputFP,
		int binaryOutput,
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
	FILE *tempOutputFP=NULL;
	char *tempOutputFileName=NULL;
	FILE **tempOutputIndexFPs=NULL;
	char **tempOutputIndexFileNames=NULL;
	int numWritten=0, numReads=0, numReadsFiltered=0;
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
	tempOutputIndexFPs = malloc(sizeof(FILE*)*numIndexes);
	if(NULL == tempOutputIndexFPs) {
		PrintError(FnName,
				"tempOutputIndexFPs",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	tempOutputIndexFileNames = malloc(sizeof(char*)*numIndexes);
	if(NULL == tempOutputIndexFileNames) {
		PrintError(FnName,
				"tempOutputIndexFileNames",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* If we are ending the search, output to the final output file.  Otherwise,
	 * output to a temporary file.
	 * */
	if(CopyForNextSearch == copyForNextSearch) {
		/* Open temporary file for the entire index search */
		tempOutputFP=OpenTmpFile(tmpDir, &tempOutputFileName);
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
			tempOutputIndexFPs[i] = OpenTmpFile(tmpDir, &tempOutputIndexFileNames[i]); 
		}
	}
	else {
		tempOutputIndexFPs[0] = tempOutputFP;
	}

	/* For each RGIndex, write temporary output */
	for(i=0;i<numIndexes;i++) { /* For each RGIndex */

		if(VERBOSE >= 0) {
			fprintf(stderr, "Searching index %d out of %d...\n", i+1, numIndexes);
		}
		numMatches = FindMatchesInIndex(indexFileNames[i],
				binaryInput,
				rg,
				offsets,
				numOffsets,
				space,
				numMismatches,
				numInsertions,
				numDeletions,
				numGapInsertions,
				numGapDeletions,
				pairedEnd,
				maxKeyMatches,
				maxNumMatches,
				forwardStrandOnly,
				numThreads,
				tempSeqFPs,
				tempOutputIndexFPs[i],
				binaryOutput,
				tmpDir,
				timing,
				totalDataStructureTime,
				totalSearchTime,
				totalOutputTime
					);
		if(VERBOSE >= 0) {
			fprintf(stderr, "Search for index %d out of %d complete.\n", i+1, numIndexes);
		}
	}

	/* Merge temporary output from each index and output to the output file. */
	if(numIndexes > 1) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "Merging the output from each index...\n");
		}

		startTime=time(NULL);
		/* Merge the temp index files into the all indexes file */
		numWritten=RGMatchesMergeFilesAndOutput(tempOutputIndexFPs,
				numIndexes,
				tempOutputFP,
				pairedEnd,
				binaryOutput,
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
			CloseTmpFile(&tempOutputIndexFPs[i],
					&tempOutputIndexFileNames[i]);
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
		assert(tempOutputFP != outputFP);
		numWritten=ReadTempReadsAndOutput(tempOutputFP,
				outputFP,
				tempSeqFP,
				pairedEnd,
				binaryOutput);
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
				NULL, /* We should have already filtered the reads */
				tempSeqFPs,
				tempSeqFileNames,
				0,
				0,
				pairedEnd,
				numThreads,
				tmpDir,
				&numReads,
				&numReadsFiltered);
		/* In this case, all the reads should be valid so we should apportion all reads */
		assert(numReadsFiltered == 0);
		assert(numReads == numWritten);

		/* Close the tempSeqFP */
		CloseTmpFile(&tempSeqFP, &tempSeqFileName);
		/* Close the temporary output file */
		CloseTmpFile(&tempOutputFP, &tempOutputFileName);
	}
	else {
		fclose(tempOutputFP);
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
		int binaryInput,
		RGBinary *rg,
		int *offsets,
		int numOffsets,
		int space,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions,
		int pairedEnd,
		int maxKeyMatches,
		int maxNumMatches,
		int forwardStrandOnly,
		int numThreads,
		FILE ***tempSeqFPs,
		FILE *indexFP,
		int binaryOutput,
		char *tmpDir,
		int timing,
		int *totalDataStructureTime,
		int *totalSearchTime,
		int *totalOutputTime)
{
	char *FnName = "FindMatchesInIndex";
	int i;
	FILE **tempOutputThreadFPs=NULL;
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
		PrintError(FnName,
				"threads",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory to pass data to threads */
	data=malloc(sizeof(ThreadIndexData)*numThreads);
	if(NULL==data) {
		PrintError(FnName,
				"data",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Allocate memory for one file pointer per thread */
	tempOutputThreadFPs=malloc(sizeof(FILE*)*numThreads); 
	if(NULL == tempOutputThreadFPs) {
		PrintError(FnName,
				"tempOutputThreadFPs",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	tempOutputThreadFileNames = malloc(sizeof(char*)*numThreads);
	if(NULL == tempOutputThreadFileNames) {
		PrintError("FindMatchesInThreades",
				"tempOutputThreadFileNames",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* If we have only one thread, output directly to the index fp */
	if(numThreads > 1) {
		/* Open files for thread output */
		for(i=0;i<numThreads;i++) {
			tempOutputThreadFPs[i] = OpenTmpFile(tmpDir, &tempOutputThreadFileNames[i]);
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
	ReadRGIndex(indexFileName, &index, binaryInput, space);
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
		data[i].binaryOutput = binaryOutput;
		data[i].index = &index;
		data[i].rg = rg;
		data[i].offsets = offsets;
		data[i].numOffsets = numOffsets;
		data[i].space = space;
		data[i].numMismatches = numMismatches;
		data[i].numInsertions = numInsertions;
		data[i].numDeletions = numDeletions;
		data[i].numGapInsertions = numGapInsertions;
		data[i].numGapDeletions = numGapDeletions;
		data[i].pairedEnd = pairedEnd;
		data[i].maxKeyMatches = maxKeyMatches;
		data[i].maxNumMatches = maxNumMatches;
		data[i].forwardStrandOnly = forwardStrandOnly;
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
			PrintError(FnName,
					"pthread_create: errCode",
					"Could not start thread",
					Exit,
					ThreadError);
		}
	}
	/* Wait for threads to return */
	for(i=0;i<numThreads;i++) {
		/* Wait for the given thread to return */
		errCode = pthread_join(threads[i],
				&status);
		/* Check the return code of the thread */
		if(0!=errCode) {
			PrintError(FnName,
					"pthread_join: errCode",
					"Thread returned an error",
					Exit,
					ThreadError);
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
		startTime = time(NULL);
		i=RGMatchesMergeThreadTempFilesIntoOutputTempFile(tempOutputThreadFPs,
				numThreads,
				indexFP,
				pairedEnd,
				binaryOutput);
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
			CloseTmpFile(&tempOutputThreadFPs[i],
					&tempOutputThreadFileNames[i]);
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
	/*
	   char *FnName="FindMatchesInIndex";
	   */
	RGMatches m;
	int numRead = 0;
	ThreadIndexData *data = (ThreadIndexData*)(arg);
	/* Function arguments */
	FILE *tempSeqFP = data->tempSeqFP;
	FILE *tempOutputFP = data->tempOutputFP;
	RGIndex *index = data->index;
	RGBinary *rg = data->rg;
	int *offsets = data->offsets;
	int numOffsets = data->numOffsets;
	int space = data->space;
	int numMismatches = data->numMismatches;
	int numInsertions = data->numInsertions;
	int numDeletions = data->numDeletions;
	int numGapInsertions = data->numGapInsertions;
	int numGapDeletions = data->numGapDeletions;
	int binaryOutput = data->binaryOutput;
	int pairedEnd = data->pairedEnd;
	int maxKeyMatches = data->maxKeyMatches;
	int maxNumMatches = data->maxNumMatches;
	int forwardStrandOnly = data->forwardStrandOnly;
	int threadID = data->threadID;
	data->numMatches = 0;

	/* Initialize match structures */
	RGMatchesInitialize(&m);

	/* For each read */
	if(VERBOSE >= 0) {
		fprintf(stderr, "\rthreadID:%d\tnumRead:[%d]",
				threadID,
				numRead);
	}
	while(EOF!=GetNextRead(tempSeqFP, 
				&m,
				pairedEnd)) {
		numRead++;

		if(VERBOSE >= 0 && numRead%FM_ROTATE_NUM==0) {
			fprintf(stderr, "\rthreadID:%d\tnumRead:[%d]",
					threadID,
					numRead);
		}

		/* Read one */
		RGReadsFindMatches(index,
				rg,
				&m.matchOne,
				offsets,
				numOffsets,
				space,
				numMismatches,
				numInsertions,
				numDeletions,
				numGapInsertions,
				numGapDeletions,
				maxKeyMatches,
				maxNumMatches,
				forwardStrandOnly);
		if(pairedEnd==1) {
			RGReadsFindMatches(index,
					rg,
					&m.matchTwo,
					offsets,
					numOffsets,
					space,
					numMismatches,
					numInsertions,
					numDeletions,
					numGapInsertions,
					numGapDeletions,
					maxKeyMatches,
					maxNumMatches,
				forwardStrandOnly);
		}

		if((0 < m.matchOne.numEntries && 1 != m.matchOne.maxReached ) ||
				(0 < m.matchTwo.numEntries && 1 != m.matchTwo.maxReached)) {
			data->numMatches++;
		}

		/* Output to file */
		RGMatchesPrint(tempOutputFP, 
				&m,
				binaryOutput); 

		/* Free matches */
		RGMatchesFree(&m);
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\rthreadID:%d\tnumRead:[%d]",
				threadID,
				numRead);
	}

	return NULL;
}
