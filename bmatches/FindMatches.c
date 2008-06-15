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
#include "Definitions.h"
#include "ReadInputFiles.h"
#include "FindMatches.h"

/* TODO */
void FindMatches(char *outputFileName,
		int binaryOutput,
		char *rgFileName,
		char *rgIndexMainListFileName,
		char *rgIndexSecondaryListFileName,
		char *sequenceFileName, 
		char *offsetsFileName,
		int binaryInput,
		int startReadNum,
		int endReadNum,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions,
		int pairedEnd,
		int maxMatches,
		int numThreads,
		char *tmpDir,
		int timing
		)
{
	int numMainIndexes=0;
	char **mainIndexFileNames=NULL;

	int numSecondaryIndexes=0;
	char **secondaryIndexFileNames=NULL;

	int *offsets=NULL;
	int numOffsets=0;

	FILE *seqFP=NULL;
	FILE **tempSeqFPs=NULL;
	char **tempSeqFileNames=NULL;
	FILE *outputFP=NULL;
	int i;

	int numMatches;
	int numReads;

	int startTime, endTime;
	int seconds, minutes, hours;
	int totalReadRGTime = 0;
	int totalDataStructureTime = 0; /* This will only give the to load and deleted the indexes and trees (excludes searching and other things) */
	int totalSearchTime = 0; /* This will only give the time searching (excludes load times and other things) */
	int totalOutputTime = 0; /* This wll only give the total time to merge and output */

	RGBinary rg;
	int startChr, startPos, endChr, endPos;

	/* Read in the main RGIndex File Names */
	numMainIndexes=ReadFileNames(rgIndexMainListFileName, &mainIndexFileNames);
	if(numMainIndexes<=0) {
		PrintError("FindMatches",
				"numMainIndexes",
				"Read zero indexes",
				Exit,
				OutOfRange);
	}

	/* Read in the RGIndex File Names */
	numSecondaryIndexes=ReadFileNames(rgIndexSecondaryListFileName, &secondaryIndexFileNames);
	if(numSecondaryIndexes<=0) {
		PrintError("FindMatches",
				"numSecondaryIndexes",
				"Read zero indexes",
				Exit,
				OutOfRange);
	}

	/* Check the indexes.
	 * 1. We want the two sets of files to have the same range.
	 * */
	CheckRGIndexes(mainIndexFileNames, 
			numMainIndexes,
			secondaryIndexFileNames,
			numSecondaryIndexes,
			binaryInput,
			&startChr,
			&startPos,
			&endChr,
			&endPos);

	/* Read in the reference genome */
	startTime = time(NULL);
	RGBinaryReadBinary(&rg,
			rgFileName);
	endTime = time(NULL);
	totalReadRGTime = endTime - startTime;

	/* Read in the offsets */
	numOffsets=ReadOffsets(offsetsFileName, &offsets);

	/* Since we may be running through many indexes and only look for a small portion
	 * of the reads (see startReadNum and endReadNum), we copy the relevant reads
	 * to a temporary file, thereby elmininating the need to iterate through the 
	 * source sequence read file for each index. 
	 * */
	/* open sequence file */
	if((seqFP=fopen(sequenceFileName, "r"))==0) {
		PrintError("FindMatches",
				sequenceFileName,
				"Could not open sequenceFileName for reading",
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
	/* Read the sequences to the thread temp files */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Reading %s into temp files.\n",
				sequenceFileName);
	}
	/* This will close the sequences file */
	numReads=ReadSequencesToTempFile(seqFP,
			&tempSeqFPs,
			&tempSeqFileNames,
			startReadNum,
			endReadNum,
			pairedEnd,
			numThreads,
			tmpDir,
			timing);
	/* Close the sequence file */
	fclose(seqFP);
	assert(numReads >= numThreads);
	if(VERBOSE >= 0) {
		fprintf(stderr, "Read %d reads from %s.\n",
				numReads,
				sequenceFileName);
	}

	/* IDEA 
	 * 		Use temp files to store the results for each index.  Once we have one
	 * 		through each index, merge the results and output to file.  Store all
	 * 		sequences that had no match in a temp file to use when searching the
	 * 		trees.
	 *
	 * 		Use temp files to store the results for each tree.  Once we have gone
	 * 		through each tree, merge the results and append to the output file.
	 * 		*/

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
		fprintf(stderr, "Processing %d reads using %d main indexes.\n",
				numReads,
				numMainIndexes);
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Do step 1: search the main indexes for all sequences */
	numMatches=FindMatchesInIndexes(mainIndexFileNames,
			binaryInput,
			&rg,
			numMainIndexes,
			offsets,
			numOffsets,
			numMismatches,
			numInsertions,
			numDeletions,
			numGapInsertions,
			numGapDeletions,
			pairedEnd,
			maxMatches,
			numThreads,
			&tempSeqFPs,
			&tempSeqFileNames,
			outputFP,
			binaryOutput,
			1,
			tmpDir,
			timing,
			&totalDataStructureTime,
			&totalSearchTime,
			&totalOutputTime
				);

	if(numReads-numMatches > 0) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Processing remaining %d reads using %d secondary indexes.\n",
					numReads-numMatches,
					numSecondaryIndexes);
			fprintf(stderr, "%s", BREAK_LINE);
		}

		/* Do step 2: search the indexes for all sequences */
		numMatches+=FindMatchesInIndexes(secondaryIndexFileNames,
				binaryInput,
				&rg,
				numSecondaryIndexes,
				offsets,
				numOffsets,
				numMismatches,
				numInsertions,
				numDeletions,
				numGapInsertions,
				numGapDeletions,
				pairedEnd,
				INT_MAX,
				numThreads,
				&tempSeqFPs,
				&tempSeqFileNames,
				outputFP,
				binaryOutput,
				0,
				tmpDir,
				timing,
				&totalDataStructureTime,
				&totalSearchTime,
				&totalOutputTime
					);
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

int FindMatchesInIndexes(char **secondaryIndexFileNames,
		int binaryInput,
		RGBinary *rg,
		int numIndexes,
		int *offsets,
		int numOffsets,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions,
		int pairedEnd,
		int maxMatches,
		int numThreads,
		FILE ***tempSeqFPs,
		char ***tempSeqFileNames,
		FILE *outputFP,
		int binaryOutput,
		int MainIndexes,
		char *tmpDir,
		int timing,
		int *totalDataStructureTime,
		int *totalSearchTime,
		int *totalOutputTime)
{
	int i, j;
	FILE *tempOutputFP=NULL;
	char *tempOutputFileName=NULL;
	FILE **tempOutputIndexFPs=NULL;
	char **tempOutputIndexFileNames=NULL;
	FILE **tempOutputThreadFPs=NULL;
	char **tempOutputThreadFileNames=NULL;
	RGIndex index;
	int numWritten, numRead;
	int startTime;
	int endTime;
	int errCode;
	ThreadIndexData *data=NULL;
	pthread_t *threads=NULL;
	FILE *tempSeqFP=NULL;
	char *tempSeqFileName=NULL;
	void *status;

	/* IDEA: for each index, split search into threads generating one output file per thread.
	 * After the threads have finished their searches, merge their output into one output file
	 * specific for each index.  After all indexes have been searched, merge the index specific
	 * output.
	 * */

	/* Allocate memory for threads */
	threads=malloc(sizeof(pthread_t)*numThreads);
	if(NULL==threads) {
		PrintError("FindMatchesInIndexes",
				"threads",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory to pass data to threads */
	data=malloc(sizeof(ThreadIndexData)*numThreads);
	if(NULL==data) {
		PrintError("FindMatchesInIndexes",
				"data",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the index specific file pointers */
	tempOutputIndexFPs = (FILE**)malloc(sizeof(FILE*)*numIndexes);
	if(NULL == tempOutputIndexFPs) {
		PrintError("FindMatchesInIndexes",
				"tempOutputIndexFPs",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	tempOutputIndexFileNames = malloc(sizeof(char*)*numIndexes);
	if(NULL == tempOutputIndexFileNames) {
		PrintError("FindMatchesInIndexes",
				"tempOutputIndexFileNames",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for one file pointer per thread */
	tempOutputThreadFPs=malloc(sizeof(FILE*)*numThreads); 
	if(NULL == tempOutputThreadFPs) {
		PrintError("FindMatchesInIndexes",
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

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Temporary files opened\n");
	}

	/* For each RGIndex, write temporary output */
	for(i=0;i<numIndexes;i++) { /* For each RGIndex */
		/* Open files for thread output */
		for(j=0;j<numThreads;j++) {
			tempOutputThreadFPs[j] = OpenTmpFile(tmpDir, &tempOutputThreadFileNames[j]);
			assert(tempOutputThreadFPs[j] != NULL);
		}

		/* Initialize index */
		index.positions = NULL;
		index.chromosomes = NULL;
		index.length = 0;
		index.numTiles = 0;
		index.tileLengths = NULL;
		index.gaps = NULL;
		index.startChr=0;
		index.startPos=0;
		index.endChr=0;
		index.endPos=0;

		/* Read in the RG Index */
		startTime = time(NULL);
		ReadRGIndex(secondaryIndexFileNames[i], &index, binaryInput);
		endTime = time(NULL);
		(*totalDataStructureTime)+=endTime - startTime;	

		if(VERBOSE >= 0) {
			fprintf(stderr, "Searching given index...\n");
		}

		/* Set position to read from the beginning of the file */
		for(j=0;j<numThreads;j++) {
			fseek((*tempSeqFPs)[j], 0, SEEK_SET);
		}

		/* Execute */
		startTime = time(NULL);
		/* Initialize arguments to threads */
		for(j=0;j<numThreads;j++) {
			data[j].tempOutputFP= tempOutputThreadFPs[j];
			data[j].binaryOutput = binaryOutput;
			data[j].index = &index;
			data[j].rg = rg;
			data[j].offsets = offsets;
			data[j].numOffsets = numOffsets;
			data[j].numMismatches = numMismatches;
			data[j].numInsertions = numInsertions;
			data[j].numDeletions = numDeletions;
			data[j].numGapInsertions = numGapInsertions;
			data[j].numGapDeletions = numGapDeletions;
			data[j].pairedEnd = pairedEnd;
			data[j].maxMatches = maxMatches;
			data[j].threadID = j;
		}

		/* Open threads */
		for(j=0;j<numThreads;j++) {
			data[j].tempSeqFP = (*tempSeqFPs)[j];
			/* Start thread */
			errCode = pthread_create(&threads[j], /* thread struct */
					NULL, /* default thread attributes */
					FindMatchesInIndex, /* start routine */
					&data[j]); /* data to routine */
			if(0!=errCode) {
				PrintError("FindMatchesInIndexes",
						"pthread_create: errCode",
						"Could not start thread",
						Exit,
						ThreadError);
			}
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "Created threadID:%d\n",
						j);
			}
		}
		/* Wait for threads to return */
		for(j=0;j<numThreads;j++) {
			/* Wait for the given thread to return */
			errCode = pthread_join(threads[j],
					&status);
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "Thread returned with errCode:%d\n",
						errCode);
			}
			/* Check the return code of the thread */
			if(0!=errCode) {
				PrintError("FindMatchesInIndexes",
						"pthread_join: errCode",
						"Thread returned an error",
						Exit,
						ThreadError);
			}
			/* Reinitialize file pointer */
			tempOutputThreadFPs[j] = data[j].tempOutputFP;
			fseek(tempOutputThreadFPs[j], 0, SEEK_SET);
			assert(tempOutputThreadFPs[j]!=NULL);
		}
		endTime = time(NULL);
		(*totalSearchTime)+=endTime - startTime;
		if(VERBOSE >= 0) {
			fprintf(stderr, "\n");
		}

		/*
		   for(j=0;j<numThreads;j++) {
		   fprintf(stderr, "Thread %d found %d num matches.\n",
		   data[j].threadID,
		   data[j].numMatches);
		   }
		   */

		/* Open a temporary file (this is reentrant) */
		tempOutputIndexFPs[i] = OpenTmpFile(tmpDir, &tempOutputIndexFileNames[i]); 

		/* Merge temp thread output into temp index output */
		/* Idea: the reads were apportioned in a given order,
		 * so merge to recover the initial order.
		 * */
		if(VERBOSE >= 0) {
			fprintf(stderr, "Merging thread temp files...\n");
		}
		j=RGMatchMergeThreadTempFilesIntoOutputTempFile(tempOutputThreadFPs,
				numThreads,
				tempOutputIndexFPs[i],
				pairedEnd,
				binaryOutput);
		if(VERBOSE >= 0) {
			fprintf(stderr, "Merged %d reads from threads.\n", j);
		}

		/* Close temp thread output */
		for(j=0;j<numThreads;j++) {
			CloseTmpFile(&tempOutputThreadFPs[j],
					&tempOutputThreadFileNames[j]);
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
			fprintf(stderr, "Searching given index complete.\n");
		}
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "Merging the output from each index...\n");
	}

	/* Merge temporary output from each tree and output to the final output file. */
	startTime=time(NULL);
	/* If we not on the main indexes, output to the final output file.  Otherwise,
	 * output to a temporary file.
	 * */
	if(MainIndexes == 1) {
		/* Open temporary file for the entire index search */
		tempOutputFP=OpenTmpFile(tmpDir, &tempOutputFileName);
	}
	else {
		/* Otherwise we are in the secondary index search so output to the final output file */
		tempOutputFP=outputFP;
	}
	/* Merge the temp index files into the all indexes file */
	numWritten=RGMatchMergeFilesAndOutput(tempOutputIndexFPs,
			numIndexes,
			tempOutputFP,
			pairedEnd,
			binaryOutput,
			maxMatches);
	endTime=time(NULL);
	(*totalOutputTime)+=endTime-startTime;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Exiting RGMatchMergeFilesAndOutput to FindMatchesInIndexes\n");
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "Found matches for %d reads.\n", numWritten);
	}

	/* Close the temporary index files */
	for(i=0;i<numIndexes;i++) {
		CloseTmpFile(&tempOutputIndexFPs[i],
				&tempOutputIndexFileNames[i]);
	}

	/* Close the temporary sequence files */
	for(i=0;i<numThreads;i++) {
		/* Close temporary file */
		CloseTmpFile(&(*tempSeqFPs)[i], &(*tempSeqFileNames)[i]);
	}

	if(MainIndexes == 1) {
		/* Go through the temporary output file and output those sequences that have 
		 * at least one match to the final output file.  For those sequences that have
		 * zero matches, output them to the temporary sequence file */

		if(VERBOSE >= 0) {
			fprintf(stderr, "Copying unmatched reads for secondary index search.\n");
		}

		/* Open a new temporary sequence file */
		tempSeqFP = OpenTmpFile(tmpDir, &tempSeqFileName);

		startTime=time(NULL);
		ReadTempSequencesAndOutput(tempOutputFP,
				outputFP,
				tempSeqFP,
				pairedEnd,
				binaryOutput);
		endTime=time(NULL);
		(*totalOutputTime)+=endTime-startTime;

		/* Move to the beginning of the sequence file */
		fseek(tempSeqFP, 0, SEEK_SET);

		if(VERBOSE >= 0) {
			fprintf(stderr, "Splitting unmatched reads into temp files.\n");
		}
		/* Now apportion the remaining sequences into temp files for the threads when 
		 * searching the secondary indexes 
		 * */
		numRead=ReadSequencesToTempFile(tempSeqFP,
				tempSeqFPs,
				tempSeqFileNames,
				0,
				0,
				pairedEnd,
				numThreads,
				tmpDir,
				timing);

		/* Close the tempSeqFP */
		CloseTmpFile(&tempSeqFP, &tempSeqFileName);
	}

	/* Close the temporary output file */
	if(MainIndexes == 1) {
		CloseTmpFile(&tempOutputFP, &tempOutputFileName);
	}
	else {
		fclose(tempOutputFP);
	}

	/* Free memory for temporary file pointers */
	free(tempOutputThreadFPs);
	free(tempOutputThreadFileNames);
	free(tempOutputIndexFPs);
	free(tempOutputIndexFileNames);
	/* Free thread data */
	free(threads);
	free(data);

	if(VERBOSE >= 0) {
		if(MainIndexes == 1) {
			fprintf(stderr, "Searching main indexes complete.\n");
		}
		else {
			fprintf(stderr, "Searching secondary indexes complete.\n");
		}
	}

	return numWritten;
}

/* TODO */
void *FindMatchesInIndex(void *arg)
{
	char *FnName="FindMatchesInIndex";
	char *sequenceName=NULL;
	char *sequence=NULL;
	char *pairedSequence=NULL;
	RGMatch sequenceMatch;
	RGMatch pairedSequenceMatch;
	int numRead = 0;
	ThreadIndexData *data = (ThreadIndexData*)(arg);
	int sequenceLength;
	int pairedSequenceLength;
	/* Function arguments */
	FILE *tempSeqFP = data->tempSeqFP;
	FILE *tempOutputFP = data->tempOutputFP;
	RGIndex *index = data->index;
	RGBinary *rg = data->rg;
	int *offsets = data->offsets;
	int numOffsets = data->numOffsets;
	int numMismatches = data->numMismatches;
	int numInsertions = data->numInsertions;
	int numDeletions = data->numDeletions;
	int numGapInsertions = data->numGapInsertions;
	int numGapDeletions = data->numGapDeletions;
	int binaryOutput = data->binaryOutput;
	int pairedEnd = data->pairedEnd;
	int maxMatches = data->maxMatches;
	int threadID = data->threadID;
	data->numMatches = 0;

	if(pairedEnd==1) {
		PrintError("FindMatchesInIndex",
				"pairedEnd",
				"Paired end not implemented",
				Exit,
				OutOfRange);
	}

	/* Initialize match structures */
	RGMatchInitialize(&sequenceMatch);
	RGMatchInitialize(&pairedSequenceMatch);
	
	/* Allocate memory for the data */
	sequenceName = malloc(sizeof(char)*SEQUENCE_NAME_LENGTH);
	sequence = malloc(sizeof(char)*SEQUENCE_LENGTH);
	pairedSequence = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL == sequenceName || NULL == sequence || NULL == pairedSequence) {
		PrintError(FnName,
				"sequenceName, sequence or pairedSequence",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* For each sequence */
	if(VERBOSE >= 0) {
		fprintf(stderr, "\rthreadID:%d\tnumRead:%d",
				threadID,
				numRead);
	}
	while(EOF!=ReadNextSequence(tempSeqFP, &sequence, &sequenceLength, &pairedSequence, &pairedSequenceLength, &sequenceName, pairedEnd)) {
		numRead++;

		if(VERBOSE >= 0 && numRead%FM_ROTATE_NUM==0) {
			fprintf(stderr, "\rthreadID:%d\tnumRead:%d",
					threadID,
					numRead);
		}

		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "\nRead: %s\n%s\n",
					sequenceName,
					sequence);
		}

		RGReadsFindMatches(index,
				rg,
				&sequenceMatch,
				sequence,
				sequenceLength,
				offsets,
				numOffsets,
				numMismatches,
				numInsertions,
				numDeletions,
				numGapInsertions,
				numGapDeletions,
				maxMatches);
		if(pairedEnd==1) {
			RGReadsFindMatches(index,
					rg,
					&pairedSequenceMatch,
					pairedSequence,
					pairedSequenceLength,
					offsets,
					numOffsets,
					numMismatches,
					numInsertions,
					numDeletions,
					numGapInsertions,
					numGapDeletions,
					maxMatches);
		}

		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Outputting to a temp file.\n");
		}

		if(sequenceMatch.numEntries > 0) {
			data->numMatches++;
		}

		/* Output to file */
		RGMatchPrint(tempOutputFP, 
				sequenceName, 
				sequence, 
				pairedSequence, 
				&sequenceMatch, 
				&pairedSequenceMatch, 
				pairedEnd, 
				binaryOutput); 

		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Freeing matches.\n");
		}
		/* Free matches */
		RGMatchFree(&sequenceMatch);
		if(pairedEnd == 1) {
			RGMatchFree(&pairedSequenceMatch);
		}
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\rthreadID:%d\tnumRead:%d",
				threadID,
				numRead);
	}

	/* Free memory */
	free(sequenceName);
	free(sequence);
	free(pairedSequence);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Returning from FindMatchesInIndex\n");
	}
	return NULL;
}
