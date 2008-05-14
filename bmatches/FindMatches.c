#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>
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
		char *rgListFileName,
		char *rgIndexMainListFileName,
		char *rgIndexListFileName,
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
		int numThreads,
		int timing
		)
{
	int numMainRGIndexes=0;
	char **rgIndexMainFileNames=NULL;

	int numRGIndexes=0;
	char **rgIndexFileNames=NULL;

	int *offsets=NULL;
	int numOffsets=0;

	FILE *seqFP=NULL;
	FILE **tempSeqFPs=NULL;
	FILE *outputFP=NULL;
	int i;

	int numMatches;
	int numReads;

	int seconds, minutes, hours;
	int totalDataStructureTime = 0; /* This will only give the to load and deleted the indexes and trees (excludes searching and other things) */
	int totalSearchTime = 0; /* This will only give the time searching (excludes load times and other things) */
	int totalOutputTime = 0; /* This wll only give the total time to merge and output */

	RGBinary rg;
	int startChr, startPos, endChr, endPos;

	/* Read in the main RGIndex File Names */
	numMainRGIndexes=ReadFileNames(rgIndexMainListFileName, &rgIndexMainFileNames);
	if(numMainRGIndexes<=0) {
		PrintError("FindMatches",
				"numMainRGIndexes",
				"Read zero indexes",
				Exit,
				OutOfRange);
	}

	/* Read in the RGIndex File Names */
	numRGIndexes=ReadFileNames(rgIndexListFileName, &rgIndexFileNames);
	if(numRGIndexes<=0) {
		PrintError("FindMatches",
				"numRGIndexes",
				"Read zero indexes",
				Exit,
				OutOfRange);
	}

	/* Check the indexes.
	 * 1. We want the two sets of files to have the same range.
	 * */
	CheckRGIndexes(rgIndexMainFileNames, 
			numMainRGIndexes,
			rgIndexFileNames,
			numRGIndexes,
			binaryInput,
			&startChr,
			&startPos,
			&endChr,
			&endPos);

	/* Read in the reference genome */
	RGBinaryRead(rgListFileName,
			&rg,
			startChr,
			startPos,
			endChr,
			endPos);

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
				"Could not allcoate memory",
				Exit,
				MallocMemory);
	}
	/* Read the sequences to the thread temp files */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Reading %s into temp files.\n",
				sequenceFileName);
	}
	numReads=ReadSequencesToTempFile(seqFP,
			&tempSeqFPs,
			startReadNum,
			endReadNum,
			pairedEnd,
			numThreads,
			timing);
	assert(numReads >= numThreads);
	if(VERBOSE >= 0) {
		fprintf(stderr, "Read %d reads from %s.\n",
				numReads,
				sequenceFileName);
	}
	/* Close the sequence file */
	fclose(seqFP);

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
		fprintf(stderr, "Processing %d reads using main indexes.\n",
				numReads);
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Do step 1: search the main indexes for all sequences */
	numMatches=FindMatchesInIndexes(rgIndexMainFileNames,
			binaryInput,
			&rg,
			numMainRGIndexes,
			offsets,
			numOffsets,
			numMismatches,
			numInsertions,
			numDeletions,
			numGapInsertions,
			numGapDeletions,
			pairedEnd,
			numThreads,
			&tempSeqFPs,
			outputFP,
			1,
			timing,
			&totalDataStructureTime,
			&totalSearchTime,
			&totalOutputTime
			);

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Procesing remaining %d reads using secondary indexes.\n",
				numReads-numMatches);
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Do step 2: search the indexes for all sequences */
	numMatches+=FindMatchesInIndexes(rgIndexFileNames,
			binaryInput,
			&rg,
			numRGIndexes,
			offsets,
			numOffsets,
			numMismatches,
			numInsertions,
			numDeletions,
			numGapInsertions,
			numGapDeletions,
			pairedEnd,
			numThreads,
			&tempSeqFPs,
			outputFP,
			0,
			timing,
			&totalDataStructureTime,
			&totalSearchTime,
			&totalOutputTime
			);

	if(VERBOSE>=0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "In total, found matches for %d out of %d reads.\n", 
				numMatches,
				numReads);
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Close output file */
	fclose(outputFP);

	/* Free main RGIndex file names */
	for(i=0;i<numMainRGIndexes;i++) {
		free(rgIndexMainFileNames[i]);
	}
	free(rgIndexMainFileNames);

	/* Free RGIndex file names */
	for(i=0;i<numRGIndexes;i++) {
		free(rgIndexFileNames[i]);
	}
	free(rgIndexFileNames);

	/* Free reference genome */
	RGBinaryDelete(&rg);

	/* Free offsets */
	free(offsets);

	/* Print timing */
	if(timing == 1) {
		/* Data structure time */
		seconds = totalDataStructureTime;
		hours = seconds/3600;
		seconds -= hours*3600;
		minutes = seconds/60;
		seconds -= minutes*60;
		fprintf(stderr, "Total time loading and deleting indexes and trees: %d hour, %d minutes and %d seconds.\n",
				hours,
				minutes,
				seconds);
		/* Search time */
		seconds = totalSearchTime;
		hours = seconds/3600;
		seconds -= hours*3600;
		minutes = seconds/60;
		seconds -= minutes*60;
		fprintf(stderr, "Total time searching indexes and trees: %d hour, %d minutes and %d seconds.\n",
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

int FindMatchesInIndexes(char **rgIndexFileNames,
		int binaryInput,
		RGBinary *rg,
		int numRGIndexes,
		int *offsets,
		int numOffsets,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions,
		int pairedEnd,
		int numThreads,
		FILE ***tempSeqFPs,
		FILE *outputFP,
		int MainIndexes,
		int timing,
		int *totalDataStructureTime,
		int *totalSearchTime,
		int *totalOutputTime)
{
	int i, j;
	FILE *tempOutputFP=NULL;
	FILE **tempOutputIndexFPs=NULL;
	FILE **tempOutputThreadFPs=NULL;
	RGIndex index;
	int numWritten, numRead;
	int startTime;
	int endTime;
	int errCode;
	ThreadIndexData *data=NULL;
	pthread_t *threads=NULL;
	FILE *tempSeqFP=NULL;
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
	tempOutputIndexFPs = (FILE**)malloc(sizeof(FILE*)*numRGIndexes);
	if(NULL == tempOutputIndexFPs) {
		PrintError("FindMatchesInIndexes",
				"tempOutputIndexFPs",
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

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Temporary files opened\n");
	}

	/* For each RGIndex, write temporary output */
	for(i=0;i<numRGIndexes;i++) { /* For each RGIndex */
		/* Open files for thread output */
		for(j=0;j<numThreads;j++) {
			tempOutputThreadFPs[j] = tmpfile();
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
		ReadRGIndex(rgIndexFileNames[i], &index, binaryInput);
		endTime = time(NULL);
		(*totalDataStructureTime)+=endTime - startTime;	

		if(VERBOSE >= 0) {
			fprintf(stderr, "\n");
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
			data[j].threadID = j;
		}

		/* Create threads */
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

		/* Open a temporary file (this is reentrant) */
		tempOutputIndexFPs[i] = tmpfile();
		assert(tempOutputIndexFPs[i]!=NULL);
		/* Merge temp thread output into temp index output */
		/* Idea: the reads were apportioned in a given order,
		 * so merge to recover the initial order.
		 * */
		RGMatchMergeThreadTempFilesIntoOutputTempFile(tempOutputThreadFPs,
				numThreads,
				tempOutputIndexFPs[i],
				pairedEnd);

		/* Close temp thread output */
		for(j=0;j<numThreads;j++) {
			fclose(tempOutputThreadFPs[j]);
		}

		/* Free memory of the RGIndex */
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
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Will enter RGMatchMergeFilesAndOutput from FindMatchesInIndexes\n");
	}
	startTime=time(NULL);
	/* If we not on the main indexes, output to the final output file.  Otherwise,
	 * output to a temporary file.
	 * */
	if(MainIndexes == 1) {
		/* Open temporary file for the entire index search */
		tempOutputFP=tmpfile();
	}
	else {
		tempOutputFP=outputFP;
	}
	/* Merge the temp index files into the all indexes file */
	numWritten=RGMatchMergeFilesAndOutput(tempOutputIndexFPs,
			numRGIndexes,
			tempOutputFP,
			pairedEnd);
	endTime=time(NULL);
	(*totalOutputTime)+=endTime-startTime;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Exiting RGMatchMergeFilesAndOutput to FindMatchesInIndexes\n");
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "Found matches for %d reads.\n", numWritten);
	}

	/* Close the temporary index files */
	for(i=0;i<numRGIndexes;i++) {
		fclose(tempOutputIndexFPs[i]);
	}

	/* Close the temporary sequence files */
	for(i=0;i<numThreads;i++) {
		/* Close temporary file */
		fclose((*tempSeqFPs)[i]);
	}

	if(MainIndexes == 1) {
		/* Go through the temporary output file and output those sequences that have 
		 * at least one match to the final output file.  For those sequences that have
		 * zero matches, output them to the temporary sequence file */

		/* Open a new temporary sequence file */
		for(i=0;i<numThreads;i++) {
			/* Open new temporary file */
			(*tempSeqFPs)[i] = tmpfile();
		}
		tempSeqFP = tmpfile();

		startTime=time(NULL);
		ReadTempSequencesAndOutput(tempOutputFP,
				outputFP,
				tempSeqFP,
				pairedEnd);
		endTime=time(NULL);
		(*totalOutputTime)+=endTime-startTime;

		/* Move to the beginning of the sequence file */
		fseek(tempSeqFP, 0, SEEK_SET);

		/* Now apportion the remaining sequences into temp files for the threads when 
		 * searching trees 
		 * */
		numRead=ReadSequencesToTempFile(tempSeqFP,
				tempSeqFPs,
				0,
				0,
				pairedEnd,
				numThreads,
				timing);

		/* Close the temporary sequence file */
		fclose(tempSeqFP);
	}

	/* Close the temporary output file */
	fclose(tempOutputFP);

	/* Free memory for temporary file pointers */
	free(tempOutputThreadFPs);
	free(tempOutputIndexFPs);
	/* Free thread data */
	free(threads);
	free(data);

	if(VERBOSE >= 0) {
		fprintf(stderr, "Searching indexes complete.\n");
	}

	return numWritten;
}

/* TODO */
void *FindMatchesInIndex(void *arg)
{
	char *sequenceName;
	char *sequence;
	char *pairedSequence;
	RGMatch sequenceMatch;
	RGMatch pairedSequenceMatch;
	int numRead = 0;
	ThreadIndexData *data = (ThreadIndexData*)(arg);
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
	int pairedEnd = data->pairedEnd;
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
	sequenceMatch.positions=NULL;
	sequenceMatch.chromosomes=NULL;
	sequenceMatch.strand=NULL;
	sequenceMatch.numEntries=0;
	pairedSequenceMatch.positions=NULL;
	pairedSequenceMatch.chromosomes=NULL;
	pairedSequenceMatch.strand=NULL;
	pairedSequenceMatch.numEntries=0;

	/* Allocate memory for the data */
	sequenceName = (char*)malloc(sizeof(char)*SEQUENCE_NAME_LENGTH);
	sequence = (char*)malloc(sizeof(char)*SEQUENCE_LENGTH);
	pairedSequence = (char*)malloc(sizeof(char)*SEQUENCE_LENGTH);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "In FindMatchesInIndex\n");
	}

	/* For each sequence */
	while(EOF!=ReadNextSequence(tempSeqFP, &sequence, &pairedSequence, &sequenceName, pairedEnd)) {
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
		/* Initialize match structures */
		sequenceMatch.positions=NULL;
		sequenceMatch.chromosomes=NULL;
		sequenceMatch.strand=NULL;
		sequenceMatch.numEntries=0;
		pairedSequenceMatch.positions=NULL;
		pairedSequenceMatch.chromosomes=NULL;
		pairedSequenceMatch.strand=NULL;
		pairedSequenceMatch.numEntries=0;

		RGReadsFindMatches(index,
				rg,
				&sequenceMatch,
				sequence,
				offsets,
				numOffsets,
				numMismatches,
				numInsertions,
				numDeletions,
				numGapInsertions,
				numGapDeletions);
		if(pairedEnd==1) {
			RGReadsFindMatches(index,
					rg,
					&pairedSequenceMatch,
					pairedSequence,
					offsets,
					numOffsets,
					numMismatches,
					numInsertions,
					numDeletions,
					numGapInsertions,
					numGapDeletions);
		}

		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Outputting to a temp file.\n");
		}

		if(sequenceMatch.numEntries > 0) {
			data->numMatches++;
		}

		/* Output to file */
		RGMatchOutputToFile(tempOutputFP, sequenceName, sequence, pairedSequence, &sequenceMatch, &pairedSequenceMatch, pairedEnd); 

		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Freeing matches.\n");
		}
		/* Free matches */
		if(sequenceMatch.numEntries > 0) {
			free(sequenceMatch.positions);
			sequenceMatch.positions=NULL;
			free(sequenceMatch.chromosomes);
			sequenceMatch.chromosomes=NULL;
			free(sequenceMatch.strand);
			sequenceMatch.strand=NULL;
			sequenceMatch.numEntries=0;
		}
		if(pairedEnd == 1 && pairedSequenceMatch.numEntries > 0) {
			free(pairedSequenceMatch.positions);
			pairedSequenceMatch.positions=NULL;
			free(pairedSequenceMatch.chromosomes);
			pairedSequenceMatch.chromosomes=NULL;
			free(pairedSequenceMatch.strand);
			pairedSequenceMatch.strand=NULL;
			pairedSequenceMatch.numEntries=0;
		}
	}
	if(VERBOSE>=0) {
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
