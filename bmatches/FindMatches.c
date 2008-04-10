#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>
#include "../blib/BError.h"
#include "../blib/RGIndex.h"
#include "../blib/RGTree.h"
#include "../blib/RGSeqPair.h"
#include "Definitions.h"
#include "ReadInputFiles.h"
#include "FindMatches.h"

/* TODO */
void RunMatches(char *outputFileName,
		int binaryOutput,
		char *rgIndexListFileName,
		char *rgTreeListFileName,
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
	int numRGIndexes=0;
	char **rgIndexFileNames=NULL;

	int numRGTrees=0;
	char **rgTreeFileNames=NULL;
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

	/* Read in the RGIndex File Names */
	numRGIndexes=ReadFileNames(rgIndexListFileName, &rgIndexFileNames);
	if(numRGIndexes<=0) {
		PrintError("RunMatches",
				"numRGIndexes",
				"Read zero indexes",
				Exit,
				OutOfRange);
	}

	/* Read in the RGTree File Names */
	numRGTrees=ReadFileNames(rgTreeListFileName, &rgTreeFileNames);
	if(numRGTrees<=0) {
		PrintError("RunMatches",
				"numRGTrees",
				"Read zero trees",
				Exit,
				OutOfRange);
	}

	/* TODO */
	/* We should probably make sure the match length of the indexes are the same */
	/* We should probably check that the chr/pos ranges of the trees is the same as the indexes */
	/* We should probably check that the number of gap insertions (matchLength) and deletions (maxGap)
	 * are within bounds. 
	 */

	/* Read in the offsets */
	numOffsets=ReadOffsets(offsetsFileName, &offsets);

	/* Since we may be running through many indexes and trees and only look for a small portion
	 * of the reads (see startReadNum and endReadNum), we copy the relevant reads
	 * to a temporary file, thereby elmininating the need to iterate through the 
	 * source sequence read file for each index or tree. 
	 * */
	/* Open sequence file */
	/* open sequence file */
	if((seqFP=fopen(sequenceFileName, "r"))==0) {
		PrintError("RunMatches",
				sequenceFileName,
				"Could not open sequenceFileName for reading",
				Exit,
				OpenFileError);
	}
	/* Allocate memory for the temp file pointers - one for each thread */
	tempSeqFPs=malloc(sizeof(FILE*)*numThreads);
	if(NULL==tempSeqFPs) {
		PrintError("RunMatches",
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
		PrintError("RunMatches",
				outputFileName,
				"Could not open outputFileName for writing",
				Exit,
				OpenFileError);
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Processing %d reads using Indexes\n",
				numReads);
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Do step 1: search the indexes for all sequences */
	numMatches=FindMatchesInIndexes(rgIndexFileNames,
			binaryInput,
			numRGIndexes,
			pairedEnd,
			numThreads,
			&tempSeqFPs,
			outputFP,
			timing,
			&totalDataStructureTime,
			&totalSearchTime,
			&totalOutputTime
			);

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Procesing remaining %d reads using Trees\n",
				numReads-numMatches);
		fprintf(stderr, "%s", BREAK_LINE);
	}
	assert(numReads-numMatches > numThreads);

	/* Do step 2: search the trees for the remaining sequences */
	numMatches+=FindMatchesInTrees(rgTreeFileNames,
			binaryInput,
			numRGTrees,
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

	/* Close temporary sequence files */
	for(i=0;i<numThreads;i++) {
		fclose(tempSeqFPs[i]);
	}

	/* Close output file */
	fclose(outputFP);

	/* Free RGIndex file names */
	for(i=0;i<numRGIndexes;i++) {
		free(rgIndexFileNames[i]);
	}
	free(rgIndexFileNames);

	/* Free RGTree file names */
	for(i=0;i<numRGTrees;i++) {
		free(rgTreeFileNames[i]);
	}
	free(rgTreeFileNames);

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
		int numRGIndexes,
		int pairedEnd,
		int numThreads,
		FILE ***tempSeqFPs,
		FILE *outputFP,
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
	int matchLength=-1;
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
		index.nodes = NULL;
		index.numNodes = 0;
		index.matchLength = 0;
		index.startChr=0;
		index.startPos=0;
		index.endChr=0;
		index.endPos=0;

		/* Read in the RG Index */
		startTime = time(NULL);
		ReadRGIndex(rgIndexFileNames[i], &index, binaryInput);
		endTime = time(NULL);
		(*totalDataStructureTime)+=endTime - startTime;	

		/* Check for the same match length across the indexes */
		if(i==0) {
			matchLength = index.matchLength;
		}
		else if(matchLength != index.matchLength) {
			PrintError("FindMatchesInIndexes",
					"index.matchLength",
					"Current index has a different matchlength than previous indexes",
					Exit,
					OutOfRange);
		}

		/* Set position to read from the beginning of the file */
		for(j=0;j<numThreads;j++) {
			fseek((*tempSeqFPs)[j], 0, SEEK_SET);
		}

		/* Execute */
		startTime = time(NULL);
		/* Initalize arguments to threads */
		for(j=0;j<numThreads;j++) {
			data[j].tempOutputFP= tempOutputThreadFPs[j];
			data[j].index = &index;
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
			fprintf(stderr, "%s", BREAK_LINE);
		}
	}

	/* Merge temporary output from each tree and output to the final output file. */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Will enter RGMatchMergeFilesAndOutput from FindMatchesInIndexes\n");
	}
	startTime=time(NULL);
	/* Open temporary file for the entire index search */
	tempOutputFP=tmpfile();
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

	/* Close the temporary sequence files and open new ones */
	for(i=0;i<numThreads;i++) {
		/* Close temporary file */
		fclose((*tempSeqFPs)[i]);
		/* Open new temporary file */
		(*tempSeqFPs)[i] = tmpfile();
	}

	/* Open a new temporary sequence file */
	tempSeqFP = tmpfile();

	/* Go through the temporary output file and output those sequences that have 
	 * at least one match to the final output file.  For those sequences that have
	 * zero matches, output them to the temporary sequence file */
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

	/* Close the temporary output file */
	fclose(tempOutputFP);
	/* Close the temporary sequence file */
	fclose(tempSeqFP);

	/* Free memory for temporary file pointers */
	free(tempOutputThreadFPs);
	free(tempOutputIndexFPs);
	/* Free thread data */
	free(threads);
	free(data);

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
	int numSkipped = 0;
	ThreadIndexData *data = (ThreadIndexData*)(arg);
	/* Function arguments */
	FILE *tempSeqFP = data->tempSeqFP;
	FILE *tempOutputFP = data->tempOutputFP;
	RGIndex *index = data->index;
	int pairedEnd = data->pairedEnd;
	int threadID = threadID;
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

		/* Find matches */
		if(strlen(sequence) == index->matchLength) {

			RGSeqPairFindMatchesInIndex(index,
					&sequenceMatch,
					sequence);
			if(pairedEnd==1) {
				assert(strlen(pairedSequence)==index->matchLength);
				RGSeqPairFindMatchesInIndex(index,
						&pairedSequenceMatch,
						pairedSequence);
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
		else {
			numSkipped++;
		}
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

/* TODO */
int FindMatchesInTrees(char **rgTreeFileNames,
		int binaryInput,
		int numRGTrees,
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
		int timing,
		int *totalDataStructureTime,
		int *totalSearchTime,
		int *totalOutputTime)
{
	int i, j;
	FILE **tempOutputTreeFPs=NULL;
	FILE **tempOutputThreadFPs=NULL;
	RGTree tree;
	int numMatches=0;
	time_t startTime;
	time_t endTime;
	int errCode;
	ThreadTreeData *data=NULL;
	pthread_t *threads=NULL;
	void *status;

	/* Allocate memory for threads */
	threads=malloc(sizeof(pthread_t)*numThreads);
	if(NULL==threads) {
		PrintError("FindMatchesInTrees",
				"threads",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory to pass data to threads */
	data=malloc(sizeof(ThreadIndexData)*numThreads);
	if(NULL==data) {
		PrintError("FindMatchesInTrees",
				"data",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for temporary RGTree output files */
	tempOutputTreeFPs=malloc(sizeof(FILE*)*numRGTrees);
	if(NULL==tempOutputTreeFPs) {
		PrintError("FindMatchesInTrees",
				"tempOutputTreeFPs",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for temporary thread output files */
	tempOutputThreadFPs=malloc(sizeof(FILE*)*numThreads);
	if(NULL==tempOutputThreadFPs) {
		PrintError("FindMatchesInTrees",
				"tempOutputTreeFPs",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* For each RGTree, write temporary output */
	for(i=0;i<numRGTrees;i++) { /* For each RGTree */

		/* Initialize tree */
		tree.nodes = NULL;
		tree.numNodes = 0;
		tree.matchLength = 0;
		tree.startChr=0;
		tree.startPos=0;
		tree.endChr=0;
		tree.endPos=0;

		/* Read in the RG Tree */
		startTime=time(NULL);
		ReadRGTree(rgTreeFileNames[i], &tree, binaryInput);
		endTime=time(NULL);
		(*totalDataStructureTime)+=endTime-startTime;

		/* Set position to read from the beginning of the file */
		for(j=0;j<numThreads;j++) {
			fseek((*tempSeqFPs)[j], 0, SEEK_SET);
		}

		/* Execute */
		startTime=time(NULL);
		/* Initialize arguments to the thread */
		for(j=0;j<numThreads;j++) {
			/* Open file for thread output */
			tempOutputThreadFPs[j] = tmpfile();
			assert(tempOutputThreadFPs[j] != NULL);
			data[j].tempSeqFP = (*tempSeqFPs)[j];
			data[j].tempOutputFP = tempOutputThreadFPs[j];
			data[j].tree = &tree;
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
			/* Start thread */
			threads[j] = 0;
			/*
			fprintf(stderr, "HERE creating thread %d (%d) (uninit:%d)\n", 
					j, 
					data[j].threadID,
					(int)threads[j]);
			fflush(stderr);
			*/
			errCode = pthread_create(&threads[j], /* thread struct */
					NULL, /* default thread attributes */
					FindMatchesInTree, /* start routine */
					&data[j]); /* data to routine */
			/*
			fprintf(stderr, "HERE created thread %d (pthread_t:%d)\n", j, (int)threads[j]);
			fflush(stderr);
			*/
			if(0!=errCode) {
				PrintError("FindMatchesInTrees",
						"pthread_create: errCode",
						"Could not start thread",
						Exit,
						ThreadError);
			}
		}
		/*
		for(j=0;j<numThreads;j++) {
			fprintf(stderr, "(2) Testing thread %d.  Value:%d\n",
					j,
					data[j].threadID);
			fflush(stderr);
		}
		*/
		/* Wait for threads to return */
		for(j=0;j<numThreads;j++) {
			/* Wait for the given thread to return */
			/*
			fprintf(stderr, "Waiting for thread %d (pthread_t:%d) to exit.\n",
					j,
					(int)threads[j]);
			fflush(stderr);
			*/
			errCode = pthread_join(threads[j],
					&status);
			/*
			fprintf(stderr, "Caught return for thread %d (pthread_t:%d) to exit.\n",
					j,
					(int)threads[j]);
			fflush(stderr);
			*/
			/* Check the return code of the thread */
			if(0!=errCode) {
				PrintError("FindMatchesInTrees",
						"pthread_join: errCode",
						"Thread returned an error",
						Exit,
						ThreadError);
			}
			assert(data[j].tempOutputFP != NULL);
			/* Reinitialize file pointer */
			tempOutputThreadFPs[j] = data[j].tempOutputFP;
			fseek(tempOutputThreadFPs[j], 0, SEEK_SET);
			assert(tempOutputThreadFPs[j]!=NULL);
		}
		endTime=time(NULL);
		(*totalSearchTime)+=endTime-startTime;

		/* Open a temporary file (this is reentrant) */
		tempOutputTreeFPs[i] = tmpfile();
		assert(tempOutputTreeFPs[i]!=NULL);
		/* Merge temp thread output into temp index output */
		/* Idea: the reads were apportioned in a given order,
		 * so merge to recover the initial order.
		 * */
		RGMatchMergeThreadTempFilesIntoOutputTempFile(tempOutputThreadFPs,
				numThreads,
				tempOutputTreeFPs[i],
				pairedEnd);

		/* Close temp thread output */
		for(j=0;j<numThreads;j++) {
			fclose(tempOutputThreadFPs[j]);
		}

		/* Free memory of the RGTree */
		startTime=time(NULL);
		RGTreeDelete(&tree);
		endTime=time(NULL);
		(*totalDataStructureTime)+=endTime-startTime;
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
		}
	}

	/* Merge temporary output from each tree and output to the final output file. */
	startTime=time(NULL);
	numMatches=RGMatchMergeFilesAndOutput(tempOutputTreeFPs,
			numRGTrees,
			outputFP,
			pairedEnd);
	endTime=time(NULL);
	(*totalOutputTime)+=endTime-startTime;

	if(VERBOSE >= 0) {
		fprintf(stderr, "Found matches for %d reads.\n", numMatches);
	}

	/* Close temporary files - THEY WILL BE DELETED */
	for(i=0;i<numRGTrees;i++) {
		fclose(tempOutputTreeFPs[i]);
	}

	/* Free memory for temporary file pointers */
	free(tempOutputThreadFPs);
	free(tempOutputTreeFPs);
	/* Free thread data */
	free(threads);
	free(data);


	return numMatches;
}

/* TODO */
void *FindMatchesInTree(void *arg)
{
	char *sequenceName;
	char *sequence;
	char *pairedSequence;
	RGMatch sequenceMatch;
	RGMatch pairedSequenceMatch;
	int numRead = 0;
	ThreadTreeData *data = (ThreadTreeData*)(arg);
	/* Function arguments */
	FILE *tempSeqFP = data->tempSeqFP;
	FILE *tempOutputFP = data->tempOutputFP;
	RGTree *tree = data->tree;
	int *offsets = data->offsets;
	int numOffsets = data->numOffsets;
	int numMismatches = data->numMismatches;
	int numInsertions = data->numInsertions;
	int numDeletions = data->numDeletions;
	int numGapInsertions = data->numGapInsertions;
	int numGapDeletions = data->numGapDeletions;
	int pairedEnd = data->pairedEnd;
	int threadID = data->threadID;
	/*
	fprintf(stderr, "HERE starting thread %d (self:%d)\n",
			threadID,
			(int)pthread_self());
	fflush(stderr);
	if(NULL==data->tempOutputFP) {
		fprintf(stderr, "Assert will catch with %d (self:%d)\n",
			threadID,
			(int)pthread_self());
	}
	assert(data->tempOutputFP!=NULL);
	*/

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
		fprintf(stderr, "In FindMatchesInTree\n");
	}

	/* For each sequence */
	while(EOF!=ReadNextSequence(tempSeqFP, &sequence, &pairedSequence, &sequenceName, pairedEnd)) 
	{
		numRead++;
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "\nRead: %s\t%s\n",
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

		/* Find matches */
		RGSeqPairFindMatchesInTree(tree,
				&sequenceMatch,
				sequence,
				&offsets,
				numOffsets,
				numMismatches,
				numInsertions,
				numDeletions,
				numGapInsertions,
				numGapDeletions);
		if(pairedEnd==1) {
			RGSeqPairFindMatchesInTree(tree,
					&pairedSequenceMatch,
					pairedSequence,
					&offsets,
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

	/* Free memory */
	free(sequenceName);
	free(sequence);
	free(pairedSequence);

	/*
	fprintf(stderr, "HERE thread %d returning\n",
			data->threadID);
	fflush(stderr);
	*/

	return NULL;
}
