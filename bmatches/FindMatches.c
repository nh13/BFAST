#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
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
		int timing
		)
{
	int numRGIndexes=0;
	char **rgIndexFileNames=NULL;

	int numRGTrees=0;
	char **rgTreeFileNames=NULL;
	int *offsets=NULL;
	int numOffsets=0;

	FILE *tempSeqFP=NULL;
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
	numReads=ReadSequencesToTempFile(sequenceFileName,
			&tempSeqFP,
			startReadNum,
			endReadNum,
			pairedEnd,
			timing);

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
			&tempSeqFP,
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
			tempSeqFP,
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

	/* Close temporary sequence file */
	fclose(tempSeqFP);

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
		FILE **tempSeqFP,
		FILE *outputFP,
		int timing,
		int *totalDataStructureTime,
		int *totalSearchTime,
		int *totalOutputTime)
{
	int i;
	FILE **tempFPs=NULL;
	FILE *tempOutputFP=NULL;
	RGIndex index;
	int numMatches=0;
	int matchLength=-1;
	int startTime;
	int endTime;

	/* Create temporary files */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Creating temporary files (one for each index)\n");
	}
	tempFPs = (FILE**)malloc(sizeof(FILE*)*numRGIndexes);
	for(i=0;i<numRGIndexes;i++) { /* For each RGIndex */
		/* Open a temporary file (this is reentrant) */
		tempFPs[i] = tmpfile();
		assert(tempFPs[i]!=NULL);
	}
	tempOutputFP=tmpfile();

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Temporary files opened\n");
	}

	/* For each RGIndex, write temporary output */
	for(i=0;i<numRGIndexes;i++) { /* For each RGIndex */

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

		/* reset pointer to temp file to the beginning of the file */
		fseek((*tempSeqFP), 0, SEEK_SET);

		/* Execute */
		startTime = time(NULL);
		FindMatchesInIndex((*tempSeqFP),
				tempFPs[i], 
				&index,
				pairedEnd);
		endTime = time(NULL);
		(*totalSearchTime)+=endTime - startTime;

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
	numMatches=RGMatchMergeFilesAndOutput(tempFPs,
			numRGIndexes,
			tempOutputFP,
			pairedEnd);
	endTime=time(NULL);
	(*totalOutputTime)+=endTime-startTime;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Exiting RGMatchMergeFilesAndOutput to FindMatchesInIndexes\n");
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "Found matches for %d reads.\n", numMatches);
	}

	/* Close the temporary sequence file */
	fclose((*tempSeqFP));

	/* Open a new temporary sequence file */
	(*tempSeqFP) = tmpfile();

	/* Go through the temporary output file and output those sequences that have 
	 * at least one match to the final output file.  For those sequences that have
	 * zero matches, output them to the temporary sequence file */
	startTime=time(NULL);
	ReadTempSequencesAndOutput(tempOutputFP,
			outputFP,
			(*tempSeqFP),
			pairedEnd);
	endTime=time(NULL);
	(*totalOutputTime)+=endTime-startTime;

	/* Move to the beginning of the sequence file */
	fseek((*tempSeqFP), 0, SEEK_SET);

	/* Close temporary files - THEY WILL BE DELETED */
	for(i=0;i<numRGIndexes;i++) {
		fclose(tempFPs[i]);
	}

	/* Close the temporary output file */
	fclose(tempOutputFP);

	/* Free memory for temporary file pointers */
	free(tempFPs);

	return numMatches;
}

/* TODO */
void FindMatchesInIndex(FILE *tempSeqFP,
		FILE *tempFP,
		RGIndex *index,
		int pairedEnd)
{
	char *sequenceName;
	char *sequence;
	char *pairedSequence;
	RGMatch sequenceMatch;
	RGMatch pairedSequenceMatch;
	int numRead = 0;
	int numSkipped = 0;

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

	if(VERBOSE >= 0) {
		fprintf(stderr, "Currently on read:\n0");
	}

	/* For each sequence */
	while(EOF!=ReadNextSequence(tempSeqFP, &sequence, &pairedSequence, &sequenceName, pairedEnd)) {
		numRead++;
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "\nRead: %s\n%s\n",
					sequenceName,
					sequence);
		}
		if(VERBOSE >= 0){
			if(numRead%FM_ROTATE_NUM == 0) {
				fprintf(stderr, "\r%d", numRead);
			}
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

			/* Output to file */
			RGMatchOutputToFile(tempFP, sequenceName, sequence, pairedSequence, &sequenceMatch, &pairedSequenceMatch, pairedEnd); 

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
	if(VERBOSE >= 0) {
		fprintf(stderr, "\r%d\nRead in %d reads with %d skipped.\n", numRead, numRead, numSkipped);
	}

	/* Free memory */
	free(sequenceName);
	free(sequence);
	free(pairedSequence);
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
		FILE *tempSeqFP,
		FILE *outputFP,
		int timing,
		int *totalDataStructureTime,
		int *totalSearchTime,
		int *totalOutputTime)
{
	int i;
	FILE **tempFPs=NULL;
	RGTree tree;
	int numMatches=0;
	time_t startTime;
	time_t endTime;

	/* Create temporary files */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Creating temporary files (one for each tree)\n");
	}
	tempFPs = (FILE**)malloc(sizeof(FILE*)*numRGTrees);
	for(i=0;i<numRGTrees;i++) { /* For each RGTree */
		/* Open a temporary file (this is reentrant) */
		tempFPs[i] = tmpfile();
		assert(tempFPs[i]!=NULL);
	}
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Temporary files opened\n");
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

		/* reset pointer to temp file to the beginning of the file */
		fseek(tempSeqFP, 0, SEEK_SET);

		/* Execute */
		startTime=time(NULL);
		FindMatchesInTree(tempSeqFP,
				tempFPs[i], 
				&tree,
				&offsets,
				numOffsets,
				numMismatches,
				numInsertions,
				numDeletions,
				numGapInsertions,
				numGapDeletions,
				pairedEnd);
		endTime=time(NULL);
		(*totalSearchTime)+=endTime-startTime;

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
	numMatches=RGMatchMergeFilesAndOutput(tempFPs,
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
		fclose(tempFPs[i]);
	}

	/* Free memory for temporary file pointers */
	free(tempFPs);

	return numMatches;
}

/* TODO */
void FindMatchesInTree(FILE* tempSeqFP,
		FILE *tempFP,
		RGTree *tree,
		int **offsets,
		int numOffsets,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions,
		int pairedEnd)
{
	char *sequenceName;
	char *sequence;
	char *pairedSequence;
	RGMatch sequenceMatch;
	RGMatch pairedSequenceMatch;
	int numRead = 0;

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

	if(VERBOSE >= 0) {
		fprintf(stderr, "Currently on read:\n0");
	}

	/* For each sequence */
	while(EOF!=ReadNextSequence(tempSeqFP, &sequence, &pairedSequence, &sequenceName, pairedEnd)) {
		numRead++;
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "\nRead: %s\t%s\n",
					sequenceName,
					sequence);
		}
		if(VERBOSE >= 0){
			if(numRead%FM_ROTATE_NUM == 0) {
				fprintf(stderr, "\r%d", numRead);
			}
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
				offsets,
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

		/* Output to file */
		RGMatchOutputToFile(tempFP, sequenceName, sequence, pairedSequence, &sequenceMatch, &pairedSequenceMatch, pairedEnd); 

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
	if(VERBOSE >= 0) {
		fprintf(stderr, "\r%d\nRead in %d reads.\n", numRead, numRead);
	}

	/* Free memory */
	free(sequenceName);
	free(sequence);
	free(pairedSequence);
}
