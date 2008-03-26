#include <stdlib.h>
#include "../blib/SRTree.h"
#include "../blib/RGTree.h"
#include "../blib/RGSeqPair.h"
#include "ReadInputFiles.h"
#include "FindMatches.h"

/* TODO */
void RunMatches(char *outputFileName,
		char *sequenceFileName, 
		char *rgTreeListFileName,
		int useSequenceTree,
		int useSequentialRGTrees,
		int startReadNum,
		int endReadNum,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int pairedEnd
		)
{
	RGTree *trees;
	int numRGTrees=0;
	char **rgTreeFileNames;
	int *offsets;
	int numOffsets=0;
	FILE **tempFPs=NULL;
	FILE *tempSeqFP=NULL;
	FILE *outputFP=NULL;
	int i, j;

	if(useSequenceTree==1) {
		fprintf(stderr, "Error.  Using the sequence tree is not implemented.  Terminating!\n");
		exit(1);
	}

	/* Read in the RGTree File Names */
	numRGTrees=ReadRGTreeFileNames(rgTreeListFileName, &rgTreeFileNames, &offsets, &numOffsets);
	if(numRGTrees<=0) {
		fprintf(stderr, "Error.  Read zero rgTrees from %s.  Terminating!\n", rgTreeListFileName);
		exit(1);
	}

	/* Allocate memory for the trees */
	trees = (RGTree*)malloc(sizeof(RGTree)*numRGTrees);
	for(i=0;i<numRGTrees;i++) {
		for(j=0;j<4;j++) {
			trees[i].root->next[j] = NULL;
		}
	}

	/* Read in all the RGTrees at once, or separately? */
	if(0==useSequentialRGTrees) {
		fprintf(stderr, "Error.  Concurrent tree reading is not implemented.  Terminating!\n");
		exit(1);

		/* IDEA
		 * 		go through each read sequentially
		 * 		*/

		/* Read in the trees all at one */
		for(i=0;i<numRGTrees;i++) { /* For each RGTree */
			/* Read the RG Tree into memory */
			ReadRGTree(rgTreeFileNames[i], &trees[i]);
		}

		/* Execute */ 

		/* Free memory from the RGTrees */
		for(i=0;i<numRGTrees;i++) { /* For each RGTree */
			RGTreeDelete(&trees[i]);
		}
	}
	else {
		/* Read in trees separately */

		/* IDEA 
		 * 		use temp files to store the results for each tree.  Once we have gone
		 * 		through each tree, merge the results and output to file.
		 * 		*/

		/* Create temporary files */
		tempFPs = (FILE**)malloc(sizeof(FILE*)*numRGTrees);
		for(i=0;i<numRGTrees;i++) { /* For each RGTree */
			/* Open a temporary file (this is reentrant) */
			tempFPs[i] = tmpfile();
		}

		/* Since we may be running through many trees and only look for a small portion
		 * of the reads (see startReadNum and endReadNum), we copy the relevant reads
		 * to a temporary file, thereby elmininating the need to iterate through the 
		 * source sequence read file for each tree 
		 * */
		ReadSequencesToTempFile(sequenceFileName,
				&tempSeqFP,
				startReadNum,
				endReadNum,
				pairedEnd);

		/* For each RGTree, write temporary output */
		for(i=0;i<numRGTrees;i++) { /* For each RGTree */
			ReadRGTree(rgTreeFileNames[i], &trees[i]);

			/* reset pointer to temp file to the beginning of the file */
			fseek(tempSeqFP, 0, SEEK_SET);

			/* Execute */
			FindMatches(tempSeqFP,
					tempFPs[i], 
					&trees[i],
					&offsets,
					numOffsets,
					numMismatches,
					numInsertions,
					numDeletions,
					pairedEnd);

			/* Free memory of the RGTree */
			RGTreeDelete(&trees[i]);
		}

		/* Open output file */
		if((outputFP=fopen(outputFileName, "w"))==0) {
			fprintf(stderr, "Error opening %s for writing.  Terminating!\n", outputFileName);
			exit(1);
		}

		/* Merge temporary output */
		RGMatchMergeFilesAndOutput(tempFPs,
				numRGTrees,
				outputFP,
				pairedEnd);

		/* Close output file */
		fclose(outputFP);

		/* Close temporary files - THEY WILL BE DELETED */
		for(i=0;i<numRGTrees;i++) {
			fclose(tempFPs[i]);
		}

		/* Free memory for temporary file pointers */
		free(tempFPs);
	}

	/* Free RGTree trees */
	free(trees);

	/* Free RGTree file names */
	for(i=0;i<numRGTrees;i++) {
		free(rgTreeFileNames[i]);
	}
	free(rgTreeFileNames);

	/* Free offsets */
	free(offsets);
}

/* TODO */
void FindMatches(FILE* tempSeqFP,
		FILE *tempFP,
		RGTree *tree,
		int **offsets,
		int numOffsets,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int pairedEnd)
{
	char *sequenceName;
	char *sequence;
	char *pairedSequence;
	RGMatch sequenceMatch;
	RGMatch pairedSequenceMatch;

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

	/* For each sequence */
	while(EOF!=ReadNextSequence(tempSeqFP, &sequence, &pairedSequence, &sequenceName, pairedEnd)) {

		/* Find matches */
		RGSeqPairFindMatches(tree,
				&sequenceMatch,
				sequence,
				offsets,
				numOffsets,
				numMismatches,
				numInsertions,
				numDeletions);
		if(pairedEnd==1) {
			RGSeqPairFindMatches(tree,
					&pairedSequenceMatch,
					pairedSequence,
					offsets,
					numOffsets,
					numMismatches,
					numInsertions,
					numDeletions);
		}
		/* Output to file */
		RGMatchOutputToFile(tempFP, sequenceName, sequence, pairedSequence, &sequenceMatch, &pairedSequenceMatch, pairedEnd); 
	}

	/* Free memory */
	free(sequenceName);
	free(sequence);
	free(pairedSequence);
}
