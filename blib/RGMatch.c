#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <string.h>
#include "BLibDefinitions.h"
#include "RGMatch.h"
#include "RGTree.h"
#include "RGMatch.h"

/* TODO */
void RGMatchRemoveDuplicates(RGMatch *s)
{
	int i;
	int prevIndex=0;

	if(s->numEntries > 0) {
		/* Merge sort the data structure */
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Merge sorting match\n");
			for(i=0;i<s->numEntries;i++) {
				fprintf(stderr, "%d\t%d\t%c\n",
						(int)s->chromosomes[i],
						s->positions[i],
						s->strand[i]);
			}
		}
		RGMatchMergeSort(s, 0, s->numEntries-1);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Sorted!\n");
			for(i=0;i<s->numEntries;i++) {
				fprintf(stderr, "%d\t%d\t%c\n",
						(int)s->chromosomes[i],
						s->positions[i],
						s->strand[i]);
			}
		}

		/* Remove duplicates */
		for(i=1;i<s->numEntries;i++) {
			if(RGMatchCompareAtIndex(s, prevIndex, s, i)==0) {
				/* ignore */
			}
			else {
				prevIndex++;
				/* Copy to prevIndex (incremented) */
				RGMatchCopyAtIndex(s, i, s, prevIndex);
			}
		}

		/* Reallocate pair */
		/* does not make sense if there are no entries */
		s->positions = (int*)realloc(s->positions, sizeof(int)*(prevIndex+1));
		s->chromosomes = (char*)realloc(s->chromosomes, sizeof(char)*(prevIndex+1));
		s->strand = (char*)realloc(s->strand, sizeof(char)*(prevIndex+1));
		s->numEntries = prevIndex+1;
	}
}

/* TODO */
/* Note: this is not in-place */
void RGMatchMergeSort(RGMatch *s, int low, int high)
{
	/* NOTE: when high-low < 20 we could use selection sort since it is faster 
	 * on smaller lengths. */

	/* MergeSort! */
	int mid = (low + high)/2;
	int start_upper = mid + 1;
	int end_upper = high;
	int start_lower = low;
	int end_lower = mid;
	int ctr, i;

	RGMatch temp;

	if(low >= high) {
		return;
	}

	/* Partition the list into two lists and then sort them recursively */
	RGMatchMergeSort(s, low, mid);
	RGMatchMergeSort(s, mid+1, high);

	/* Allocate temporary memory */
	temp.positions= (int*)malloc((high-low+1)*sizeof(int));
	temp.chromosomes = (char*)malloc((high-low+1)*sizeof(char));
	temp.strand = (char*)malloc((high-low+1)*sizeof(char));
	temp.numEntries = high-low+1;

	/* Merge the two lists */
	ctr = 0;
	while( (start_lower<=end_lower) && (start_upper<=end_upper) )
	{
		if(RGMatchCompareAtIndex(s, start_lower, s, start_upper) <= 0) {
			RGMatchCopyAtIndex(s, start_lower, &temp, ctr);
			start_lower++;
		}
		else {
			RGMatchCopyAtIndex(s, start_upper, &temp, ctr);
			start_upper++;
		}
		ctr++;
		assert(ctr<high-low+1);
	}
	if(start_lower<=end_lower) {
		while(start_lower<=end_lower) {
			RGMatchCopyAtIndex(s, start_lower, &temp, ctr);
			ctr++;
			start_lower++;
		}
	}
	else {
		while(start_upper<=end_upper) {
			RGMatchCopyAtIndex(s, start_upper, &temp, ctr);
			ctr++;
			start_upper++;
		}
	}
	for(i=low, ctr=0;i<=high;i++, ctr++) {
		RGMatchCopyAtIndex(&temp, ctr, s, i);
	}
	free(temp.positions);
	free(temp.chromosomes);
	free(temp.strand);
}

/* TODO */
void RGMatchOutputToFile(FILE *fp,
		char *sequenceName,
		char *sequence,
		char *pairedSequence,
		RGMatch *sequenceMatch,
		RGMatch *pairedSequenceMatch,
		int pairedEnd)
{
	int i;
	/* Print the matches to the output file */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "In RGMatchOutputToFile.\n");
	}

	/* Print sequence name */
	fprintf(fp, "%s\n", sequenceName);
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "%s\n", sequenceName);
	}

	/* Print first sequence */
	fprintf(fp, "%s", sequence);
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "%s", sequence);
	}

	/* Print first sequence matches */
	fprintf(fp, "\t%d", sequenceMatch->numEntries);
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "\t%d", sequenceMatch->numEntries);
	}

	for(i=0;i<sequenceMatch->numEntries;i++) {
		fprintf(fp, "\t%d\t%d\t%c", 
				(int)sequenceMatch->chromosomes[i],
				sequenceMatch->positions[i],
				sequenceMatch->strand[i]);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "\t%d\t%d\t%c", 
					(int)sequenceMatch->chromosomes[i],
					sequenceMatch->positions[i],
					sequenceMatch->strand[i]);
		}
	}
	fprintf(fp, "\n");
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "\n");
	}

	/* Print Paired end if necessary */
	if(pairedEnd == 1) {
		/* Print paired sequence */
		fprintf(fp, "%s", pairedSequence);

		/* Print first pairedSequence matches */
		fprintf(fp, "\t%d", pairedSequenceMatch->numEntries);
		for(i=0;i<pairedSequenceMatch->numEntries;i++) {
			fprintf(fp, "\t%d\t%d\t%c", 
					(int)pairedSequenceMatch->chromosomes[i],
					pairedSequenceMatch->positions[i],
					pairedSequenceMatch->strand[i]);
		}
		fprintf(fp, "\n");
	}
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "\nExiting RGMatchOutputToFile.\n");
	}
}

/* TODO */
void RGMatchMergeFilesAndOutput(FILE **tempFPs,
		int numFiles,
		FILE *outputFP,
		int pairedEnd)
{
	int i;
	RGMatch match;
	RGMatch pairedMatch;
	int continueReading = 1;
	char **sequenceNames;
	char **sequences;
	char **pairedSequences;

	/* Allocate memory for the sequenceNames, sequences and pairedSequences */
	sequenceNames = (char**)malloc(sizeof(char*)*numFiles);
	sequences= (char**)malloc(sizeof(char*)*numFiles);
	pairedSequences = (char**)malloc(sizeof(char*)*numFiles);
	for(i=0;i<numFiles;i++) {
		sequenceNames[i] = (char*)malloc(sizeof(char)*SEQUENCE_NAME_LENGTH);
		sequences[i] = (char*)malloc(sizeof(char)*SEQUENCE_LENGTH);
		pairedSequences[i] = (char*)malloc(sizeof(char)*SEQUENCE_LENGTH);
	}

	/* Seek to the beginning of the files */
	for(i=0;i<numFiles;i++) {
		fseek(tempFPs[i], 0, SEEK_SET);
	}

	/* Read in each sequence/match one at a time */
	while(continueReading == 1) {

		/* Initialize match */
		match.positions=NULL;
		match.chromosomes=NULL;
		match.strand=NULL;
		match.numEntries=0;
		pairedMatch.positions=NULL;
		pairedMatch.chromosomes=NULL;
		pairedMatch.strand=NULL;
		pairedMatch.numEntries=0;

		/* Read one sequence from each file */ 
		for(i=0;continueReading==1 && i<numFiles;i++) {
			if(RGMatchGetNextFromFile(tempFPs[i],
						sequenceNames[i],
						sequences[i],
						pairedSequences[i],
						&match,
						&pairedMatch,
						pairedEnd)==EOF) {
				continueReading=0;
			}
		}

		if(continueReading==1) {

			/* Error checking */
			for(i=1;i<numFiles;i++) {
				/* Make sure we are reading the same sequence */
				assert(strcmp(sequenceNames[i], sequenceNames[0])==0);
				assert(strcmp(sequences[i], pairedSequences[i])==0);
			}

			/* Remove duplicates */
			RGMatchRemoveDuplicates(&match);
			if(pairedEnd==1) {
				RGMatchRemoveDuplicates(&pairedMatch);
			}

			/* Print to output file */
			RGMatchOutputToFile(outputFP,
					sequenceNames[0],
					sequences[0],
					pairedSequences[0],
					&match,
					&pairedMatch,
					pairedEnd);

			/* Free memory */
			free(match.positions);
			free(match.chromosomes);
			free(match.strand);
			free(pairedMatch.positions);
			free(pairedMatch.chromosomes);
			free(pairedMatch.strand);
		}
	}

	/* Free memory */
	for(i=0;i<numFiles;i++) {
		free(sequenceNames[i]);
		free(sequences[i]);
		free(pairedSequences[i]);
	}
	free(sequenceNames);
	free(sequences);
	free(pairedSequences);
}

/* TODO */
int RGMatchGetNextFromFile(FILE *fp,
		char *sequenceName,
		char *sequence,
		char *pairedSequence,
		RGMatch *sequenceMatch,
		RGMatch *pairedSequenceMatch,
		int pairedEnd)
{
	int i;
	int tempInt;
	/* Read the matches from the input file */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "In RGMatchGetNextFromFile.\n");
	}

	/* Read sequence name */
	if(fscanf(fp, "%s", sequenceName)==EOF) {
		return EOF;
	}

	/* Read first sequence */
	if(fscanf(fp, "%s", sequence)==EOF) {
		fprintf(stderr, "Error.  Could not read in sequence.  Terminating!\n");
		exit(1);
	}

	/* Read in the number of matches */
	if(fscanf(fp, "%d", &sequenceMatch->numEntries)==EOF) {
		fprintf(stderr, "Error.  Could not read in the number of matches.  Terminating!\n");
		exit(1);
	}

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "%s\n%s\n%d\n",
				sequenceName,
				sequence,
				sequenceMatch->numEntries);
	}

	/* Allocate memory for the matches */
	sequenceMatch->positions = (int*)malloc(sizeof(int)*(sequenceMatch->numEntries));
	sequenceMatch->chromosomes = (char*)malloc(sizeof(char)*(sequenceMatch->numEntries));
	sequenceMatch->strand = (char*)malloc(sizeof(char)*(sequenceMatch->numEntries));

	/* Read first sequence matches */
	for(i=0;i<sequenceMatch->numEntries;i++) {
		if(fscanf(fp, "%d %d %c", 
					&tempInt,
					&sequenceMatch->positions[i],
					&sequenceMatch->strand[i])==EOF) {
			fprintf(stderr, "Error.  Could not readin the %dth match.  Terminating!\n", i);
			exit(1);
		}
		sequenceMatch->chromosomes[i] = tempInt;
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "chr%d:%d:%c\t",
					(int)sequenceMatch->chromosomes[i],
					sequenceMatch->positions[i],
					sequenceMatch->strand[i]);
		}
	}
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "\n");
	}

	/* Read Paired end if necessary */
	if(pairedEnd == 1) {
		/* Read in the number of matches */
		if(fscanf(fp, "%d", &pairedSequenceMatch->numEntries)==EOF) {
			fprintf(stderr, "Error.  Could not read in the number of paired matches.  Terminating!\n");
			exit(1);
		}

		/* Allocate memory for the matches */
		pairedSequenceMatch->positions = (int*)malloc(sizeof(int)*(pairedSequenceMatch->numEntries));
		pairedSequenceMatch->chromosomes = (char*)malloc(sizeof(char)*(pairedSequenceMatch->numEntries));
		pairedSequenceMatch->strand = (char*)malloc(sizeof(char)*(pairedSequenceMatch->numEntries));

		/* Read first pairedSequence matches */
		for(i=0;i<pairedSequenceMatch->numEntries;i++) {
			if(fscanf(fp, "%d %d %c", 
						&tempInt,
						&pairedSequenceMatch->positions[i],
						&pairedSequenceMatch->strand[i])==EOF) {
				fprintf(stderr, "Error.  Could not readin the %dth paired match.  Terminating!\n", i);
				exit(1);
			}
			pairedSequenceMatch->chromosomes[i] = tempInt;
		}
	}

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Exiting RGMatchGetNextFromFile.\n");
	}
	return 1;
}

/* TODO */
int RGMatchCompareAtIndex(RGMatch *mOne, int indexOne, RGMatch *mTwo, int indexTwo) 
{
	assert(indexOne >= 0 && indexOne < mOne->numEntries);
	assert(indexTwo >= 0 && indexTwo < mTwo->numEntries);
	if( ((int)mOne->chromosomes[indexOne]) < ((int)mTwo->chromosomes[indexTwo]) ||
			(((int)mOne->chromosomes[indexOne]) == ((int)mTwo->chromosomes[indexTwo]) && mOne->positions[indexOne] < mTwo->positions[indexTwo]) ||
			(((int)mOne->chromosomes[indexOne]) == ((int)mTwo->chromosomes[indexTwo]) && mOne->positions[indexOne] == mTwo->positions[indexTwo] && mOne->strand[indexOne] <= mTwo->strand[indexTwo])) {
		return -1;
	}
	else if(((int)mOne->chromosomes[indexOne]) == ((int)mTwo->chromosomes[indexTwo]) && mOne->positions[indexOne] == mTwo->positions[indexTwo] && mOne->strand[indexOne] == mTwo->strand[indexTwo]) {
		return 0;
	}
	else {
		return 1;
	}
}

/* TODO */
void RGMatchCopyAtIndex(RGMatch *src, int srcIndex, RGMatch *dest, int destIndex)
{
	if(!(srcIndex >= 0 && srcIndex < src->numEntries)) {
		fprintf(stderr, "Error. srcIndex:%d\tnumEntries:%d\n",
				srcIndex,
				src->numEntries);
	}
	if(!(destIndex >= 0 && destIndex < dest->numEntries)) {
		fprintf(stderr, "Error. destIndex:%d\tnumEntries:%d\n",
				destIndex,
				dest->numEntries);
	}

	assert(srcIndex >= 0 && srcIndex < src->numEntries);
	assert(destIndex >= 0 && destIndex < dest->numEntries);
	dest->positions[destIndex] = src->positions[srcIndex];
	dest->chromosomes[destIndex] = src->chromosomes[srcIndex];
	dest->strand[destIndex] = src->strand[srcIndex];
}
