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

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "In GMatchRemoveDuplicates\n");
	}

	if(s->numEntries > 0) {
		/* Quick sort the data structure */
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Quick sorting match\n");
			for(i=0;i<s->numEntries;i++) {
				fprintf(stderr, "%d\t%d\t%c\n",
						s->chromosomes[i],
						s->positions[i],
						s->strand[i]);
			}
		}
		RGMatchQuickSort(s, 0, s->numEntries-1);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Sorted!\n");
			for(i=0;i<s->numEntries;i++) {
				fprintf(stderr, "%d\t%d\t%c\n",
						s->chromosomes[i],
						s->positions[i],
						s->strand[i]);
			}
		}

		/* Remove duplicates */
		prevIndex=0;
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
		s->chromosomes = (unsigned char*)realloc(s->chromosomes, sizeof(unsigned char)*(prevIndex+1));
		s->strand = (char*)realloc(s->strand, sizeof(char)*(prevIndex+1));
		s->numEntries = prevIndex+1;
	}
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Exiting GMatchRemoveDuplicates\n");
	}
}

/* TODO */
void RGMatchQuickSort(RGMatch *s, int low, int high)
{
	int i;
	int pivot=-1;
	RGMatch *temp;

	if(low < high) {
		/* Allocate memory for the temp used for swapping */
		temp=(RGMatch*)malloc(sizeof(RGMatch));
		temp->numEntries=1;
		temp->chromosomes=(unsigned char*)malloc(sizeof(unsigned char));
		temp->positions=(int*)malloc(sizeof(int));
		temp->strand=(char*)malloc(sizeof(char));

		pivot = (low+high)/2;

		RGMatchCopyAtIndex(s, pivot, temp, 0);
		RGMatchCopyAtIndex(s, high, s, pivot);
		RGMatchCopyAtIndex(temp, 0, s, high);

		pivot = low;

		for(i=low;i<high;i++) {
			if(RGMatchCompareAtIndex(s, i, s, high) <= 0) {
				RGMatchCopyAtIndex(s, i, temp, 0);
				RGMatchCopyAtIndex(s, pivot, s, i);
				RGMatchCopyAtIndex(temp, 0, s, pivot);
				pivot++;
			}
		}
		RGMatchCopyAtIndex(s, pivot, temp, 0);
		RGMatchCopyAtIndex(s, high, s, pivot);
		RGMatchCopyAtIndex(temp, 0, s, high);

		/* Free temp before the recursive call, otherwise we have a worst
		 * case of O(n) space (NOT IN PLACE) 
		 * */
		free(temp->chromosomes);
		free(temp->positions);
		free(temp->strand);
		free(temp);

		RGMatchQuickSort(s, low, pivot-1);
		RGMatchQuickSort(s, pivot+1, high);
	}
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
	assert(fp!=NULL);
	/* Print the matches to the output file */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "In RGMatchOutputToFile.\n");
	}

	/* Print sequence name */
	fprintf(fp, "%s\n", sequenceName);
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "sequenceName:%s\n", sequenceName);
	}

	/* Print first sequence */
	fprintf(fp, "%s", sequence);
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "sequence:%s", sequence);
	}

	/* Print first sequence matches */
	fprintf(fp, "\t%d", sequenceMatch->numEntries);
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "\t%d", sequenceMatch->numEntries);
	}

	for(i=0;i<sequenceMatch->numEntries;i++) {
		fprintf(fp, "\t%d\t%d\t%c", 
				sequenceMatch->chromosomes[i],
				sequenceMatch->positions[i],
				sequenceMatch->strand[i]);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "\t%d\t%d\t%c", 
					sequenceMatch->chromosomes[i],
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
					pairedSequenceMatch->chromosomes[i],
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
int RGMatchMergeFilesAndOutput(FILE **tempFPs,
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
	int numMatches=0;

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
				if(pairedEnd==1) {
					assert(strcmp(sequences[i], pairedSequences[i])==0);
				}
			}

			/* Remove duplicates */
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "numEntries:%d\n",
						match.numEntries);
				for(i=0;i<match.numEntries;i++) {
					fprintf(stderr, "%d\t%d\n",
							match.chromosomes[i],
							match.positions[i]);
				}
				fprintf(stderr, "Removing duplicates\n");
			}
			RGMatchRemoveDuplicates(&match);
			if(pairedEnd==1) {
				RGMatchRemoveDuplicates(&pairedMatch);
			}

			/* Print to output file */
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "Printing to the output file\n");
			}
			if(match.numEntries > 0) {
				numMatches++;
			}
			RGMatchOutputToFile(outputFP,
					sequenceNames[0],
					sequences[0],
					pairedSequences[0],
					&match,
					&pairedMatch,
					pairedEnd);

			/* Free memory */
			if(match.numEntries > 0) {
				free(match.positions);
				free(match.chromosomes);
				free(match.strand);
			}
			if(pairedMatch.numEntries > 0) {
				free(pairedMatch.positions);
				free(pairedMatch.chromosomes);
				free(pairedMatch.strand);
			}
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

	return numMatches;
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
		fprintf(stderr, "sequenceName:%s\nsequence:%s\nnumEntries:%d\n",
				sequenceName,
				sequence,
				sequenceMatch->numEntries);
	}

	/* Allocate memory for the matches */
	sequenceMatch->positions = (int*)malloc(sizeof(int)*(sequenceMatch->numEntries));
	sequenceMatch->chromosomes = (unsigned char*)malloc(sizeof(unsigned char)*(sequenceMatch->numEntries));
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
					sequenceMatch->chromosomes[i],
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
		pairedSequenceMatch->chromosomes = (unsigned char*)malloc(sizeof(unsigned char)*(pairedSequenceMatch->numEntries));
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
	if(mOne->chromosomes[indexOne] < mTwo->chromosomes[indexTwo] ||
			(mOne->chromosomes[indexOne] == mTwo->chromosomes[indexTwo] && mOne->positions[indexOne] < mTwo->positions[indexTwo]) ||
			(mOne->chromosomes[indexOne] == mTwo->chromosomes[indexTwo] && mOne->positions[indexOne] == mTwo->positions[indexTwo] && mOne->strand[indexOne] < mTwo->strand[indexTwo])) {
		return -1;
	}
	else if(mOne->chromosomes[indexOne] ==  mTwo->chromosomes[indexTwo] && mOne->positions[indexOne] == mTwo->positions[indexTwo] && mOne->strand[indexOne] == mTwo->strand[indexTwo]) {
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
