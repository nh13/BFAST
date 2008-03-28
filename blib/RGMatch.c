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
	int prevPosition=-1;
	char prevChromosome=100;
	char prevStrand='z'; 
	RGMatch t;
	int curIndex=0;

	/* Merge sort the data structure */
	RGMatchMergeSort(s, 0, s->numEntries-1);

	/* Allocate temporary */
	t.positions = (int*)malloc(sizeof(int)*(s->numEntries));
	t.chromosomes = (char*)malloc(sizeof(char)*(s->numEntries));
	t.strand = (char*)malloc(sizeof(char)*(s->numEntries));

	/* Remove duplicates */
	for(i=0;i<s->numEntries;i++) {
		if(s->chromosomes[i] == prevChromosome &&
				s->positions[i] == prevPosition &&
				s->strand[i] == prevStrand) {
			/* Ignore */
		}
		else {
			/* Copy over to temporary pair */
			t.positions[curIndex] = s->positions[i];
			t.chromosomes[curIndex] = s->chromosomes[i];
			t.strand[curIndex] = s->strand[i];
			curIndex++;

			/* Save previous */
			prevPosition = s->positions[i];
			prevChromosome = s->chromosomes[i];
			prevStrand = s->strand[i];
		}
	}

	/* Reallocate pair */
	s->positions = (int*)realloc(s->positions, sizeof(int)*curIndex);
	s->chromosomes = (char*)realloc(s->chromosomes, sizeof(char)*curIndex);
	s->strand = (char*)realloc(s->strand, sizeof(char)*curIndex);

	/* Copy over */
	for(i=0;i<curIndex;i++) {
		s->positions[i] = t.positions[i];
		s->chromosomes[i] = t.chromosomes[i];
		s->strand[i] = t.strand[i];
	}

	/* Free temporary memory */
	free(t.positions);
	free(t.chromosomes);
	free(t.strand);
}

/* TODO */
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

	int *tempPositions;
	char *tempChromosomes;
	char *tempStrand;

	if(low >= high) {
		return;
	}

	/* Partition the list into two lists and then sort them recursively */
	RGMatchMergeSort(s, low, mid);
	RGMatchMergeSort(s, mid+1, high);

	/* Allocate temporary memory */
	tempPositions = (int*)malloc((high-low+1)*sizeof(int));
	tempChromosomes = (char*)malloc((high-low+1)*sizeof(char));
	tempStrand = (char*)malloc((high-low+1)*sizeof(char));

	/* Merge the two lists */
	ctr = 0;
	while( (start_lower<=end_lower) && (start_upper<=end_upper) )
	{
		if( ((int)s->chromosomes[start_lower]) < ((int)s->chromosomes[start_upper]) ||
				(((int)s->chromosomes[start_lower]) == ((int)s->chromosomes[start_upper]) && s->positions[start_lower] < s->positions[start_upper]) ||
				(((int)s->chromosomes[start_lower]) == ((int)s->chromosomes[start_upper]) && s->positions[start_lower] == s->positions[start_upper] && s->strand[start_lower] <= s->strand[start_upper])) {
			tempPositions[ctr] = s->positions[start_lower];
			tempChromosomes[ctr] = s->chromosomes[start_lower];
			tempStrand[ctr] = s->strand[start_lower];
			start_lower++;
		}
		else {
			tempPositions[ctr] = s->positions[start_upper];
			tempChromosomes[ctr] = s->chromosomes[start_upper];
			tempStrand[ctr] = s->strand[start_upper];
			start_upper++;
		}
		ctr++;
		assert(ctr<high-low+1);
	}
	if(start_lower<=end_lower) {
		while(start_lower<=end_lower) {
			tempPositions[ctr] = s->positions[start_lower];
			tempChromosomes[ctr] = s->chromosomes[start_lower];
			tempStrand[ctr] = s->strand[start_lower];
			ctr++;
			start_lower++;
		}
	}
	else {
		while(start_upper<=end_upper) {
			tempPositions[ctr] = s->positions[start_upper];
			tempChromosomes[ctr] = s->chromosomes[start_upper];
			tempStrand[ctr] = s->strand[start_upper];
			ctr++;
			start_upper++;
		}
	}
	for(i=low, ctr=0;i<=high;i++, ctr++) {
		s->positions[i] = tempPositions[ctr];
		s->chromosomes[i] = tempChromosomes[ctr];
		s->strand[i] = tempStrand[ctr];
	}
	free(tempPositions);
	free(tempChromosomes);
	free(tempStrand);
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

	/* Print sequence name */
	fprintf(fp, "%s\n", sequenceName);

	/* Print first sequence */
	fprintf(fp, "%s", sequence);

	/* Print first sequence matches */
	fprintf(fp, "\t%d", sequenceMatch->numEntries);
	for(i=0;i<sequenceMatch->numEntries;i++) {
		fprintf(fp, "\t%d\t%d\t%c", 
				(int)sequenceMatch->chromosomes[i],
				sequenceMatch->positions[i],
				sequenceMatch->strand[i]);
	}
	fprintf(fp, "\n");

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

	return 1;
}
