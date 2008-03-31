#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "ReadInputFiles.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/RGMatch.h"
#include "Definitions.h"
#include "RunAligner.h"

/* TODO */
void RunAligner(RGBinary *rgList,
		char *matchesFileName,
		char *scoringMatrixFileName,
		int algorithm,
		int offsetLength,
		int pairedEnd,
		char *outputID,
		char *outputDir)
{
	int i;
	FILE *outputFP=NULL;
	FILE *matchesFP=NULL;
	FILE *matchFP=NULL;
	char **matchFileNames=NULL;
	int numMatchFileNames=0;
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	char tempFileName[MAX_FILENAME_LENGTH]="\0";

	/* Open matches file */
	if((matchesFP=fopen(matchesFileName, "r"))==0) {
		fprintf(stderr, "Error.  Could not open %s for reading.  Terminating!\n", matchesFileName);
		exit(1);
	}

	/* Read in the list of file names*/
	while(fscanf(matchesFP, "%s", tempFileName)!=EOF) {
		numMatchFileNames++;
		/* Reallocate memory */
		matchFileNames = (char**)realloc(matchFileNames, sizeof(char*)*numMatchFileNames);
		/* Allocate memory */
		matchFileNames[numMatchFileNames-1] = (char*)malloc(sizeof(char)*(strlen(tempFileName+1)));
		/* Copy over file name */
		strcpy(matchFileNames[numMatchFileNames-1], tempFileName);
	}

	/* Create output file name */
	sprintf(outputFileName, "%sblatter.aligned.file.%s.%d.%s",
			outputDir,
			outputID,
			algorithm,
			BLATTER_ALIGN_FILE_EXTENSION);

	/* Open output file */
	if((outputFP=fopen(outputFileName, "w"))==0) {
		fprintf(stderr, "Error.  Could not open %s for writing.  Terminating!\n",
				outputFileName);
		exit(1);
	}

	for(i=0;i<numMatchFileNames;i++) { /* For each match file name */

		/* Open current match file */
		if((matchFP=fopen(matchFileNames[i], "r"))==0) {
			fprintf(stderr, "Error.  Could not open %s for reading.  Terminating!\n",
					matchFileNames[i]);
			exit(1);
		}

		/* Run selected algorithm */
		switch(algorithm) {
			case 0:
				RunDynamicProgramming(matchFP,
						rgList,
						scoringMatrixFileName,
						algorithm,
						offsetLength,
						pairedEnd,
						outputFP);
				break;
			default:
				fprintf(stderr, "Error.  Could not understand algorithm option %d.  Terminating!\n",
						algorithm);
				exit(1);
				break;
		}
		/* Close the match file */
		fclose(matchFP);
		matchFP=NULL;
	}

	/* Close output file */
	fclose(outputFP);

	/* Free memory */
	for(i=0;i<numMatchFileNames;i++) {
		free(matchFileNames[i]);
	}
	free(matchFileNames);
}

/* TODO */
void RunDynamicProgramming(FILE *matchFP,
		RGBinary *rgList,
		char *scoringMatrixFileName,
		int algorithm,
		int offsetLength,
		int pairedEnd,
		FILE *outputFP)
{
	ScoringMatrix sm;
	char sequenceName[SEQUENCE_NAME_LENGTH]="\0";
	char sequence[SEQUENCE_LENGTH]="\0";
	char pairedSequence[SEQUENCE_LENGTH]="\0";
	RGMatch sequenceMatch;
	RGMatch pairedSequenceMatch;
	int matchLength;
	char *reference=NULL;
	int i;
	
	/* Allocate memory for the reference */
	reference = (char*)malloc(sizeof(char)*(2*offsetLength+SEQUENCE_LENGTH));

	/* Read in scoring matrix */
	ReadScoringMatrix(scoringMatrixFileName, &sm); 

	/* Initialize match */
	sequenceMatch.positions=NULL;
	sequenceMatch.chromosomes=NULL;
	sequenceMatch.strand=NULL;
	sequenceMatch.numEntries=0;
	pairedSequenceMatch.positions=NULL;
	pairedSequenceMatch.chromosomes=NULL;
	pairedSequenceMatch.strand=NULL;
	pairedSequenceMatch.numEntries=0;

	/* Go through each sequence in the match file */
	while(EOF!=RGMatchGetNextFromFile(matchFP, 
				sequenceName, 
				sequence, 
				pairedSequence, 
				&sequenceMatch,
				&pairedSequenceMatch,
				pairedEnd)) {


		/* Get the sequence length */
		matchLength = strlen(sequence);
		assert(matchLength+2*offsetLength < SEQUENCE_LENGTH);

		/* Run the aligner */
		assert(pairedEnd==0);
		for(i=0;i<sequenceMatch.numEntries;i++) { /* For each match */
			/* Get the appropriate reference sequence */
			GetSequence(rgList, 
					sequenceMatch.chromosomes[i], 
					sequenceMatch.positions[i],
					sequenceMatch.strand[i],
					matchLength,
					offsetLength,
					reference);
			
			/* TODO */
			/* - Write GetSequence */
			/* - need a way to store the best match so far */
			/* - decide output format */

			/* Free memory */
			if(sequenceMatch.numEntries > 0) {
				free(sequenceMatch.positions);
				free(sequenceMatch.chromosomes);
				free(sequenceMatch.strand);
			}
			if(pairedSequenceMatch.numEntries > 0) {
				free(pairedSequenceMatch.positions);
				free(pairedSequenceMatch.chromosomes);
				free(pairedSequenceMatch.strand);
			}
			/* Initialize match */
			sequenceMatch.positions=NULL;
			sequenceMatch.chromosomes=NULL;
			sequenceMatch.strand=NULL;
			sequenceMatch.numEntries=0;
			pairedSequenceMatch.positions=NULL;
			pairedSequenceMatch.chromosomes=NULL;
			pairedSequenceMatch.strand=NULL;
			pairedSequenceMatch.numEntries=0;
		}
	}

	/* Free memory */
	free(reference);
}

/* TODO */
void GetSequence(RGBinary *rgList,
		int chromosome,
		int position,
		char strand,
		int matchLength,
		int offsetLength,
		char *reference)
{
}
