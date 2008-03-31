#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "ReadInputFiles.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/RGMatch.h" /* To read in the matches */
#include "../blib/RGSeqPair.h" /* For reverse compliment */
#include "Definitions.h"
#include "RunAligner.h"

/* TODO */
void RunAligner(RGBinary *rgBinary,
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
						rgBinary,
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
		RGBinary *rgBinary,
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
	int referenceLength=0;
	int i;

	/* Allocate memory for the reference */
	referenceLength = 2*offsetLength + SEQUENCE_LENGTH;
	reference = (char*)malloc(sizeof(char)*(referenceLength+1));
	reference[referenceLength] = '\0'; /* Add null terminator */

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
			GetSequenceFromReferenceGenome(rgBinary, 
					sequenceMatch.chromosomes[i], 
					sequenceMatch.positions[i],
					sequenceMatch.strand[i],
					reference,
					referenceLength);

			/* TODO */
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
void GetSequenceFromReferenceGenome(RGBinary *rgBinary,
		int chromosome,
		int position,
		char strand,
		char *reference,
		int referenceLength)
{
	int i;
	int chrIndex;
	int posIndex;
	int byteIndex;
	int curByteIndex;
	unsigned char curChar;
	unsigned char repeat;
	char *reverseCompliment;
	int numCharsPerByte;
	/* We assume that we can hold 2 [acgt] (nts) in each byte */
	assert(ALPHABET_SIZE==4);
	numCharsPerByte=ALPHABET_SIZE/2;


	int startPos = (strand == FORWARD)?(position):(position-referenceLength-1);
	int endPos = (strand == FORWARD)?(position+referenceLength-1):(position);

	/* Get chromosome index in rgBinary */
	chrIndex = chromosome - rgBinary->startChr;

	/* Check chromosome bounds */
	if(chromosome < rgBinary->startChr || chromosome > rgBinary->endChr) {
		fprintf(stderr, "Error.  In GetSequenceFromReferenceGenome: chromosome [%d] out of range [%d,%d].  Terminating!\n",
				chromosome,
				rgBinary->startChr,
				rgBinary->endChr);
		exit(1);
	}


	/* Check position bounds */
	if(startPos < rgBinary->chromosomes[chrIndex].startPos || endPos > rgBinary->chromosomes[chrIndex].endPos) {
		fprintf(stderr, "Error.  In GetSequenceFromReferenceGenome: start position [%d] out of range [%d,%d] for chromosome [%d].  Terminating!\n",
				startPos,
				rgBinary->chromosomes[chrIndex].startPos,
				rgBinary->chromosomes[chrIndex].endPos,
				chromosome);
		exit(1);
	}

	if(position < rgBinary->chromosomes[chrIndex].startPos || position > rgBinary->chromosomes[chrIndex].endPos) {
		fprintf(stderr, "Error.  In GetSequenceFromReferenceGenome: end position [%d] out of range [%d,%d] for chromosome [%d].  Terminating!\n",
				endPos,
				rgBinary->chromosomes[chrIndex].startPos,
				rgBinary->chromosomes[chrIndex].endPos,
				chromosome);
		exit(1);
	}

	/* Get the reference sequence */
	assert(startPos <= endPos);
	for(i=startPos;i<=endPos;i++) {
		/* Get position index in rgBinary */
		posIndex = i - rgBinary->chromosomes[chrIndex].startPos;
		byteIndex = i%numCharsPerByte; /* Which bits in the byte */
		curByteIndex = (posIndex - byteIndex)/numCharsPerByte; /* Get which byte */

		/* Get the sequence */
		curChar = rgBinary->chromosomes[chrIndex].sequence[curByteIndex];
		repeat = 0;
		switch(byteIndex) {
			case 0:
				/* left-most 2-bits */
				repeat = curChar & 0xC0; /* zero out the irrelevant bits */
				switch(repeat) {
					case 0x00:
						repeat = 0; 
						break;
					case 0x40:
						repeat = 1;
						break;
					case 0x80:
						repeat = 2;
						break;
					default:
						fprintf(stderr, "Error.  In GetSequenceFromReferenceGenome, could not understand repeat [%c].  Terminating!\n",
								repeat);
						exit(1);
						break;
				}
				break;
				/* third and fourth bits from the left */
				curChar = curChar & 0x30; /* zero out the irrelevant bits */
				switch(curChar) {
					case 0x00:
						curChar = 'a';
						break;
					case 0x10:
						curChar = 'c';
						break;
					case 0x20:
						curChar = 'g';
						break;
					case 0x30:
						curChar = 't';
						break;
					default:
						break;
				}
				break;
			case 1:
				/* third and fourth bits from the right */
				repeat = curChar & 0x0C; /* zero out the irrelevant bits */
				switch(repeat) {
					case 0x00:
						repeat = 0;
						break;
					case 0x04:
						repeat = 1;
						break;
					case 0x08:
						repeat = 2;
						break;
					default:
						fprintf(stderr, "Error.  In GetSequenceFromReferenceGenome, could not understand repeat [%c].  Terminating!\n",
								repeat);
						exit(1);
						break;
				}
				break;
			case 3:
				/* right-most 2-bits */
				curChar = curChar & 0x03; /* zero out the irrelevant bits */
				switch(curChar) {
					case 0x00:
						curChar = 'a';
						break;
					case 0x01:
						curChar = 'c';
						break;
					case 0x02:
						curChar = 'g';
						break;
					case 0x03:
						curChar = 't';
						break;
					default:
						break;
				}
				break;
			default:
				fprintf(stderr, "Error.  In GetSequenceFromReferenceGenome, could not understand byteIndex [%d].  Terminating!\n",
						byteIndex);
				exit(1);
		}
		/* Update based on repeat */
		switch(repeat) {
			case 0:
				/* ignore, not a repeat */
				break;
			case 1:
				/* repeat, convert char to upper */
				curChar=ToUpper(curChar);
				break;
			case 2:
				/* N character */
				curChar='N';
				break;
			default:
				fprintf(stderr, "Error.  In GetSequenceFromReferenceGenome, could not understand repeat indexed [%d].  Terminating!\n",
						repeat);
				break;
		}
		/* Update reference */
		reference[i=startPos] = curChar;
	}
	/* Get the reverse compliment if necessary */
	if(strand == FORWARD) {
		/* ignore */
	}
	else if(strand == REVERSE) {
		/* Get the reverse compliment */
		reverseCompliment = (char*)malloc(sizeof(char)*(referenceLength+1));
		GetReverseCompliment(reference, reverseCompliment, referenceLength);
		strcpy(reference, reverseCompliment);
		free(reverseCompliment);
	}
	else {
		fprintf(stderr, "Error.  Could not understand strandedness [%c].  Terminating!\n", strand);
		exit(1);
	}
}

