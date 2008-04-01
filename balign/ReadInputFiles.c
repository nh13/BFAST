#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "Definitions.h"
#include "../blib/BLibDefinitions.h"
#include "ReadInputFiles.h"

/* TODO */
/* What about Ns in the genome ? */
void ReadReferenceGenome(char *rgFileName, 
		RGBinary *rg,
		int startChr,
		int startPos,
		int endChr,
		int endPos)
{
	FILE *fpRG=NULL;
	char c;
	char repeat;
	int curChr;
	int curPos;
	int numChrs=0;
	int numPosRead=0;
	int continueReading=1;
	int byteIndex;
	int numCharsPerByte;
	int sequenceIndex;

	char **chrFileNames=NULL;
	int numChrFileNames=0;
	char defaultFileName[MAX_FILENAME_LENGTH]="\0";
	char header[MAX_FILENAME_LENGTH]="\0";

	/* We assume that we can hold 2 [acgt] (nts) in each byte */
	assert(ALPHABET_SIZE==4);
	numCharsPerByte=ALPHABET_SIZE/2;

	rg->startChr=startChr;
	rg->startPos=startPos;
	rg->endChr=endChr;
	rg->endPos=endPos;
	rg->chromosomes=NULL;
	rg->numChrs=endChr-startChr+1;

	/* open file */
	if((fpRG=fopen(rgFileName, "r"))==0) {
		fprintf(stderr, "Error opening %s for reading.  Terminating!\n", rgFileName);
		exit(1);
	}

	/* Read in file names */
	if(VERBOSE>=0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Reading in chromosomes from %s.\n", rgFileName);
	}
	while(fscanf(fpRG, "%s", defaultFileName)!=EOF) {
		numChrFileNames++;
		chrFileNames = realloc(chrFileNames, sizeof(char*)*(numChrFileNames));
		chrFileNames[numChrFileNames-1] = (char*)malloc(sizeof(char)*MAX_FILENAME_LENGTH);
		strcpy(chrFileNames[numChrFileNames-1], defaultFileName);
		if(VERBOSE>1) {
			fprintf(stderr, "%d:%s\n", numChrFileNames, chrFileNames[numChrFileNames-1]);
		}
	}
	if(VERBOSE>=0) {
		fprintf(stderr, "Read in %d chromosomes from %s.\n", numChrFileNames, rgFileName);
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* close file */
	fclose(fpRG);

	if(endChr < 1) {
		startChr=1;
		startPos=1;
		endChr=numChrFileNames;
		endPos=INT_MAX;
	}
	else {
		if(endPos < 1) {
			endPos=INT_MAX;
		}
		if(startChr < 1) {
			startChr=1;
			startPos=1;
		}
		else if(startPos < 1) {
			startPos = 1;
		}
	}
	assert(startChr>=1 && startChr<=numChrFileNames);
	assert(startPos>=1);
	assert(endChr>=1 && endChr<=numChrFileNames);
	assert(endPos >= 1);
	assert(startPos <= endPos);

	if(VERBOSE>=0) {
		fprintf(stderr, "%s", BREAK_LINE);
	}
	if(VERBOSE>=0) {
		fprintf(stderr, "Will read from chr%d:%d to chr%d:%d.\n",
				startChr,
				startPos,
				endChr,
				endPos);
	}

	/* Read in the the sequence */
	for(curChr=startChr;curChr<=endChr;curChr++) {

		/* Initialize the number of positions read */
		numPosRead=0;

		if(VERBOSE>=0) {
			fprintf(stderr, "Reading in chromosome %d from %s.\n", 
					curChr,
					chrFileNames[curChr-1]); 
		}

		/* open file */
		assert(curChr <= numChrFileNames);
		if((fpRG=fopen(chrFileNames[curChr-1], "r"))==0) {
			fprintf(stderr, "Error opening %s for reading.  Terminating!\n", chrFileNames[curChr-1]);
			exit(1);
		}

		/* Read in header */
		if(fscanf(fpRG, "%s", header) == EOF) {
			fprintf(stderr, "Error.  Could not read header from %s.  Terminating!\n",
					chrFileNames[curChr-1]);
			exit(1);
		}

		/* Update the number of chromosomes */
		numChrs++;

		/* Reallocate memory */
		rg->chromosomes = (RGBinaryChr*)realloc(rg->chromosomes, numChrs*sizeof(RGBinaryChr));
		rg->chromosomes[numChrs-1].sequence = NULL;

		/* Read in Chromosome */
		curPos=1;
		continueReading=1;
		while(continueReading==1 && fscanf(fpRG, "%c", &c) > 0) {
			repeat = c;
			c=ToLower(c);

			/* Reallocate memory in increments.  This allows us to avoid having
			 * to reallocate memory every iteration through the loop, thus speeding
			 * up computation.  
			 * */
			if(c == '\n') {
				/* ignore */
			}
			else {
				if(startChr == curChr && startPos > curPos) {
					/* ignore, not within bounds yet */
				}
				else if(endChr == curChr && curPos > endPos) {
					/* exit the loop */
					continueReading=0;
					/* Decrement to avoid the increment below */
					curPos--;
				}
				else {
					byteIndex = numPosRead%numCharsPerByte;
					sequenceIndex = (numPosRead - byteIndex)/2;
					if(byteIndex==0) {
						/* Allocate once we have filled up the byte */
						rg->chromosomes[numChrs-1].sequence = (unsigned char*)realloc(rg->chromosomes[numChrs-1].sequence, sizeof(unsigned char)*(sequenceIndex+1));
						/* Initialize byte */
						rg->chromosomes[numChrs-1].sequence[sequenceIndex] = 0;
					}
					/* Insert the sequence correctly (as opposed to incorrectly) */
					InsertSequenceLetterIntoByte(&rg->chromosomes[numChrs-1].sequence[sequenceIndex], byteIndex, c, repeat);
					numPosRead++;
				}
				curPos++;
			}
		}
		/* Decrement position */
		curPos--;
		/* Update our our output */

		/* Add metadata */
		if(curChr == startChr) {
			rg->chromosomes[numChrs-1].startPos = startPos;
		}
		else {
			rg->chromosomes[numChrs-1].startPos = 1;
		}
		rg->chromosomes[numChrs-1].endPos = curPos;
		rg->chromosomes[numChrs-1].chromosome = curChr;

		if(VERBOSE>=0) {
			fprintf(stderr, "Read %d bases on chr%d:%d-%d from %s.\n",
					numPosRead,
					rg->chromosomes[numChrs-1].chromosome ,
					rg->chromosomes[numChrs-1].startPos,
					rg->chromosomes[numChrs-1].endPos,
					chrFileNames[curChr-1]);
		}

		/* Reallocate to reduce memory (fit exactly) */
		rg->chromosomes[numChrs-1].sequence = realloc(rg->chromosomes[numChrs-1].sequence, sizeof(char)*(numPosRead));

		/* Close file */
		fclose(fpRG);

	}
	assert(numChrs == rg->numChrs);
	/* Add metadata */
	rg->startChr = rg->chromosomes[0].chromosome;
	rg->startPos = rg->chromosomes[0].startPos;
	rg->endChr = rg->chromosomes[rg->numChrs-1].chromosome;
	rg->endPos = rg->chromosomes[rg->numChrs-1].endPos;

	if(VERBOSE>=0) {
		fprintf(stderr, "In total read from chr%d:%d to chr%d:%d.\n",
				rg->startChr,
				rg->startPos,
				rg->endChr,
				rg->endPos);
		fprintf(stderr, "%s", BREAK_LINE);
	}

}

/* TODO */
char ToLower(char a) 
{
	switch(a) {
		case 'A':
			return 'a';
			break;
		case 'C':
			return 'c';
			break;
		case 'G':
			return 'g';
			break;
		case 'T':
			return 't';
			break;
		default:
			return a;
	}
}

/* TODO */
char ToUpper(char a) 
{
	switch(a) {
		case 'a':
			return 'A';
			break;
		case 'c':
			return 'C';
			break;
		case 'g':
			return 'G';
			break;
		case 't':
			return 'T';
			break;
		default:
			return a;
	}
}

/* TODO */
void InsertSequenceLetterIntoByte(unsigned char *dest,
		int byteIndex,
		char src,
		char repeat)
{
	int numCharsPerByte;
	/* We assume that we can hold 2 [acgt] (nts) in each byte */
	assert(ALPHABET_SIZE==4);
	numCharsPerByte=ALPHABET_SIZE/2;

	switch(byteIndex%numCharsPerByte) {
		case 0:
			(*dest)=0;
			/* left-most 2-bits will hold the repeat*/
			switch(repeat) {
				case 'a':
				case 'c':
				case 'g':
				case 't':
					/* zero */
					(*dest) = (*dest) | 0x00;
					break;
				case 'A':
				case 'C':
				case 'G':
				case 'T':
					/* one */
					(*dest) = (*dest) | 0x40;
					break;
				case 'N':
					/* two */
					(*dest) = (*dest) | 0x80;
					break;
				default:
					fprintf(stderr, "Error.  In InsertSequenceLetterIntoByte, could not undertsand repeat [%c].  Terminating!\n",
							repeat);
					exit(1);
			}
			/* third and fourth bits from the left will hold the sequence */
			switch(src) {
				case 'a':
					(*dest) = (*dest) | 0x00;
					break;
				case 'c':
					(*dest) = (*dest) | 0x10;
					break;
				case 'g':
					(*dest) = (*dest) | 0x20;
					break;
				case 't':
					(*dest) = (*dest) | 0x30;
					break;
				default:
					break;
			}
			break;
		case 1:
			/* third and fourth bits from the right will hold the repeat*/
			switch(repeat) {
				case 'a':
				case 'c':
				case 'g':
				case 't':
					/* zero */
					(*dest) = (*dest) | 0x00;
					break;
				case 'A':
				case 'C':
				case 'G':
				case 'T':
					/* one */
					(*dest) = (*dest) | 0x04;
					break;
				case 'N':
					/* two */
					(*dest) = (*dest) | 0x08;
					break;
				default:
					fprintf(stderr, "Error.  In InsertSequenceLetterIntoByte, could not undertsand repeat [%c].  Terminating!\n",
							repeat);
					exit(1);
			}
			/* right most 2-bits will hold the sequence */
			switch(src) {
				case 'a':
					(*dest) = (*dest) | 0x00;
					break;
				case 'c':
					(*dest) = (*dest) | 0x01;
					break;
				case 'g':
					(*dest) = (*dest) | 0x02;
					break;
				case 't':
					(*dest) = (*dest) | 0x03;
					break;
				default:
					break;
			}
			break;
		default:
			fprintf(stderr, "Error.  Incorrect byteIndex [%d].  Terminating!\n", byteIndex);
			exit(1);
	}
}

/* TODO */
int ReadScoringMatrix(char *scoringMatrixFileName, ScoringMatrix *sm)
{
	int i, j;
	FILE *fp;

	/* Open the scoring matrix file */
	if((fp=fopen(scoringMatrixFileName, "r"))==0) {
		fprintf(stderr, "Error.  Could not open %s for reading.  Terminating!\n",
				scoringMatrixFileName);
		exit(1);
	}

	/* Read in the gap open penalty */
	if(fscanf(fp, "%lf", &sm->gapOpenPenalty)==EOF) {
		fprintf(stderr, "Error.  Could not read gap open penalty from %s.  Terminating!\n",
				scoringMatrixFileName);
		exit(1);
	}

	/* Read in the gap close penalty */
	if(fscanf(fp, "%lf", &sm->gapExtensionPenalty)==EOF) {
		fprintf(stderr, "Error.  Could not read gap extension penalty from %s.  Terminating!\n",
				scoringMatrixFileName);
		exit(1);
	}

	/* Assume the key is acgt */
	assert(ALPHABET_SIZE==4);
	/* Allocate memory for the key */
	sm->key = (char*)malloc(sizeof(char)*ALPHABET_SIZE);
	/* Read in the key */
	/* Assume the key is acgt */
	sm->key[0] = 'a';
	sm->key[1] = 'c';
	sm->key[2] = 'g';
	sm->key[3] = 't';

	/* Allocate memory for the scores */
	sm->scores = (double**)malloc(sizeof(double*)*ALPHABET_SIZE);
	for(i=0;i<ALPHABET_SIZE;i++) {
		sm->scores[i] = (double*)malloc(sizeof(double)*ALPHABET_SIZE);
	}
	/* Read in the score matrix */
	for(i=0;i<ALPHABET_SIZE;i++) { /* Read row */
		for(j=0;j<ALPHABET_SIZE;j++) { /* Read column */
			if(fscanf(fp, "%lf", &sm->scores[i][j])==EOF) {
				fprintf(stderr, "Error.  Could not read scoring matrix cell [%d,%d].  Terminating!\n", i, j);
				exit(1);
			}
		}
	}

	/* Close the file */
	fclose(fp);

	/* TODO */
	/* Check that the matrix is symmetric */

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Read scoring matrix\n");
		fprintf(stderr, "gapOpenPenalty:%lf\ngapExtensionPenalty:%lf\n",
				sm->gapOpenPenalty,
				sm->gapExtensionPenalty);
		for(i=0;i<ALPHABET_SIZE;i++) { 
			for(j=0;j<ALPHABET_SIZE;j++) { 
				fprintf(stderr, "%lf\t", sm->scores[i][j]);
			}
			fprintf(stderr, "\n");
		}
	}

	return 1;
}
