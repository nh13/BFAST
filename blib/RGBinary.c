#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "BError.h"
#include "BLib.h"
#include "BLibDefinitions.h"
#include "RGBinary.h"

/* TODO */
void RGBinaryRead(char *rgFileName, 
		RGBinary *rg,
		int startChr,
		int startPos,
		int endChr,
		int endPos)
{
	int i;
	FILE *fpRG=NULL;
	char c;
	char original;
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

	/* Initialize the data structure for holding the rg */
	rg->startChr=startChr;
	rg->startPos=startPos;
	rg->endChr=endChr;
	rg->endPos=endPos;
	rg->chromosomes=NULL;
	rg->numChrs=endChr-startChr+1;

	/*****/
	/* Read in the file names for each chromosome */
	/*****/

	/* open file */
	if((fpRG=fopen(rgFileName, "r"))==0) {
		PrintError("RGBinaryRead",
				rgFileName,
				"Could not open rgFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Read in file names */
	if(VERBOSE>=0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Reading in chromosomes from %s.\n", rgFileName);
	}
	while(fscanf(fpRG, "%s", defaultFileName)!=EOF) {
		numChrFileNames++;
		chrFileNames = realloc(chrFileNames, sizeof(char*)*(numChrFileNames));
		if(NULL==chrFileNames) {
			PrintError("RGBinaryRead",
					"chrFileNames",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		chrFileNames[numChrFileNames-1] = (char*)malloc(sizeof(char)*MAX_FILENAME_LENGTH);
		if(NULL==chrFileNames[numChrFileNames-1]) {
			PrintError("RGBinaryRead",
					"chrFileNames[numChrFileNames-1]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
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

	/* Update the start and end locations based on the files read in. 
	 * This is only based off of the number of files, not the file names
	 * themselves. */
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

	/*****/
	/* Read in the sequence for each chromosome. */
	/*****/

	/* Read in the the sequence */
	for(curChr=startChr;curChr<=endChr;curChr++) {

		/* Initialize the number of positions read for the chromosome. */
		numPosRead=0;

		if(VERBOSE>=0) {
			fprintf(stderr, "Reading in chromosome %d from %s.\n", 
					curChr,
					chrFileNames[curChr-1]); 
		}

		/* open file */
		assert(curChr <= numChrFileNames);
		if((fpRG=fopen(chrFileNames[curChr-1], "r"))==0) {
			PrintError("RGBinaryRead",
					chrFileNames[curChr-1],
					"Could not open chrFileNames[] for reading",
					Exit,
					OpenFileError);
		}

		/* Read in header */
		if(fscanf(fpRG, "%s", header) == EOF) {
			PrintError("RGBinaryRead",
					chrFileNames[curChr-1],
					"Could not read header from the current file",
					Exit,
					EndOfFile);
		}

		/* Update the number of chromosomes */
		numChrs++;

		/* Reallocate memory to store one more chromosome. */
		rg->chromosomes = (RGBinaryChr*)realloc(rg->chromosomes, numChrs*sizeof(RGBinaryChr));
		if(NULL == rg->chromosomes) {
			PrintError("RGBinaryRead",
					"rg->chromosomes",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		rg->chromosomes[numChrs-1].sequence = NULL;

		/* Read in chromosome. */
		if(VERBOSE >= 0) {
			fprintf(stderr, "Reading in [chr,pos]:\n[-1,-1]");
		}
		curPos=1;
		continueReading=1;
		while(continueReading==1 && fscanf(fpRG, "%c", &original) > 0) {
			/* original - will be the original sequence.  Possibilities include:
			 * Non-repeat sequence: a,c,g,t
			 * Repeat sequence: A,C,G,T
			 * Null Character: N,n
			 * */

			original=TransformFromIUPAC(original);

			/* Save the character as lower */
			c=ToLower(original);

			if(VERBOSE >= 0) {
				if(curPos%READ_ROTATE_NUM==0) {
					fprintf(stderr, "\r[%d,%d]",
							curChr,
							curPos);
				}
			}

			/* Reallocate memory in increments.  This allows us to avoid having
			 * to reallocate memory every iteration through the loop, thus speeding
			 * up computation.  
			 * */
			if(c == '\n') {
				/* ignore */
			}
			else {
				/* Validate base pair */
				if(ValidateBasePair(original)==0) {
					fprintf(stderr, "Base:[%c]\n", original);
					PrintError("RGBinaryRead",
							"original",
							"Not a valid base pair",
							Exit,
							OutOfRange);
				}
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
						rg->chromosomes[numChrs-1].sequence = realloc(rg->chromosomes[numChrs-1].sequence, sizeof(unsigned char)*(sequenceIndex+1));
						if(NULL == rg->chromosomes[numChrs-1].sequence) {
							PrintError("RGBinaryRead",
									"rg->chromosomes[numChrs-1].sequence",
									"Could not reallocate memory",
									Exit,
									ReallocMemory);
						}
						/* Initialize byte */
						rg->chromosomes[numChrs-1].sequence[sequenceIndex] = 0;
					}
					/* Insert the sequence correctly (as opposed to incorrectly) */
					RGBinaryInsertBase(&rg->chromosomes[numChrs-1].sequence[sequenceIndex], byteIndex, c, original);
					numPosRead++;
				}
				curPos++;
			}
		}
		/* Decrement position */
		curPos--;
		/* Update our our output */
		if(VERBOSE >= 0) {
			fprintf(stderr, "\r[%d,%d]\n", curChr, curPos);
		}

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
		if(NULL == rg->chromosomes[numChrs-1].sequence) {
			PrintError("RGBinaryRead",
					"rg->chromosomes[numChrs-1].sequence",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}

		/* Close file */
		fclose(fpRG);

	}
	assert(numChrs == rg->numChrs);

	/* Add final metadata */
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

	/* Free memory */
	for(i=0;i<numChrFileNames;i++) {
		free(chrFileNames[i]);
	}
	if(numChrFileNames>0) {
		free(chrFileNames);
	}
}

/* TODO */
void RGBinaryDelete(RGBinary *rg)
{
	int i;

	/* Free each chromosome */
	for(i=0;i<rg->numChrs;i++) {
		free(rg->chromosomes[i].sequence);
	}
	/* Free the chromosomes */
	free(rg->chromosomes);
	rg->chromosomes = NULL;

	/* Initialize structure */
	rg->numChrs = 0;
	rg->startChr = 0;
	rg->startPos = 0;
	rg->endChr = 0;
	rg->endPos = 0;
}

/* TODO */
void RGBinaryInsertBase(unsigned char *dest,
		int byteIndex,
		char src,
		char repeat)
{
	/*********************************
	 * In four bits we hold two no:
	 *
	 * left two bits:
	 * 		0 - No repat and [acgt]
	 * 		1 - Repeat and [acgt]
	 * 		2 - [nN]
	 * 		3 - undefined
	 *
	 * right-most:
	 * 		0 - a
	 * 		1 - c
	 * 		2 - g
	 * 		3 - t
	 *********************************
	 * */
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
				case 'n':
					/* two */
					(*dest) = (*dest) | 0x80;
					break;
				default:
					PrintError("RGBinaryInsertSequenceLetterIntoByte",
							NULL,
							"Could not understand case 0 repeat",
							Exit,
							OutOfRange);
			}
			/* third and fourth bits from the left will hold the sequence */
			switch(src) {
				case 'n': /* This does not matter if the base is an n */
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
					PrintError("RGBinaryInsertSequenceLetterIntoByte",
							NULL,
							"Could not understand case 0 base",
							Exit,
							OutOfRange);
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
				case 'n':
					/* two */
					(*dest) = (*dest) | 0x08;
					break;
				default:
					PrintError("RGBinaryInsertSequenceLetterIntoByte",
							NULL,
							"Could not understand case 1 repeat",
							Exit,
							OutOfRange);
			}
			/* right most 2-bits will hold the sequence */
			switch(src) {
				case 'n': /* This does not matter if the base is an n */
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
					PrintError("RGBinaryInsertSequenceLetterIntoByte",
							NULL,
							"Could not understand case 1 base",
							Exit,
							OutOfRange);
			}
			break;
		default:
			PrintError("RGBinaryInsertSequenceLetterIntoByte",
					NULL,
					"Could not understand byteIndex",
					Exit,
					OutOfRange);
	}
}

/* TODO */
void RGBinaryGetSequence(RGBinary *rgBinary,
		int chromosome,
		int position,
		char strand,
		int offsetLength,
		char *reference,
		int matchLength,
		int *returnReferenceLength,
		int *returnPosition)
{
	int chrIndex;
	char *reverseCompliment;
	int numCharsPerByte;
	int startPos, endPos;
	int referenceLength = 2*offsetLength + matchLength;
	int curPos;
	/* We assume that we can hold 2 [acgt] (nts) in each byte */
	assert(ALPHABET_SIZE==4);
	numCharsPerByte=ALPHABET_SIZE/2;
	startPos=-1;
	endPos=-1;
	if(FORWARD == strand) {
		startPos = position - offsetLength;
		endPos = position + matchLength - 1 + offsetLength;
	}
	else if(REVERSE == strand) {
		startPos = position - matchLength + 1 - offsetLength;
		endPos = position + offsetLength;
	}
	else {
		PrintError("RGBinaryGetSequence",
				NULL,
				"Could not recognize strand",
				Exit,
				OutOfRange);
	}

	/* Get chromosome index in rgBinary */
	chrIndex = chromosome - (int)rgBinary->startChr;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "In RGBinaryGetSequence: user chromosome [%d] position [%d] strand [%c].\n",
				chromosome,
				position,
				strand);
		fprintf(stderr, "In RGBinaryGetSequence: chromosome [%d] with range [%d,%d].\n",
				chromosome,
				rgBinary->startChr,
				endPos == position + offsetLength);
	}

	/* Check chromosome bounds */
	if(chromosome < (int)rgBinary->startChr || chromosome > (int)rgBinary->endChr) {
		PrintError("RGBinaryGetSequence",
				NULL,
				"Chromosome is out of range",
				Exit,
				OutOfRange);
	}

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "In RGBinaryGetSequence: start position [%d] range [%d,%d] and adjusted range [%d,%d]\n",
				position,
				rgBinary->chromosomes[chrIndex].startPos,
				rgBinary->chromosomes[chrIndex].endPos,
				startPos,
				endPos);
	}

	/* Check position bounds */
	if(startPos < (int)rgBinary->chromosomes[chrIndex].startPos || endPos > (int)rgBinary->chromosomes[chrIndex].endPos) {
		/* Adjust position bounds if possible */
		if(startPos < (int)rgBinary->chromosomes[chrIndex].startPos) {
			/* Check that we have enough sequence */
			startPos = rgBinary->chromosomes[chrIndex].startPos;
		}
		if(endPos > (int)rgBinary->chromosomes[chrIndex].endPos) {
			/* Check that we have enough sequence */
			endPos = rgBinary->chromosomes[chrIndex].endPos;
		}
	}

	/* Get the reference sequence */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "startPos:%d\tendPos:%d\n",
				startPos,
				endPos);
	}
	/* Update the reference length */
	referenceLength = endPos - startPos + 1;
	assert(startPos <= endPos);
	assert(startPos >= 1);
	for(curPos=startPos;curPos<=endPos;curPos++) {
		reference[curPos-startPos] = RGBinaryGetBase(rgBinary, chromosome, curPos);
	}
	reference[curPos-startPos] = '\0';

	/* Get the reverse compliment if necessary */
	if(strand == FORWARD) {
		/* ignore */
	}
	else if(strand == REVERSE) {
		/* Get the reverse compliment */
		reverseCompliment = (char*)malloc(sizeof(char)*(referenceLength+1));
		if(NULL == reverseCompliment) {
			PrintError("RGBinaryGetSequence",
					"reverseCompliment",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		GetReverseComplimentAnyCase(reference, reverseCompliment, referenceLength);
		strcpy(reference, reverseCompliment);
		/* Free the memory */
		free(reverseCompliment);
	}
	else {
		PrintError("RGBinaryGetSequence",
				"strand",
				"Could not understand strand",
				Exit,
				OutOfRange);
	}
	/* Update start pos and reference length */
	(*returnReferenceLength) = referenceLength;
	(*returnPosition) = startPos;
	if(VERBOSE>=DEBUG) {
		fprintf(stderr, "Exiting RGBinaryGetSequence:[%s] length [%d] referenceLength [%d] and startPos [%d].\n",
				reference,
				(int)strlen(reference),
				referenceLength,
				startPos);
	}
	assert(referenceLength==strlen(reference));
}

char RGBinaryGetBase(RGBinary *rg,
		int chromosome,
		int position) 
{
	assert(chromosome >= rg->startChr && chromosome <= rg->endChr);

	int numCharsPerByte=ALPHABET_SIZE/2;
	unsigned char curByte, curChar;
	int repeat;
	int chrIndex = chromosome - rg->startChr;

	assert(position >= rg->chromosomes[chrIndex].startPos && position <= rg->chromosomes[chrIndex].endPos);

	/* For DNA */
	assert(numCharsPerByte == 2);

	/* The index in the sequence for the given position */
	int posIndex = position - rg->chromosomes[chrIndex].startPos;
	int byteIndex = posIndex%numCharsPerByte; /* Which bits in the byte */
	posIndex = (posIndex - byteIndex)/numCharsPerByte; /* Get which byte */

	/* Get the current byte */
	curByte= rg->chromosomes[chrIndex].sequence[posIndex];

	/* Extract base */
	repeat = 0;
	curChar = 'E';
	switch(byteIndex) {
		case 0:
			/* left-most 2-bits */
			repeat = curByte & 0xC0; /* zero out the irrelevant bits */
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
					PrintError("RGBinaryGetSequence",
							NULL,
							"Could not understand case 0 repeat",
							Exit,
							OutOfRange);
					break;
			}                /* third and fourth bits from the left */
			curByte = curByte & 0x30; /* zero out the irrelevant bits */
			switch(curByte) {
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
					PrintError("RGBinaryGetSequence",
							NULL,
							"Could not understand case 0 base",
							Exit,
							OutOfRange);
					break;
			}
			break;
		case 1:
			/* third and fourth bits from the right */
			repeat = curByte & 0x0C; /* zero out the irrelevant bits */
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
					PrintError("RGBinaryGetSequence",
							NULL,
							"Could not understand case 1 repeat",
							Exit,
							OutOfRange);
					break;
			}
			/* right-most 2-bits */
			curByte = curByte & 0x03; /* zero out the irrelevant bits */
			switch(curByte) {
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
					PrintError("RGBinaryGetSequence",
							NULL,
							"Could not understand case 1 base",
							Exit,
							OutOfRange);
					break;
			}
			break;
		default:
			PrintError("RGBinaryGetSequence",
					"byteIndex",
					"Could not understand byteIndex",
					Exit,
					OutOfRange);
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
			fprintf(stderr, "Error.  In RGBinaryGetSequence, could not understand repeat indexed [%d].  Terminating!\n",
					repeat);
			break;
	}
	/* Error check */
	switch(curChar) {
		case 'a':
		case 'c':
		case 'g':
		case 't':
		case 'A':
		case 'C':
		case 'G':
		case 'T':
		case 'N':
			break;
		default:
			PrintError("RGBinaryGetSequence",
					NULL,
					"Could not understand base",
					Exit,
					OutOfRange);
	}
	return curChar;
}

int RGBinaryIsRepeat(RGBinary *rg,
		int chromosome,
		int position)
{
	char curBase = RGBinaryGetBase(rg,
			chromosome,
			position);

	switch(curBase) {
		case 'A':
		case 'C':
		case 'G':
		case 'T':
			return 1;
			break;
		default:
			return 0;
			break;
	}
}

int RGBinaryIsN(RGBinary *rg,
		int chromosome, 
		int position) 
{
	char curBase = RGBinaryGetBase(rg,
			chromosome,
			position);

	return ( (curBase == 'n' || curBase == 'N')?1:0);
}
