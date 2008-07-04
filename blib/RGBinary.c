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
		int32_t startChr,
		int32_t startPos,
		int32_t endChr,
		int32_t endPos)
{
	int32_t i;
	FILE *fpRG=NULL;
	int8_t c;
	int8_t original;
	int32_t curChr;
	int32_t curPos;
	int32_t numChrs=0;
	int32_t numPosRead=0;
	int32_t continueReading=1;
	int32_t byteIndex;
	int32_t numCharsPerByte;

	char **chrFileNames=NULL;
	int32_t numChrFileNames=0;
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
		chrFileNames[numChrFileNames-1] = malloc(sizeof(char)*MAX_FILENAME_LENGTH);
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
	assert(startChr!=endChr || startPos <= endPos);

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
		rg->chromosomes = realloc(rg->chromosomes, numChrs*sizeof(RGBinaryChr));
		if(NULL == rg->chromosomes) {
			PrintError("RGBinaryRead",
					"rg->chromosomes",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		rg->chromosomes[numChrs-1].sequence = NULL;
		rg->chromosomes[numChrs-1].numBytes = 0;

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
					if(byteIndex==0) {
						rg->chromosomes[numChrs-1].numBytes++;
						/* Allocate once we have filled up the byte */
						rg->chromosomes[numChrs-1].sequence = realloc(rg->chromosomes[numChrs-1].sequence, sizeof(uint8_t)*(rg->chromosomes[numChrs-1].numBytes));
						if(NULL == rg->chromosomes[numChrs-1].sequence) {
							PrintError("RGBinaryRead",
									"rg->chromosomes[numChrs-1].sequence",
									"Could not reallocate memory",
									Exit,
									ReallocMemory);
						}
						/* Initialize byte */
						rg->chromosomes[numChrs-1].sequence[rg->chromosomes[numChrs-1].numBytes-1] = 0;
					}
					/* Insert the sequence correctly (as opposed to incorrectly) */
					RGBinaryInsertBase(&rg->chromosomes[numChrs-1].sequence[rg->chromosomes[numChrs-1].numBytes-1], byteIndex, c, original);
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
		assert(numCharsPerByte==2);
		if(rg->chromosomes[numChrs-1].numBytes != (1+rg->chromosomes[numChrs-1].endPos - rg->chromosomes[numChrs-1].startPos + 1)/numCharsPerByte) {
			fprintf(stderr, "%d\t%d\t%d\t%d\n",
					rg->chromosomes[numChrs-1].numBytes,
					rg->chromosomes[numChrs-1].startPos,
					rg->chromosomes[numChrs-1].endPos,
					rg->chromosomes[numChrs-1].endPos-rg->chromosomes[numChrs-1].startPos+1);
		}
		/* we must add one since there could be an odd number of positions */
		assert(rg->chromosomes[numChrs-1].numBytes == (1+rg->chromosomes[numChrs-1].endPos - rg->chromosomes[numChrs-1].startPos + 1)/numCharsPerByte);
		rg->chromosomes[numChrs-1].sequence = realloc(rg->chromosomes[numChrs-1].sequence, sizeof(int8_t)*(rg->chromosomes[numChrs-1].numBytes));
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
	/* Add final metadata */
	rg->numChrs = numChrs;
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
	if(numChrFileNames>0) {
		for(i=0;i<numChrFileNames;i++) {
			free(chrFileNames[i]);
		}
		free(chrFileNames);
	}
}

void RGBinaryReadBinary(RGBinary *rg,
		char *rgFileName)
{
	char *FnName="RGBinaryReadBinary";
	FILE *fpRG;
	int i;
	int32_t numCharsPerByte;
	/* We assume that we can hold 2 [acgt] (nts) in each byte */
	assert(ALPHABET_SIZE==4);
	numCharsPerByte=ALPHABET_SIZE/2;

	if(VERBOSE>=0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Reading in reference genome from %s.\n", rgFileName);
	}

	/* Open output file */
	if((fpRG=fopen(rgFileName, "rb"))==0) {
		PrintError(FnName,
				rgFileName,
				"Could not open rgFileName for writing",
				Exit,
				OpenFileError);
	}

	/* Read RGBinary information */
	if( fread(&rg->numChrs, sizeof(int32_t), 1, fpRG)!=1 ||
			fread(&rg->startChr, sizeof(int32_t), 1, fpRG)!=1 ||
			fread(&rg->startPos, sizeof(int32_t), 1, fpRG)!=1 ||
			fread(&rg->endChr, sizeof(int32_t), 1, fpRG)!=1 ||
			fread(&rg->endPos, sizeof(int32_t), 1, fpRG)!=1) {
		PrintError(FnName,
				NULL,
				"Could not read RGBinary information",
				Exit,
				ReadFileError);
	}

	/* Allocate memory for the chromosomes */
	rg->chromosomes = malloc(sizeof(RGBinaryChr)*rg->numChrs);
	if(NULL==rg->chromosomes) {
		PrintError(FnName,
				"rg->chromosomes",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Read each chromosome */
	for(i=0;i<rg->numChrs;i++) {
		/* Read RGChr information */
		if(fread(&rg->chromosomes[i].chromosome, sizeof(int32_t), 1, fpRG)!=1 ||
				fread(&rg->chromosomes[i].startPos, sizeof(int32_t), 1, fpRG)!=1 ||
				fread(&rg->chromosomes[i].endPos, sizeof(int32_t), 1, fpRG)!=1 ||
				fread(&rg->chromosomes[i].numBytes, sizeof(uint32_t), 1, fpRG)!=1) {
			PrintError(FnName,
					NULL,
					"Could not read RGChr information",
					Exit,
					ReadFileError);
		}
		/* Allocate memory for the sequence */
		rg->chromosomes[i].sequence = malloc(sizeof(uint8_t)*rg->chromosomes[i].numBytes);
		if(NULL==rg->chromosomes[i].sequence) {
			PrintError(FnName,
					"rg->chromosomes[i].sequence",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Read sequence */
		if(fread(rg->chromosomes[i].sequence, sizeof(uint8_t), rg->chromosomes[i].numBytes, fpRG)!=rg->chromosomes[i].numBytes) {
			PrintError(FnName,
					NULL,
					"Could not read sequence",
					Exit,
					ReadFileError);
		}
	}

	/* Close the output file */
	fclose(fpRG);

	if(VERBOSE>=0) {
		fprintf(stderr, "In total read from chr%d:%d to chr%d:%d.\n",
				rg->startChr,
				rg->startPos,
				rg->endChr,
				rg->endPos);
		fprintf(stderr, "%s", BREAK_LINE);
	}
}

void RGBinaryWriteBinary(RGBinary *rg,
		char *rgFileName) 
{
	char *FnName="RGBinaryWriteBinary";
	FILE *fpRG;
	int i;
	int32_t numCharsPerByte;
	/* We assume that we can hold 2 [acgt] (nts) in each byte */
	assert(ALPHABET_SIZE==4);
	numCharsPerByte=ALPHABET_SIZE/2;

	fprintf(stderr, "%s", BREAK_LINE);
	fprintf(stderr, "Outputting to %s\n", rgFileName);

	/* Open output file */
	if((fpRG=fopen(rgFileName, "wb"))==0) {
		PrintError(FnName,
				rgFileName,
				"Could not open rgFileName for writing",
				Exit,
				OpenFileError);
	}

	/* Output RGBinary information */
	if(fwrite(&rg->numChrs, sizeof(int32_t), 1, fpRG) != 1 ||
			fwrite(&rg->startChr, sizeof(int32_t), 1, fpRG) != 1 || 
			fwrite(&rg->startPos, sizeof(int32_t), 1, fpRG) != 1 ||
			fwrite(&rg->endChr, sizeof(int32_t), 1, fpRG) != 1 ||
			fwrite(&rg->endPos, sizeof(int32_t), 1, fpRG) != 1) {
		PrintError(FnName,
				NULL,
				"Could not output rg header",
				Exit,
				WriteFileError);
	}

	/* Output each chromosome */
	for(i=0;i<rg->numChrs;i++) {
		/* Output RGChr information */
		if(fwrite(&rg->chromosomes[i].chromosome, sizeof(int32_t), 1, fpRG) != 1 || 
				fwrite(&rg->chromosomes[i].startPos, sizeof(int32_t), 1, fpRG) != 1 ||
				fwrite(&rg->chromosomes[i].endPos, sizeof(int32_t), 1, fpRG) != 1 ||
				fwrite(&rg->chromosomes[i].numBytes, sizeof(uint32_t), 1, fpRG) != 1 ||
				/* Output sequence */
				fwrite(rg->chromosomes[i].sequence, sizeof(uint8_t), rg->chromosomes[i].numBytes, fpRG) != rg->chromosomes[i].numBytes) {
			PrintError(FnName,
					NULL,
					"Could not output rg chromosome",
					Exit,
					WriteFileError);
		}
	}

	fclose(fpRG);

	fprintf(stderr, "Output complete.\n");
	fprintf(stderr, "%s", BREAK_LINE);
}

/* TODO */
void RGBinaryDelete(RGBinary *rg)
{
	int32_t i;

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
void RGBinaryInsertBase(uint8_t *dest,
		int32_t byteIndex,
		int8_t src,
		int8_t repeat)
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
	int32_t numCharsPerByte;
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
		int32_t chromosome,
		int32_t position,
		int8_t strand,
		int32_t offsetLength,
		char *reference,
		int32_t matchLength,
		int32_t *returnReferenceLength,
		int32_t *returnPosition)
{
	int32_t chrIndex;
	char *reverseCompliment;
	int32_t numCharsPerByte;
	int32_t startPos, endPos;
	int32_t referenceLength = 2*offsetLength + matchLength;
	int32_t curPos;
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
		reverseCompliment = malloc(sizeof(char)*(referenceLength+1));
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

int8_t RGBinaryGetBase(RGBinary *rg,
		int32_t chromosome,
		int32_t position) 
{
	assert(chromosome >= rg->startChr && chromosome <= rg->endChr);

	int32_t numCharsPerByte=ALPHABET_SIZE/2;
	uint8_t curByte, curChar;
	int32_t repeat;
	int32_t chrIndex = chromosome - rg->startChr;

	assert(position >= rg->chromosomes[chrIndex].startPos && position <= rg->chromosomes[chrIndex].endPos);

	/* For DNA */
	assert(numCharsPerByte == 2);

	/* The index in the sequence for the given position */
	int32_t posIndex = position - rg->chromosomes[chrIndex].startPos;
	int32_t byteIndex = posIndex%numCharsPerByte; /* Which bits in the byte */
	posIndex = (posIndex - byteIndex)/numCharsPerByte; /* Get which byte */

	/* Get the current byte */
	assert(posIndex >= 0 && posIndex < rg->chromosomes[chrIndex].numBytes);
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
			/* repeat, convert int8_t to upper */
			curChar=ToUpper(curChar);
			break;
		case 2:
			/* N int8_tacter */
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

int32_t RGBinaryIsRepeat(RGBinary *rg,
		int32_t chromosome,
		int32_t position)
{
	int8_t curBase = RGBinaryGetBase(rg,
			chromosome,
			position);

	return RGBinaryIsBaseRepeat(curBase);
}

int32_t RGBinaryIsBaseRepeat(int8_t curBase)
{
	switch(curBase) {
		/* Lower case is repat */
		case 'a':
		case 'c':
		case 'g':
		case 't':
			return 1;
			break;
		case 'A':
		case 'G':
		case 'C':
		case 'T':
		default:
			return 0;
			break;
	}
}

int32_t RGBinaryIsN(RGBinary *rg,
		int32_t chromosome, 
		int32_t position) 
{
	int8_t curBase = RGBinaryGetBase(rg,
			chromosome,
			position);

	return RGBinaryIsBaseN(curBase);
}

int32_t RGBinaryIsBaseN(int8_t curBase)
{
	return ( (curBase == 'n' || curBase == 'N')?1:0);
}
