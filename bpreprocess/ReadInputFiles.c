#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include "Definitions.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "ReadInputFiles.h"

/* TODO */
/* Note: to save memory, we could store the genome in binary format (we reduce the size
 * required to store the genome by four */
int ReadReferenceGenome(char *rgListFileName, 
		int binaryInput,
		RGList *rgList,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int endPosOverhang)
{
	FILE *fpRG=NULL;
	char c;
	int curChr;
	int curPos;
	int numChrs=0;
	int numPosRead=0;
	int continueReading=1;

	char **chrFileNames=NULL;
	int numChrFileNames=0;
	char defaultFileName[MAX_FILENAME_LENGTH]="\0";
	char header[MAX_FILENAME_LENGTH]="\0";

	assert(binaryInput==0);

	rgList->startChr=startChr;
	rgList->startPos=startPos;
	rgList->endChr=endChr;
	rgList->endPos=endPos;
	rgList->chromosomes=NULL;
	rgList->numChrs=endChr-startChr+1;

	/* open file */
	if((fpRG=fopen(rgListFileName, "r"))==0) {
		PrintError("ReadReferenceGenome",
				rgListFileName,
				"Could not open rgListFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Read in file names */
	if(VERBOSE>=0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Reading in chromosomes from %s.\n", rgListFileName);
	}
	while(fscanf(fpRG, "%s", defaultFileName)!=EOF) {
		numChrFileNames++;
		chrFileNames = realloc(chrFileNames, sizeof(char*)*(numChrFileNames));
		if(NULL==chrFileNames) {
			PrintError("ReadReferenceGenome",
					"chrFileNames",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		chrFileNames[numChrFileNames-1] = (char*)malloc(sizeof(char)*MAX_FILENAME_LENGTH);
		if(NULL==chrFileNames[numChrFileNames-1]) {
			PrintError("ReadReferenceGenome",
					"chrFileNames[]",
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
		fprintf(stderr, "Read in %d chromosomes from %s.\n", numChrFileNames, rgListFileName);
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
				endPos+endPosOverhang);
	}

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
			PrintError("ReadReferenceGenome",
					chrFileNames[curChr-1],
					"Could not open chromosome file name for reading",
					Exit,
					OpenFileError);
		}

		/* Read in header */
		if(fscanf(fpRG, "%s", header) == EOF) {
			PrintError("ReadReferenceGenome",
					chrFileNames[curChr-1],
					"Could not read in header from chromsome file",
					Exit,
					EndOfFile);
		}

		/* Update the number of chromosomes */
		numChrs++;

		/* Reallocate memory */
		rgList->chromosomes = (RGChr*)realloc(rgList->chromosomes, numChrs*sizeof(RGChr));
		if(NULL == rgList->chromosomes) {
			PrintError("ReadReferenceGenome",
					"rgList->chromosomes",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		rgList->chromosomes[numChrs-1].sequence = NULL;

		/* Read in Chromosome */
		curPos=1;
		continueReading=1;
		while(continueReading==1 && fscanf(fpRG, "%c", &c) > 0) {
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
				else if(endChr == curChr && curPos > endPos + endPosOverhang) {
					/* exit the loop */
					continueReading=0;
					/* Decrement to avoid the increment below */
					curPos--;
				}
				else {
					if(numPosRead%RIF_REALLOCATE_INCREMENT==0) {
						rgList->chromosomes[numChrs-1].sequence = realloc(rgList->chromosomes[numChrs-1].sequence, sizeof(char)*(numPosRead+RIF_REALLOCATE_INCREMENT));
						if(NULL==rgList->chromosomes[numChrs-1].sequence) {
							PrintError("ReadReferenceGenome",
									"rgList->chromosomes[].sequence",
									"Could not reallocate memory",
									Exit,
									ReallocMemory);
						}
						rgList->chromosomes[numChrs-1].sequence[numPosRead] = c;
						numPosRead++;
					}
					curPos++;
				}
			}
		}
		/* Decrement position */
		curPos--;
		/* Update our our output */

		/* Add metadata */
		if(curChr == startChr) {
			rgList->chromosomes[numChrs-1].startPos = startPos;
		}
		else {
			rgList->chromosomes[numChrs-1].startPos = 1;
		}
		rgList->chromosomes[numChrs-1].endPos = curPos;
		rgList->chromosomes[numChrs-1].chromosome = curChr;

		if(VERBOSE>=0) {
			fprintf(stderr, "Read %d bases on chr%d:%d-%d from %s.\n",
					numPosRead,
					rgList->chromosomes[numChrs-1].chromosome ,
					rgList->chromosomes[numChrs-1].startPos,
					rgList->chromosomes[numChrs-1].endPos,
					chrFileNames[curChr-1]);
		}

		/* Reallocate to reduce memory (fit exactly) */
		rgList->chromosomes[numChrs-1].sequence = realloc(rgList->chromosomes[numChrs-1].sequence, sizeof(char)*(numPosRead));
		if(NULL == rgList->chromosomes[numChrs-1].sequence) {
			PrintError("ReadReferenceGenome",
					"rgList->chromosomes[].sequence",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}

		/* Close file */
		fclose(fpRG);

	}
	assert(numChrs == rgList->numChrs);
	/* Add metadata */
	rgList->startChr = rgList->chromosomes[0].chromosome;
	rgList->startPos = rgList->chromosomes[0].startPos;
	rgList->endChr = rgList->chromosomes[rgList->numChrs-1].chromosome;
	rgList->endPos = rgList->chromosomes[rgList->numChrs-1].endPos;

	if(VERBOSE>=0) {
		fprintf(stderr, "In total read from chr%d:%d to chr%d:%d.\n",
				rgList->startChr,
				rgList->startPos,
				rgList->endChr,
				rgList->endPos);
		fprintf(stderr, "%s", BREAK_LINE);
	}

	return continueReading;
}

void ReadGaps(char *gapFileName,
		int **gaps,
		int *numGaps,
		int *maxGap,
		int matchLength)
{
	FILE *fp;
	int tempGap;

	if((fp=fopen(gapFileName, "r"))==0) {
		PrintError("ReadGaps",
				gapFileName,
				"Could not open gapFileName for reading",
				Exit,
				OpenFileError);
	}
	(*maxGap)=0;
	(*numGaps)=0;
	while(fscanf(fp, "%d", &tempGap)!=EOF) {
		(*numGaps)++;
		(*gaps)=realloc((*gaps), sizeof(int)*(*numGaps));
		if(NULL == (*gaps)) {
			PrintError("ReadGaps",
					"(*gap)",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		(*gaps)[(*numGaps)-1] = tempGap;
		if((*numGaps)>1 && (*gaps)[(*numGaps)-2] >= (*gaps)[(*numGaps)-1]) {
			PrintError("ReadGaps",
					gapFileName,
					"Gaps in gapFileName must be in ascending order",
					Exit,
					OutOfRange);
		}
		/* Make sure that we don't go backwards */
		assert((*gaps)[(*numGaps)-1] > -1.0*fabs(matchLength));
		if((*gaps)[(*numGaps)-1]>(*maxGap)) {
			(*maxGap) = (*gaps)[(*numGaps)-1];
		}
	}
	if((*numGaps) <= 0) {
		PrintError("ReadGaps",
				gapFileName,
				"Could not read any gaps from gapFileName",
				Exit,
				OutOfRange);
	}
	fclose(fp);

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

int ValidateSequence(char *sequence, int length)
{
	int i;
	for(i=0;i<length;i++) {
		switch(sequence[i]) {
			case 'a':
				break;
			case 'c':
				break;
			case 'g':
				break;
			case 't':
				break;
			default:
				return 0;
		}
	}
	return 1;
}

