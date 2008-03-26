#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "Definitions.h"
#include "../blib/BLibDefinitions.h"
#include "ReadInputFiles.h"

/* TODO */
void ReadReferenceGenome(char *rgListFileName, 
		RGList *rgList,
		int startChr,
		int startPos,
		int endChr,
		int endPos)
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

	fprintf(stderr, "HERE:%d,%d,%d,%d\n",
			startChr,
			startPos,
			endChr,
			endPos);

	rgList->startChr=startChr;
	rgList->startPos=startPos;
	rgList->endChr=endChr;
	rgList->endPos=endPos;
	rgList->chromosomes=NULL;
	rgList->numChrs=endChr-startChr+1;

	/* open file */
	if((fpRG=fopen(rgListFileName, "r"))==0) {
		fprintf(stderr, "Error opening %s for reading.  Terminating!\n", rgListFileName);
		exit(1);
	}

	/* Read in file names */
	if(VERBOSE>0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Reading in chromosomes from %s.\n", rgListFileName);
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
	if(VERBOSE>0) {
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

	if(VERBOSE>0) {
		fprintf(stderr, "%s", BREAK_LINE);
	}
	if(VERBOSE>1) {
		fprintf(stderr, "Will read from chr%d:%d to chr%d:%d.\n",
				startChr,
				startPos,
				endChr,
				endPos);
	}

	for(curChr=startChr;curChr<=endChr;curChr++) {

		if(VERBOSE>0) {
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
		if(VERBOSE>1) {
			fprintf(stderr, "Header:%s\n", header);
		}

		/* Update the number of chromosomes */
		numChrs++;

		/* Reallocate memory */
		rgList->chromosomes = (RGChr*)realloc(rgList->chromosomes, numChrs*sizeof(RGChr));
		rgList->chromosomes[numChrs-1].sequence = NULL;

		/* Read in Chromosome */
		curPos=1;
		continueReading=1;
		while(continueReading==1 && fscanf(fpRG, "%c", &c) > 0) {
			if(VERBOSE >= 0 && 0 == (curPos%RIF_ROTATE_NUM)) {
				fprintf(stderr, "\r%d\t%d", curChr, curPos);
			}
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
					if(numPosRead%RIF_REALLOCATE_INCREMENT==0) {
						rgList->chromosomes[numChrs-1].sequence = realloc(rgList->chromosomes[numChrs-1].sequence, sizeof(char)*(numPosRead+RIF_REALLOCATE_INCREMENT));
					}
					rgList->chromosomes[numChrs-1].sequence[numPosRead] = c;
					numPosRead++;
				}
				curPos++;
			}
		}
		/* Decrement position */
		curPos--;
		/* Update our our output */
		if(VERBOSE >= 0) {
			fprintf(stderr, "\r%d\t%d", curChr, curPos);
		}

		/* Add metadata */
		if(curChr == startChr) {
			rgList->chromosomes[numChrs-1].startPos = startPos;
		}
		else {
			rgList->chromosomes[numChrs-1].startPos = 1;
		}
		rgList->chromosomes[numChrs-1].endPos = curPos;
		rgList->chromosomes[numChrs-1].chromosome = curChr;

		if(VERBOSE>1) {
			fprintf(stderr, "\n");
			fprintf(stderr, "Read %d bases on chr%d:%d-%d from %s.\n",
					numPosRead,
					rgList->chromosomes[numChrs-1].chromosome ,
					rgList->chromosomes[numChrs-1].startPos,
					rgList->chromosomes[numChrs-1].endPos,
					chrFileNames[curChr-1]);
		}

		/* Reallocate to reduce memory (fit exactly) */
		rgList->chromosomes[numChrs-1].sequence = realloc(rgList->chromosomes[numChrs-1].sequence, sizeof(char)*(numPosRead));

		/* Close file */
		fclose(fpRG);

	}
	assert(numChrs == rgList->numChrs);
	if(VERBOSE == 0) {
		fprintf(stderr, "\n");
	}
	if(VERBOSE>0) {
		fprintf(stderr, "In total read from chr%d:%d to chr%d:%d.\n",
				rgList->startChr,
				rgList->startPos,
				rgList->endChr,
				rgList->endPos);
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
