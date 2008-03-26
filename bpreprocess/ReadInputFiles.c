#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "Definitions.h"
#include "../blib/BLibDefinitions.h"
#include "ReadInputFiles.h"

/* TODO */
void ReadReferenceGenome(char *rgFileName, 
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
	int i;

	char **chrFileNames=NULL;
	int numChrFileNames=0;
	char defaultFileName[MAX_FILENAME_LENGTH]="\0";
	char header[MAX_FILENAME_LENGTH]="\0";

	rgList->startChr=startChr;
	rgList->startPos=startPos;
	rgList->endChr=endChr;
	rgList->endPos=endPos;

	/* open file */
	if((fpRG=fopen(rgFileName, "r"))==0) {
		fprintf(stderr, "Error opening %s for reading.  Terminating!\n", rgFileName);
		exit(1);
	}

	/* Read in file names */
	while(fscanf(fpRG, "%s", defaultFileName)!=EOF) {
		numChrFileNames++;
		chrFileNames = realloc(chrFileNames, sizeof(char*)*(numChrFileNames));
		chrFileNames[numChrFileNames-1] = (char*)malloc(sizeof(char)*MAX_FILENAME_LENGTH);
		strcpy(chrFileNames[numChrFileNames-1], defaultFileName);
	}

	/* close file */
	fclose(fpRG);

	for(curChr=startChr;curChr<=endChr;curChr++) {

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
		rgList->chromosomes = realloc(rgList->chromosomes, numChrs*sizeof(RGChr));

		/* Read in Chromosome */
		curPos=1;
		while(fscanf(fpRG, "%c", &c) > 0) {
			if(0 == (curPos%RIF_ROTATE_NUM)) {
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
				if( (startChr == curChr && startPos <= curChr) ||
						(endChr == curChr && curChr <= endPos)) {
					/* ignore, not within bounds */
				}
				else {
					if(numPosRead%RIF_REALLOCATE_INCREMENT==0) {
						rgList->chromosomes[numChrs-1].sequence = realloc(rgList->chromosomes[numChrs-1].sequence, sizeof(char)*(numPosRead+RIF_REALLOCATE_INCREMENT));
					}
					rgList->chromosomes[numChrs-1].sequence[curPos] = c;
					numPosRead++;
				}
				curPos++;
			}
		}

		/* Add metadata */
		if(curChr == startChr) {
			rgList->chromosomes[numChrs-1].startPos = startPos;
		}
		else {
			rgList->chromosomes[numChrs-1].startPos = 0;
		}
		if(curChr == endChr) {
			rgList->chromosomes[numChrs-1].endPos = endPos;
		}
		else {
			rgList->chromosomes[numChrs-1].endPos = 0;
		}
		rgList->chromosomes[numChrs-1].chromosome = curChr;;

		/* Reallocate to reduce memory (fit exactly) */
		rgList->chromosomes[numChrs-1].sequence = realloc(rgList->chromosomes[numChrs-1].sequence, sizeof(char)*(numPosRead));

		/* Close file */
		fclose(fpRG);

	}

	fprintf(stderr, "\n");

	for(i=0;i<numChrs;i++) {
		fprintf(stderr, "Read in chromsome %d starting at %d ending at %d.\n",
				rgList->chromosomes[i].chromosome,
				rgList->chromosomes[i].startPos,
				rgList->chromosomes[i].endPos);
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
