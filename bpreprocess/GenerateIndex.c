#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "math.h"
#include "../blib/RGIndex.h"
#include "../blib/BLibDefinitions.h"
#include "Definitions.h"
#include "GenerateTree.h"
#include "GenerateIndex.h"

/* TODO */
void GenerateIndex(RGList *rgList,
		int matchLength, 
		char *outputID, 
		char *outputDir)
{
	int i, j, k;
	char *curSequence=NULL;
	char outputFileName[1024];
	FILE *fp;

	RGIndex index;

	/* Allocate memory for holding the temporary sequences */
	curSequence = (char*)malloc(sizeof(char)*matchLength);

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Initialize the index */
	index.nodes = NULL;
	index.numNodes = 0;
	index.matchLength = matchLength;
	index.startChr = rgList->startChr;
	index.startPos = rgList->startPos;
	index.endChr = rgList->endChr;
	index.endPos = rgList->endPos;

	/* Create the output file name */
	sprintf(outputFileName, "%sblatter.index.file.%s.%d.%d.%d.%d.%d.%s", 
			outputDir,
			outputID,
			index.startChr,
			index.startPos,
			index.endChr,
			index.endPos,
			matchLength,
			BLATTER_INDEX_FILE_EXTENSION);

	if(VERBOSE >=0) {
		fprintf(stderr, "Generating an index file, outputing to %s.\n",
				outputFileName);
		fprintf(stderr, "Currently on [chr, pos, size]:\n%d\t%d\t%d",
				-1, -1, -1);
	}

	/* For every chromosome in the reference genome */
	for(i=0;i<rgList->numChrs;i++) {
		/* For every starting position of the l-mer */
		for(j=rgList->chromosomes[i].startPos;j <= rgList->chromosomes[i].endPos - matchLength + 1;j++) {
			if(VERBOSE >=0 && j%GT_ROTATE_NUM==0) {
				fprintf(stderr, "\r%2d\t%12d\t%10.2lfMB", 
						rgList->chromosomes[i].chromosome,
						j,
						RGIndexGetSize(&index, RGT_MEGABYTES));
			}
			/* Copy over sequences */
			for(k=j;k<j+matchLength;k++) {
				curSequence[k-j] = rgList->chromosomes[i].sequence[k-rgList->chromosomes[i].startPos];
			}
			if(VERBOSE>=DEBUG) {
				fprintf(stderr, "Inserting [%s] at chr%d:%d.\n",
						curSequence,
						rgList->chromosomes[i].chromosome,
						j);
			}
			/* Only insert if it is a valid sequence */
			if(ValidateSequence(curSequence, matchLength)==1) {
				/* Insert sequence into the index */
				RGIndexInsert(&index, curSequence, matchLength, rgList->chromosomes[i].chromosome, j);
			}
		}
		if(VERBOSE >=0) {
			fprintf(stderr, "\r%2d\t%12d\t%10.2lfMB", 
					rgList->chromosomes[i].chromosome,
					j-1,
					RGIndexGetSize(&index, RGT_MEGABYTES));
		}
	}
	/* Clean up the index */
	if(VERBOSE >= 0) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Cleaning up the index.\n");
	}
	RGIndexCleanUpIndex(&index);

	if(VERBOSE >= 0) {
		fprintf(stderr, "Outputting index to %s\n", outputFileName);
	}


	/* Open the output file */
	if(!(fp=fopen(outputFileName, "wb"))) {
		fprintf(stderr, "Error. Could not open %s for writing.  Terminating!\n",
				outputFileName);
		exit(1);
	}   

	RGIndexPrintIndex(fp, &index);

	/* Close the output file */
	fclose(fp);

	/* Free memory */
	RGIndexDelete(&index);

	if(VERBOSE >= 0) {
		fprintf(stderr, "Generated index successfully!\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Free memory */
	free(curSequence);
}
