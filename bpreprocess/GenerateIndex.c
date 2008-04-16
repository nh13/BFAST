#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "sys/malloc.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/RGIndex.h"
#include "../blib/RGIndexLayout.h"
#include "../blib/RGBinary.h"
#include "Definitions.h"
#include "GenerateIndex.h"

/* TODO */
void GenerateIndex(RGBinary *rg,
		RGIndexLayout *rgLayout,
		int numThreads,
		char *outputID,
		char *outputDir,
		int binaryOutput)
{
	int i, j;
	char outputFileName[ MAX_FILENAME_LENGTH]="\0"; 
	FILE *fp=NULL;

	RGIndex index;

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
	}

	for(i=0;i<rgLayout->numIndexes;i++) { /* For each index to create */

		/* Initialize the index */
		index.positions=NULL;
		index.chromosomes=NULL;
		index.length=0;
		index.startChr = rg->startChr;
		index.startPos = rg->startPos;
		index.endChr = rg->endChr;
		index.endPos = rg->endPos;

		/* Copy over index information */
		index.totalLength = 0;
		index.numTiles = (unsigned int)rgLayout->numTiles[i];
		/* Allocate memory and copy over tile lengths */
		index.tileLengths = malloc(sizeof(unsigned int)*rgLayout->numTiles[i]);
		if(NULL == index.tileLengths) {
			PrintError("GenerateIndex",
					"index.tileLengths",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		for(j=0;j<rgLayout->numTiles[i];j++) {
			index.tileLengths[j] = rgLayout->tileLengths[i][j];
			index.totalLength += rgLayout->tileLengths[i][j];
		}
		/* Allocate memory and copy over gaps */
		index.gaps = malloc(sizeof(unsigned int)*(rgLayout->numTiles[i]-1));
		if(NULL == index.gaps) {
			PrintError("GenerateIndex",
					"index.gaps",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		for(j=0;j<rgLayout->numTiles[i]-1;j++) {
			index.gaps[j] = rgLayout->gaps[i][j];
			index.totalLength += rgLayout->gaps[i][j];
		}

		/* Create the index */
		if(VERBOSE >=0) {
			fprintf(stderr, "Creating the index...\n");
		}
		RGIndexCreate(&index, rg, 1, 0);
		if(VERBOSE >= 0) {
			fprintf(stderr, "Insertions complete.\n");
		}

		/* Clean up the index */
		if(VERBOSE >= 0) {
			fprintf(stderr, "Cleaning up...\n");
		}
		RGIndexCleanUpIndex(&index, rg, numThreads);
		if(VERBOSE >= 0) {
			fprintf(stderr, "Finishing cleaning up.\n");
			fprintf(stderr, "Index size is %.3lfGB.\n",
					RGIndexGetSize(&index, GIGABYTES));
		}

		/* Create the file name */
		sprintf(outputFileName, "%sblatter.index.file.%s.%d.%d.%d.%d.%d.%s",
				outputDir,
				outputID,
				index.startChr,
				index.startPos,
				index.endChr,
				index.endPos,
				i+1,
				BLATTER_INDEX_FILE_EXTENSION);

		/* Open the output file */
		if(!(fp=fopen(outputFileName, "wb"))) {
			PrintError("GenerateIndex",
					outputFileName,
					"Could not open outputFileName for writing",
					Exit,
					OpenFileError);
		}   

		/* Output the index to file */
		if(VERBOSE >= 0) {
			fprintf(stderr, "Outputting index to %s\n", outputFileName);
		}
		RGIndexPrint(fp, &index, binaryOutput);
		/* Close the output file */
		fclose(fp);

		/* Free memory */
		if(VERBOSE >= 0) {
			fprintf(stderr, "Deleting index...\n");
		}
		RGIndexDelete(&index);

		if(VERBOSE >= 0) {
			fprintf(stderr, "Generated index successfully!\n");
			fprintf(stderr, "%s", BREAK_LINE);
		}
	}
}
