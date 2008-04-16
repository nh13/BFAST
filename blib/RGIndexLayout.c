#include <stdio.h>
#include <stdlib.h>
#include "../blib/BError.h"
#include "../blib/BLibDefinitions.h"
#include "RGIndexLayout.h"

/* TODO */
void RGIndexLayoutRead(char *layoutFileName, RGIndexLayout *rgLayout)
{
	int i, j;
	FILE *fp=NULL;

	/* Open the file */
	if((fp=fopen(layoutFileName, "r"))==0) {
		PrintError("RGIndexLayoutRead",
				NULL,
				"Could not open index layout file reading",
				Exit,
				OpenFileError);
	}

	/* Initialize the layout data structure */
	rgLayout->numIndexes=0;
	rgLayout->numTiles=NULL;
	rgLayout->tileLengths=NULL;
	rgLayout->gaps=NULL;

	/* Read in the number of indexes */
	if(EOF==fscanf(fp, "%d", &rgLayout->numIndexes)) {
		PrintError("RGIndexLayoutRead",
				NULL,
				"Could not read the number of indexes",
				Exit,
				EndOfFile);
	}

	/* Allocate memory for the number of tiles */
	rgLayout->numTiles = realloc(rgLayout->numTiles, sizeof(int)*rgLayout->numIndexes);
	if(NULL==rgLayout->numTiles) {
		PrintError("RGIndexLayoutRead",
				"rgLayout->numTiles",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for tileLengths */
	rgLayout->tileLengths = realloc(rgLayout->tileLengths, sizeof(int*)*rgLayout->numIndexes);
	if(NULL==rgLayout->tileLengths) {
		PrintError("RGIndexLayoutRead",
				"rgLayout->tileLengths",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for gaps */
	rgLayout->gaps = realloc(rgLayout->gaps, sizeof(int*)*rgLayout->numIndexes);
	if(NULL==rgLayout->gaps) {
		PrintError("RGIndexLayoutRead",
				"rgLayout->gaps",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Read the indexes */
	for(i=0;i<rgLayout->numIndexes;i++) {
		/* Read in the number of tiles */
		if(EOF==fscanf(fp, "%d", &rgLayout->numTiles[i])) {
			PrintError("RGIndexLayoutRead",
					NULL,
					"Could not read the number of indexes",
					Exit,
					EndOfFile);
		}

		/* Allocate memory for the individual tile lengths and gaps */
		rgLayout->tileLengths[i] = malloc(sizeof(int)*rgLayout->numTiles[i]);
		if(NULL==rgLayout->tileLengths[i]) {
			PrintError("RGIndexLayoutRead",
					"rgLayout->tileLengths[i]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		rgLayout->gaps[i] = malloc(sizeof(int)*(rgLayout->numTiles[i]-1));
		if(NULL==rgLayout->gaps[i]) {
			PrintError("RGIndexLayoutRead",
					"rgLayout->gaps[i]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}

		/* Read in the tiles and gaps */
		for(j=0;j<rgLayout->numTiles[i];j++) {
			if(EOF==fscanf(fp, "%d", &rgLayout->tileLengths[i][j])) {
				PrintError("RGIndexLayoutRead",
						"rgLayout->tileLengths[i][j]",
						"Could not read in tile length",
						Exit,
						EndOfFile);
			}
			if(j<rgLayout->numTiles[i]-1) {
				if(EOF==fscanf(fp, "%d", &rgLayout->gaps[i][j])) {
					PrintError("RGIndexLayoutRead",
							"rgLayout->gaps[i][j]",
							"Could not read in tile length",
							Exit,
							EndOfFile);
				}
			}
		}
	}

	/* Close the file */
	fclose(fp);
}

void RGIndexLayoutDelete(RGIndexLayout *rgLayout)
{
	int i;

	/* Free the tile lengths and gaps */
	for(i=0;i<rgLayout->numIndexes;i++) {
		free(rgLayout->tileLengths[i]);
		free(rgLayout->gaps[i]);
	}
	free(rgLayout->tileLengths);
	free(rgLayout->gaps);
	rgLayout->tileLengths = NULL;
	rgLayout->gaps = NULL;
	/* Free the numTiles */
	free(rgLayout->numTiles);
	rgLayout->numTiles = NULL;
	/* Reinitialize the number of indexes */
	rgLayout->numIndexes = 0;
}
