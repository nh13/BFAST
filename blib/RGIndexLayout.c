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
	rgLayout->hashLengths=0;
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

	/* Allocate memory for the hashLengths */
	rgLayout->hashLengths = malloc(sizeof(int)*rgLayout->numIndexes);
	if(NULL==rgLayout->hashLengths) {
		PrintError("RGIndexLayoutRead",
				"rgLayout->hashLengths",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the number of tiles */
	rgLayout->numTiles = malloc(sizeof(int)*rgLayout->numIndexes);
	if(NULL==rgLayout->numTiles) {
		PrintError("RGIndexLayoutRead",
				"rgLayout->numTiles",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for tileLengths */
	rgLayout->tileLengths = malloc(sizeof(int*)*rgLayout->numIndexes);
	if(NULL==rgLayout->tileLengths) {
		PrintError("RGIndexLayoutRead",
				"rgLayout->tileLengths",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for gaps */
	rgLayout->gaps = malloc(sizeof(int*)*rgLayout->numIndexes);
	if(NULL==rgLayout->gaps) {
		PrintError("RGIndexLayoutRead",
				"rgLayout->gaps",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Read the indexes */
	for(i=0;i<rgLayout->numIndexes;i++) {
		/* Read in the hash length */
		if(EOF==fscanf(fp, "%u", &rgLayout->hashLengths[i])) {
			PrintError("RGIndexLayoutRead",
					NULL,
					"Could not read the hashLengths",
					Exit,
					EndOfFile);
		}

		/* Read in the number of tiles */
		if(EOF==fscanf(fp, "%d", &rgLayout->numTiles[i])) {
			PrintError("RGIndexLayoutRead",
					NULL,
					"Could not read the number of tiles",
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

	/* Free the tile lengths for each layout*/
	for(i=0;i<rgLayout->numIndexes;i++) {
		free(rgLayout->tileLengths[i]);
		rgLayout->tileLengths[i]=NULL;
	}
	/* Free the gaps */
	for(i=0;i<rgLayout->numIndexes;i++) {
		free(rgLayout->gaps[i]);
		rgLayout->gaps[i]=NULL;
	}
	free(rgLayout->hashLengths);
	rgLayout->hashLengths = NULL;
	free(rgLayout->tileLengths);
	rgLayout->tileLengths = NULL;
	free(rgLayout->gaps);
	rgLayout->gaps = NULL;
	/* Free the numTiles */
	free(rgLayout->numTiles);
	rgLayout->numTiles = NULL;
	/* Reinitialize the number of indexes */
	rgLayout->numIndexes = 0;
}
