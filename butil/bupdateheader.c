#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/RGIndex.h"
#include "bupdateheader.h"

#define Name "bupdateheader"

/* Updates the header of a bfast index file.  This is useful for 
 * converting old index files, with obsolete headers, as to work
 * with new bfast index files. */

int main(int argc, char *argv[]) 
{
	FILE *fp=NULL;
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	int missing;
	RGIndex index;

	if(argc == 3) {

		strcpy(inputFileName, argv[1]);
		missing = atoi(argv[2]);
		sprintf(outputFileName, "%s.%d.fixed",
				inputFileName,
				missing);

		/* Initialize */
		RGIndexInitialize(&index);

		fprintf(stderr, "Fixing %s.\n", inputFileName);
		if(!(fp=fopen(inputFileName, "rb"))) {
			PrintError(Name,
					inputFileName,
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}
		switch(missing) {
			case 0:
				/* Read in without the color space flag */
				Zero(fp, &index); 
				break;
			default:
				PrintError(Name,
						"missing",
						"Could not understand what is missing",
						Exit,
						OutOfRange);
				break;
		}
		fclose(fp);

		/* Print */
		fprintf(stderr, "Outputting to %s.\n", outputFileName);
		if(!(fp=fopen(outputFileName, "wb"))) {
			PrintError(Name,
					outputFileName,
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}
		RGIndexPrint(fp, &index, 1);
		fclose(fp);

		/* Free index */
		RGIndexDelete(&index);
	}
	else {
		fprintf(stderr, "%s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast index file>\n");
		fprintf(stderr, "\t<missing:\n\t\t0: color space flag>\n");
	}

	return 0;
}

/* TODO */
void Zero(FILE *fp, RGIndex *index)
{
	/* Read in header */
	if(fread(&index->length, sizeof(uint32_t), 1, fp)!=1
			|| fread(&index->hashWidth, sizeof(uint32_t), 1, fp)!=1
			|| fread(&index->hashLength, sizeof(int64_t), 1, fp)!=1
			|| fread(&index->totalLength, sizeof(int32_t), 1, fp)!=1
			|| fread(&index->numTiles, sizeof(int32_t), 1, fp)!=1
			|| fread(&index->repeatMasker, sizeof(int32_t), 1, fp)!=1
			|| fread(&index->startChr, sizeof(int32_t), 1, fp)!=1
			|| fread(&index->startPos, sizeof(int32_t), 1, fp)!=1
			|| fread(&index->endChr, sizeof(int32_t), 1, fp)!=1
			|| fread(&index->endPos, sizeof(int32_t), 1, fp)!=1) {
		PrintError("RGIndexReadHeader",
				NULL,
				"Could not read header",
				Exit,
				ReadFileError);
	}

	/* Add color space */
	index->colorSpace = 0;

	/* Error checking */
	assert(index->length > 0);
	assert(index->hashWidth > 0);
	assert(index->hashLength > 0);
	assert(index->totalLength > 0);
	assert(index->numTiles > 0);
	assert(index->repeatMasker == 0 || index->repeatMasker == 1);
	assert(index->colorSpace == 0 || index->colorSpace == 1);
	assert(index->startChr > 0);
	assert(index->startPos > 0);
	assert(index->endChr > 0);
	assert(index->endPos > 0);

	assert(index->length > 0);

	/* Allocate memory for the positions */
	index->positions = malloc(sizeof(uint32_t)*index->length);
	if(NULL == index->positions) {
		PrintError("RGIndexRead",
				"index->positions",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the chromosomes */
	index->chromosomes = malloc(sizeof(uint8_t)*index->length);
	if(NULL == index->chromosomes) {
		PrintError("RGIndexRead",
				"index->chromosomes",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the starts */
	index->starts = malloc(sizeof(uint32_t)*index->hashLength);
	if(NULL == index->starts) {
		PrintError("RGIndexRead",
				"index->starts",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Allocate memory for the ends */
	index->ends = malloc(sizeof(uint32_t)*index->hashLength);
	if(NULL == index->ends) {
		PrintError("RGIndexRead",
				"index->ends",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the tile lengths */
	index->tileLengths = malloc(sizeof(int32_t)*index->numTiles);
	if(NULL == index->tileLengths) {
		PrintError("RGIndexRead",
				"index->tileLengths",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the gaps */
	index->gaps = malloc(sizeof(int32_t)*(index->numTiles-1));
	if(NULL == index->gaps) {
		PrintError("RGIndexRead",
				"index->gaps",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Read in positions */
	if(fread(index->positions, sizeof(uint32_t), index->length, fp)!=index->length) {
		PrintError("RGIndexRead",
				NULL,
				"Could not read in positions",
				Exit,
				ReadFileError);
	}

	/* Read in the chromosomes */
	if(fread(index->chromosomes, sizeof(uint8_t), index->length, fp)!=index->length) {
		PrintError("RGIndexRead",
				NULL,
				"Could not read in chromosomes",
				Exit,
				ReadFileError);
	}

	/* Read in starts */
	if(fread(index->starts, sizeof(uint32_t), index->hashLength, fp)!=index->hashLength) {
		PrintError("RGIndexRead",
				NULL,
				"Could not read in starts",
				Exit,
				ReadFileError);
	}

	/* Read in ends */
	if(fread(index->ends, sizeof(uint32_t), index->hashLength, fp)!=index->hashLength) {
		PrintError("RGIndexRead",
				NULL,
				"Could not read in ends",
				Exit,
				ReadFileError);
	}

	/* Read the tileLengths */
	if(fread(index->tileLengths, sizeof(int32_t), index->numTiles, fp)!=index->numTiles) {
		PrintError("RGIndexRead",
				NULL,
				"Could not read in tile lengths",
				Exit,
				ReadFileError);
	}

	/* Read the gaps */
	if(fread(index->gaps, sizeof(int32_t), index->numTiles-1, fp)!= (index->numTiles-1)) {
		PrintError("RGIndexRead",
				NULL,
				"Could not read in gaps",
				Exit,
				ReadFileError);
	}

}

