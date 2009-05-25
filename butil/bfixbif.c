#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <zlib.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BLib.h"
#include "../blib/BError.h"
#include "../blib/RGIndex.h"
#include "../blib/AlignedRead.h"
#include "../blib/AlignedEnd.h"
#include "../blib/AlignedEntry.h"
#include "../blib/ScoringMatrix.h"
#include "bfixbif.h"

#define Name "bfixbif"

/* Updates a pre 0.5.0 bif to 0.5.0 */
int main(int argc, char *argv[]) 
{
	char inputFileName[MAX_FILENAME_LENGTH]="\0";

	if(2 <= argc) {
		int i;

		for(i=1;i<argc;i++) {
			strcpy(inputFileName, argv[i]);
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Updating %s.\n", inputFileName);
			if(NULL!=strstr(inputFileName,  BFAST_INDEX_FILE_EXTENSION)) {
				ConvertBIF(inputFileName);
			}
			else {
				PrintError(Name,
						"input file",
						"Could not recognize input file extension",
						Warn,
						OutOfRange);
			}
			fprintf(stderr, "%s", BREAK_LINE);
		}
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<input file>\n");
	}

	return 0;
}

void ConvertBIF(char *inputFileName) 
{
	char *FnName="ConvertBIF";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	RGIndex index;

	strcpy(outputFileName, inputFileName);
	strcat(outputFileName, ".new");
	assert(0!=strcmp(inputFileName, outputFileName));

	RGIndexInitialize(&index);
	RGIndexReadOld(&index, inputFileName);

	/* Fix package version */
	free(index.packageVersion);
	index.packageVersionLength = (int)strlen(PACKAGE_VERSION);
	index.packageVersion = malloc(sizeof(char)*(index.packageVersionLength+1));
	if(NULL==index.packageVersion) {
		PrintError(FnName,
				"index->packageVersion",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	strcpy(index.packageVersion, PACKAGE_VERSION);

	RGIndexPrint(outputFileName, &index);
	RGIndexDelete(&index);

	fprintf(stderr, "Outputted to %s.\n",
			outputFileName);
}

/* TODO */
void RGIndexReadOld(RGIndex *index, char *rgIndexFileName)
{
	char *FnName="RGIndexReadOld";

	FILE *fp;

	if(VERBOSE >= 0) {
		fprintf(stderr, "Reading index from %s.\n",
				rgIndexFileName);
	}

	/* open file */
	if((fp=fopen(rgIndexFileName, "r"))==0) {
		PrintError(FnName,
				rgIndexFileName,
				"Could not open rgIndexFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Read in the header */
	RGIndexReadHeaderOld(fp, index);

	assert(index->length > 0);

	/* Allocate memory for the positions */
	index->positions = malloc(sizeof(uint32_t)*index->length);
	if(NULL == index->positions) {
		PrintError(FnName,
				"index->positions",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the contigs */
	if(index->contigType == Contig_8) {
		index->contigs_8 = malloc(sizeof(uint8_t)*index->length);
		if(NULL == index->contigs_8) {
			PrintError(FnName,
					"index->contigs",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	else {
		index->contigs_32 = malloc(sizeof(uint32_t)*index->length);
		if(NULL == index->contigs_32) {
			PrintError(FnName,
					"index->contigs",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	/* Allocate memory for the starts */
	index->starts = malloc(sizeof(uint32_t)*index->hashLength);
	if(NULL == index->starts) {
		PrintError(FnName,
				"index->starts",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Read in positions */
	if(fread(index->positions, sizeof(uint32_t), index->length, fp)!=index->length) {
		PrintError(FnName,
				NULL,
				"Could not read in positions",
				Exit,
				ReadFileError);
	}

	/* Read in the contigs */
	if(index->contigType == Contig_8) {
		if(fread(index->contigs_8, sizeof(uint8_t), index->length, fp)!=index->length) {
			PrintError(FnName,
					NULL,
					"Could not read in contigs_8",
					Exit,
					ReadFileError);
		}
	}
	else {
		if(fread(index->contigs_32, sizeof(uint32_t), index->length, fp)!=index->length) {
			PrintError(FnName,
					NULL,
					"Could not read in contigs_32",
					Exit,
					ReadFileError);
		}
	}

	/* Read in starts */
	if(fread(index->starts, sizeof(uint32_t), index->hashLength, fp)!=index->hashLength) {
		PrintError(FnName,
				NULL,
				"Could not read in starts",
				Exit,
				ReadFileError);
	}

	/* close file */
	fclose(fp);

	if(VERBOSE >= 0) {
		fprintf(stderr, "Read index from %s.\n",
				rgIndexFileName);
	}
}

/* TODO */
void RGIndexReadHeaderOld(FILE *fp, RGIndex *index)
{
	char *FnName = "RGIndexReadHeaderOld";
	/* Read in header */
	if(fread(&index->id, sizeof(int32_t), 1, fp) != 1 ||
			fread(&index->packageVersionLength, sizeof(int32_t), 1, fp) != 1) {
		PrintError(FnName,
				NULL,
				"Could not read header",
				Exit,
				ReadFileError);
	}
	index->packageVersion = malloc(sizeof(char)*(index->packageVersionLength+1));
	if(NULL==index->packageVersion) {
		PrintError(FnName,
				"index->packageVersion",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	if(fread(index->packageVersion, sizeof(char), index->packageVersionLength, fp) != index->packageVersionLength ||
			fread(&index->length, sizeof(int64_t), 1, fp) != 1 ||
			fread(&index->contigType, sizeof(int32_t), 1, fp) != 1 ||
			fread(&index->startContig, sizeof(int32_t), 1, fp) != 1 ||
			fread(&index->startPos, sizeof(int32_t), 1, fp) != 1 ||
			fread(&index->endContig, sizeof(int32_t), 1, fp) != 1 ||
			fread(&index->endPos, sizeof(int32_t), 1, fp) != 1 ||
			fread(&index->width, sizeof(int32_t), 1, fp) != 1 ||
			fread(&index->keysize, sizeof(int32_t), 1, fp) != 1 ||
			fread(&index->repeatMasker, sizeof(int32_t), 1, fp) != 1 ||
			fread(&index->space, sizeof(int32_t), 1, fp) != 1 ||
			fread(&index->hashWidth, sizeof(uint32_t), 1, fp) != 1 ||
			fread(&index->hashLength, sizeof(int64_t), 1, fp) != 1) {
		PrintError(FnName,
				NULL,
				"Could not read header",
				Exit,
				ReadFileError);
	}
	index->packageVersion[index->packageVersionLength]='\0';
	/* Allocate memory for the mask */
	index->mask = malloc(sizeof(int32_t)*index->width);
	if(NULL==index->mask) {
		PrintError(FnName,
				"index->mask",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Read the mask */
	if(fread(index->mask, sizeof(int32_t), index->width, fp) != index->width) {
		PrintError(FnName,
				NULL,
				"Could not read header",
				Exit,
				ReadFileError);
	}

	/* Error checking */
	assert(index->id == (int)BFAST_ID);
	assert(index->length > 0);
	assert(index->contigType == Contig_8 || index->contigType == Contig_32);
	assert(index->startContig > 0);
	assert(index->startPos > 0);
	assert(index->endContig > 0);
	assert(index->endPos > 0);
	assert(index->width > 0);
	assert(index->keysize > 0);
	assert(index->repeatMasker == 0 || index->repeatMasker == 1);
	assert(index->space == NTSpace || index->space == ColorSpace);
	assert(index->hashWidth > 0);
	assert(index->hashLength > 0);
}
