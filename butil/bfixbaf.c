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
#include "bfixbaf.h"

#define Name "bfixbaf"

/* Updates a 0.4.1 to 0.4.7 BAF file to 0.5.0 */

/* 0.4.1 - 0.4.5 will update MQ */

int main(int argc, char *argv[]) 
{
	char scoringMatrixFileName[MAX_FILENAME_LENGTH]="\0";
	int avgMismatchQuality;
	int space;
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	char oldVersionString[MAX_FILENAME_LENGTH]="\0";
	int oldVersion[3]={0,0,0};

	if(6 <= argc) {
		int i;
		strcpy(oldVersionString, argv[1]);
		strcpy(scoringMatrixFileName, argv[2]);
		avgMismatchQuality = atoi(argv[3]);
		space = atoi(argv[4]);

		/* Check version */
		ParsePackageVersion(oldVersionString, 
				&oldVersion[0],
				&oldVersion[1],
				&oldVersion[2]);

		if(0 != oldVersion[0] ||
				4 != oldVersion[1] ||
				oldVersion[2] < 1 ||
				7 < oldVersion[2]) {
			PrintError(Name,
					oldVersionString,
					"Version was out of range",
					Exit,
					OutOfRange);
		}

		for(i=5;i<argc;i++) {
			strcpy(inputFileName, argv[i]);
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Updating %s.\n", inputFileName);
			if(NULL!=strstr(inputFileName,  BFAST_ALIGNED_FILE_EXTENSION)) {
				ConvertBAF(inputFileName,
						scoringMatrixFileName,
						avgMismatchQuality,
						space,
						oldVersion[2]);
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
		fprintf(stderr, "\t<old version number (ex. 0.4.1)>\n");
		fprintf(stderr, "\t<scoring matrix file name>\n");
		fprintf(stderr, "\t<average mismatch quality>\n");
		fprintf(stderr, "\t<space 0: NT Space 1: Color Space>\n");
		fprintf(stderr, "\t<input file>\n");
	}

	return 0;
}

void ConvertBAF(char *inputFileName,
		char*scoringMatrixFileName,
		int32_t avgMismatchQuality,
		int32_t space,
		int32_t minorVersionNumber)
{
	char *FnName="ConvertBAF";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fpIn = NULL;
	gzFile fpOut = NULL;
	AlignedRead a;
	double mismatchScore;
	ScoringMatrix sm;

	/* Read in scoring matrix */
	ScoringMatrixInitialize(&sm);
	ScoringMatrixRead(scoringMatrixFileName, &sm, space);
	/* Calculate mismatch score */
	/* Assumes all match scores are the same and all substitution scores are the same */
	if(space == NTSpace) {
		mismatchScore = sm.ntMatch - sm.ntMismatch;
	}
	else {
		mismatchScore = sm.colorMatch - sm.colorMismatch;
	}

	strcpy(outputFileName, inputFileName);
	strcat(outputFileName, ".new");
	assert(0!=strcmp(inputFileName, outputFileName));

	if(!(fpIn = fopen(inputFileName, "rb"))) {
		PrintError(FnName,
				inputFileName,
				"Could not open file for reading",
				Exit,
				OpenFileError);
	}
	if(!(fpOut = gzopen(outputFileName, "wb"))) {
		PrintError(FnName,
				outputFileName,
				"Could not open file for writing",
				Exit,
				OpenFileError);
	}

	AlignedReadInitialize(&a);
	if(1 <= minorVersionNumber && minorVersionNumber <= 5) {
		while(0 < AlignedReadReadOld(&a, fpIn, 0)) {
			AlignedReadUpdateMappingQuality(&a, 
					mismatchScore,
					avgMismatchQuality);
			AlignedReadPrint(&a, fpOut);
			AlignedReadFree(&a);
		}
	}
	else {
		while(0 < AlignedReadReadOld(&a, fpIn, 1)) {
			AlignedReadPrint(&a, fpOut);
		}
	}

	fclose(fpIn);
	gzclose(fpOut);

	fprintf(stderr, "Outputted to %s.\n",
			outputFileName);
}

/* TODO */
int32_t AlignedReadReadOld(AlignedRead *a,
		FILE *inputFP,
		int32_t withMQ)
{
	char *FnName = "AlignedReadReadOld";
	int32_t i;

	assert(a != NULL);

	/* Allocate memory for the read name */
	a->readName = malloc(sizeof(char)*SEQUENCE_NAME_LENGTH);
	if(a->readName == NULL) {
		if(NULL == a->readName) {
			PrintError(FnName,
					"a->readName",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}

	/* Read the read name, paired end flag, space flag, and the number of entries for both entries */
	if(fread(&a->readNameLength, sizeof(int32_t), 1, inputFP) != 1) {
		/* Free read name before leaving */
		free(a->readName);
		a->readName=NULL;
		return EOF;
	}
	if(fread(a->readName, sizeof(char), a->readNameLength, inputFP) != a->readNameLength ||
			fread(&a->space, sizeof(int32_t), 1, inputFP) != 1 ||
			fread(&a->numEnds, sizeof(int32_t), 1, inputFP) != 1) {
		PrintError(FnName,
				NULL,
				"Could not read from file",
				Exit,
				ReadFileError);
	}
	/* Add the null terminator */
	a->readName[a->readNameLength]='\0';
	/* Reallocate to conserve memory */
	if(0 < a->readNameLength) {
		a->readName = realloc(a->readName, sizeof(char)*(a->readNameLength+1));
		if(NULL == a->readName) {
			PrintError(FnName,
					"a->readName",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
	}
	else {
		free(a->readName);
		a->readName=NULL;
	}

	/* Allocate memory for the ends */
	a->ends = malloc(sizeof(AlignedEnd)*a->numEnds);
	if(NULL==a->ends) {
		PrintError(FnName,
				"a->ends",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Read the alignment */
	for(i=0;i<a->numEnds;i++) {
		AlignedEndInitialize(&a->ends[i]);
		if(EOF==AlignedEndReadOld(&a->ends[i],
					inputFP,
					a->space,
					withMQ)) {
			PrintError(FnName,
					NULL,
					"Could not read a->ends[i]",
					Exit,
					EndOfFile);
		}
	}

	return 1;
}

/* TODO */
int32_t AlignedEndReadOld(AlignedEnd *a,
		FILE *inputFP,
		int32_t space,
		int32_t withMQ)
{
	char *FnName = "AlignedEndReadOld";
	int32_t i;

	if(fread(&a->readLength, sizeof(int32_t), 1, inputFP) != 1 ||
			fread(&a->qualLength, sizeof(int32_t), 1, inputFP) != 1) {
		return EOF;
	}
	/* Allocate memory for the alignment */
	if(a->read == NULL) {
		a->read = malloc(sizeof(char)*(1+a->readLength));
		if(NULL == a->read) {
			PrintError(FnName,
					"a->read",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	if(a->qual == NULL) {
		a->qual = malloc(sizeof(char)*(1+a->qualLength));
		if(NULL == a->qual) {
			PrintError(FnName,
					"a->qual",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}

	if(fread(a->read, sizeof(char), a->readLength, inputFP) != a->readLength ||
			fread(a->qual, sizeof(char), a->qualLength, inputFP) != a->qualLength ||
			fread(&a->numEntries, sizeof(int32_t), 1, inputFP) != 1) {
		PrintError(FnName,
				"a->reads, a->qual, and a->numEntries",
				"Could not read from file",
				Exit,
				ReadFileError);
	}
	/* Add the null terminator to strings */
	a->read[a->readLength]='\0';
	a->qual[a->qualLength]='\0';
	AlignedEndReallocate(a,
			a->numEntries);

	for(i=0;i<a->numEntries;i++) {
		AlignedEntryInitialize(&a->entries[i]);
		if(EOF == AlignedEntryReadOld(&a->entries[i],
					inputFP,
					space,
					withMQ)) {
			PrintError(FnName,
					"a->entries[i]",
					"Could not read from file",
					Exit,
					ReadFileError);
		}
	}

	return 1;
}

/* TODO */
int32_t AlignedEntryReadOld(AlignedEntry *a,
		FILE *inputFP,
		int32_t space,
		int32_t withMQ)
{
	char *FnName = "AlignedEntryReadOld";

	assert(NULL != a);

	/* Allocate memory for the alignment */
	if(a->read == NULL) {
		a->read = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == a->read) {
			PrintError(FnName,
					"a->read",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	if(a->reference == NULL) {
		a->reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == a->reference) {
			PrintError(FnName,
					"a->reference",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	if(space == ColorSpace) {
		if(a->colorError == NULL) {
			a->colorError = malloc(sizeof(char)*SEQUENCE_LENGTH);
			if(NULL == a->colorError) {
				PrintError(FnName,
						"a->colorError",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}
	}

	if(fread(&a->contigNameLength, sizeof(int32_t), 1, inputFP) != 1) {
		return EOF;
	}
	/* Copy over contig name */
	a->contigName = malloc(sizeof(char)*(a->contigNameLength+1));
	if(NULL==a->contigName) {
		PrintError(FnName,
				"a->contigName",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	if(fread(a->contigName, sizeof(char), a->contigNameLength, inputFP) != a->contigNameLength ||
			fread(&a->contig, sizeof(uint32_t), 1, inputFP) != 1 ||
			fread(&a->position, sizeof(uint32_t), 1, inputFP) != 1 ||
			fread(&a->strand, sizeof(char), 1, inputFP) != 1 ||
			fread(&a->score, sizeof(double), 1, inputFP) != 1) {
		return EOF;
	}

	/* This is the only difference between 0.4.1-0.4.5 and 0.4.6-0.4.7 */
	if(1 == withMQ &&
			fread(&a->mappingQuality, sizeof(int32_t), 1, inputFP) != 1) {
		return EOF;
	}

	if(fread(&a->referenceLength, sizeof(uint32_t), 1, inputFP) != 1 ||
			fread(&a->length, sizeof(uint32_t), 1, inputFP) != 1 ||
			fread(a->read, sizeof(char), a->length, inputFP) != a->length ||
			fread(a->reference, sizeof(char), a->length, inputFP) != a->length) {
		return EOF;
	}
	/* Add the null terminator to strings */
	a->contigName[a->contigNameLength]='\0';
	a->read[a->length]='\0';
	a->reference[a->length]='\0';
	if(ColorSpace==space) {
		if(fread(a->colorError, sizeof(char), a->length, inputFP) != a->length) {
			return EOF;
		}
		a->colorError[a->length]='\0';
	} /* Reallocate to conserve memory */
	assert(a->length > 0);
	a->read = realloc(a->read, sizeof(char)*(a->length+1));
	if(NULL == a->read) {
		PrintError(FnName,
				"a->read",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}
	/* Reference */
	a->reference = realloc(a->reference, sizeof(char)*(a->length+1));
	if(NULL == a->reference) {
		PrintError(FnName,
				"a->reference",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}
	/* Color error, if necessary */
	if(space == ColorSpace) {
		a->colorError = realloc(a->colorError, sizeof(char)*(a->length+1));
		if(NULL == a->colorError) {
			PrintError(FnName,
					"a->colorError",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
	}
	else {
		assert(NULL == a->colorError);
	}

	return 1;
}
