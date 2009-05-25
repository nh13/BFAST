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
#include "../blib/RGMatches.h"
#include "../blib/RGMatch.h"
#include "bfixbmf.h"

#define Name "bfixbmf"

/* Updates a pre 0.5.0 bmf to 0.5.0 */
int main(int argc, char *argv[]) 
{
	char inputFileName[MAX_FILENAME_LENGTH]="\0";

	if(2 <= argc) {
		int i;

		for(i=1;i<argc;i++) {
			strcpy(inputFileName, argv[i]);
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Updating %s.\n", inputFileName);
			if(NULL!=strstr(inputFileName, BFAST_MATCHES_FILE_EXTENSION)) {
				ConvertBMF(inputFileName);
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

void ConvertBMF(char *inputFileName) 
{
	char *FnName="ConvertBMF";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fpIn = NULL;
	gzFile fpOut = NULL;
	RGMatches m;

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

	RGMatchesInitialize(&m);
	while(EOF != RGMatchesReadOld(fpIn, &m)) {
		RGMatchesPrint(fpOut, &m);
		RGMatchesFree(&m);
	}

	fclose(fpIn);
	gzclose(fpOut);

	fprintf(stderr, "Outputted to %s.\n",
			outputFileName);
}

/* TODO */
int32_t RGMatchesReadOld(FILE *fp,
		RGMatches *m)
{
	char *FnName = "RGMatchesReadOld";
	int32_t i;

	/* Read the matches from the input file */
	/* Read read name length */
	if(fread(&m->readNameLength, sizeof(int32_t), 1, fp) != 1) {
		return EOF;
	}
	assert(m->readNameLength < SEQUENCE_NAME_LENGTH);
	assert(m->readNameLength > 0);

	/* Allocate memory for the read name */
	m->readName = malloc(sizeof(char)*(m->readNameLength + 1));
	if(NULL == m->readName) {
		PrintError(FnName,
				"m->readName",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Read in read name */
	if(fread(m->readName, sizeof(char), m->readNameLength, fp)!=m->readNameLength) {
		PrintError(FnName,
				"m->readName",
				"Could not read in read name",
				Exit,
				ReadFileError);
	}
	m->readName[m->readNameLength]='\0';
	/* Read numEnds */
	if(fread(&m->numEnds, sizeof(int32_t), 1, fp)!=1) {
		PrintError(FnName,
				"numEnds",
				"Could not read in numEnds",
				Exit,
				ReadFileError);
	}

	/* Allocate the ends */
	m->ends = malloc(sizeof(RGMatch)*m->numEnds);
	if(NULL == m->ends) {
		PrintError(FnName,
				"m->ends",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Read each end */
	for(i=0;i<m->numEnds;i++) {
		/* Initialize */
		RGMatchInitialize(&m->ends[i]);
		/* Read */
		RGMatchReadOld(fp,
				&m->ends[i]);
	}

	return 1;
}

/* TODO */
int32_t RGMatchReadOld(FILE *fp,
		RGMatch *m)
{

	char *FnName = "RGMatchReadOld";

	/* Read in the read length */
	if(fread(&m->readLength, sizeof(int32_t), 1, fp) != 1 ||
			fread(&m->qualLength, sizeof(int32_t), 1, fp) != 1) {
		if(feof(fp) != 0) {
			return EOF;
		}
		else {
			PrintError(FnName,
					"m->readLength",
					"Could not read in read length",
					Exit,
					ReadFileError);
		}
	}
	assert(m->readLength < SEQUENCE_LENGTH);
	assert(m->readLength > 0);

	/* Allocate memory for the read */
	m->read = malloc(sizeof(char)*(m->readLength+1));
	if(NULL==m->read) {
		PrintError(FnName,
				"read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	m->qual = malloc(sizeof(char)*(m->qualLength+1));
	if(NULL==m->qual) {
		PrintError(FnName,
				"qual",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Read in the read */
	if(fread(m->read, sizeof(char), m->readLength, fp)!=m->readLength ||
			fread(m->qual, sizeof(char), m->qualLength, fp)!=m->qualLength) {
		PrintError(FnName,
				"m->read",
				"Could not read in the read and qual",
				Exit,
				ReadFileError);
	}
	m->read[m->readLength]='\0';
	m->qual[m->qualLength]='\0';

	/* Read in if we have reached the maximum number of matches */
	if(fread(&m->maxReached, sizeof(int32_t), 1, fp)!=1) {
		PrintError(FnName,
				"m->maxReached",
				"Could not read in m->maxReached",
				Exit,
				ReadFileError);
	}
	assert(0 == m->maxReached || 1 == m->maxReached);

	/* Read in the number of matches */
	if(fread(&m->numEntries, sizeof(int32_t), 1, fp)!=1) {
		PrintError(FnName,
				"m->numEntries",
				"Could not read in m->numEntries",
				Exit,
				ReadFileError);
	}
	assert(m->numEntries >= 0);

	/* Allocate memory for the matches */
	RGMatchReallocate(m, m->numEntries);

	/* Read first sequence matches */
	if(fread(m->contigs, sizeof(uint32_t), m->numEntries, fp)!=m->numEntries) {
		PrintError(FnName,
				"m->contigs",
				"Could not read in contigs",
				Exit,
				ReadFileError);
	}
	if(fread(m->positions, sizeof(uint32_t), m->numEntries, fp)!=m->numEntries) {
		PrintError(FnName,
				"m->positions",
				"Could not read in positions",
				Exit,
				ReadFileError);
	}
	if(fread(m->strands, sizeof(char), m->numEntries, fp)!=m->numEntries) {
		PrintError(FnName,
				"m->strands",
				"Could not read in strand",
				Exit,
				ReadFileError);
	}

	return 1;
}
