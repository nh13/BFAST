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
#include "../blib/RGBinary.h"
#include "../blib/AlignedRead.h"
#include "../blib/AlignedEnd.h"
#include "../blib/AlignedEntry.h"
#include "../blib/ScoringMatrix.h"
#include "bfixbrg.h"

#define Name "bfixbrg"

/* Updates a pre 0.5.0 brg to 0.5.0 */
int main(int argc, char *argv[]) 
{
	char inputFileName[MAX_FILENAME_LENGTH]="\0";

	if(2 <= argc) {
		int i;

		for(i=1;i<argc;i++) {
			strcpy(inputFileName, argv[i]);
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Updating %s.\n", inputFileName);
			if(NULL!=strstr(inputFileName,  BFAST_RG_FILE_EXTENSION)) {
				ConvertBRG(inputFileName);
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

void ConvertBRG(char *inputFileName) 
{
	char *FnName="ConvertBRG";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	RGBinary rg;

	strcpy(outputFileName, inputFileName);
	strcat(outputFileName, ".new");
	assert(0!=strcmp(inputFileName, outputFileName));

	RGBinaryReadBinaryOld(&rg, inputFileName);

	/* Fix package version */
	free(rg.packageVersion);
	rg.packageVersionLength = (int)strlen(PACKAGE_VERSION);
	rg.packageVersion = malloc(sizeof(char)*(rg.packageVersionLength+1));
	if(NULL==rg.packageVersion) {
		PrintError(FnName,
				"rg->packageVersion",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	strcpy(rg.packageVersion, PACKAGE_VERSION);

	RGBinaryWriteBinary(&rg, outputFileName);
	RGBinaryDelete(&rg);

	fprintf(stderr, "Outputted to %s.\n",
			outputFileName);
}

/* Read from file */
void RGBinaryReadBinaryOld(RGBinary *rg,
		char *rgFileName)
{
	char *FnName="RGBinaryReadBinaryOld";
	FILE *fpRG;
	int32_t i;
	int32_t numCharsPerByte;
	int64_t numPosRead=0;
	/* We assume that we can hold 2 [acgt] (nts) in each byte */
	assert(ALPHABET_SIZE==4);
	numCharsPerByte=ALPHABET_SIZE/2;

	if(VERBOSE>=0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Reading in reference genome from %s.\n", rgFileName);
	}

	/* Open output file */
	if((fpRG=fopen(rgFileName, "rb"))==0) {
		PrintError(FnName,
				rgFileName,
				"Could not open rgFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Read RGBinary information */
	if( fread(&rg->id, sizeof(int32_t), 1, fpRG)!=1 ||
			fread(&rg->packageVersionLength, sizeof(int32_t), 1, fpRG)!=1) {
		PrintError(FnName,
				NULL,
				"Could not read RGBinary information",
				Exit,
				ReadFileError);
	}
	assert(0<rg->packageVersionLength);
	rg->packageVersion = malloc(sizeof(char)*(rg->packageVersionLength+1));
	if(NULL==rg->packageVersion) {
		PrintError(FnName,
				"rg->packageVersion",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	if(fread(rg->packageVersion, sizeof(char), rg->packageVersionLength, fpRG)!=rg->packageVersionLength ||
			fread(&rg->numContigs, sizeof(int32_t), 1, fpRG)!=1 ||
			fread(&rg->space, sizeof(int32_t), 1, fpRG)!=1) {
		PrintError(FnName,
				NULL,
				"Could not read RGBinary information",
				Exit,
				ReadFileError);
	}
	rg->packageVersion[rg->packageVersionLength]='\0';
	/* Check id */
	if(BFAST_ID != rg->id) {
		PrintError(FnName,
				"rg->id",
				"The id did not match",
				Exit,
				OutOfRange);
	}
	/*
	CheckPackageCompatibility(rg->packageVersion,
			BFASTReferenceGenomeFile);
			*/

	assert(rg->numContigs > 0);
	assert(rg->space == NTSpace|| rg->space == ColorSpace);
	rg->packed = RGBinaryPacked;

	/* Allocate memory for the contigs */
	rg->contigs = malloc(sizeof(RGBinaryContig)*rg->numContigs);
	if(NULL==rg->contigs) {
		PrintError(FnName,
				"rg->contigs",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Read each contig */
	for(i=0;i<rg->numContigs;i++) {
		/* Read contig name length */
		if(fread(&rg->contigs[i].contigNameLength, sizeof(int32_t), 1, fpRG)!=1) {
			PrintError(FnName,
					NULL,
					"Could not read contig name length",
					Exit,
					ReadFileError);
		}
		assert(rg->contigs[i].contigNameLength > 0);
		/* Allocate memory */
		rg->contigs[i].contigName = malloc(sizeof(char)*(rg->contigs[i].contigNameLength+1));
		if(NULL==rg->contigs[i].contigName) {
			PrintError(FnName,
					"contigName",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Read RGContig information */
		if(fread(rg->contigs[i].contigName, sizeof(char), rg->contigs[i].contigNameLength, fpRG) != rg->contigs[i].contigNameLength ||
				fread(&rg->contigs[i].sequenceLength, sizeof(int32_t), 1, fpRG)!=1 ||
				fread(&rg->contigs[i].numBytes, sizeof(uint32_t), 1, fpRG)!=1) {
			PrintError(FnName,
					NULL,
					"Could not read RGContig information",
					Exit,
					ReadFileError);
		}

		/* It should be packed */
		assert(numCharsPerByte == (rg->contigs[i].sequenceLength + (rg->contigs[i].sequenceLength % 2))/rg->contigs[i].numBytes);
		/* Add null terminator */
		rg->contigs[i].contigName[rg->contigs[i].contigNameLength]='\0';
		/* Allocate memory for the sequence */
		rg->contigs[i].sequence = malloc(sizeof(char)*rg->contigs[i].numBytes);
		if(NULL==rg->contigs[i].sequence) {
			PrintError(FnName,
					"rg->contigs[i].sequence",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Read sequence */
		if(fread(rg->contigs[i].sequence, sizeof(char), rg->contigs[i].numBytes, fpRG)!=rg->contigs[i].numBytes) {
			PrintError(FnName,
					NULL,
					"Could not read sequence",
					Exit,
					ReadFileError);
		}

		numPosRead += rg->contigs[i].sequenceLength;
	}


	/* Close the output file */
	fclose(fpRG);

	if(VERBOSE>=0) {
		fprintf(stderr, "In total read %d contigs for a total of %lld bases\n",
				rg->numContigs,
				(long long int)numPosRead);
		fprintf(stderr, "%s", BREAK_LINE);
	}
}
