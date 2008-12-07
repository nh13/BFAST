#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/RGIndex.h"
#include "bfix.h"

#define Name "bfix"

/* Updates an old bfast reference genome file or bfast index file
 * to the newest format. */

int main(int argc, char *argv[]) 
{
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	int version=V_0_1_13;

	if(3 == argc) {
		strcpy(inputFileName, argv[1]);
		version = atoi(argv[2]);
		assert(version == V_0_1_13); /* Only one supported */

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Updating %s.\n", inputFileName);
		if(NULL!=strstr(inputFileName, BFAST_RG_FILE_EXTENSION)) {
			ConvertRGBinary(inputFileName, version);
		}
		else if(NULL!=strstr(inputFileName, BFAST_INDEX_FILE_EXTENSION)) {
			ConvertRGIndex(inputFileName, version);
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
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<input file>\n");
		fprintf(stderr, "\t<version 0: <= 0.1.13>\n");
	}

	return 0;
}

void ConvertRGBinary(char *inputFileName,
		int version)
{
	char *FnName="ConvertRGBinary";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";

	strcpy(outputFileName, inputFileName);
	strcat(outputFileName, ".new");
	assert(0!=strcmp(inputFileName, outputFileName));
	switch(version) {
		case V_0_1_13:
			ConvertRGBinaryFrom_0_1_13(inputFileName, outputFileName);
			break;
		default:
			PrintError(FnName,
					"version",
					"Could not recognize version",
					Exit,
					OutOfRange);
	}
	fprintf(stderr, "Outputted to %s.\n",
			outputFileName);
}

void ConvertRGIndex(char *inputFileName,
		int version)
{
	char *FnName="ConvertRGIndex";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";

	strcpy(outputFileName, inputFileName);
	strcat(outputFileName, ".new");
	assert(0!=strcmp(inputFileName, outputFileName));
	switch(version) {
		case V_0_1_13:
			ConvertRGIndexFrom_0_1_13(inputFileName, outputFileName);
			break;
		default:
			PrintError(FnName,
					"version",
					"Could not recognize version",
					Exit,
					OutOfRange);
	}
	fprintf(stderr, "Outputted to %s.\n",
			outputFileName);
}

/* 0.1.13 Functions */
void ConvertRGBinaryFrom_0_1_13(char *inputFileName,
		char *outputFileName) 
{
	/* Changes:
	 * - added package version
	 *   */
	char *FnName="ConvertRGBinaryFrom_0_1_13";
	RGBinary_0_1_13 rg_0_1_13;
	RGBinary rg;
	/* Read old */
	RGBinaryReadBinary_0_1_13(&rg_0_1_13,
			inputFileName);

	/* Free unconserved fields */

	/* Copy over conserved fields */
	rg.id = rg_0_1_13.id;
	rg.contigs = rg_0_1_13.contigs;
	rg.numContigs = rg_0_1_13.numContigs;
	rg.space = rg_0_1_13.space;

	/* Modify changed fields */
	rg.packed = RGBinaryPacked;

	/* Add new fields */
	rg.packageVersionLength = (int)strlen(PACKAGE_VERSION);
	rg.packageVersion = malloc(sizeof(int8_t)*(rg.packageVersionLength+1));
	if(NULL==rg.packageVersion) {
		PrintError(FnName,
				"rg.packageVersion",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	strcpy((char*)rg.packageVersion, PACKAGE_VERSION);

	/* Print new */
	RGBinaryWriteBinary(&rg, outputFileName);

	/* Free memory */
	RGBinaryDelete(&rg);
}

void ConvertRGIndexFrom_0_1_13(char *inputFileName,
		char *outputFileName) 
{
	/* Changes:
	 * - added package version
	 * - new hash 
	 *   */
	char *FnName="ConvertRGIndexFrom_0_1_13";
	FILE *fp=NULL;
	RGIndex_0_1_13 index_0_1_13;
	RGIndex index;
	int64_t i;
	uint32_t prevStart=0;

	/* Read old */
	if(!(fp=fopen(inputFileName, "rb"))) {
		PrintError(FnName,
				inputFileName,
				"Could not open file for reading",
				Exit,
				OpenFileError);
	}
	RGIndexRead_0_1_13(fp, &index_0_1_13, BinaryInput);
	fclose(fp);
	fp=NULL;

	/* Free unconserved fields */
	free(index_0_1_13.ends);
	index_0_1_13.ends=NULL;

	/* Copy over conserved fields */
	index.id = index_0_1_13.id;
	index.contigs_8 = index_0_1_13.contigs_8;
	index.contigs_32 = index_0_1_13.contigs_32;
	index.positions = index_0_1_13.positions;
	index.length = index_0_1_13.length;
	index.contigType = index_0_1_13.contigType;
	index.startContig = index_0_1_13.startContig;
	index.startPos = index_0_1_13.startPos;
	index.endContig = index_0_1_13.endContig;
	index.endPos = index_0_1_13.endPos;
	index.width = index_0_1_13.width;
	index.keysize = index_0_1_13.keysize;
	index.mask = index_0_1_13.mask;
	index.repeatMasker = index_0_1_13.repeatMasker;
	index.space = index_0_1_13.space;
	index.hashWidth = index_0_1_13.hashWidth;
	index.hashLength = index_0_1_13.hashLength;
	index.starts = index_0_1_13.starts;

	/* Modify changed fields */
	/* Update starts */
	/* Go through hash and reset all UINT_MAX starts */
	for(i=index.hashLength-1, prevStart=UINT_MAX;
			0<=i;
			i--) {
		if(UINT_MAX == index.starts[i]) {
			index.starts[i] = prevStart;
		}
		else {
			prevStart = index.starts[i];
		}
	}


	/* Add new fields */
	index.packageVersionLength = (int)strlen(PACKAGE_VERSION);
	index.packageVersion = malloc(sizeof(int8_t)*(index.packageVersionLength+1));
	if(NULL==index.packageVersion) {
		PrintError(FnName,
				"index.packageVersion",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	strcpy((char*)index.packageVersion, PACKAGE_VERSION);

	/* Print new */
	if(!(fp=fopen(outputFileName, "wb"))) {
		PrintError(FnName,
				outputFileName,
				"Could not open file for writing",
				Exit,
				OpenFileError);
	}
	RGIndexPrint(fp, &index);
	fclose(fp);

	/* Free memory */
	RGIndexDelete(&index);
}

void RGBinaryReadBinary_0_1_13(RGBinary_0_1_13 *rg,
		char *rgFileName)
{       
	char *FnName="RGBinaryReadBinary";
	FILE *fpRG;
	int i;
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
			fread(&rg->numContigs, sizeof(int32_t), 1, fpRG)!=1 ||
			fread(&rg->space, sizeof(int32_t), 1, fpRG)!=1) {
		PrintError(FnName,
				NULL,
				"Could not read RGBinary information",
				Exit,
				ReadFileError);
	}

	/* Check id */
	if(BFAST_ID != rg->id) {
		PrintError(FnName,
				"rg->id",
				"The id did not match",
				Exit,
				OutOfRange);
	}

	assert(rg->numContigs > 0);
	assert(rg->space == NTSpace|| rg->space == ColorSpace);

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
		rg->contigs[i].contigName = malloc(sizeof(int8_t)*(rg->contigs[i].contigNameLength+1));
		if(NULL==rg->contigs[i].contigName) {
			PrintError(FnName,
					"contigName",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Read RGContig information */
		if(fread(rg->contigs[i].contigName, sizeof(int8_t), rg->contigs[i].contigNameLength, fpRG) != rg->contigs[i].contigNameLength ||
				fread(&rg->contigs[i].sequenceLength, sizeof(int32_t), 1, fpRG)!=1 ||
				fread(&rg->contigs[i].numBytes, sizeof(uint32_t), 1, fpRG)!=1) {
			PrintError(FnName,
					NULL, 
					"Could not read RGContig information",
					Exit,
					ReadFileError);
		}
		/* Add null terminator */
		rg->contigs[i].contigName[rg->contigs[i].contigNameLength]='\0';
		/* Allocate memory for the sequence */
		rg->contigs[i].sequence = malloc(sizeof(uint8_t)*rg->contigs[i].numBytes);
		if(NULL==rg->contigs[i].sequence) {
			PrintError(FnName,
					"rg->contigs[i].sequence",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Read sequence */
		if(fread(rg->contigs[i].sequence, sizeof(uint8_t), rg->contigs[i].numBytes, fpRG)!=rg->contigs[i].numBytes) {
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

void RGIndexRead_0_1_13(FILE *fp, RGIndex_0_1_13 *index, int32_t binaryInput)
{
	char *FnName="RGIndexRead";
	int64_t i;
	uint32_t tempInt;

	/* Read in the header */
	RGIndexReadHeader_0_1_13(fp, index, binaryInput);

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
	/* Allocate memory for the ends */
	index->ends = malloc(sizeof(uint32_t)*index->hashLength);
	if(NULL == index->ends) {
		PrintError(FnName,
				"index->ends",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	if(binaryInput == TextInput) {
		/* Read the positions and contigs */
		for(i=0;i<index->length;i++) {
			if(fscanf(fp, "%u\t%u\n",
						&index->positions[i],
						&tempInt)==EOF) {
				PrintError(FnName,
						NULL,
						"Could not read in contig and position",
						Exit,
						EndOfFile);
			}
			if(index->contigType == Contig_8) {
				index->contigs_8[i] = (uint8_t)tempInt;
			}
			else {
				index->contigs_32[i] = tempInt;
			}
		}

		/* Read the positions and contigs */
		for(i=0;i<index->hashLength;i++) {
			if(fscanf(fp, "%u\t%u\n",
						&index->starts[i],
						&index->ends[i])==EOF) {
				PrintError(FnName,
						NULL,
						"Could not read in starts and ends",
						Exit,
						EndOfFile);
			}
		}
	}
	else {
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

		/* Read in ends */
		if(fread(index->ends, sizeof(uint32_t), index->hashLength, fp)!=index->hashLength) {
			PrintError(FnName,
					NULL,
					"Could not read in ends",
					Exit,
					ReadFileError);
		}           
	}   
}

void RGIndexReadHeader_0_1_13(FILE *fp, RGIndex_0_1_13 *index, int32_t binaryInput)
{
	char *FnName = "RGIndexReadHeader";
	int i;
	char tempChar;
	long long int tempLongLongInt[2];
	/* Read in header */
	if(binaryInput == 0) {
		if(fscanf(fp, "%d %lld %d %d %d %d %d %d %d %d %d %u %lld",
					&index->id,
					&tempLongLongInt[0],
					&index->contigType,
					&index->startContig,
					&index->startPos,
					&index->endContig,
					&index->endPos,
					&index->width,
					&index->keysize,
					&index->repeatMasker,
					&index->space,
					&index->hashWidth,
					&tempLongLongInt[1])==EOF) {
			PrintError(FnName,
					NULL,
					"Could not read header",
					Exit,
					EndOfFile);
		}
		index->length = tempLongLongInt[0];
		index->hashLength = tempLongLongInt[1];
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
		for(i=0;i<index->width;i++) {
			if(fscanf(fp, "%c", &tempChar)==EOF) {
				PrintError(FnName,
						NULL,
						"Could not read header",
						Exit,
						EndOfFile);
			}
			switch(tempChar) {
				case '0':
					index->mask[i] = 0;
					break;
				case '1':
					index->mask[i] = 1;
					break;
				default:
					PrintError(FnName,
							NULL,
							"Could not read mask",
							Exit,
							OutOfRange);
			}
		}
	}
	else {
		if(fread(&index->id, sizeof(int32_t), 1, fp) != 1 ||
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

