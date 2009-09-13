#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <zlib.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/AlignedRead.h"
#include "../blib/AlignedReadConvert.h"
#include "../blib/ScoringMatrix.h"
#include "Definitions.h"
#include "Filter.h"
#include "InputOutputToFiles.h"

/* TODO */
void ReadInputFilterAndOutput(RGBinary *rg,
		char *inputFileName,
		int algorithm,
		int queueLength,
		char *outputID,
		char *outputDir,
		int outputFormat)
{
	char *FnName="ReadInputFilterAndOutput";
	gzFile fp=NULL;
	int32_t i, j;
	int64_t counter, foundType, numNotReported, numReported;
	AlignedRead *aBuffer;
	int32_t aBufferLength=0;
	int32_t numRead, aBufferIndex;
	char reportedFileName[MAX_FILENAME_LENGTH]="\0";
	gzFile fpReportedGZ=NULL;
	FILE *fpReported=NULL;
	char notReportedFileName[MAX_FILENAME_LENGTH]="\0";
	gzFile fpNotReported=NULL;
	char fileExtension[256]="\0";

	/* Open the input file */
	if(!(fp=gzopen(inputFileName, "rb"))) {
		PrintError(FnName,
				inputFileName,
				"Could not open inputFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Get file extension for the output files */
	switch(outputFormat) {
		case BAF:
			strcpy(fileExtension, BFAST_ALIGNED_FILE_EXTENSION);
			break;
		case MAF:
			strcpy(fileExtension, BFAST_MAF_FILE_EXTENSION);
			break;
		case GFF:
			strcpy(fileExtension, BFAST_GFF_FILE_EXTENSION);
			break;
		case SAM:
			strcpy(fileExtension, BFAST_SAM_FILE_EXTENSION);
			break;
		default:
			PrintError(FnName,
					"outputFormat",
					"Could not understand output format",
					Exit,
					OutOfRange);
	}
	/* Create output file names */
	sprintf(notReportedFileName, "%s%s.not.reported.file.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			BFAST_ALIGNED_FILE_EXTENSION);
	sprintf(reportedFileName, "%s%s.reported.file.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			fileExtension);

	/* Open output files, if necessary */
	if(!(fpNotReported=gzopen(notReportedFileName, "wb"))) {
		PrintError(FnName,
				notReportedFileName,
				"Could not open notReportedFileName for writing",
				Exit,
				OpenFileError);
	}
	if(BAF == outputFormat) {
		if(!(fpReportedGZ=gzopen(reportedFileName, "wb"))) {
			PrintError(FnName,
					reportedFileName,
					"Could not open reportedFileName for writing",
					Exit,
					OpenFileError);
		}
	}
	else {
		if(!(fpReported=fopen(reportedFileName, "wb"))) {
			PrintError(FnName,
					reportedFileName,
					"Could not open reportedFileName for writing",
					Exit,
					OpenFileError);
		}
	}

	AlignedReadConvertPrintHeader(fpReported, rg, outputFormat);

	aBufferLength=queueLength;
	aBuffer=malloc(sizeof(AlignedRead)*aBufferLength);
	if(NULL == aBuffer) {
		PrintError(FnName,
				"aBuffer",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Go through each read */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Processing reads, currently on:\n0");
	}
	counter = numReported = numNotReported = numRead = 0;

	while(0 != (numRead = GetAlignedReads(fp, aBuffer, aBufferLength))) {

		for(aBufferIndex=0;aBufferIndex<numRead;aBufferIndex++) {

			if(VERBOSE >= 0 && counter%ALIGNENTRIES_READ_ROTATE_NUM==0) {
				fprintf(stderr, "\r%lld",
						(long long int)counter);
			}

			/* Filter */
			foundType=FilterAlignedRead(&aBuffer[aBufferIndex],
					algorithm);

			if(NoneFound == foundType) {
				/* Print to Not Reported file */
				AlignedReadPrint(&aBuffer[aBufferIndex], fpNotReported);
				numNotReported++;

				/* Free the alignments for output */
				for(i=0;i<aBuffer[aBufferIndex].numEnds;i++) {
					for(j=0;j<aBuffer[aBufferIndex].ends[i].numEntries;j++) {
						AlignedEntryFree(&aBuffer[aBufferIndex].ends[i].entries[j]);
					}
					aBuffer[aBufferIndex].ends[i].numEntries=0;
				}
			}

			/* Print to Output file */
			AlignedReadConvertPrintOutputFormat(&aBuffer[aBufferIndex], rg, fpReported, fpReportedGZ, outputID, outputFormat, BinaryOutput);
			numReported++;

			/* Free memory */
			AlignedReadFree(&aBuffer[aBufferIndex]);
			/* Increment counter */
			counter++;
		}
	}
	if(VERBOSE>=0) {
		fprintf(stderr, "\r%lld\n",
				(long long int)counter);
	}

	/* Close output files, if necessary */
	if(BAF == outputFormat) {
		gzclose(fpReportedGZ);
	}
	else {
		fclose(fpReported);
	}
	gzclose(fpNotReported);
	/* Close the input file */
	gzclose(fp);

	free(aBuffer);
}

int32_t GetAlignedReads(gzFile fp, AlignedRead *aBuffer, int32_t maxToRead)
{
	int32_t numRead=0;
	while(numRead < maxToRead) {
		AlignedReadInitialize(&aBuffer[numRead]);
		if(EOF == AlignedReadRead(&aBuffer[numRead], fp)) {
			return numRead;
		}
		numRead++;
	}
	return numRead;
}
