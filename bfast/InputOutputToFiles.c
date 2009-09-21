#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <zlib.h>
#include "BLibDefinitions.h"
#include "BError.h"
#include "AlignedRead.h"
#include "AlignedReadConvert.h"
#include "ScoringMatrix.h"
#include "FilterAlignments.h"
#include "InputOutputToFiles.h"

/* TODO */
void ReadInputFilterAndOutput(RGBinary *rg,
		char *inputFileName,
		int algorithm,
		int queueLength,
		int outputFormat,
		char *outputID,
		char *unmappedFileName)
{
	char *FnName="ReadInputFilterAndOutput";
	gzFile fp=NULL;
	int32_t i, j;
	int64_t counter, foundType, numUnmapped, numReported;
	AlignedRead *aBuffer;
	int32_t aBufferLength=0;
	int32_t numRead, aBufferIndex;
	gzFile fpReportedGZ=NULL;
	FILE *fpReported=NULL;
	gzFile fpUnmapped=NULL;

	/* Open the input file */
	if(NULL == inputFileName) {
		if(!(fp=gzdopen(fileno(stdin), "rb"))) {
			PrintError(FnName, "stdin", "Could not open stdin for reading", Exit, OpenFileError);
		}
	}
	else {
		if(!(fp=gzopen(inputFileName, "rb"))) {
			PrintError(FnName, inputFileName, "Could not open inputFileName for reading", Exit, OpenFileError);
		}
	}

	/* Open output files, if necessary */
	if(BAF == outputFormat) {
		if(!(fpReportedGZ=gzdopen(fileno(stdout), "wb"))) {
			PrintError(FnName, "stdout", "Could not open stdout for writing", Exit, OpenFileError);
		}
	}
	else {
		if(!(fpReported=fdopen(fileno(stdout), "wb"))) {
			PrintError(FnName, "stdout", "Could not open stdout for writing", Exit, OpenFileError);
		}
	}
	if(NULL != unmappedFileName) {
		if(!(fpUnmapped=gzopen(unmappedFileName, "wb"))) {
			PrintError(FnName, unmappedFileName, "Could not open unmappedFileName for writing", Exit, OpenFileError);
		}
	}

	AlignedReadConvertPrintHeader(fpReported, rg, outputFormat);

	aBufferLength=queueLength;
	aBuffer=malloc(sizeof(AlignedRead)*aBufferLength);
	if(NULL == aBuffer) {
		PrintError(FnName, "aBuffer", "Could not allocate memory", Exit, MallocMemory);
	}

	/* Go through each read */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Processing reads, currently on:\n0");
	}
	counter = numReported = numUnmapped = numRead = 0;

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
				if(NULL != unmappedFileName) {
					AlignedReadPrint(&aBuffer[aBufferIndex], fpUnmapped);
				}
				numUnmapped++;

				/* Free the alignments for output */
				for(i=0;i<aBuffer[aBufferIndex].numEnds;i++) {
					for(j=0;j<aBuffer[aBufferIndex].ends[i].numEntries;j++) {
						AlignedEntryFree(&aBuffer[aBufferIndex].ends[i].entries[j]);
					}
					aBuffer[aBufferIndex].ends[i].numEntries=0;
				}
			}
			else {
				numReported++;
			}

			/* Increment counter */
			counter++;
		}

		/* Print to Output file */
		for(aBufferIndex=0;aBufferIndex<numRead;aBufferIndex++) {
			AlignedReadConvertPrintOutputFormat(&aBuffer[aBufferIndex], rg, fpReported, fpReportedGZ, (NULL == outputID) ? "" : outputID, outputFormat, BinaryOutput);

			/* Free memory */
			AlignedReadFree(&aBuffer[aBufferIndex]);
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
	if(NULL != unmappedFileName) {
		gzclose(fpUnmapped);
	}
	/* Close the input file */
	gzclose(fp);
	free(aBuffer);

	if(VERBOSE>=0) {
		fprintf(stderr, "Found %lld mappings after filtering.\n%lld reads were left unmapped\n",
				(long long int)numReported,
				(long long int)numUnmapped);
	}
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
