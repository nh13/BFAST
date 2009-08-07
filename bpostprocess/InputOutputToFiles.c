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
		char *outputID,
		char *outputDir,
		int outputFormat)
{
	char *FnName="ReadInputFilterAndOutput";
	gzFile fp=NULL;
	int32_t i, j;
	int64_t counter, foundType, numNotReported, numReported;
	AlignedRead a;
	char reportedFileName[MAX_FILENAME_LENGTH]="\0";
	gzFile fpReportedGZ=NULL;
	FILE *fpReported=NULL;
	char notReportedFileName[MAX_FILENAME_LENGTH]="\0";
	gzFile fpNotReportedGZ=NULL;
	FILE *fpNotReported=NULL;
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
			strcpy(fileExtension,  BFAST_ALIGNED_FILE_EXTENSION);
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
			fileExtension);
	sprintf(reportedFileName, "%s%s.reported.file.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			fileExtension);

	/* Open output files, if necessary */
	if(BAF == outputFormat) { 
		if(!(fpNotReportedGZ=gzopen(notReportedFileName, "wb"))) {
			PrintError(FnName,
					notReportedFileName,
					"Could not open notReportedFileName for writing",
					Exit,
					OpenFileError);
		}
	}
	else {
		if(!(fpNotReported=fopen(notReportedFileName, "wb"))) {
			PrintError(FnName,
					notReportedFileName,
					"Could not open notReportedFileName for writing",
					Exit,
					OpenFileError);
		}
	}
	AlignedReadConvertPrintHeader(fpNotReported, rg, outputFormat);
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

	/* Initialize */
	AlignedReadInitialize(&a);

	/* Go through each read */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Processing reads, currently on:\n0");
	}
	counter = numReported = numNotReported = 0;
	while(EOF != AlignedReadRead(&a, fp)) {
		if(VERBOSE >= 0 && counter%ALIGNENTRIES_READ_ROTATE_NUM==0) {
			fprintf(stderr, "\r%lld",
					(long long int)counter);
		}

		/* Filter */
		foundType=FilterAlignedRead(&a,
				algorithm);

		if(NoneFound == foundType) {
			/* Print to Not Reported file */
			AlignedReadConvertPrintOutputFormat(&a, rg, fpNotReported, fpNotReportedGZ, outputID, outputFormat, BinaryOutput);
			numNotReported++;

			/* Free the alignments for output */
			for(i=0;i<a.numEnds;i++) {
				for(j=0;j<a.ends[i].numEntries;j++) {
					AlignedEntryFree(&a.ends[i].entries[j]);
				}
				a.ends[i].numEntries=0;
			}
		}

		/* Print to Output file */
		AlignedReadConvertPrintOutputFormat(&a, rg, fpReported, fpReportedGZ, outputID, outputFormat, BinaryOutput);
		numReported++;

		/* Free memory */
		AlignedReadFree(&a);
		/* Increment counter */
		counter++;
	}
	if(VERBOSE>=0) {
		fprintf(stderr, "\r%lld\n",
				(long long int)counter);
	}

	/* Close output files, if necessary */
	if(BAF == outputFormat) {
		gzclose(fpReportedGZ);
		gzclose(fpNotReportedGZ);
	}
	else {
		fclose(fpReported);
		fclose(fpNotReported);
	}
	/* Close the input file */
	gzclose(fp);
}
