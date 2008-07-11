#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/AlignEntries.h"
#include "Definitions.h"
#include "Filter.h"
#include "InputOutputToFiles.h"

/* TODO */
void ReadInputFilterAndOutput(char *inputFileName,
		int binaryInput,
		int pairedEnd,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int algorithmReads,
		int minScoreReads,
		int algorithmReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int meanDistancePaired,
		int chrAbPaired,
		int inversionsPaired,
		char *outputID,
		char *outputDir,
		char *tmpDir,
		int outputFormat)
{
	char *FnName="ReadInputFilterAndOutput";
	FILE *fp=NULL;
	int64_t counter, foundType, numChrAb, numInversions, numNotReported, numReported;
	AlignEntries a;
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fpOut=NULL;
	char chrAbFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fpChrAb=NULL;
	char inversionsFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fpInversions=NULL;
	char notReportedFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fpNotReported=NULL;

	/* Open the input file */
	if(!(fp=fopen(inputFileName, "rb"))) {
		PrintError(FnName,
				inputFileName,
				"Could not open inputFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Create output file names */
	sprintf(chrAbFileName, "%s%s.chrab.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			BFAST_MAF_FILE_EXTENSION);
	sprintf(inversionsFileName, "%s%s.inversion.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			BFAST_MAF_FILE_EXTENSION);
	sprintf(notReportedFileName, "%s%s.not.reported.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			BFAST_MAF_FILE_EXTENSION);
	sprintf(outputFileName, "%s%s.report.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			BFAST_MAF_FILE_EXTENSION);
	/* Open output files, if necessary */
	if(chrAbPaired == 1) {
		assert(1==pairedEnd);
		if(!(fpChrAb=fopen(chrAbFileName, "wb"))) {
			PrintError(FnName,
					chrAbFileName,
					"Could not open chrAbFileName for writing",
					Exit,
					OpenFileError);
		}
	}
	if(inversionsPaired == 1) {
		assert(1==pairedEnd);
		if(!(fpInversions=fopen(inversionsFileName, "wb"))) {
			PrintError(FnName,
					inversionsFileName,
					"Could not open inversionsFileName for writing",
					Exit,
					OpenFileError);
		}
	}
	if(!(fpNotReported=fopen(notReportedFileName, "wb"))) {
		PrintError(FnName,
				notReportedFileName,
				"Could not open notReportedFileName for writing",
				Exit,
				OpenFileError);
	}
	if(!(fpOut=fopen(outputFileName, "wb"))) {
		PrintError(FnName,
				outputFileName,
				"Could not open outputFileName for writing",
				Exit,
				OpenFileError);
	}

	/* Initialize */
	AlignEntriesInitialize(&a);

	/* Go through each read */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Processing reads, currently on:\n0");
	}
	counter = numReported = numNotReported = numChrAb = numInversions = 0;
	while(EOF != AlignEntriesRead(&a, fp)) {
		if(VERBOSE >= 0 && counter%ALIGNENTRIES_READ_ROTATE_NUM==0) {
			fprintf(stderr, "\r%lld",
					(long long int)counter);
		}
		/* Filter */
		foundType=FilterAlignEntries(&a,
				algorithmReads,
				minScoreReads,
				startChr,        
				startPos,
				endChr,        
				endPos,
				pairedEnd,        
				algorithmReadsPaired,
				minScoreReadsPaired,
				minDistancePaired,          
				maxDistancePaired,
				meanDistancePaired);

		/* Print the apporiate files based on the return type */
		switch(foundType) {
			case NoneFound:
				/* Print to Not Reported file */
				PrintAlignEntriesToOutputFormat(&a, fpNotReported, outputFormat);
				numNotReported++;
				break;
			case Found:
				/* Print to Output file */
				PrintAlignEntriesToOutputFormat(&a, fpOut, outputFormat);
				numReported++;
				break;
			case ChrAb:
				assert(pairedEnd == 1);
				/* Print to Chromosomal Abnormalities file */
				PrintAlignEntriesToOutputFormat(&a, fpChrAb, outputFormat);
				AlignEntriesPrint(&a, fpChrAb);
				numChrAb++;
				break;
			case Inversion:
				assert(pairedEnd == 1);
				/* Print to Inversions file */
				PrintAlignEntriesToOutputFormat(&a, fpInversions, outputFormat);
				numInversions++;
				break;
			default:
				PrintError(FnName,
						"foundType",
						"Could not understand foundType",
						Exit,
						OutOfRange);
				break;
		}

		/* Free memory */
		AlignEntriesFree(&a);
		/* Increment counter */
		counter++;
	}
	if(VERBOSE>=0) {
		fprintf(stderr, "\r%lld\n",
				(long long int)counter);
		fprintf(stderr, "Out of %lld reads:\n", (long long int)counter);
		fprintf(stderr, "Found alignments for %lld reads.\n", (long long int)numReported);
		fprintf(stderr, "Could not unambiguously align %lld reads.\n", (long long int)numNotReported);
		if(1==chrAbPaired) {
			fprintf(stderr, "Found %lld paired end reads with chromosomal abnormalities.\n", (long long int)numChrAb);
		}
		if(1==inversionsPaired) {
			fprintf(stderr, "Found %lld inverted paired end reads.\n", (long long int)numInversions); 
		}
	}

	/* Close output files, if necessary */
	fclose(fpOut);
	fclose(fpNotReported);
	if(inversionsPaired == 1) {
		assert(1==pairedEnd);
		fclose(fpInversions);
	}
	if(chrAbPaired == 1) {
		assert(1==pairedEnd);
		fclose(fpChrAb);
	}

	/* Close the input file */
	fclose(fp);
}

/* TODO */
void ReadInputAndOutput(char *inputFileName,
		int binaryInput,
		int pairedEnd,
		char *outputID,
		char *outputDir,
		char *tmpDir,
		int outputFormat)
{
	char *FnName="ReadInputAndOutput";
	FILE *fp=NULL;
	int64_t counter;
	AlignEntries a;
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fpOut=NULL;

	/* Open the input file */
	if(!(fp=fopen(inputFileName, "rb"))) {
		PrintError(FnName,
				inputFileName,
				"Could not open inputFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Create output file names */
	sprintf(outputFileName, "%s%s.report.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			BFAST_MAF_FILE_EXTENSION);
	/* Open output file */
	if(!(fpOut=fopen(outputFileName, "wb"))) {
		PrintError(FnName,
				outputFileName,
				"Could not open outputFileName for writing",
				Exit,
				OpenFileError);
	}

	/* Initialize */
	AlignEntriesInitialize(&a);

	/* Go through each read */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Processing reads, currently on:\n0");
	}
	counter = 0;
	while(EOF != AlignEntriesRead(&a, fp)) {
		if(VERBOSE >= 0 && counter%ALIGNENTRIES_READ_ROTATE_NUM==0) {
			fprintf(stderr, "\r%lld",
					(long long int)counter);
		}

		assert(a.pairedEnd == pairedEnd);

		/* Print to Output file */
		PrintAlignEntriesToOutputFormat(&a, fpOut, outputFormat);

		/* Free memory */
		AlignEntriesFree(&a);
		/* Increment counter */
		counter++;
	}
	if(VERBOSE>=0) {
		fprintf(stderr, "\r%lld\n",
				(long long int)counter);
	}

	/* Close output files, if necessary */
	fclose(fpOut);

	/* Close the input file */
	fclose(fp);
}

/* TODO */
void PrintAlignEntriesToOutputFormat(AlignEntries *a, 
		FILE *fp,
		int outputFormat)
{
	char *FnName = "PrintAlignEntriesToOutputFormat";
	switch(outputFormat) {
		case MAF:
			break;
		default:
			PrintError(FnName,
					"outputFormat",
					"Could not understand outputFormat",
					Exit,
					OutOfRange);
			break;
	}
}
