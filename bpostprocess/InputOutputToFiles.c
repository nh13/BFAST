#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
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
	char fileExtension[256]="\0";

	/* Open the input file */
	if(!(fp=fopen(inputFileName, "rb"))) {
		PrintError(FnName,
				inputFileName,
				"Could not open inputFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Get file extension for the output files */
	switch(outputFormat) {
		case BAF:
			strcpy(fileExtension,  BFAST_ALIGN_FILE_EXTENSION);
			break;
		case MAF:
			strcpy(fileExtension, BFAST_MAF_FILE_EXTENSION);
			break;
		default:
				PrintError(FnName,
						"outputFormat",
						"Could not understand output format",
						Exit,
						OutOfRange);
	}
	/* Create output file names */
	sprintf(chrAbFileName, "%s%s.chrab.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			fileExtension);
	sprintf(inversionsFileName, "%s%s.inversion.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			fileExtension);
	sprintf(notReportedFileName, "%s%s.not.reported.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			fileExtension);
	sprintf(outputFileName, "%s%s.report.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			fileExtension);

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
		PrintHeader(fpChrAb, outputFormat);
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
		PrintHeader(fpInversions, outputFormat);
	}
	if(!(fpNotReported=fopen(notReportedFileName, "wb"))) {
		PrintError(FnName,
				notReportedFileName,
				"Could not open notReportedFileName for writing",
				Exit,
				OpenFileError);
	}
	PrintHeader(fpNotReported, outputFormat);
	if(!(fpOut=fopen(outputFileName, "wb"))) {
		PrintError(FnName,
				outputFileName,
				"Could not open outputFileName for writing",
				Exit,
				OpenFileError);
	}
	PrintHeader(fpOut, outputFormat);

	/* Initialize */
	AlignEntriesInitialize(&a);

	/* Go through each read */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Processing reads, currently on:\n0");
	}
	counter = numReported = numNotReported = numChrAb = numInversions = 0;
	while(EOF != AlignEntriesRead(&a, fp, pairedEnd, SpaceDoesNotMatter)) {
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
void PrintHeader(FILE *fp,
		int outputFormat) 
{
	char *FnName = "PrintHeader";
	switch(outputFormat) {
		case BAF:
			/* Do nothing */
			break;
		case MAF:
			if(0>fprintf(fp, "##maf version=1 scoring=%s\n",
						PROGRAM_NAME)) {
				PrintError(FnName,
						"header",
						"Could not write to file",
						Exit,
						WriteFileError);
			}
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

/* TODO */
void PrintAlignEntriesToOutputFormat(AlignEntries *a, 
		FILE *fp,
		int outputFormat)
{
	char *FnName = "PrintAlignEntriesToOutputFormat";
	switch(outputFormat) {
		case BAF:
			AlignEntriesPrint(a, fp);
			break;
		case MAF:
			PrintAlignEntriesToMAF(a, fp);
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

/* TODO */
void PrintAlignEntriesToMAF(AlignEntries *a,
		FILE *fp)
{
	/*
	   char *FnName="PrintAlignEntriesToMAF";
	   */
	int i;

	/* Get Data */
	if(0==a->pairedEnd) {
		for(i=0;i<a->numEntriesOne;i++) {
			PrintAlignEntryToMAF(&a->entriesOne[i], a->readName, a->pairedEnd, a->colorSpace, 1, fp); 
		}
	}
	else {
		for(i=0;i<a->numEntriesOne;i++) {
			PrintAlignEntryToMAF(&a->entriesOne[i], a->readName, a->pairedEnd, a->colorSpace, 1, fp); 
		}
		for(i=0;i<a->numEntriesTwo;i++) {
			PrintAlignEntryToMAF(&a->entriesTwo[i], a->readName, a->pairedEnd, a->colorSpace, 2, fp); 
		}
	}

}

/* TODO */
void PrintAlignEntryToMAF(AlignEntry *a,
		char *readName,
		int pairedEnd,
		int colorSpace,
		int readNum,
		FILE *fp)
{
	char *FnName="PrintAlignEntryToMAF";
	int i;
	int originalReferenceLength=0;
	int originalReadLength=0; 

	/* Recover original lengths */
	for(i=0;i<a->length;i++) {
		if(a->reference[i] != GAP) {
			originalReferenceLength++;
		}
		if(a->read[i] != GAP) {
			originalReadLength++;
		}
	}

	/* Print the score */
	if(colorSpace == 1) {
		if(0>fprintf(fp, "a score=%lf paired-end=%d read=%d color-errors=%s\n",
					a->score,
					pairedEnd,
					readNum,
					a->colorError)) {
			PrintError(FnName,
					NULL,
					"Could not write to file",
					Exit,
					WriteFileError);
		}
	}
	else {
		if(0>fprintf(fp, "a score=%lf paired-end=%d read=%d\n",
					a->score,
					pairedEnd,
					readNum)) {
			PrintError(FnName,
					NULL,
					"Could not write to file",
					Exit,
					WriteFileError);
		}
	}

	/* Print the reference */
	if(0>fprintf(fp, "s chr%d %u %d %c %d %s\n",
				a->chromosome,
				a->position-1, /* zero based */
				a->length,
				a->strand,
				originalReferenceLength,
				a->reference)) {
		PrintError(FnName,
				NULL,
				"Could not write to file",
				Exit,
				WriteFileError);
	}
	/* Print the read */
	if(0>fprintf(fp, "s %s %u %d %c %d %s\n\n", /* Include a blank line */
				readName,
				0,
				a->length,
				a->strand,
				originalReadLength,
				a->read)) {
		PrintError(FnName,
				NULL,
				"Could not write to file",
				Exit,
				WriteFileError);
	}
}
