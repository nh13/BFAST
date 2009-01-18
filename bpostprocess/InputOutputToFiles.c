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
void ReadInputFilterAndOutput(RGBinary *rg,
		char *inputFileName,
		int binaryInput,
		int pairedEnd,
		int startContig,
		int startPos,
		int endContig,
		int endPos,
		int algorithmReads,
		int minScoreReads,
		int maxMismatches,
		int maxColorErrors,
		int algorithmReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int maxMismatchesPaired,
		int maxColorErrorsPaired,
		int contigAbPaired,
		int inversionsPaired,
		int unpaired,
		char *outputID,
		char *outputDir,
		int outputFormat)
{
	char *FnName="ReadInputFilterAndOutput";
	FILE *fp=NULL;
	int64_t counter, foundType, numContigAb, numUnpaired, numInversions, numNotReported, numReported;
	AlignEntries a;
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fpOut=NULL;
	char contigAbFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fpContigAb=NULL;
	char inversionsFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fpInversions=NULL;
	char unpairedFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fpUnpaired=NULL;
	char notReportedFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fpNotReported=NULL;
	char fileExtension[256]="\0";

	assert(binaryInput == BinaryInput ||
			binaryInput == TextInput);

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
			strcpy(fileExtension,  BFAST_ALIGNED_FILE_EXTENSION);
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
	sprintf(contigAbFileName, "%s%s.contigab.file.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			fileExtension);
	sprintf(inversionsFileName, "%s%s.inversion.file.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			fileExtension);
	sprintf(unpairedFileName, "%s%s.unpaired.file.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			fileExtension);
	sprintf(notReportedFileName, "%s%s.not.reported.file.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			fileExtension);
	sprintf(outputFileName, "%s%s.reported.file.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			fileExtension);

	/* Open output files, if necessary */
	if(contigAbPaired == 1) {
		assert(1==pairedEnd);
		if(!(fpContigAb=fopen(contigAbFileName, "wb"))) {
			PrintError(FnName,
					contigAbFileName,
					"Could not open contigAbFileName for writing",
					Exit,
					OpenFileError);
		}
		PrintHeader(fpContigAb, outputFormat);
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
	if(unpaired == 1) {
		assert(1==pairedEnd);
		if(!(fpUnpaired=fopen(unpairedFileName, "wb"))) {
			PrintError(FnName,
					unpairedFileName,
					"Could not open unpairedFileName for writing",
					Exit,
					OpenFileError);
		}
		PrintHeader(fpUnpaired, outputFormat);
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
	counter = numReported = numNotReported = numContigAb = numUnpaired = numInversions = 0;
	while(EOF != AlignEntriesRead(&a, fp, pairedEnd, SpaceDoesNotMatter, binaryInput)) {
		if(VERBOSE >= 0 && counter%ALIGNENTRIES_READ_ROTATE_NUM==0) {
			fprintf(stderr, "\r%lld",
					(long long int)counter);
		}
		/* Filter */
		foundType=FilterAlignEntries(&a,
				algorithmReads,
				minScoreReads,
				startContig,        
				startPos,
				endContig,        
				endPos,
				maxMismatches,
				maxColorErrors,
				pairedEnd,        
				algorithmReadsPaired,
				unpaired,
				minScoreReadsPaired,
				minDistancePaired,          
				maxDistancePaired,
				maxMismatchesPaired,
				maxColorErrorsPaired);

		/* Print the apporiate files based on the return type */
		switch(foundType) {
			case NoneFound:
				/* Print to Not Reported file */
				PrintAlignEntriesToOutputFormat(&a, rg, fpNotReported, outputFormat, binaryInput);
				numNotReported++;
				break;
			case Found:
				/* Print to Output file */
				PrintAlignEntriesToOutputFormat(&a, rg, fpOut, outputFormat, binaryInput);
				numReported++;
				break;
			case ContigAb:
				assert(pairedEnd == 1);
				if(contigAbPaired == 1) {
					/* Print to Contig Abnormalities file */
					PrintAlignEntriesToOutputFormat(&a, rg, fpContigAb, outputFormat, binaryInput);
					numContigAb++;
				}
				break;
			case Unpaired:
				assert(pairedEnd == 1);
				if(unpaired == 1) {
					/* Print to Unpaired file */
					PrintAlignEntriesToOutputFormat(&a, rg, fpUnpaired, outputFormat, binaryInput);
					numUnpaired++;
				}
				break;
			case Inversion:
				assert(pairedEnd == 1);
				if(inversionsPaired == 1) {
					/* Print to Inversions file */
					PrintAlignEntriesToOutputFormat(&a, rg, fpInversions, outputFormat, binaryInput);
					numInversions++;
				}
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
		if(1==contigAbPaired) {
			fprintf(stderr, "Found %lld paired end reads with contigomosomal abnormalities.\n", (long long int)numContigAb);
		}
		if(1==inversionsPaired) {
			fprintf(stderr, "Found %lld inverted paired end reads.\n", (long long int)numInversions); 
		}
		if(1==unpaired) {
			fprintf(stderr, "Found %lld unpaired paired end reads.\n", (long long int)numUnpaired);
		}
	}

	/* Close output files, if necessary */
	fclose(fpOut);
	fclose(fpNotReported);
	if(inversionsPaired == 1) {
		assert(1==pairedEnd);
		fclose(fpInversions);
	}
	if(contigAbPaired == 1) {
		assert(1==pairedEnd);
		fclose(fpContigAb);
	}
	if(unpaired == 1) {
		assert(1==pairedEnd);
		fclose(fpUnpaired);
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
		RGBinary *rg,
		FILE *fp,
		int outputFormat,
		int binaryInput)
{
	char *FnName = "PrintAlignEntriesToOutputFormat";
	switch(outputFormat) {
		case BAF:
			AlignEntriesPrint(a, fp, binaryInput);
			break;
		case MAF:
			PrintAlignEntriesToMAF(a, rg, fp);
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
		RGBinary *rg,
		FILE *fp)
{
	/*
	   char *FnName="PrintAlignEntriesToMAF";
	   */
	int i;

	/* Get Data */
	if(0==a->pairedEnd) {
		for(i=0;i<a->numEntriesOne;i++) {
			PrintAlignEntryToMAF(&a->entriesOne[i], rg, a->readName, a->pairedEnd, a->space, 1, fp); 
		}
	}
	else {
		for(i=0;i<a->numEntriesOne;i++) {
			PrintAlignEntryToMAF(&a->entriesOne[i], rg, a->readName, a->pairedEnd, a->space, 1, fp); 
		}
		for(i=0;i<a->numEntriesTwo;i++) {
			PrintAlignEntryToMAF(&a->entriesTwo[i], rg, a->readName, a->pairedEnd, a->space, 2, fp); 
		}
	}

}

/* TODO */
void PrintAlignEntryToMAF(AlignEntry *a,
		RGBinary *rg,
		char *readName,
		int pairedEnd,
		int space,
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
	assert(originalReferenceLength == a->referenceLength);

	/* Print the score */
	if(space == ColorSpace) {
		if(0>fprintf(fp, "a score=%lf paired-end=%d read=%d color-errors=%s contig-index=%d\n",
					a->score,
					pairedEnd,
					readNum,
					a->colorError,
					a->contig)) {
			PrintError(FnName,
					NULL,
					"Could not write to file",
					Exit,
					WriteFileError);
		}
	}
	else {
		assert(space == NTSpace);
		if(0>fprintf(fp, "a score=%lf paired-end=%d read=%d contig-index=%d\n",
					a->score,
					pairedEnd,
					readNum,
					a->contig)) {
			PrintError(FnName,
					NULL,
					"Could not write to file",
					Exit,
					WriteFileError);
		}
	}

	/* Make sure the contig reported is within bounds of the reference genome */
	assert(1 <= a->contig && a->contig <= rg->numContigs);

	/* Print the reference */
	if(0>fprintf(fp, "s %s %u %d %c %d %s\n",
				a->contigName,
				a->position-1, /* zero based */
				originalReferenceLength,
				a->strand,
				rg->contigs[a->contig-1].sequenceLength, /* original contig length */
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
				originalReadLength, /* We align the full read */
				a->strand,
				originalReadLength, /* original read length */
				a->read)) {
		PrintError(FnName,
				NULL,
				"Could not write to file",
				Exit,
				WriteFileError);
	}
}
