#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/AlignedRead.h"
#include "../blib/AlignedReadConvert.h"
#include "Definitions.h"
#include "Filter.h"
#include "InputOutputToFiles.h"

/* TODO */
void ReadInputFilterAndOutput(RGBinary *rg,
		char *inputFileName,
		int binaryInput,
		int startContig,
		int startPos,
		int endContig,
		int endPos,
		int algorithmReads,
		int minScores,
		int maxMismatches,
		int maxColorErrors,
		int minDistancePaired,
		int maxDistancePaired,
		int useDistancePaired,
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
	AlignedRead a;
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
		case GFF:
			strcpy(fileExtension, BFAST_GFF_FILE_EXTENSION);
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
		if(!(fpContigAb=fopen(contigAbFileName, "wb"))) {
			PrintError(FnName,
					contigAbFileName,
					"Could not open contigAbFileName for writing",
					Exit,
					OpenFileError);
		}
		AlignedReadConvertPrintHeader(fpContigAb, outputFormat);
	}
	if(inversionsPaired == 1) {
		if(!(fpInversions=fopen(inversionsFileName, "wb"))) {
			PrintError(FnName,
					inversionsFileName,
					"Could not open inversionsFileName for writing",
					Exit,
					OpenFileError);
		}
		AlignedReadConvertPrintHeader(fpInversions, outputFormat);
	}
	if(unpaired == 1) {
		if(!(fpUnpaired=fopen(unpairedFileName, "wb"))) {
			PrintError(FnName,
					unpairedFileName,
					"Could not open unpairedFileName for writing",
					Exit,
					OpenFileError);
		}
		AlignedReadConvertPrintHeader(fpUnpaired, outputFormat);
	}
	if(!(fpNotReported=fopen(notReportedFileName, "wb"))) {
		PrintError(FnName,
				notReportedFileName,
				"Could not open notReportedFileName for writing",
				Exit,
				OpenFileError);
	}
	AlignedReadConvertPrintHeader(fpNotReported, outputFormat);
	if(!(fpOut=fopen(outputFileName, "wb"))) {
		PrintError(FnName,
				outputFileName,
				"Could not open outputFileName for writing",
				Exit,
				OpenFileError);
	}
	AlignedReadConvertPrintHeader(fpOut, outputFormat);

	/* Initialize */
	AlignedReadInitialize(&a);

	/* Go through each read */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Processing reads, currently on:\n0");
	}
	counter = numReported = numNotReported = numContigAb = numUnpaired = numInversions = 0;
	while(EOF != AlignedReadRead(&a, fp, binaryInput)) {
		if(VERBOSE >= 0 && counter%ALIGNENTRIES_READ_ROTATE_NUM==0) {
			fprintf(stderr, "\r%lld",
					(long long int)counter);
		}
		/* Filter */
		foundType=FilterAlignedRead(&a,
				algorithmReads,
				minScores,
				startContig,        
				startPos,
				endContig,        
				endPos,
				maxMismatches,
				maxColorErrors,
				minDistancePaired,          
				maxDistancePaired,
				useDistancePaired,
				contigAbPaired,
				inversionsPaired,
				unpaired);

		/* Print the apporiate files based on the return type */
		switch(foundType) {
			case NoneFound:
				/* Print to Not Reported file */
				AlignedReadConvertPrintOutputFormat(&a, rg, fpNotReported, outputFormat, binaryInput);
				numNotReported++;
				break;
			case Found:
				/* Print to Output file */
				AlignedReadConvertPrintOutputFormat(&a, rg, fpOut, outputFormat, binaryInput);
				numReported++;
				break;
			case ContigAb:
				if(contigAbPaired == 1) {
					/* Print to Contig Abnormalities file */
					AlignedReadConvertPrintOutputFormat(&a, rg, fpContigAb, outputFormat, binaryInput);
					numContigAb++;
				}
				break;
			case Unpaired:
				if(unpaired == 1) {
					/* Print to Unpaired file */
					AlignedReadConvertPrintOutputFormat(&a, rg, fpUnpaired, outputFormat, binaryInput);
					numUnpaired++;
				}
				break;
			case Inversion:
				if(inversionsPaired == 1) {
					/* Print to Inversions file */
					AlignedReadConvertPrintOutputFormat(&a, rg, fpInversions, outputFormat, binaryInput);
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
		AlignedReadFree(&a);
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
		fclose(fpInversions);
	}
	if(contigAbPaired == 1) {
		fclose(fpContigAb);
	}
	if(unpaired == 1) {
		fclose(fpUnpaired);
	}

	/* Close the input file */
	fclose(fp);
}
