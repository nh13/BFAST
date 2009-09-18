#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include <config.h>
#include "BLibDefinitions.h"
#include "AlignedRead.h"
#include "AlignedReadConvert.h"
#include "BError.h"
#include "BLib.h"

#define Name "bfast bafconvert"
#define BAFCONVERT_ROTATE_NUM 100000

/* Converts a balign file to the specified output format.
 * */

int BfastBAFConvert(int argc, char *argv[])
{
	FILE *fpIn=NULL, *fpOut=NULL;
	gzFile fpInGZ=NULL, fpOutGZ=NULL;
	long long int counter;
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	char brgFileName[MAX_FILENAME_LENGTH]="\0";
	char outputID[MAX_FILENAME_LENGTH]="\0";
	char *last;
	int outputType=BAF;
	int outputSubType=TextOutput;
	int inputType=BinaryInput;
	int c;
	AlignedRead a;
	RGBinary rg;

	// Get parameters
	while((c = getopt(argc, argv, "i:o:r:O:")) >= 0) {
		switch(c) {
			case 'i': strcpy(inputFileName, optarg); break;
			case 'O': outputType = atoi(optarg); break;
			case 'r': strcpy(brgFileName, optarg); break;
			case 'o': strcpy(outputID, optarg); break;
			default: fprintf(stderr, "Unrecognized option: -%c\n", c); return 1;
		}
	}

	if(argc == optind) {
		fprintf(stderr, "%s %s\n", "bfast", PACKAGE_VERSION);
		fprintf(stderr, "\nUsage:%s [options]\n", Name);
		fprintf(stderr, "\t-i\t\tbfast aligned file name\n");
		fprintf(stderr, "\t-O\t\toutput type:\n"
				"\t\t\t\t0-BAF text to BAF binary\n"
				"\t\t\t\t1-BAF binary to BAF text\n"
				"\t\t\t\t2-BAF binary to MAF\n"
				"\t\t\t\t3-BAF binary to GFF (v2)\n"
				"\t\t\t\t4-BAF binary to SAM (v.%s)\n",
				BFAST_SAM_VERSION
			   );
		fprintf(stderr, "\t-r\t\tbfast reference genome file (not required for BAF output)>\n");
		fprintf(stderr, "\t-o\t\toutput ID (required only for SAM output)\n");
		fprintf(stderr, "\nsend bugs to %s\n",
				PACKAGE_BUGREPORT);
		return 1;
	}

	/* Create output file name */
	last = StrStrGetLast(inputFileName,
			BFAST_ALIGNED_FILE_EXTENSION);
	if(NULL == last) {
		PrintError(Name,
				inputFileName,
				"Could not recognize file extension",
				Exit,
				OutOfRange);
	}

	strncpy(outputFileName, inputFileName, (last - inputFileName));
	switch(outputType) {
		case 0:
			outputType=BAF;
			inputType=TextInput;
			outputSubType=BinaryOutput;
			strcat(outputFileName, "binary.");
			strcat(outputFileName, BFAST_ALIGNED_FILE_EXTENSION);
			break;
		case 1:
			outputType=BAF;
			inputType=BinaryInput;
			outputSubType=TextOutput;
			strcat(outputFileName, "text.");
			strcat(outputFileName, BFAST_ALIGNED_FILE_EXTENSION);
			break;
		case 2:
			outputType=MAF;
			inputType=BinaryInput;
			outputSubType=TextOutput;
			strcat(outputFileName, BFAST_MAF_FILE_EXTENSION);
			RGBinaryReadBinary(&rg,
					brgFileName);
			break;
		case 3:
			outputType=GFF;
			inputType=BinaryInput;
			outputSubType=TextOutput;
			strcat(outputFileName, BFAST_GFF_FILE_EXTENSION);
			RGBinaryReadBinary(&rg,
					brgFileName);
			break;
		case 4:
			outputType=SAM;
			inputType=BinaryInput;
			outputSubType=TextOutput;
			strcat(outputFileName, BFAST_SAM_FILE_EXTENSION);
			RGBinaryReadBinary(&rg,
					brgFileName);
			break;
		default:
			PrintError(Name,
					NULL,
					"Could not understand output type",
					Exit,
					OutOfRange);
	}

	/* Open the input file */
	if(BinaryInput == inputType) {
		if(!(fpInGZ=gzopen(inputFileName, "rb"))) {
			PrintError(Name,
					inputFileName,
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}
	}
	else {
		if(!(fpIn=fopen(inputFileName, "rb"))) {
			PrintError(Name,
					inputFileName,
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}
	}
	/* Open the output file */
	if(BinaryOutput == outputSubType) {
		if(!(fpOutGZ=gzopen(outputFileName, "wb"))) {
			PrintError(Name,
					outputFileName,
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}
	}
	else {
		if(!(fpOut=fopen(outputFileName, "wb"))) {
			PrintError(Name,
					outputFileName,
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}
	}

	/* Print Header */
	AlignedReadConvertPrintHeader(fpOut, &rg, outputType);
	/* Initialize */
	AlignedReadInitialize(&a);
	counter = 0;
	fprintf(stderr, "Currently on:\n0");
	/* Read in each match */
	while((TextInput == inputType && EOF != AlignedReadReadText(&a, fpIn)) ||
			(BinaryInput == inputType && EOF != AlignedReadRead(&a, fpInGZ))) {
		if(counter%BAFCONVERT_ROTATE_NUM==0) {
			fprintf(stderr, "\r%lld",
					counter);
		}
		counter++;
		/* Print each match */
		AlignedReadConvertPrintOutputFormat(&a,
				&rg,
				fpOut,
				fpOutGZ,
				outputID,
				outputType,
				outputSubType);
		AlignedReadFree(&a);
	}
	fprintf(stderr, "\r%lld\n",
			counter);
	/* Close the input file */
	if(TextInput == inputType) {
		fclose(fpIn);
	}
	else {
		gzclose(fpInGZ);
	}
	/* Close the output file */
	if(TextOutput == outputSubType) {
		fclose(fpOut);
	}
	else {
		gzclose(fpOutGZ);
	}
	if(MAF == outputType ||
			SAM == outputType) {
		RGBinaryDelete(&rg);
	}

	fprintf(stderr, "Terminating successfully!\n");
	return 0;
}
