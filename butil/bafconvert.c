#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/AlignedRead.h"
#include "../blib/AlignedReadConvert.h"
#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "bafconvert.h"

#define Name "bafconvert"
#define BAFCONVERT_ROTATE_NUM 100000

/* Converts a balign file to the specified output format.
 * */

int main(int argc, char *argv[])
{
	FILE *fpIn, *fpOut;
	long long int counter;
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	char rgFileName[MAX_FILENAME_LENGTH]="\0";
	char outputID[MAX_FILENAME_LENGTH]="\0";
	char *last;
	int outputType;
	int outputSubType=TextOutput;
	int inputType=BinaryInput;
	AlignedRead a;
	RGBinary rg;

	if(3 <= argc &&
			argc <= 5) {
		strcpy(inputFileName, argv[1]);
		outputType = atoi(argv[2]);

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
				if(4 != argc) {
					PrintError(Name,
							NULL,
							"A bfast reference genome file must be given",
							Exit,
							OutOfRange);
				}
				if(argc < 4) {
					PrintError(Name,
							"bfast reference genome file",
							"Command line argument",
							Exit,
							OutOfRange);
				}
				strcpy(rgFileName, argv[3]);
				RGBinaryReadBinary(&rg,
						rgFileName);
				break;
			case 3:
				outputType=GFF;
				inputType=BinaryInput;
				outputSubType=TextOutput;
				strcat(outputFileName, BFAST_GFF_FILE_EXTENSION);
				if(argc < 4) {
					PrintError(Name,
							"bfast reference genome file",
							"Command line argument",
							Exit,
							OutOfRange);
				}
				strcpy(rgFileName, argv[3]);
				RGBinaryReadBinary(&rg,
						rgFileName);
				break;
			case 4:
				outputType=SAM;
				inputType=BinaryInput;
				outputSubType=TextOutput;
				strcat(outputFileName, BFAST_SAM_FILE_EXTENSION);
				if(argc < 4) {
					PrintError(Name,
							"bfast reference genome file",
							"Command line argument",
							Exit,
							OutOfRange);
				}
				strcpy(rgFileName, argv[3]);
				RGBinaryReadBinary(&rg,
						rgFileName);
				if(argc < 5) {
					PrintError(Name,
							"output ID",
							"Command line argument",
							Exit,
							OutOfRange);
				}
				strcpy(outputID, argv[4]);
				break;
			default:
				PrintError(Name,
						NULL,
						"Could not understand output type",
						Exit,
						OutOfRange);
		}

		/* Open the input file */
		if(!(fpIn=fopen(inputFileName, "rb"))) {
			PrintError(Name,
					inputFileName,
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}
		/* Open the output file */
		if(!(fpOut=fopen(outputFileName, "wb"))) {
			PrintError(Name,
					outputFileName,
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}
		
		/* Print Header */
		AlignedReadConvertPrintHeader(fpOut, &rg, outputType);
		/* Initialize */
		AlignedReadInitialize(&a);
		counter = 0;
		fprintf(stderr, "Currently on:\n0");
		/* Read in each match */
		while(EOF != AlignedReadRead(&a, fpIn, inputType)) {
			if(counter%BAFCONVERT_ROTATE_NUM==0) {
				fprintf(stderr, "\r%lld",
						counter);
			}
			counter++;
			/* Print each match */
			AlignedReadConvertPrintOutputFormat(&a,
					&rg,
					fpOut,
					outputID,
					outputType,
					outputSubType);
			AlignedReadFree(&a);
		}
		fprintf(stderr, "\r%lld\n",
				counter);
		/* Close the input file */
		fclose(fpIn);
		/* Close the output file */
		fclose(fpOut);
		if(MAF == outputType ||
				SAM == outputType) {
			RGBinaryDelete(&rg);
		}

		fprintf(stderr, "Terminating successfully!\n");
	}
	else {
		fprintf(stderr, "Usage:%s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast aligned file name>\n");
		fprintf(stderr, "\t<output type:\n"
				"\t\t0-BAF text to BAF binary\n"
				"\t\t1-BAF binary to BAF text\n"
				"\t\t2-BAF binary to MAF\n"
				"\t\t3-BAF binary to GFF (v2)\n"
				"\t\t4-BAF binary to SAM (v.%s)\n",
				BFAST_SAM_VERSION
			   );
		fprintf(stderr, "\t<bfast reference genome file (not required for BAF output)\n");
		fprintf(stderr, "\t<output ID (required only for SAM output)>\n");
	}
	return 0;
}
