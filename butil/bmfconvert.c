#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "../blib/RGMatches.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "bmfconvert.h"

#define Name "bmfconvert"
#define BMFCONVERT_ROTATE_NUM 100000
#define BMFCONVERT_FASTQ 2

/* Converts a bmatches file from binary to plaintext or vice versa.
 * */

int main(int argc, char *argv[])
{

	FILE *fpIn, *fpOut;
	int binaryInput = 0, binaryOutput = 0;
	long long int counter;
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	int outputType = 0;
	char *last;
	RGMatches m;

	if(argc == 3) {
		strcpy(inputFileName, argv[1]);
		outputType = atoi(argv[2]);

		last = StrStrGetLast(inputFileName,
				BFAST_MATCHES_FILE_EXTENSION);
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
				binaryInput = TextInput;
				binaryOutput = BinaryOutput;
				strcat(outputFileName, "text.");
				strcat(outputFileName, BFAST_MATCHES_FILE_EXTENSION);
				break;
			case 1:
				binaryInput = BinaryInput;
				binaryOutput = TextOutput;
				strncpy(outputFileName, inputFileName, (last - inputFileName));
				strcat(outputFileName, "binary.");
				strcat(outputFileName, BFAST_MATCHES_FILE_EXTENSION);
				break;
			case 2:
				binaryInput = BinaryInput;
				strcat(outputFileName, BFAST_MATCHES_READS_FILTERED_FILE_EXTENSION);
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
		/* Initialize */
		RGMatchesInitialize(&m);
		counter = 0;
		fprintf(stderr, "Currently on:\n0");
		/* Read in each match */
		while(EOF != RGMatchesRead(fpIn, &m, binaryInput)) {
			if(counter%BMFCONVERT_ROTATE_NUM==0) {
				fprintf(stderr, "\r%lld",
						counter);
			}
			counter++;
			/* Print each match */
			switch(outputType) {
				case 0:
				case 1:
					RGMatchesPrint(fpOut, &m, binaryOutput);
					break;
				case 2:
					RGMatchesPrintFastq(fpOut, &m);
					break;
				default:
					PrintError(Name,
							NULL,
							"Could not understand output type",
							Exit,
							OutOfRange);
			}
			RGMatchesFree(&m);
		}
		fprintf(stderr, "\r%lld\n",
				counter);
		/* Close the input file */
		fclose(fpIn);
		/* Close the output file */
		fclose(fpOut);

		fprintf(stderr, "Terminating successfully!\n");
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast matches file name>\n");
		fprintf(stderr, "\t\t0-BMF text to BMF binary\n");
		fprintf(stderr, "\t\t1-BMF binary to BMF text\n");
		fprintf(stderr, "\t\t2-BMF binary to FASTQ\n");
	}
	return 0;
}
