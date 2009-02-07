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

/* Converts a bmatches file from binary to plaintext or vice versa.
 * */

int main(int argc, char *argv[])
{

	FILE *fpIn, *fpOut;
	int binaryInput, binaryOutput;
	long long int counter;
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	char *last;
	RGMatches m;

	if(argc == 3) {
		strcpy(inputFileName, argv[1]);
		binaryInput = atoi(argv[2]);

		assert(TextInput==binaryInput || BinaryInput==binaryInput);

		/* Set binary for output */
		binaryOutput = (binaryInput==TextInput)?BinaryOutput:TextOutput;
		assert(binaryInput!=binaryOutput);

		/* Create output file name */
		last = StrStrGetLast(inputFileName,
				BFAST_MATCHES_FILE_EXTENSION);
		if(NULL == last) {
			PrintError(Name,
					inputFileName,
					"Could not recognize file extension",
					Exit,
					OutOfRange);
		} 
		if(TextOutput == binaryOutput) {
			strncpy(outputFileName, inputFileName, (last - inputFileName));
			strcat(outputFileName, "text.");
			strcat(outputFileName, BFAST_MATCHES_FILE_EXTENSION);
		}
		else {
			strncpy(outputFileName, inputFileName, (last - inputFileName));
			strcat(outputFileName, "binary.");
			strcat(outputFileName, BFAST_MATCHES_FILE_EXTENSION);
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
		while(EOF != RGMatchesRead(fpIn, &m, PairedEndDoesNotMatter, binaryInput)) {
			if(counter%BMFCONVERT_ROTATE_NUM==0) {
				fprintf(stderr, "\r%lld",
						counter);
			}
			counter++;
			/* Print each match */
			RGMatchesPrint(fpOut, &m, binaryOutput);
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
		fprintf(stderr, "\t<input type: 0-text 1-binary>\n");
	}
	return 0;
}
