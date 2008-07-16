#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "../blib/RGMatches.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "bmfconvert.h"

#define Name "bmfconvert"
#define BMFCONVERT_ROTATE_NUM 100000

/* Converts a bmatches file from binary to plaintext or vice versa.
 * */

int main(int argc, char *argv[])
{

	FILE *fpIn, *fpOut;
	int binaryInput, binaryOutput;
	int pairedEnd;
	long long int counter;
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	RGMatches m;

	if(argc == 4) {
		strcpy(inputFileName, argv[1]);
		binaryInput = atoi(argv[2]);
		pairedEnd = atoi(argv[3]);

		assert(0==binaryInput || 1==binaryInput);
		assert(0==pairedEnd || 1==pairedEnd);

		/* Creat output file name */
		sprintf(outputFileName, "%s.converted",
				inputFileName);

		/* Set binary for output */
		binaryOutput = (binaryInput+1)%2;
		assert(0==binaryOutput || 1==binaryOutput);
		assert(binaryInput!=binaryOutput);

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
		while(EOF != RGMatchesRead(fpIn, &m, pairedEnd, binaryInput)) {
			if(counter%BMFCONVERT_ROTATE_NUM==0) {
				fprintf(stderr, "\r%lld",
						counter);
			}
			counter++;
			/* Print each match */
			RGMatchesPrint(fpOut, &m, pairedEnd, binaryOutput);
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
		fprintf(stderr, "%s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<matches file name>\n");
		fprintf(stderr, "\t<input type: 0-text 1-binary>\n");
		fprintf(stderr, "\t<0-single-end 1-paired end>\n");
	}
	return 0;
}
