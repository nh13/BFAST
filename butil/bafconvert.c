#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/AlignEntries.h"
#include "../blib/BError.h"
#include "bafconvert.h"

#define Name "bafconvert"
#define BAFCONVERT_ROTATE_NUM 100000

/* Converts a bmatches file from binary to plaintext or vice versa.
 * */

int main(int argc, char *argv[])
{

	FILE *fpIn, *fpOut;
	int binaryInput, binaryOutput;
	long long int counter;
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	AlignEntries a;

	if(argc == 4) {
		strcpy(inputFileName, argv[1]);
		binaryInput = atoi(argv[2]);

		assert(TextInput==binaryInput || BinaryInput==binaryInput);

		/* Creat output file name */
		sprintf(outputFileName, "%s.converted",
				inputFileName);

		/* Set binary for output */
		binaryOutput = (binaryInput==TextInput)?BinaryOutput:TextOutput;
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
		AlignEntriesInitialize(&a);
		counter = 0;
		fprintf(stderr, "Currently on:\n0");
		/* Read in each match */
		while(EOF != AlignEntriesRead(&a, fpIn, PairedEndDoesNotMatter, SpaceDoesNotMatter, binaryInput)) {
			if(counter%BAFCONVERT_ROTATE_NUM==0) {
				fprintf(stderr, "\r%lld",
						counter);
			}
			counter++;
			/* Print each match */
			AlignEntriesPrint(&a, fpOut, binaryOutput);
			AlignEntriesFree(&a);
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
		fprintf(stderr, "\t<bfast matches file name>\n");
		fprintf(stderr, "\t<input type: 0-text 1-binary>\n");
	}
	return 0;
}
