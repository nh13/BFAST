#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include "../blib/AlignedEntry.h"
#include "../blib/AlignedRead.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BLib.h"
#include "../blib/BError.h"
#include "bmergesorted.h"

#define Name "bmergesorted"
#define BMERGESORTED_ROTATE_NUM 100000
#define BSORT_MAX_LINE_LENGTH 100

/* Merges two sorted baf files
 * */

void MergeFiles(FILE *fpOne,
		FILE *fpTwo,
		FILE *fpOut)
{
	int64_t ctr=0;
	fpos_t posOne;
	fpos_t posTwo;
	AlignedRead aOne, aTwo;


	fprintf(stderr, "Currently on:\n0");
	AlignedReadInitialize(&aOne);
	AlignedReadInitialize(&aTwo);
	while(0 == feof(fpOne) &&
			0 == feof(fpTwo)) {

		assert(0 == fgetpos(fpOne, &posOne));
		assert(0 == fgetpos(fpTwo, &posTwo));

		if(EOF == AlignedReadRead(&aOne, fpOne, BinaryInput) ||
				EOF == AlignedReadRead(&aTwo, fpTwo, BinaryInput)) {
			if(0 == feof(fpOne)) {
				assert(0 == fsetpos(fpOne, &posOne));
			}
			if(0 == feof(fpTwo)) {
				assert(0 == fsetpos(fpTwo, &posTwo));
			}
		}
		else {
			ctr++;
			if(0 == ctr % BMERGESORTED_ROTATE_NUM) {
				fprintf(stderr, "\r%lld",
						(long long int)ctr);
			}
			if(AlignedReadCompareAll(&aOne, &aTwo) <= 0) {
				AlignedReadPrint(&aOne, fpOut, BinaryOutput);
				fsetpos(fpTwo, &posTwo);
			}
			else {
				AlignedReadPrint(&aTwo, fpOut, BinaryOutput);
				fsetpos(fpOne, &posOne);
			}
		}
		AlignedReadFree(&aOne);
		AlignedReadFree(&aTwo);
	}
	while(EOF != AlignedReadRead(&aOne, fpOne, BinaryInput)) {
		ctr++;
		if(0 == ctr % BMERGESORTED_ROTATE_NUM) {
			fprintf(stderr, "\r%lld",
					(long long int)ctr);
		}
		AlignedReadPrint(&aOne, fpOut, BinaryOutput);
		AlignedReadFree(&aOne);
	}
	while(EOF != AlignedReadRead(&aTwo, fpTwo, BinaryInput)) {
		ctr++;
		if(0 == ctr % BMERGESORTED_ROTATE_NUM) {
			fprintf(stderr, "\r%lld",
					(long long int)ctr);
		}
		AlignedReadPrint(&aTwo, fpOut, BinaryOutput);
		AlignedReadFree(&aTwo);
	}
	fprintf(stderr, "\r%lld\n",
			(long long int)ctr);
}


int main(int argc, char *argv[])
{
	char inputFileName[2][MAX_FILENAME_LENGTH]={"\0","\0"};
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *inputFP[2]={NULL,NULL};
	FILE *outputFP=NULL;
	int i;

	if(argc == 4) {
		strcpy(inputFileName[0], argv[1]);
		strcpy(inputFileName[1], argv[2]);
		sprintf(outputFileName, "bfast.%s.baf",
				argv[3]);

		/* Open input */
		for(i=0;i<2;i++) {
			if(!(inputFP[i] = fopen(inputFileName[i], "rb"))) {
				PrintError(Name,
						inputFileName[i],
						"Could not open file for reading",
						Exit,
						OpenFileError);
			}
		}
		/* Open output */
		if(!(outputFP = fopen(outputFileName, "wb"))) {
			PrintError(Name,
					outputFileName,
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}

		/* Split entries and print */
		MergeFiles(inputFP[0],
				inputFP[1],
				outputFP);

		/* Close files */
		fclose(outputFP);
		/* Close input */
		for(i=0;i<2;i++) {
			fclose(inputFP[i]);
		}

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully!\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast *sorted* report file name>\n");
		fprintf(stderr, "\t<bfast *sorted* report file name>\n");
		fprintf(stderr, "\t<output ID>\n");
	}
	return 0;
}
