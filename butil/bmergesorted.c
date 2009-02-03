#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include "../blib/AlignEntry.h"
#include "../blib/AlignEntries.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BLib.h"
#include "../blib/BError.h"
#include "bsort.h"

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
	AlignEntries aOne, aTwo;


	fprintf(stderr, "Currently on:\n0");
	AlignEntriesInitialize(&aOne);
	AlignEntriesInitialize(&aTwo);
	while(0 == feof(fpOne) &&
			0 == feof(fpTwo)) {

		assert(0 == fgetpos(fpOne, &posOne));
		assert(0 == fgetpos(fpTwo, &posTwo));

		if(EOF == AlignEntriesRead(&aOne, fpOne, PairedEndDoesNotMatter, SpaceDoesNotMatter, BinaryInput) ||
				EOF == AlignEntriesRead(&aTwo, fpTwo, PairedEndDoesNotMatter, SpaceDoesNotMatter, BinaryInput)) {
			assert(0 == fsetpos(fpOne, &posOne));
			assert(0 == fsetpos(fpTwo, &posTwo));
		}
		else {
			ctr++;
			if(0 == ctr % BMERGESORTED_ROTATE_NUM) {
				fprintf(stderr, "\r%lld",
						(long long int)ctr);
			}
			if(AlignEntriesCompareAll(&aOne, &aTwo) <= 0) {
				AlignEntriesPrint(&aOne, fpOut, BinaryOutput);
				fsetpos(fpTwo, &posTwo);
			}
			else {
				AlignEntriesPrint(&aTwo, fpOut, BinaryOutput);
				fsetpos(fpOne, &posOne);
			}
		}
		AlignEntriesFree(&aOne);
		AlignEntriesFree(&aTwo);
	}
	while(EOF != AlignEntriesRead(&aOne, fpOne, PairedEndDoesNotMatter, SpaceDoesNotMatter, BinaryInput)) {
		ctr++;
		if(0 == ctr % BMERGESORTED_ROTATE_NUM) {
			fprintf(stderr, "\r%lld",
					(long long int)ctr);
		}
		AlignEntriesPrint(&aOne, fpOut, BinaryOutput);
		AlignEntriesFree(&aOne);
	}
	while(EOF != AlignEntriesRead(&aTwo, fpTwo, PairedEndDoesNotMatter, SpaceDoesNotMatter, BinaryInput)) {
		ctr++;
		if(0 == ctr % BMERGESORTED_ROTATE_NUM) {
			fprintf(stderr, "\r%lld",
					(long long int)ctr);
		}
		AlignEntriesPrint(&aTwo, fpOut, BinaryOutput);
		AlignEntriesFree(&aTwo);
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
