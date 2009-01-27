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
#include "btranslocations.h"

#define Name "btranslocations"
#define BSORT_ROTATE_NUM 100000

/* Outputs unique paired end alignments for which each end is
 * on a different contig.
 * */

int main(int argc, char *argv[])
{
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	char outputID[MAX_FILENAME_LENGTH]="\0";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *inputFP=NULL;
	FILE *outputFP=NULL;
	AlignEntries a;
	int64_t numRead, numPrinted;

	if(argc == 3) {
		strcpy(inputFileName, argv[1]);
		strcpy(outputID, argv[2]);

		if(!(inputFP = fopen(inputFileName, "rb"))) {
			PrintError(Name,
					inputFileName,
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}
		
		/* Create output file name 
		 * TODO */
		sprintf(outputFileName, "bfast.translocations.%s.baf",
				outputID);
		if(!(outputFP = fopen(outputFileName, "wb"))) {
			PrintError(Name,
					outputFileName,
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}

		numRead = numPrinted = 0;
		AlignEntriesInitialize(&a);
		fprintf(stderr, "Reading in from %s.\nOutputting to %s.\nCurrently on:\n0",
				inputFileName,
				outputFileName);
		while(EOF != AlignEntriesRead(&a,
					inputFP,
					PairedEnd,
					SpaceDoesNotMatter,
					BinaryInput)) {
			numRead++;
			if(0 == numRead%BTRANSLOCATIONS_ROTATE_NUM) {
				fprintf(stderr, "\r%lld",
						(long long int)numRead);
			}
			if(a.numEntriesOne == 1 &&
					a.numEntriesTwo == 1 &&
					a.entriesOne[0].contig != a.entriesTwo[0].contig) {
				AlignEntriesPrint(&a,
						outputFP,
						BinaryOutput);
				numPrinted++;
			}
		}
				fprintf(stderr, "\r%lld\n",
						(long long int)numRead);

		/* Close files */
		fclose(outputFP);
		fclose(inputFP);

		fprintf(stderr, "Read in %lld and outputted %lld paired end alignments.\n",
				(long long int)numRead,
				(long long int)numPrinted);
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully!\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast report file name>\n");
		fprintf(stderr, "\t<output ID>\n");
	}
	return 0;
}
