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
	char outputRange[2][MAX_FILENAME_LENGTH]={"\0","\0"};
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *inputFP=NULL;
	FILE *outputFP=NULL;
	AlignEntries a;
	int64_t numRead, numPrinted, i, numToSatisfy;
	Range start, end;

	if(argc == 6) {
		strcpy(outputRange[0], argv[1]);
		strcpy(outputRange[1], argv[2]);
		numToSatisfy=atoi(argv[3]);
		assert(0 == numToSatisfy || 1 == numToSatisfy);
		strcpy(outputID, argv[4]);

		ParseRange(&start, outputRange[0]);
		ParseRange(&end, outputRange[1]);

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
		fprintf(stderr, "Outputting to %s.\n",
				outputFileName);

		numRead = numPrinted = 0;
		for(i=5;i<argc;i++) {
			strcpy(inputFileName, argv[i]);
			if(!(inputFP = fopen(inputFileName, "rb"))) {
				PrintError(Name,
						inputFileName,
						"Could not open file for reading",
						Exit,
						OpenFileError);
			}

			AlignEntriesInitialize(&a);
			fprintf(stderr, "Reading in from %s.\nCurrently on:\n0",
					inputFileName);
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
					if(numToSatisfy < CheckRange(&start, a.entriesOne[0].contig, a.entriesOne[0].position) + 
							CheckRange(&end, a.entriesTwo[0].contig, a.entriesTwo[0].position)) {
						AlignEntriesPrint(&a,
								outputFP,
								BinaryOutput);
						numPrinted++;
					}
				}
			}
			fprintf(stderr, "\r%lld\n",
					(long long int)numRead);
			fclose(inputFP);
		}

		/* Close files */
		fclose(outputFP);

		fprintf(stderr, "Read in %lld and outputted %lld paired end alignments.\n",
				(long long int)numRead,
				(long long int)numPrinted);
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully!\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<end one range (contig1-contig2:pos1-pos2)\n");
		fprintf(stderr, "\t<end two range (contig1-contig2:pos1-pos2)\n");
		fprintf(stderr, "\t<require both 0: require one range to be satisfied 1: required both ranges to be satisfied>\n");
		fprintf(stderr, "\t<output ID>\n");
		fprintf(stderr, "\t<bfast report file names>\n");
	}
	return 0;
}

void ParseRange(Range *r,
		char *string)
{
	char *FnName="ParseRange";
	if(4 != sscanf(string, "%d-%d:%d-%d\n",
				&r->contigStart,
				&r->contigEnd,
				&r->positionStart,
				&r->positionEnd)) {
		PrintError(FnName,
				string,
				"Could not parse string.  Should be in %d-%d:%d-%d format",
				Exit,
				OutOfRange);
	}
}

int32_t CheckRange(Range *r,
		int32_t contig,
		int32_t position)
{
	if(r->contigStart < contig ||
			(r->contigStart == contig && r->positionStart <= position)) {
		if(contig < r->contigEnd ||
				(contig == r->contigEnd && position <= r->positionEnd)) {
			return 1;
		}
	}
	return 0;
}
