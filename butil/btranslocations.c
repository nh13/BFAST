#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <zlib.h>

#include "../blib/AlignedEntry.h"
#include "../blib/AlignedRead.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BLib.h"
#include "../blib/BError.h"
#include "btranslocations.h"

#define Name "btranslocations"
#define BSORT_ROTATE_NUM 100000

/* Outputs unique paired two alignments for which each two is
 * on a different contig.
 * */

int main(int argc, char *argv[])
{
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	char outputID[MAX_FILENAME_LENGTH]="\0";
	char outputRange[2][MAX_FILENAME_LENGTH]={"\0","\0"};
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	gzFile inputFP=NULL;
	gzFile outputFP=NULL;
	AlignedRead a;
	int64_t numRead, numPrinted, i, numToSatisfy;
	Range one, two;

	if(6 <= argc) {
		strcpy(outputRange[0], argv[1]);
		strcpy(outputRange[1], argv[2]);
		numToSatisfy=atoi(argv[3]);
		assert(0 == numToSatisfy || 1 == numToSatisfy);
		strcpy(outputID, argv[4]);

		ParseRange(&one, outputRange[0]);
		ParseRange(&two, outputRange[1]);

		/* Create output file name 
		 * TODO */
		sprintf(outputFileName, "bfast.translocation.file.%s.%s.%s.baf",
				outputID,
				outputRange[0],
				outputRange[1]);
		if(!(outputFP = gzopen(outputFileName, "wb"))) {
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
			if(!(inputFP = gzopen(inputFileName, "rb"))) {
				PrintError(Name,
						inputFileName,
						"Could not open file for reading",
						Exit,
						OpenFileError);
			}

			AlignedReadInitialize(&a);
			fprintf(stderr, "Reading in from %s.\nCurrently on:\n0",
					inputFileName);
			while(EOF != AlignedReadRead(&a,
						inputFP)) {
				numRead++;
				if(0 == numRead%BTRANSLOCATIONS_ROTATE_NUM) {
					fprintf(stderr, "\r%lld",
							(long long int)numRead);
				}
				if(2 == a.numEnds &&
						1 == a.ends[0].numEntries &&
						1 == a.ends[1].numEntries) {
					if(numToSatisfy < CheckRange(&one, a.ends[0].entries[0].contig, a.ends[0].entries[0].position) + 
							CheckRange(&two, a.ends[1].entries[0].contig, a.ends[1].entries[0].position)) {
						AlignedReadPrint(&a,
								outputFP);
						fflush(outputFP);
						numPrinted++;
					}
					else if(numToSatisfy < CheckRange(&two, a.ends[0].entries[0].contig, a.ends[0].entries[0].position) + 
							CheckRange(&one, a.ends[1].entries[0].contig, a.ends[1].entries[0].position)) {
						AlignedReadPrint(&a,
								outputFP);
						fflush(outputFP);
						numPrinted++;
					}
				}
				AlignedReadFree(&a);
			}
			fprintf(stderr, "\r%lld\n",
					(long long int)numRead);
			gzclose(inputFP);
		}

		/* Close files */
		gzclose(outputFP);

		fprintf(stderr, "Read in %lld and outputted %lld paired two alignments.\n",
				(long long int)numRead,
				(long long int)numPrinted);
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully!\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<range one (contig1-contig2:pos1-pos2)\n");
		fprintf(stderr, "\t<range two (contig1-contig2:pos1-pos2)\n");
		fprintf(stderr, "\t<0: require one range to be satisfied 1: required both ranges to be satisfied>\n");
		fprintf(stderr, "\t<output ID>\n");
		fprintf(stderr, "\t<bfast report file names>\n");
	}
	return 0;
}
