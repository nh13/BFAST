#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/RGIndex.h"
#include "bheader.h"

#define Name "bheader"

/* Prints the header of a bfast index file; the header completely 
 * defines the index. */

int main(int argc, char *argv[]) 
{
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	int i;

	if(2 <= argc) {
		for(i=1;i<argc;i++) {
			strcpy(inputFileName, argv[i]);

			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Getting info for %s.\n", inputFileName);
			if(NULL!=strstr(inputFileName, BFAST_RG_FILE_EXTENSION)) {
				RGBinaryPrintInfo(inputFileName);
			}
			else if(NULL!=strstr(inputFileName, BFAST_INDEX_FILE_EXTENSION)) {
				RGIndexPrintInfo(inputFileName);
			}
			else {
				PrintError(Name,
						"input file",
						"Could not recognize input file extension",
						Warn,
						OutOfRange);
			}
			fprintf(stderr, "%s", BREAK_LINE);
		}
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<input file>\n");
	}

	return 0;
}


