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
	FILE *fp=NULL;
	char inputFileName[MAX_FILENAME_LENGTH]="\0";

	if(argc == 2) {

		strcpy(inputFileName, argv[1]);

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
					Exit,
					OutOfRange);
		}
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<input file>\n");
	}

	return 0;
}


