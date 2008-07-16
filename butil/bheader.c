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
		if(!(fp=fopen(inputFileName, "rb"))) {
			PrintError("bheader",
					inputFileName,
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}

		if(NULL!=strstr(inputFileName, BFAST_RG_FILE_EXTENSION)) {
			fprintf(stdout, "Warning.  Not implemented.\n");
		}
		else if(NULL!=strstr(inputFileName, BFAST_INDEX_FILE_EXTENSION)) {
			RGIndexPrintInfo(fp, 1);
		}

		fclose(fp);
	}
	else {
		fprintf(stderr, "%s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast index file>\n");
	}

	return 0;
}


