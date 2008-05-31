#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/RGIndex.h"
#include "bheader.h"

int main(int argc, char *argv[]) 
{
	FILE *fp=NULL;
	char inputFileName[MAX_FILENAME_LENGTH]="\0";

	if(argc == 2) {

		strcpy(inputFileName, argv[1]);

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
			RGIndex index;
			RGIndexReadHeader(fp, &index, 1);
			RGIndexPrintHeader(stdout, &index, 0);

		}

		fclose(fp);
	}
	else {
		fprintf(stdout, "Please give a input file name.  Terminating!\n");
	}

	return 0;
}


