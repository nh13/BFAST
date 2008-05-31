#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/RGIndex.h"
#include "bfixhash.h"

int main(int argc, char *argv[]) 
{
	int64_t i;
	FILE *fp=NULL;
	char indexFileName[MAX_FILENAME_LENGTH]="\0";
	char rgFileName[MAX_FILENAME_LENGTH]="\0";

	if(argc == 3) {
		RGBinary rg;
		RGIndex index;

		strcpy(rgFileName, argv[1]);
		strcpy(indexFileName, argv[2]);

		/* Read in the rg binary file */
		RGBinaryReadBinary(&rg, rgFileName);

		if(!(fp=fopen(indexFileName, "rb"))) {
			PrintError("bfixhash",
					indexFileName,
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}

		/* Read the index */
		RGIndexRead(fp, &index, 1);
		fclose(fp);

		/* Test that it is sorted correctly */
		for(i=1;i<index.length;i++) {
			assert(RGIndexCompareAt(&index, &rg, i-1, i) <= 0);
		}


		/* Fix the hash by recreating it */
		RGIndexCreateHash(&index, &rg);

		/* Print the new index */ 
		if(!(fp=fopen(indexFileName, "wb"))) {
			PrintError("bfixhash",
					indexFileName,
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}
		RGIndexPrint(fp, &index, 1);
		fclose(fp);

		/* Delete the index */
		RGIndexDelete(&index);

		/* Delete the rg */
		RGBinaryDelete(&rg);

	}
	else {
		fprintf(stdout, "Please give a reference genome file name then an index file name.  Terminating!\n");
	}

	return 0;
}


