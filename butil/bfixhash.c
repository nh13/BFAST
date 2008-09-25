#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/RGIndex.h"
#include "bfixhash.h"

#define Name "bfixhash"

/* Regenerates the hash lookup table for an index given 
 * a new hash width. 
 * */

int main(int argc, char *argv[]) 
{
	FILE *fp=NULL;
	char indexFileName[MAX_FILENAME_LENGTH]="\0";
	char rgFileName[MAX_FILENAME_LENGTH]="\0";
	int hashWidth;

	if(argc == 4) {
		RGBinary rg;
		RGIndex index;

		strcpy(rgFileName, argv[1]);
		strcpy(indexFileName, argv[2]);
		hashWidth = atoi(argv[3]);

		/* Read in the rg binary file */
		RGBinaryReadBinary(&rg, rgFileName);

		/* Read the index */
		fprintf(stderr, "Reading in index from %s.\n",
				indexFileName);
		if(!(fp=fopen(indexFileName, "rb"))) {
			PrintError("bfixhash",
					indexFileName,
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}
		RGIndexRead(fp, &index, 1);
		fclose(fp);

		/* Free hash */
		free(index.starts);
		index.starts=NULL;
		free(index.ends);
		index.ends=NULL;

		/* Update to new width */
		index.hashWidth = hashWidth;
		index.hashLength = pow(4, index.hashWidth);

		/* Fix the hash by recreating it */
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Fixing hash.\n");
		RGIndexCreateHash(&index, &rg);

		/* Print the new index */ 
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Outputting to %s.\n", indexFileName);
		if(!(fp=fopen(indexFileName, "wb"))) {
			PrintError("bfixhash",
					indexFileName,
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}
		RGIndexPrint(fp, &index, 1);
		fclose(fp);

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Cleaning up.\n");
		/* Delete the index */
		RGIndexDelete(&index);
		/* Delete the rg */
		RGBinaryDelete(&rg);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast reference genome file>\n");
		fprintf(stderr, "\t<bfast index file>\n");
		fprintf(stderr, "\t<new hash width>\n");
	}

	return 0;
}


