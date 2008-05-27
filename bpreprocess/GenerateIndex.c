#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/RGIndex.h"
#include "../blib/RGIndexLayout.h"
#include "../blib/RGBinary.h"
#include "Definitions.h"
#include "GenerateIndex.h"

/* TODO */
void GenerateIndex(RGBinary *rg,
		RGIndexLayout *rgLayout,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int numThreads,
		char *outputID,
		char *outputDir,
		int binaryOutput)
{
	int32_t i;
	char outputFileName[ MAX_FILENAME_LENGTH]="\0"; 
	FILE *fp=NULL;

	RGIndex index;

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
	}

	for(i=0;i<rgLayout->numIndexes;i++) { /* For each index to create */

		/* Create the index */
		if(VERBOSE >=0) {
			fprintf(stderr, "Creating the index...\n");
		}
		RGIndexCreate(&index, 
				rg, 
				rgLayout,
				startChr,
				startPos,
				endChr,
				endPos,
				i,
				numThreads,
				1, /* Include repeat masker */
				0); /* Do not include Ns */

		if(VERBOSE >= 0) {
			fprintf(stderr, "Index created.\n");
			fprintf(stderr, "Index size is %.3lfGB.\n",
					RGIndexGetSize(&index, GIGABYTES));
		}

		/* Create the file name */
		sprintf(outputFileName, "%sbfast.index.file.%s.%d.%d.%d.%d.%d.%s",
				outputDir,
				outputID,
				index.startChr,
				index.startPos,
				index.endChr,
				index.endPos,
				i+1,
				BFAST_INDEX_FILE_EXTENSION);

		/* Open the output file */
		if(!(fp=fopen(outputFileName, "wb"))) {
			PrintError("GenerateIndex",
					outputFileName,
					"Could not open outputFileName for writing",
					Exit,
					OpenFileError);
		}   

		/* Output the index to file */
		if(VERBOSE >= 0) {
			fprintf(stderr, "Outputting index to %s\n", outputFileName);
		}
		RGIndexPrint(fp, &index, binaryOutput);
		/* Close the output file */
		fclose(fp);

		/* Free memory */
		if(VERBOSE >= 0) {
			fprintf(stderr, "Deleting index...\n");
		}
		RGIndexDelete(&index);

		if(VERBOSE >= 0) {
			fprintf(stderr, "Generated index successfully!\n");
			fprintf(stderr, "%s", BREAK_LINE);
		}
	}
}
