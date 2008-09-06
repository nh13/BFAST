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
		int colorSpace,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int repeatMasker,
		int numThreads,
		char *outputID,
		char *outputDir,
		char *tmpDir,
		int binaryOutput)
{
	char *FnName="GenerateIndex";
	int32_t i;
	char outputFileName[ MAX_FILENAME_LENGTH]="\0"; 
	FILE *fp=NULL;
	RGIndex index;

	/* Adjust start and end based on reference genome */
	/* Adjust start */
	if(startChr < rg->startChr) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Warning: startChr was less than reference genome's start chromosome.\n");
			fprintf(stderr, "Defaulting to reference genome's start chromosome and position: chr%d:%d.\n",
					rg->startChr,
					rg->startPos);
		}
		startChr = rg->startChr;
		startPos = rg->startPos;
	}
	else if(startChr == rg->startChr &&
			startPos < rg->startPos) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Warning: startPos was less than reference genome's start position.\n");
			fprintf(stderr, "Defaulting to reference genome's start position: %d.\n",
					rg->startPos);
		}
		startPos = rg->startPos;
	}
	/* Adjust end */
	if(endChr > rg->endChr) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Warning: endChr was greater than reference genome's end chromosome.\n");
			fprintf(stderr, "Defaulting to reference genome's end chromosome and position: chr%d:%d.\n",
					rg->endChr,
					rg->endPos);
		}
		endChr = rg->endChr;
		endPos = rg->endPos;
	}
	else if(endChr == rg->endChr &&
			endPos > rg->endPos) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Warning: endPos was greater than reference genome's end position.\n");
			fprintf(stderr, "Defaulting to reference genome's end position: %d.\n",
					rg->endPos);
		}
		endPos = rg->endPos;
	}

	/* Check that the start and end bounds are ok */
	if(startChr > endChr) {
		fprintf(stderr, "startChr:%d\nendChr:%d\n",
				startChr,
				endChr);
		PrintError(FnName,
				NULL,
				"The start chromosome is greater than the end chromosome",
				Exit,
				OutOfRange);
	}
	else if(startChr == endChr &&
			startPos > endPos) {
		PrintError(FnName,
				NULL,
				"The start position is greater than the end position on the same chromosome",
				Exit,
				OutOfRange);
	}

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
				colorSpace,
				startChr,
				startPos,
				endChr,
				endPos,
				i,
				numThreads,
				repeatMasker,
				0, /* Do not include Ns */
				tmpDir);

		if(VERBOSE >= 0) {
			fprintf(stderr, "Index created.\n");
			fprintf(stderr, "Index size is %.3lfGB.\n",
					RGIndexGetSize(&index, GIGABYTES));
		}

		/* Create the file name */
		sprintf(outputFileName, "%s%s.index.file.%s.%d.%d.%d.%d.%d.%d.%s",
				outputDir,
				PROGRAM_NAME,
				outputID,
				index.startChr,
				index.startPos,
				index.endChr,
				index.endPos,
				repeatMasker,
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
