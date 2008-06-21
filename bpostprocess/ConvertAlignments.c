#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include "../blib/BError.h"
#include "../blib/AlignEntry.h"
#include "../blib/BLib.h"
#include "Definitions.h"
#include "PrintOutputFiles.h"
#include "ConvertAlignments.h"

/* TODO */
void ConvertAlignments(char *inputFileName,
		int inputFormat,
		int uniqueMatches,
		int bestScore,
		int minScore,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int regionLength,
		char *outputID,
		char *outputDir,
		char *tmpDir,
		int outputFormat)
{
	char *FnName="ConvertAlignments";
	int i, j;
	FILE *inputFP;
	ChrFiles chrFiles;
	RGFiles rgFiles;
	int numReads, tmpNumReads;

	/* Open the input file */
	if(0==(inputFP=fopen(inputFileName, "r"))) {
		PrintError(FnName,
				inputFileName,
				"Could not open inputFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Split the input file by chromosome and filter while we're at it */
	numReads=PrintAlignEntriesToTempFilesByChr(inputFP,
			&rgFiles,
			uniqueMatches,
			bestScore,
			minScore,
			startChr,
			startPos,
			endChr,
			endPos,
			tmpDir);

	/* Close file */
	fclose(inputFP);

	/* Go through each chromosome */
	tmpNumReads = 0;
	for(i=startChr;i<=endChr;i++) {

		/* Initialize chrFiles */
		chrFiles.files=NULL;
		chrFiles.numFiles=0;

		/* Move to the beginning of the tmp file */
		fseek(rgFiles.chrFiles[i-rgFiles.startChr], 0, SEEK_SET);

		/* Print the input file into tmp files for each 
		 * region.
		 * */
		assert(inputFormat == BAlignFile);
		tmpNumReads += PrintAlignEntriesToTempFilesWithinChr(rgFiles.chrFiles[i-rgFiles.startChr],
				uniqueMatches,
				bestScore,
				minScore,
				i,
				(i==startChr)?startPos:0,
				(i==endChr)?endPos:INT_MAX,
				regionLength,
				tmpDir,
				&chrFiles);

		/* Move tmp files to the start */
		for(j=0;j<chrFiles.numFiles;j++) {
			fseek(chrFiles.files[j], 0, SEEK_SET);
		}

		PrintAlignEntries(&chrFiles,
				regionLength,
				outputFormat,
				outputDir,
				outputID,
				i);

		/* Free memory */
		free(chrFiles.files);
		chrFiles.files = NULL;
		free(chrFiles.fileNames);
		chrFiles.fileNames = NULL;

		/* Close tmp file for the chromosome */
		CloseTmpFile(&rgFiles.chrFiles[i-rgFiles.startChr], &rgFiles.chrFileNames[i-rgFiles.startChr]);
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Free memory */
	free(rgFiles.chrFiles);
	rgFiles.chrFiles = NULL;
	free(rgFiles.chrFileNames);
	rgFiles.chrFileNames = NULL;

	assert(tmpNumReads == numReads);
}
