#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include "../blib/BError.h"
#include "../blib/AlignEntry.h"
#include "../blib/BLib.h"
#include "Definitions.h"
#include "PrintOutputFiles.h"
#include "FilterAlignments.h"

/* TODO */
void FilterAlignments(char *inputFileName,
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
	char *FnName="FilterAlignment";
	int i, j;
	FILE *inputFP;
	ChrFiles chrFiles;
	RGFiles rgFiles;

	/* Open the input file */
	if(0==(inputFP=fopen(inputFileName, "r"))) {
		PrintError(FnName,
				inputFileName,
				"Could not open inputFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Split the input file by chromosome and filter while we're at it */
	PrintAlignEntriesToTempFilesByChr(inputFP,
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
		PrintAlignEntriesToTempFilesWithinChr(rgFiles.chrFiles[i-rgFiles.startChr],
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
}

int FilterEntries(AlignEntry **entries,
		int numEntries,
		int uniqueMatches,
		int bestScore,
		int minScore,
		int startChr,
		int startPos,
		int endChr,
		int endPos)
{
	int i;
	char *FnName = "FilterEntries";
	int minIndex = -1;
	int numMinIndexes = 0; 
	assert(uniqueMatches != 1 || bestScore != 1);
	assert(uniqueMatches != 0 || bestScore != 0);
	assert(numEntries > 0);

	/* Filter all entries with score < minScore and do not align within bounds */
	for(i=0;i<numEntries;i++) {
		/* Check filter conditions */
		if((*entries)[i].score < minScore ||
				(*entries)[i].chromosome < startChr ||
				(*entries)[i].chromosome > endChr ||
				((*entries)[i].chromosome == startChr && ((*entries)[i].position + (*entries)[i].length) < startPos) ||  
				((*entries)[i].chromosome == endPos && (*entries)[i].position > endPos)) {
			/* Filter */

			/* Copy entry at end to here */
			if(i<numEntries-1) {
				AlignEntryCopy(&(*entries)[numEntries-1], &(*entries)[i]);
			}
			/* Free memory and reallocate */
			AlignEntryFree(&(*entries)[numEntries-1]);
			numEntries--;
			(*entries) = realloc((*entries), sizeof(AlignEntry)*numEntries);
			if(numEntries != 0 && NULL == (*entries)) {
				PrintError(FnName,
						"(*entries)",
						"Could not reallocate memory",
						Exit,
						ReallocMemory);
			}
			/* Decrement i since we want to reconsider the current position (we just swapped) */
			if(i<numEntries-1) {
				i--;
			}
		}
		else if(minIndex == -1 || (*entries)[minIndex].score <= (*entries)[i].score) {
			if(minIndex >= 0 && (*entries)[minIndex].score == (*entries)[i].score) {
				/* We should store if there are more than one place with the same score */
				numMinIndexes++;
			}
			else {
				/* Get the entry with the best score */
				minIndex = i;
				numMinIndexes=1;
			}
		}
	}

	/* Apply uniqueMatches and bestScore */
	if(uniqueMatches == 1) {
		if(1==numEntries) {
			return 1;
		}
		else {
			/* Free all memory and return 0 */
			/* Free memory and reallocate */
			for(i=0;i<numEntries;i++) {
				AlignEntryFree(&(*entries)[i]);
			}
			free((*entries));
			(*entries)=NULL;
			return 0;
		}
	}
	else if(bestScore == 1) {
		if(numMinIndexes != 1) {
			/* We have multiple entries with the same best score */
			/* Free all memory and return 0 */
			/* Free memory and reallocate */
			for(i=0;i<numEntries;i++) {
				AlignEntryFree(&(*entries)[i]);
			}
			free((*entries));
			(*entries)=NULL;
			return 0;
		}
		else {
			/* Should have only one entry with the best score */
			assert(minIndex >=0);
			assert(numMinIndexes == 1);
			/* Copy min index to the front */
			AlignEntryCopy(&(*entries)[minIndex], &(*entries)[0]);
			/* Free memory and reallocate */
			for(i=1;i<numEntries;i++) {
				AlignEntryFree(&(*entries)[i]);
			}
			numEntries=1;
			(*entries) = realloc((*entries), sizeof(AlignEntry)*numEntries);
			if(NULL == (*entries)) {
				PrintError(FnName,
						"(*entries)",
						"Could not reallocate memory",
						Exit,
						ReallocMemory);
			}
			return 1;
		}
	}
	PrintError(FnName,
			NULL,
			"Control reached unintended point",
			Exit,
			OutOfRange);
	return -1; 
}
