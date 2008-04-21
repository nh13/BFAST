#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include "../blib/BError.h"
#include "../blib/AlignEntry.h"
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
		char *outputID,
		char *outputDir,
		int outputFormat)
{
	char *FnName="FilterAlignment";
	int i, j;
	FILE *inputFP;
	ChrFiles chrFiles;

	/* Open the input file */
	if(0==(inputFP=fopen(inputFileName, "r"))) {
		PrintError(FnName,
				inputFileName,
				"Could not open inputFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Go through each chromosome */
	for(i=startChr;i<=endChr;i++) {

		/* Go to the beginning of the file */
		fseek(inputFP, 0, SEEK_SET);

		/* Initialize chrFiles */
		chrFiles.files=NULL;
		chrFiles.numFiles=0;

		/* Print the input file into tmp files for each 
		 * region.
		 * */
		switch(inputFormat) {
			case BAlignFile:
				PrintAlignEntriesToTempFiles(inputFP,
						uniqueMatches,
						bestScore,
						minScore,
						i,
						(i==startChr)?startPos:0,
						i,
						(i==endChr)?endPos:INT_MAX,
						&chrFiles);
				break;
			case WigFile:
			case BedFile:
				PrintError(FnName,
						"inputFormat",
						"Input format not supported",
						Exit,
						OutOfRange);
				break;
			default:
				PrintError(FnName,
						"inputFormat",
						"Input format not recognized",
						Exit,
						OutOfRange);
				break;
		}

		/* Move tmp files to the start */
		for(j=0;j<chrFiles.numFiles;j++) {
			fseek(chrFiles.files[j], 0, SEEK_SET);
		}

		/* Input format = ? */
		switch(inputFormat) {
			case BAlignFile:
				/* BAlignFile input, output format = ? */
				PrintAlignEntries(&chrFiles,
						outputFormat,
						outputDir,
						outputID,
						i);
				break;
			case WigFile:
				PrintError(FnName,
						"inputFormat",
						"Input format not supported",
						Exit,
						OutOfRange);
				break;
			case BedFile:
				PrintError(FnName,
						"inputFormat",
						"Input format not supported",
						Exit,
						OutOfRange);
				break;
			default:
				PrintError(FnName,
						"inputFormat",
						"Input format not recognized",
						Exit,
						OutOfRange);
				break;
		}

		/* close the files */
		for(j=0;j<chrFiles.numFiles;j++) {
			fclose(chrFiles.files[j]);
		}

		/* Free memory */
		free(chrFiles.files);
		chrFiles.files = NULL;
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\n");
	}

	/* Close file */
	fclose(inputFP);
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

	/* Filter all entries with score < minSCore and do not align within bounds */
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
			free((*entries)[numEntries-1].read);
			free((*entries)[numEntries-1].reference);
			numEntries--;
			(*entries) = realloc((*entries), sizeof(AlignEntry)*numEntries);
			if(NULL == (*entries)) {
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
			if((*entries)[minIndex].score == (*entries)[i].score) {
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
				free((*entries)[i].read);
				free((*entries)[i].reference);
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
				free((*entries)[i].read);
				free((*entries)[i].reference);
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
				free((*entries)[i].read);
				free((*entries)[i].reference);
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
