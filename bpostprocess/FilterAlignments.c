#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
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
	int i;
	FILE **filteredFPs=NULL;
	FILE *fp;
	int numChrs = endChr - startChr + 1;
	AlignEntry *entries=NULL;
	int numEntries=0;

	/* Allocate memory for the tmp files */
	filteredFPs = malloc(sizeof(FILE*)*numChrs);
	if(NULL == filteredFPs) {
		PrintError(FnName,
				"filteredFPs",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Initialize */
	for(i=0;i<numChrs;i++) {
		filteredFPs[i]=tmpfile();
	}

	/* Open the file */
	if(0==(fp=fopen(inputFileName, "r"))) {
		PrintError(FnName,
				inputFileName,
				"Could not open inputFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Read in all alignments for one read */
	switch(inputFormat) {
		/* TODO : separate each case into its own function */
		case BAlignFile:
			for(numEntries = AlignEntryGetOneRead(&entries, fp);
					numEntries > 0;
					numEntries = AlignEntryGetOneRead(&entries, fp)) {
				/* Apply filtering */
				numEntries = FilterEntries(&entries, 
						numEntries,
						uniqueMatches,
						bestScore,
						minScore,
						startChr,
						startPos,
						endChr,
						endPos);

				/* Output to correct tmp file */
				assert(numEntries == 1 || numEntries == 0);
				if(numEntries == 1) {
					/* Print entry to tmp file */
					AlignEntryPrint(&entries[0], filteredFPs[entries[0].chromosome - startChr]);
					/* Free entries */
					free(entries[0].read);
					free(entries[0].reference);
					free(entries);
					entries=NULL;
				}
			}
			/* Move filtered files to the start */
			for(i=0;i<numChrs;i++) {
				fseek(filteredFPs[i], 0, SEEK_SET);
			}
			/* Sort each file by chr pos */
			for(i=0;i<numChrs;i++) {
				/* Get all the entries from the temp file */
				numEntries=AlignEntryGetAll(&entries, filteredFPs[i]);
				/* Sort the entries */
				AlignEntryQuickSort(&entries, 0, numEntries-1, AlignEntrySortByChrPos);
				/* Close the temp file */
				fclose(filteredFPs[i]);
				filteredFPs[i] = NULL;
				/* Print the entries to file */
				ConvertAndPrint((void*)entries,
						numEntries,
						inputFormat,
						i+startChr,
						outputID,
						outputDir,
						outputFormat);
				/* Free the entries */
			}
			break;
		default:
			PrintError(FnName,
					"inputFormat",
					"Input format not supported",
					Exit,
					OutOfRange);
	}

	/* Close file */
	fclose(fp);

	/* Free the filtered FPs */
	for(i=0;i<numChrs;i++) {
		free(filteredFPs[i]);
	}
	free(filteredFPs);

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
		else if(minIndex == -1 || (*entries)[minIndex].score < (*entries)[i].score) {
			/* Get the entry with the best score */
			minIndex = i;
			numMinIndexes = 1;
		}
		else if((*entries)[minIndex].score < (*entries)[i].score) {
			/* We should store if there are more than one place with the same score */
			numMinIndexes++;
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
			return 0;
		}
	}
	else if(bestScore == 1) {
		if(numMinIndexes > 1) {
			/* We have multiple entries with the same best score */
			/* Free all memory and return 0 */
			/* Free memory and reallocate */
			for(i=0;i<numEntries;i++) {
				free((*entries)[i].read);
				free((*entries)[i].reference);
			}
			free((*entries));
			return 0;
		}
		else {
			/* Should have only one entry with the best score */
			assert(minIndex >=0 && numMinIndexes == 1);
			/* Copy min index to the front */
			AlignEntryCopy(&(*entries)[0], &(*entries)[minIndex]);
			/* Free memory and reallocate */
			for(i=1;i<numEntries;i++) {
				free((*entries)[i].read);
				free((*entries)[i].reference);
			}
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
