#include <stdlib.h>
#include <stdio.h>
#include "../blib/AlignEntries.h"
#include "Filter.h"

/* TODO */
int FilterAlignEntries(AlignEntries *a,
		int algorithmReads,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int pairedEnd,
		int algorithmReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int meanDistancePaired,
		int chrAbPaired,
		int inversionsPaired,
		TmpFP **tmpFPs,
		int numTmpFPs,
		FILE *fpChrAb,
		FILE *fpInversions,
		FILE *fpNotReported) 
{
	assert(a->pairedEnd == pairedEnd);
	int printedPaired=0;

	/* Filter for reads that are not paired end */
	if(pairedEnd == 0) {
		FilterMultipleAlignEntry(&a->entriesOne,
				&numEntriesOne,
				algorithmReads,
				minScoreReads,
				startChr,
				startPos,
				endChr,
				endPos);
	}
	else {
		/* We should return whether we printed as an chrAb or inversion */
		printedPaired=FilterOneAlignEntries(a,
				algorithmReads,
				minScoreReads,
				startChr,
				startPos,
				endChr,
				endPos,
				algorithmReadsPaired,
				minScoreReadsPaired,
				minDistancePaired,
				maxDistancePaired,
				meanDistancePaired,
				chrAbPaired,
				inversionsPaired,
				fpChrAb,
				fpInversions);
	}
	/* TODO */
}

/* TODO */
void FilterMultipleAlignEntry(AlignEntry **a,
		int *numEntries,
		int algorithmReads,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos)
{
	char *FnName = "FilterMultipleAlignEntry";
	int i;
	int freeEntries=0;
	int maxScore = INT_MAX;
	double numMaxScore = 0;
	int maxScoreIndex = -1;

	/* Filter each entry for this read */
	for(i=0;i<(*numEntries) && 0==freeEntries;) {
		/* Check if we should filter */
		if(1==FilterAlignEntry(&(*a)[i],
					minScoreReads,
					startChr,
					startPos,
					endChr,
					endPos)) {
			/* Copy the one on the end to here */
			if(i<(*numEntries)-1) { /* Make sure we are not at the end */
				AlignEntryCopy(&(*a)[(*numEntries)-1], &(*a)[i]); 
			}
			/* Free memory of the entry at the end */
			AlignEntryFree(&(*a)[(*numEntries)-1]);
			/* Reallocate the number of entries */
			(*numEntries)--; /* Decrement by one */
			(*a) = realloc((*a), sizeof(AlignEntry)*(*numEntries));
			if(NULL == (*a) && (*numEntries) > 0) {
				PrintError(FnName,
						"(*a)",
						"Could not reallocate memory",
						Exit,
						ReallocMemory);
			}
			/* Do not increment i */
		}
		else if(0==algorithmReads && i > 0) {
			/* This mean that two ore more entries passed the filters.  Therefore
			 * we will not have a unique alignment for this read.  Free memory and 
			 * return */
			freeEntries = 1;
			/* Do not increment i and exit the loop */
		}
		else {
			if(1==algorithmReads) {
				if((*a)[i].score == maxScore) {
					assert(numMaxScore > 0);
					numMaxScore++;
				}
				else if((*a)[i].score > maxScore) {
					maxScore = (*a)[i].score;
					numMaxScore = 1;
					maxScoreIndex = i;
				}
			}
			/* Move to the next */
			i++;
		}
	}
	/* Check if there is a best score (unique best score) */
	if(1==algorithmReads) {
		if(numMaxScore == 1) {
			if(maxScoreIndex > 0) {
				/* Copy from index to front */
				AlignEntryCopy(&(*a)[maxScoreIndex], &(*a)[0]); 
			}
			/* Free memory of the entry at the other entries */
			for(i=1;i<(*numEntries);i++) {
				AlignEntryFree(&(*a)[i]);
			}
			/* Reallocate the number of entries */
			(*numEntries)=1; /* Decrement by one */
			(*a) = realloc((*a), sizeof(AlignEntry)*(*numEntries));
			if(NULL == (*a) && (*numEntries) > 0) {
				PrintError(FnName,
						"(*a)",
						"Could not reallocate memory",
						Exit,
						ReallocMemory);
			}
		}
		else {
			freeEntries = 1;
		}
	}
	/* Check if we should free the entries */
	if(1==freeEntries) {
		/* Free all the entries */
		for(i=0;i<(*numEntries);i++) {
			AlignEntryFree(&(*a)[i]);
		}
		free((*a));
		(*a) = NULL;
		/* Set the number of entries to zero */
		(*numEntries)=0; 
	}
}

/* TODO */
int FilterAlignEntry(AlignEntry *a,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos)
{
	/* Check if the alignment has at least the minimum score */
	if(a->score < minScoreReads) {
		return 1;
	}
	/* Check the genomic location */
	if(a->chromosome < startChr || 
			(a->chromosome == startChr && a->position < startPos) ||
			(a->chromosome == endChr && a->position > endPos) ||
			(a->chromosome > endChr)) {
		return 1;
	}
	return 0;
}

/* TODO */
int FilterOneAlignEntries(AlignEntries *a,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int algorithmReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int meanDistancePaired,
		int chrAbPaired,
		int inversionsPaired,
		FILE *fpChrAb,
		FILE *fpInversions,
		FILE *fpNotReported) 
{
	int maxScoreOne, maxScoreTwo;
	int maxScoreOutsideOne, maxScoreOutsideTwo;
	int uniqueOne, uniqueTwo;
	int uniqueOutsideOne, uniqueOutsideTwo;
	int freeMemory;
	int i, j;

	/* Go through all possible pairs of alignments */
	for(i=0, freeMemory=0;i<a->numEntriesOne && 0==freeMemory;i++) { /* First read */
		for(j=0;j<a->numEntriesTwo && 0==freeMemory;j++) { /* Second read */
		}
	}
}
