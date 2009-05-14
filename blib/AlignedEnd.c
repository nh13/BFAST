#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <float.h>

#include "BLibDefinitions.h"
#include "ScoringMatrix.h"
#include "BError.h"
#include "BLib.h"
#include "AlignedEntry.h"
#include "AlignedEnd.h"

/* TODO */
int32_t AlignedEndPrint(AlignedEnd *a,
		FILE *outputFP,
		int32_t space,
		int32_t binaryOutput)
{
	int32_t i;
	assert(NULL != a->read);
	if(binaryOutput == TextOutput) {
		if(fprintf(outputFP, "%s\t%s\t%d\n",
					a->read,
					a->qual,
					a->numEntries) < 0) {
			return EOF;
		}
	}
	else {
		if(fwrite(&a->readLength, sizeof(int32_t), 1, outputFP) != 1 ||
				fwrite(&a->qualLength, sizeof(int32_t), 1, outputFP) != 1 ||
				fwrite(a->read, sizeof(char), a->readLength, outputFP) != a->readLength ||
				fwrite(a->qual, sizeof(char), a->qualLength, outputFP) != a->qualLength ||
				fwrite(&a->numEntries, sizeof(int32_t), 1, outputFP) != 1) {
			return EOF;
		}
	}

	for(i=0;i<a->numEntries;i++) {
		if(EOF == AlignedEntryPrint(&a->entries[i],
					outputFP,
					space,
					binaryOutput)) {
			return EOF;
		}
	}

	return 1;
}

/* TODO */
int32_t AlignedEndRead(AlignedEnd *a,
		FILE *inputFP,
		int32_t space,
		int32_t binaryInput)
{
	char *FnName = "AlignedEndRead";
	char read[SEQUENCE_LENGTH]="\0";
	char qual[SEQUENCE_LENGTH]="\0";
	int32_t i;

	if(binaryInput == TextInput) {
		if(fscanf(inputFP, "%s %s %d",
					read,
					qual,
					&a->numEntries) < 3) {
			return EOF;
		}

		a->readLength = strlen(read);
		a->qualLength = strlen(qual);

		/* Allocate memory for the alignment */
		if(a->read == NULL) {
			a->read = malloc(sizeof(char)*(1+a->readLength));
			if(NULL == a->read) {
				PrintError(FnName,
						"a->read",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}
		if(a->qual == NULL) {
			a->qual = malloc(sizeof(char)*(1+a->qualLength));
			if(NULL == a->qual) {
				PrintError(FnName,
						"a->qual",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}
		/* Copy over */
		strcpy(a->read, read);
		strcpy(a->qual, qual);
	}
	else {
		if(fread(&a->readLength, sizeof(int32_t), 1, inputFP) != 1 ||
				fread(&a->qualLength, sizeof(int32_t), 1, inputFP) != 1) {
			return EOF;
		}
		/* Allocate memory for the alignment */
		if(a->read == NULL) {
			a->read = malloc(sizeof(char)*(1+a->readLength));
			if(NULL == a->read) {
				PrintError(FnName,
						"a->read",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}
		if(a->qual == NULL) {
			a->qual = malloc(sizeof(char)*(1+a->qualLength));
			if(NULL == a->qual) {
				PrintError(FnName,
						"a->qual",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}

		if(fread(a->read, sizeof(char), a->readLength, inputFP) != a->readLength ||
				fread(a->qual, sizeof(char), a->qualLength, inputFP) != a->qualLength ||
				fread(&a->numEntries, sizeof(int32_t), 1, inputFP) != 1) {
			PrintError(FnName,
					"a->reads, a->qual, and a->numEntries",
					"Could not read from file",
					Exit,
					ReadFileError);
		}
		/* Add the null terminator to strings */
		a->read[a->readLength]='\0';
		a->qual[a->qualLength]='\0';
	}

	/* Allocate room for the the entries */
	AlignedEndReallocate(a,
			a->numEntries);

	for(i=0;i<a->numEntries;i++) {
		AlignedEntryInitialize(&a->entries[i]);
		if(EOF == AlignedEntryRead(&a->entries[i],
					inputFP,
					space,
					binaryInput)) {
			PrintError(FnName,
					"a->entries[i]",
					"Could not read from file",
					Exit,
					ReadFileError);
		}
	}

	return 1;
}

/* TODO */
int32_t AlignedEndRemoveDuplicates(AlignedEnd *end,
		int32_t sortOrder)
{
	/*
	   char *FnName="AlignedEndRemoveDuplicates";
	   */
	int32_t i, prevIndex;

	if(end->numEntries > 1) {
		/* Sort the entries */
		AlignedEndQuickSort(end, sortOrder, 0);
		/*
		   AlignedEndMergeSort(end, sortOrder, 0);
		   */

		/* Remove duplicates */
		prevIndex=0;
		for(i=1;i<end->numEntries;i++) {
			if(AlignedEntryCompareAtIndex(end->entries, prevIndex, end->entries, i, sortOrder)==0) {
				/* Do nothing */
			}
			else {
				/* Increment prevIndex */
				prevIndex++;
				/* Copy to prevIndex (incremented) */
				AlignedEntryCopyAtIndex(end->entries, prevIndex, end->entries, i);
			}
		}

		/* Reallocate */
		AlignedEndReallocate(end, prevIndex+1);
	}
	return end->numEntries;
}

/* TODO */
/* Log-n space */
/* Do not use, since it is buggy and has not been updated lately */  
void AlignedEndQuickSort(AlignedEnd *a,
		int32_t sortOrder,
		int32_t showPercentComplete)
{

	double curPercent = 0.0;
	AlignedEntryQuickSort(&a->entries,
			0,
			a->numEntries-1,
			sortOrder,
			showPercentComplete,
			&curPercent,
			a->numEntries);
}

/* TODO */
/* O(n) space, but really double */
void AlignedEndMergeSort(AlignedEnd *a,
		int32_t sortOrder,
		int32_t showPercentComplete)
{
	double curPercent = 0.0;
	AlignedEntryMergeSort(&a->entries,
			0,
			a->numEntries-1,
			sortOrder,
			showPercentComplete,
			&curPercent,
			a->numEntries);
}

/* TODO */
int32_t AlignedEndCompare(AlignedEnd *a,
		AlignedEnd *b, 
		int32_t sortOrder)
{
	assert(1 == a->numEntries);
	assert(1 == b->numEntries);

	return AlignedEntryCompare(&a->entries[0],
			&b->entries[0],
			sortOrder);
}

/* TODO */
void AlignedEndCopyAtIndex(AlignedEnd *dest, int32_t destIndex, AlignedEnd *src, int32_t srcIndex)
{
	if(dest != src || srcIndex != destIndex) {
		AlignedEndCopy(&(dest[destIndex]), &(src[srcIndex]));
	}
}

/* TODO */
void AlignedEndCopy(AlignedEnd *dest, AlignedEnd *src)
{
	char *FnName = "AlignedEndCopy";
	int32_t i;
	if(src != dest) {
		assert(src->read != NULL);
		/* read */
		dest->readLength = src->readLength;
		dest->read = realloc(dest->read, sizeof(char)*(src->readLength+1));
		if(NULL == dest->read) {
			PrintError(FnName,
					"dest->read",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		assert(src->read != NULL);
		strcpy(dest->read, src->read);
		/* qual */
		dest->qualLength = src->qualLength;
		dest->qual = realloc(dest->qual, sizeof(char)*(src->qualLength+1));
		if(NULL == dest->qual) {
			PrintError(FnName,
					"dest->qual",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		assert(src->qual != NULL);
		strcpy(dest->qual, src->qual);
		/* Reallocate */
		AlignedEndReallocate(dest,
				src->numEntries);
		/* Copy entries */
		for(i=0;i<dest->numEntries;i++) {
			AlignedEntryCopy(&dest->entries[i], 
					&src->entries[i]);
		}
	}
}

void AlignedEndAllocate(AlignedEnd *a,
		char *read,
		char *qual,
		int32_t numEntries)
{
	char *FnName="AlignedEndAllocate";
	int32_t i;

	/* Allocate */
	assert(NULL != read);
	a->readLength = strlen(read);
	a->read = malloc(sizeof(char)*(1+a->readLength));
	if(NULL == a->read) {
		PrintError(FnName,
				"a->read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(NULL != qual);
	a->qualLength = strlen(qual);
	a->qual = malloc(sizeof(char)*(1+a->qualLength));
	if(NULL == a->qual) {
		PrintError(FnName,
				"a->qual",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	a->numEntries = numEntries;
	a->entries = malloc(sizeof(AlignedEntry)*a->numEntries);
	if(NULL == a->entries && 0 < numEntries) {
		PrintError(FnName,
				"a->entries",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}

	/* Copy over */
	strcpy(a->read, read);
	strcpy(a->qual, qual);

	/* Initialize */
	for(i=0;i<a->numEntries;i++) {
		AlignedEntryInitialize(&a->entries[i]);
	}

}

void AlignedEndReallocate(AlignedEnd *a,
		int32_t numEntries)
{
	char *FnName="AlignedEndReallocate";
	int32_t i;

	if(numEntries < a->numEntries) {
		for(i=numEntries;i<a->numEntries;i++) {
			AlignedEntryFree(&a->entries[i]);
		}
	}

	/* Reallocate */
	a->entries = realloc(a->entries, sizeof(AlignedEntry)*numEntries);
	if(NULL == a->entries && 0 < numEntries) {
		PrintError(FnName,
				"a->entries",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}

	for(i=a->numEntries;i<numEntries;i++) {
		AlignedEntryInitialize(&a->entries[i]);
	}

	a->numEntries = numEntries;

}

void AlignedEndFree(AlignedEnd *a)
{
	int32_t i;

	free(a->read);
	free(a->qual);
	for(i=0;i<a->numEntries;i++) {
		AlignedEntryFree(&a->entries[i]);
	}
	free(a->entries);
	AlignedEndInitialize(a);
}

void AlignedEndInitialize(AlignedEnd *a) 
{
	a->read=NULL;
	a->readLength=0;
	a->qual=NULL;
	a->qualLength=0;
	a->numEntries=0;
	a->entries=NULL;
}

void AlignedEndUpdateMappingQuality(AlignedEnd *a,
		double mismatchScore,
		int avgMismatchQuality)
{
	double bestScore=INT_MIN, nextBestScore=INT_MIN;
	double bestMappingQuality = 0.0;
	int32_t i;

	/* Get best and next best score */
	for(i=0;i<a->numEntries;i++) {
		if(bestScore < a->entries[i].score) {
			bestScore = a->entries[i].score;
		}
		else if(nextBestScore < a->entries[i].score) {
			nextBestScore = a->entries[i].score;
		}
	}

	if(bestScore < 0) {
		bestScore = 0;
	}
	if(nextBestScore < 0) {
		nextBestScore = 0;
	}
	assert(nextBestScore <= bestScore);

	bestMappingQuality = ( (bestScore - nextBestScore)/mismatchScore )*avgMismatchQuality;

	for(i=0;i<a->numEntries;i++) {
		if(a->entries[i].score < bestScore) {
			a->entries[i].mappingQuality = 0;
		}
		else {
			a->entries[i].mappingQuality = bestMappingQuality;
		}
	}
}
