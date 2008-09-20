#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "BError.h"
#include "BLib.h"
#include "AlignEntry.h"
#include "AlignEntries.h"

/* TODO */
void AlignEntriesPrint(AlignEntries *a,
		FILE *outputFP)
{
	char *FnName = "AlignEntriesPrint";
	int i;

	assert(a->pairedEnd == 1 || a->numEntriesTwo == 0);

	/* Print the read name and paired end flag */
	if(fprintf(outputFP, "%s\t%d\t%d\t%d\t%d\n",
				a->readName,
				a->pairedEnd,
				a->colorSpace,
				a->numEntriesOne,
				a->numEntriesTwo) < 0) {
		PrintError(FnName,
				NULL,
				"Could not write to file",
				Exit,
				WriteFileError);
	}

	for(i=0;i<a->numEntriesOne;i++) {
		if(EOF == AlignEntryPrint(&a->entriesOne[i],
					outputFP,
					a->colorSpace)) {
			PrintError(FnName,
					"entriesOne",
					"Could not write to file",
					Exit,
					WriteFileError);
		}
	}

	if(a->pairedEnd==1) {
		for(i=0;i<a->numEntriesTwo;i++) {
			if(EOF == AlignEntryPrint(&a->entriesTwo[i],
						outputFP,
						a->colorSpace)) {
				PrintError(FnName,
						"entriesTwo",
						"Could not write to file",
						Exit,
						WriteFileError);
			}
		}
	}
}

/* TODO */
int AlignEntriesRead(AlignEntries *a,
		FILE *inputFP,
		int pairedEnd,
		int colorSpace)
{
	char *FnName = "AlignEntriesRead";
	int i;

	/* Allocate memory for the read name */
	a->readName = malloc(sizeof(char)*SEQUENCE_NAME_LENGTH);
	if(a->readName == NULL) {
		if(NULL == a->readName) {
			PrintError(FnName,
					"a->readName",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}

	/* Read the read name, paired end flag, color space flag, and the number of entries for both entries */
	if(fscanf(inputFP, "%s %d %d %d %d",
				a->readName,
				&a->pairedEnd,
				&a->colorSpace,
				&a->numEntriesOne,
				&a->numEntriesTwo)==EOF) {
		return EOF;
	}

	assert(a->pairedEnd == 1 || a->numEntriesTwo == 0);
	if(a->pairedEnd != pairedEnd) {
		PrintError(FnName,
				"a->pairedEnd != pairedEnd",
				"Paired end does not match",
				Exit,
				OutOfRange);
	}
	if(colorSpace != SpaceDoesNotMatter &&
			a->colorSpace != colorSpace) {
		PrintError(FnName,
				"a->colorSpace != colorSpace",
				"Color space does not match",
				Exit,
				OutOfRange);
	}

	/* Allocate memory for the first entry */ 
	a->entriesOne = malloc(sizeof(AlignEntry)*a->numEntriesOne);
	if(a->numEntriesOne > 0 && NULL==a->entriesOne) {
		if(NULL == a->entriesOne) {
			PrintError(FnName,
					"a->entriesOne",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	/* Allocate memory for the second entry */
	a->entriesTwo = malloc(sizeof(AlignEntry)*a->numEntriesTwo);
	if(a->numEntriesTwo > 0 && NULL==a->entriesTwo) {
		assert(a->pairedEnd == 1);
		if(NULL == a->entriesTwo) {
			PrintError(FnName,
					"a->entriesTwo",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}

	/* Read the alignment */
	for(i=0;i<a->numEntriesOne;i++) {
		AlignEntryInitialize(&a->entriesOne[i]);
		if(EOF==AlignEntryRead(&a->entriesOne[i],
					inputFP,
					colorSpace)) {
			PrintError(FnName, 
					NULL, 
					"Could not read entriesOne",
					Exit,
					EndOfFile);
		}
	}
	if(a->pairedEnd == 1) {
		for(i=0;i<a->numEntriesTwo;i++) {
			AlignEntryInitialize(&a->entriesTwo[i]);
			if(EOF==AlignEntryRead(&a->entriesTwo[i],
						inputFP,
						colorSpace)) {
				PrintError(FnName, 
						NULL, 
						"Could not read entriesTwo",
						Exit,
						EndOfFile);
			}
		}
	}

	return 1;
}

/* TODO */
void AlignEntriesRemoveDuplicates(AlignEntries *a,
		int sortOrder)
{
	/* First entry */
	a->numEntriesOne = AlignEntryRemoveDuplicates(&a->entriesOne,
			a->numEntriesOne,
			sortOrder);
	/* Second entry */
	a->numEntriesTwo = AlignEntryRemoveDuplicates(&a->entriesTwo,
			a->numEntriesTwo,
			sortOrder);
}

/* TODO */
/* Log-n space */
/* Do not use, since it is buggy and has not been updated lately */  
void AlignEntriesQuickSort(AlignEntries *a,
		int sortOrder,
		int showPercentComplete)
{
	double percentComplete;
	/* Sort the first entry */
	percentComplete=0.0;
	AlignEntryQuickSort(&a->entriesOne,
			0,
			a->numEntriesOne-1,
			sortOrder,
			showPercentComplete,
			&percentComplete,
			a->numEntriesOne-1);
	/* Sort the second entry */
	percentComplete=0.0;
	AlignEntryQuickSort(&a->entriesTwo,
			0,
			a->numEntriesTwo-1,
			sortOrder,
			showPercentComplete,
			&percentComplete,
			a->numEntriesTwo-1);
}

/* TODO */
/* O(n) space, but really double */
void AlignEntriesMergeSort(AlignEntries *a,
		int sortOrder,
		int showPercentComplete)
{
	double percentComplete;
	/* Sort the first entry */
	percentComplete=0.0;
	AlignEntryMergeSort(&a->entriesOne,
			0,
			a->numEntriesOne-1,
			sortOrder,
			showPercentComplete,
			&percentComplete,
			a->numEntriesOne-1);
	/* Sort the second entry */
	percentComplete=0.0;
	AlignEntryMergeSort(&a->entriesTwo,
			0,
			a->numEntriesTwo-1,
			sortOrder,
			showPercentComplete,
			&percentComplete,
			a->numEntriesTwo-1);
}

/* TODO */
void AlignEntriesReallocate(AlignEntries *a,
		int numEntriesOne,
		int numEntriesTwo,
		int pairedEnd)
{
	char *FnName = "AlignEntriesReallocate";
	int i;

	/* we have to free if we are reducing the number of entries */
	if(numEntriesOne < a->numEntriesOne) {
		for(i=numEntriesOne;i<a->numEntriesOne;i++) {
			AlignEntryFree(&a->entriesOne[i]);
		}
	}
	if(numEntriesTwo < a->numEntriesTwo) {
		for(i=numEntriesTwo;i<a->numEntriesTwo;i++) {
			AlignEntryFree(&a->entriesTwo[i]);
		}
	}

	a->numEntriesOne = numEntriesOne;
	a->numEntriesTwo = numEntriesTwo;
	a->pairedEnd = pairedEnd;

	a->readName = realloc(a->readName, sizeof(char)*SEQUENCE_NAME_LENGTH);
	if(a->readName == NULL) {
		if(NULL == a->readName) {
			PrintError(FnName,
					"a->readName",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}

	assert(a->pairedEnd == 1 || a->numEntriesTwo == 0);

	/* Allocate memory for the entries */ 
	a->entriesOne = realloc(a->entriesOne, sizeof(AlignEntry)*a->numEntriesOne);
	if(a->numEntriesOne > 0 && NULL==a->entriesOne) {
		if(NULL == a->entriesOne) {
			PrintError(FnName,
					"a->entriesOne",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	a->entriesTwo = realloc(a->entriesTwo, sizeof(AlignEntry)*a->numEntriesTwo);
	if(a->numEntriesTwo > 0 && NULL==a->entriesTwo) {
		assert(a->pairedEnd == 1);
		if(NULL == a->entriesTwo) {
			PrintError(FnName,
					"a->entriesTwo",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
}

/* TODO */
void AlignEntriesAllocate(AlignEntries *a,
		int numEntriesOne,
		int numEntriesTwo,
		int pairedEnd)
{
	char *FnName = "AlignEntriesAllocate";
	int i;

	a->numEntriesOne = numEntriesOne;
	a->numEntriesTwo = numEntriesTwo;
	a->pairedEnd = pairedEnd;

	a->readName = malloc(sizeof(char)*SEQUENCE_NAME_LENGTH);
	if(a->readName == NULL) {
		if(NULL == a->readName) {
			PrintError(FnName,
					"a->readName",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}

	assert(a->pairedEnd == 1 || a->numEntriesTwo == 0);

	/* Allocate memory for the entries */ 
	a->entriesOne = malloc(sizeof(AlignEntry)*a->numEntriesOne);
	if(a->numEntriesOne > 0 && NULL==a->entriesOne) {
		if(NULL == a->entriesOne) {
			PrintError(FnName,
					"a->entriesOne",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	/* Initialize */
	for(i=0;i<a->numEntriesOne;i++) {
		AlignEntryInitialize(&a->entriesOne[i]);
	}
	a->entriesTwo = malloc(sizeof(AlignEntry)*a->numEntriesTwo);
	if(a->numEntriesTwo > 0 && NULL==a->entriesTwo) {
		assert(a->pairedEnd == 1);
		if(NULL == a->entriesTwo) {
			PrintError(FnName,
					"a->entriesTwo",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	/* Initialize */
	for(i=0;i<a->numEntriesTwo;i++) {
		AlignEntryInitialize(&a->entriesTwo[i]);
	}
}

/* TODO */
void AlignEntriesFree(AlignEntries *a)
{
	int i;
	for(i=0;i<a->numEntriesOne;i++) {
		AlignEntryFree(&a->entriesOne[i]);
	}
	for(i=0;i<a->numEntriesTwo;i++) {
		AlignEntryFree(&a->entriesTwo[i]);
	}
	free(a->readName);
	free(a->entriesOne);
	free(a->entriesTwo);
	AlignEntriesInitialize(a);
}

/* TODO */
void AlignEntriesInitialize(AlignEntries *a) 
{
	a->readName=NULL;
	a->entriesOne=NULL;
	a->numEntriesOne=0;
	a->entriesTwo=NULL;
	a->numEntriesTwo=0;
	a->pairedEnd=0;
}

/* TODO */
void AlignEntriesKeepOnly(AlignEntries *a,
		int indexOne,
		int indexTwo,
		int pairedEnd)
{
	assert(pairedEnd == a->pairedEnd);
	/* First read */
	assert(0 <= indexOne && indexOne < a->numEntriesOne);
	/* Copy to the front */
	if(indexOne > 0) {
		AlignEntryCopy(&a->entriesOne[indexOne], &a->entriesOne[0]);
	}

	/* Only for paired end */
	if(1==pairedEnd) {
		assert(0 <= indexTwo && indexTwo < a->numEntriesTwo);
		if(indexTwo > 0) {
			AlignEntryCopy(&a->entriesTwo[indexTwo], &a->entriesTwo[0]);
		}
		/* Reallocate */
		AlignEntriesReallocate(a,
				1,
				1,
				1);
	}
	else {
		AlignEntriesReallocate(a,
				1,
				0,
				0);
	}
}

void AlignEntriesCopy(AlignEntries *src, AlignEntries *dst) 
{
	int i;

	/* Free and Allocate destination */
	AlignEntriesFree(dst);
	AlignEntriesAllocate(dst,
			src->numEntriesOne,
			src->numEntriesTwo,
			src->pairedEnd);
	/* Copy over */
	for(i=0;i<src->numEntriesOne;i++) {
		AlignEntryCopy(&src->entriesOne[i], &dst->entriesOne[i]);
	}
	for(i=0;i<src->numEntriesTwo;i++) {
		AlignEntryCopy(&src->entriesTwo[i], &dst->entriesTwo[i]);
	}
}
