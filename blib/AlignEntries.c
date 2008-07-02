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
void AlignEntriesPrint(AlignEntries *aEntries,
		FILE *outputFP)
{
	char *FnName = "AlignEntriesPrint";
	int i;

	assert(aEntries->pairedEnd == 1 || aEntries->numEntriesTwo == 0);

	/* Print the read name and paired end flag */
	if(fprintf(outputFP, "%s\t%d\t%d\t%d\n",
				aEntries->readName,
				aEntries->pairedEnd,
				aEntries->numEntriesOne,
				aEntries->numEntriesTwo) < 0) {
		PrintError(FnName,
				NULL,
				"Could not write to file",
				Exit,
				WriteFileError);
	}

	for(i=0;i<aEntries->numEntriesOne;i++) {
		if(EOF == AlignEntryPrint(&aEntries->entriesOne[i],
					outputFP)) {
			PrintError(FnName,
					"entriesOne",
					"Could not write to file",
					Exit,
					WriteFileError);
		}
	}

	if(aEntries->pairedEnd==1) {
		for(i=0;i<aEntries->numEntriesTwo;i++) {
			if(EOF == AlignEntryPrint(&aEntries->entriesTwo[i],
						outputFP)) {
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
int AlignEntriesRead(AlignEntries *aEntries,
		FILE *inputFP)
{
	char *FnName = "AlignEntriesRead";
	int i;

	/* Allocate memory for the read name */
	aEntries->readName = malloc(sizeof(char)*SEQUENCE_NAME_LENGTH);
	if(aEntries->readName == NULL) {
		if(NULL == aEntries->readName) {
			PrintError(FnName,
					"aEntries->readName",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}

	/* Read the read name, paired end flag, and the number of entries for both entries */
	if(fscanf(inputFP, "%s %d %d %d",
				aEntries->readName,
				&aEntries->pairedEnd,
				&aEntries->numEntriesOne,
				&aEntries->numEntriesTwo)==EOF) {
		return EOF;
	}

	assert(aEntries->pairedEnd == 1 || aEntries->numEntriesTwo == 0);

	/* Allocate memory for the first entry */ 
	aEntries->entriesOne = malloc(sizeof(AlignEntry)*aEntries->numEntriesOne);
	if(aEntries->numEntriesOne > 0 && NULL==aEntries->entriesOne) {
		if(NULL == aEntries->entriesOne) {
			PrintError(FnName,
					"aEntries->entriesOne",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	/* Allocate memory for the second entry */
	aEntries->entriesTwo = malloc(sizeof(AlignEntry)*aEntries->numEntriesTwo);
	if(aEntries->numEntriesTwo > 0 && NULL==aEntries->entriesTwo) {
		assert(aEntries->pairedEnd == 1);
		if(NULL == aEntries->entriesTwo) {
			PrintError(FnName,
					"aEntries->entriesTwo",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}

	/* Read the alignment */
	for(i=0;i<aEntries->numEntriesOne;i++) {
		AlignEntryInitialize(&aEntries->entriesOne[i]);
		if(EOF==AlignEntryRead(&aEntries->entriesOne[i],
					inputFP)) {
			PrintError(FnName, 
					NULL, 
					"Could not read entriesOne",
					Exit,
					EndOfFile);
		}
	}
	if(aEntries->pairedEnd == 1) {
		for(i=0;i<aEntries->numEntriesTwo;i++) {
			AlignEntryInitialize(&aEntries->entriesTwo[i]);
			if(EOF==AlignEntryRead(&aEntries->entriesTwo[i],
						inputFP)) {
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
void AlignEntriesRemoveDuplicates(AlignEntries *aEntries,
		int sortOrder)
{
	/* First entry */
	aEntries->numEntriesOne = AlignEntryRemoveDuplicates(&aEntries->entriesOne,
			aEntries->numEntriesOne,
			sortOrder);
	/* Second entry */
	aEntries->numEntriesTwo = AlignEntryRemoveDuplicates(&aEntries->entriesTwo,
			aEntries->numEntriesTwo,
			sortOrder);
}

/* TODO */
/* Log-n space */
/* Do not use, since it is buggy and has not been updated lately */  
void AlignEntriesQuickSort(AlignEntries *aEntries,
		int sortOrder,
		int showPercentComplete)
{
	double percentComplete;
	/* Sort the first entry */
	percentComplete=0.0;
	AlignEntryQuickSort(&aEntries->entriesOne,
			0,
			aEntries->numEntriesOne-1,
			sortOrder,
			showPercentComplete,
			&percentComplete,
			aEntries->numEntriesOne-1);
	/* Sort the second entry */
	percentComplete=0.0;
	AlignEntryQuickSort(&aEntries->entriesTwo,
			0,
			aEntries->numEntriesTwo-1,
			sortOrder,
			showPercentComplete,
			&percentComplete,
			aEntries->numEntriesTwo-1);
}

/* TODO */
/* O(n) space, but really double */
void AlignEntriesMergeSort(AlignEntries *aEntries,
		int sortOrder,
		int showPercentComplete)
{
	double percentComplete;
	/* Sort the first entry */
	percentComplete=0.0;
	AlignEntryMergeSort(&aEntries->entriesOne,
			0,
			aEntries->numEntriesOne-1,
			sortOrder,
			showPercentComplete,
			&percentComplete,
			aEntries->numEntriesOne-1);
	/* Sort the second entry */
	percentComplete=0.0;
	AlignEntryMergeSort(&aEntries->entriesTwo,
			0,
			aEntries->numEntriesTwo-1,
			sortOrder,
			showPercentComplete,
			&percentComplete,
			aEntries->numEntriesTwo-1);
}

void AlignEntriesReallocate(AlignEntries *aEntries,
		int numEntriesOne,
		int numEntriesTwo,
		int pairedEnd)
{
	char *FnName = "AlignEntriesReallocate";
	int i;

	aEntries->numEntriesOne = numEntriesOne;
	aEntries->numEntriesTwo = numEntriesTwo;
	aEntries->pairedEnd = pairedEnd;

	aEntries->readName = realloc(aEntries->readName, sizeof(char)*SEQUENCE_NAME_LENGTH);
	if(aEntries->readName == NULL) {
		if(NULL == aEntries->readName) {
			PrintError(FnName,
					"aEntries->readName",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}

	assert(aEntries->pairedEnd == 1 || aEntries->numEntriesTwo == 0);

	/* Allocate memory for the entries */ 
	aEntries->entriesOne = realloc(aEntries->entriesOne, sizeof(AlignEntry)*aEntries->numEntriesOne);
	if(aEntries->numEntriesOne > 0 && NULL==aEntries->entriesOne) {
		if(NULL == aEntries->entriesOne) {
			PrintError(FnName,
					"aEntries->entriesOne",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	/* Initialize */
	for(i=0;i<aEntries->numEntriesOne;i++) {
		AlignEntryInitialize(&aEntries->entriesOne[i]);
	}
	aEntries->entriesTwo = realloc(aEntries->entriesTwo, sizeof(AlignEntry)*aEntries->numEntriesTwo);
	if(aEntries->numEntriesTwo > 0 && NULL==aEntries->entriesTwo) {
		assert(aEntries->pairedEnd == 1);
		if(NULL == aEntries->entriesTwo) {
			PrintError(FnName,
					"aEntries->entriesTwo",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	/* Initialize */
	for(i=0;i<aEntries->numEntriesTwo;i++) {
		AlignEntryInitialize(&aEntries->entriesTwo[i]);
	}

}

void AlignEntriesAllocate(AlignEntries *aEntries,
		int numEntriesOne,
		int numEntriesTwo,
		int pairedEnd)
{
	char *FnName = "AlignEntriesAllocate";
	int i;

	aEntries->numEntriesOne = numEntriesOne;
	aEntries->numEntriesTwo = numEntriesTwo;
	aEntries->pairedEnd = pairedEnd;

	aEntries->readName = malloc(sizeof(char)*SEQUENCE_NAME_LENGTH);
	if(aEntries->readName == NULL) {
		if(NULL == aEntries->readName) {
			PrintError(FnName,
					"aEntries->readName",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}

	assert(aEntries->pairedEnd == 1 || aEntries->numEntriesTwo == 0);

	/* Allocate memory for the entries */ 
	aEntries->entriesOne = malloc(sizeof(AlignEntry)*aEntries->numEntriesOne);
	if(aEntries->numEntriesOne > 0 && NULL==aEntries->entriesOne) {
		if(NULL == aEntries->entriesOne) {
			PrintError(FnName,
					"aEntries->entriesOne",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	/* Initialize */
	for(i=0;i<aEntries->numEntriesOne;i++) {
		AlignEntryInitialize(&aEntries->entriesOne[i]);
	}
	aEntries->entriesTwo = malloc(sizeof(AlignEntry)*aEntries->numEntriesTwo);
	if(aEntries->numEntriesTwo > 0 && NULL==aEntries->entriesTwo) {
		assert(aEntries->pairedEnd == 1);
		if(NULL == aEntries->entriesTwo) {
			PrintError(FnName,
					"aEntries->entriesTwo",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	/* Initialize */
	for(i=0;i<aEntries->numEntriesTwo;i++) {
		AlignEntryInitialize(&aEntries->entriesTwo[i]);
	}
}

void AlignEntriesFree(AlignEntries *aEntries)
{
	int i;
	for(i=0;i<aEntries->numEntriesOne;i++) {
		AlignEntryFree(&aEntries->entriesOne[i]);
	}
	for(i=0;i<aEntries->numEntriesTwo;i++) {
		AlignEntryFree(&aEntries->entriesTwo[i]);
	}
	free(aEntries->readName);
	free(aEntries->entriesOne);
	free(aEntries->entriesTwo);
	AlignEntriesInitialize(aEntries);
}

void AlignEntriesInitialize(AlignEntries *aEntries) 
{
	aEntries->readName=NULL;
	aEntries->entriesOne=NULL;
	aEntries->numEntriesOne=0;
	aEntries->entriesTwo=NULL;
	aEntries->numEntriesTwo=0;
	aEntries->pairedEnd=0;
}
