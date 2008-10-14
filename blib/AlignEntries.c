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
		FILE *outputFP,
		int binaryOutput)
{
	char *FnName = "AlignEntriesPrint";
	int i;
	int32_t tempReadNameLength;

	assert(a->pairedEnd == 1 || a->numEntriesTwo == 0);

	if(binaryOutput == TextOutput) {
		/* Print the read name and paired end flag */
		if(fprintf(outputFP, "%s\t%d\t%d\t%d\t%d\n",
					a->readName,
					a->pairedEnd,
					a->space,
					a->numEntriesOne,
					a->numEntriesTwo) < 0) {
			PrintError(FnName,
					NULL,
					"Could not write to file",
					Exit,
					WriteFileError);
		}
	}
	else {
		assert(a!=NULL);
		tempReadNameLength = (int)strlen(a->readName);
		if(fwrite(&tempReadNameLength, sizeof(int32_t), 1, outputFP) != 1 ||
				fwrite(a->readName, sizeof(char), tempReadNameLength, outputFP) != tempReadNameLength ||
				fwrite(&a->pairedEnd, sizeof(int32_t), 1, outputFP) != 1 ||
				fwrite(&a->space, sizeof(int32_t), 1, outputFP) != 1 ||
				fwrite(&a->numEntriesOne, sizeof(int32_t), 1, outputFP) != 1 ||
				fwrite(&a->numEntriesTwo, sizeof(int32_t), 1, outputFP) != 1) {
			PrintError(FnName,
					NULL,
					"Could not write to file",
					Exit,
					WriteFileError);
		}
	}

	for(i=0;i<a->numEntriesOne;i++) {
		if(EOF == AlignEntryPrint(&a->entriesOne[i],
					outputFP,
					a->space,
					binaryOutput)) {
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
						a->space,
						binaryOutput)) {
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
		int space,
		int binaryInput)
{
	char *FnName = "AlignEntriesRead";
	int i;
	int32_t tempReadNameLength;

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

	/* Read the read name, paired end flag, space flag, and the number of entries for both entries */
	if(binaryInput == TextInput) {
		if(fscanf(inputFP, "%s %d %d %d %d",
					a->readName,
					&a->pairedEnd,
					&a->space,
					&a->numEntriesOne,
					&a->numEntriesTwo)==EOF) {
			/* Free read name before leaving */
			free(a->readName);
			a->readName=NULL;
			return EOF;
		}
	}
	else {
		if(fread(&tempReadNameLength, sizeof(int32_t), 1, inputFP) != 1) {
			/* Free read name before leaving */
			free(a->readName);
			a->readName=NULL;
			return EOF;
		}
		if(fread(a->readName, sizeof(char), tempReadNameLength, inputFP) != tempReadNameLength ||
				fread(&a->pairedEnd, sizeof(int32_t), 1, inputFP) != 1 ||
				fread(&a->space, sizeof(int32_t), 1, inputFP) != 1 ||
				fread(&a->numEntriesOne, sizeof(int32_t), 1, inputFP) != 1 ||
				fread(&a->numEntriesTwo, sizeof(int32_t), 1, inputFP) != 1) {
			PrintError(FnName,
					NULL,
					"Could not read from file",
					Exit,
					ReadFileError);
		}
		/* Add the null terminator */
		a->readName[tempReadNameLength]='\0';
	}

	if(a->pairedEnd == 0 && a->numEntriesTwo > 0) {
		PrintError(FnName,
				"a->pairedEnd == 0 && a->numEntriesTwo > 0",
				"Expecting single end data",
				Exit,
				OutOfRange);
	}
	if(pairedEnd != PairedEndDoesNotMatter &&
			a->pairedEnd != pairedEnd) {
		PrintError(FnName,
				"a->pairedEnd != pairedEnd",
				"Paired end does not match",
				Exit,
				OutOfRange);
	}
	if(space != SpaceDoesNotMatter &&
			a->space != space) {
		PrintError(FnName,
				"a->space != space",
				"Space does not match",
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
					a->space,
					binaryInput)) {
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
						a->space,
						binaryInput)) {
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
		int pairedEnd,
		int space)
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
	a->space = space;

	/* Don't change read name */
	/*
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
	*/

	assert(a->pairedEnd == SingleEnd || a->numEntriesTwo == PairedEnd);

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
		char *readName,
		int numEntriesOne,
		int numEntriesTwo,
		int pairedEnd,
		int space)
{
	char *FnName = "AlignEntriesAllocate";
	int i;

	a->numEntriesOne = numEntriesOne;
	a->numEntriesTwo = numEntriesTwo;
	a->pairedEnd = pairedEnd;
	a->space = space;

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
	/* Copy over */
	strcpy(a->readName, readName);

	assert(a->pairedEnd == 1 || a->numEntriesTwo == 0);

	/* Allocate memory for the entries */ 
	if(a->numEntriesOne > 0) {
		a->entriesOne = malloc(sizeof(AlignEntry)*a->numEntriesOne);
		if(NULL==a->entriesOne) {
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
	}
	else {
		a->entriesOne = NULL;
	}
	if(a->numEntriesTwo > 0) {
		a->entriesTwo = malloc(sizeof(AlignEntry)*a->numEntriesTwo);
		if(NULL==a->entriesTwo) {
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
	else {
		a->entriesTwo = NULL;
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
	a->space=NTSpace;
}

/* TODO */
void AlignEntriesKeepOnly(AlignEntries *a,
		int indexOne,
		int indexTwo,
		int pairedEnd,
		int space)
{
	assert(pairedEnd == a->pairedEnd);
	
	/* First read */
	assert(0 <= indexOne && indexOne < a->numEntriesOne);
	/* Copy to the front */
	if(indexOne > 0) {
		AlignEntryCopy(&a->entriesOne[indexOne], &a->entriesOne[0]);
	}

	/* Only for paired end */
	if(PairedEnd==pairedEnd) {
		/* Second read */
		assert(0 <= indexTwo && indexTwo < a->numEntriesTwo);
		if(indexTwo > 0) {
			AlignEntryCopy(&a->entriesTwo[indexTwo], &a->entriesTwo[0]);
		}
		/* Reallocate */
		AlignEntriesReallocate(a,
				1,
				1,
				PairedEnd,
				space);
	}
	else {
		AlignEntriesReallocate(a,
				1,
				0,
				SingleEnd,
				space);
	}
}

void AlignEntriesCopy(AlignEntries *src, AlignEntries *dst) 
{
	int i;

	/* Free and Allocate destination */
	AlignEntriesFree(dst);
	AlignEntriesAllocate(dst,
			src->readName,
			src->numEntriesOne,
			src->numEntriesTwo,
			src->pairedEnd,
			src->space);
	/* Copy over */
	for(i=0;i<src->numEntriesOne;i++) {
		AlignEntryCopy(&src->entriesOne[i], &dst->entriesOne[i]);
	}
	for(i=0;i<src->numEntriesTwo;i++) {
		AlignEntryCopy(&src->entriesTwo[i], &dst->entriesTwo[i]);
	}
}
