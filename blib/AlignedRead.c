#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "BError.h"
#include "BLib.h"
#include "AlignedEnd.h"
#include "AlignedRead.h"

/* TODO */
void AlignedReadPrint(AlignedRead *a,
		FILE *outputFP,
		int32_t binaryOutput)
{
	char *FnName = "AlignedReadPrint";
	int32_t i;

	if(binaryOutput == TextOutput) {
		/* Print32_t the read name and paired end flag */
		if(fprintf(outputFP, "%s\t%d\t%d\n",
					a->readName,
					a->space,
					a->numEnds) < 0) {
			PrintError(FnName,
					NULL,
					"Could not write to file",
					Exit,
					WriteFileError);
		}
	}
	else {
		assert(a!=NULL);
		a->readNameLength = (int)strlen(a->readName);
		if(fwrite(&a->readNameLength, sizeof(int32_t), 1, outputFP) != 1 ||
				fwrite(a->readName, sizeof(char), a->readNameLength, outputFP) != a->readNameLength ||
				fwrite(&a->space, sizeof(int32_t), 1, outputFP) != 1 ||
				fwrite(&a->numEnds, sizeof(int32_t), 1, outputFP) != 1) {
			PrintError(FnName,
					NULL,
					"Could not write to file",
					Exit,
					WriteFileError);
		}
	}

	for(i=0;i<a->numEnds;i++) {
		if(EOF == AlignedEndPrint(&a->ends[i],
					outputFP,
					a->space,
					binaryOutput)) {
			PrintError(FnName,
					"a->ends[i]",
					"Could not write to file",
					Exit,
					WriteFileError);
		}
	}
}

/* TODO */
int32_t AlignedReadRead(AlignedRead *a,
		FILE *inputFP,
		int32_t binaryInput)
{
	char *FnName = "AlignedReadRead";
	int32_t i;

	assert(a != NULL);

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
		if(fscanf(inputFP, "%s %d %d",
					a->readName,
					&a->space,
					&a->numEnds)==EOF) {
			/* Free read name before leaving */
			free(a->readName);
			a->readName=NULL;
			return EOF;
		}
		a->readNameLength = (int)strlen(a->readName);
	}
	else {
		if(fread(&a->readNameLength, sizeof(int32_t), 1, inputFP) != 1) {
			/* Free read name before leaving */
			free(a->readName);
			a->readName=NULL;
			return EOF;
		}
		if(fread(a->readName, sizeof(char), a->readNameLength, inputFP) != a->readNameLength ||
				fread(&a->space, sizeof(int32_t), 1, inputFP) != 1 ||
				fread(&a->numEnds, sizeof(int32_t), 1, inputFP) != 1) {
			PrintError(FnName,
					NULL,
					"Could not read from file",
					Exit,
					ReadFileError);
		}
		/* Add the null terminator */
		a->readName[a->readNameLength]='\0';
	}

	/* Reallocate to conserve memory */
	if(0 < a->readNameLength) {
		a->readName = realloc(a->readName, sizeof(char)*(a->readNameLength+1));
		if(NULL == a->readName) {
			PrintError(FnName,
					"a->readName",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
	}
	else {
		free(a->readName);
		a->readName=NULL;
	}

	/* Allocate memory for the first entry */ 
	a->ends = malloc(sizeof(AlignedEnd)*a->numEnds);
	if(NULL==a->ends) {
		PrintError(FnName,
				"a->ends",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Read the alignment */
	for(i=0;i<a->numEnds;i++) {
		AlignedEndInitialize(&a->ends[i]);
		if(EOF==AlignedEndRead(&a->ends[i],
					inputFP,
					a->space,
					binaryInput)) {
			PrintError(FnName, 
					NULL, 
					"Could not read a->ends[i]",
					Exit,
					EndOfFile);
		}
	}

	return 1;
}

/* TODO */
void AlignedReadRemoveDuplicates(AlignedRead *a,
		int32_t sortOrder)
{
	int32_t i;
	/* First entry */
	for(i=0;i<a->numEnds;i++) {
		AlignedEndRemoveDuplicates(&a->ends[i],
				sortOrder);
	}
}

/* TODO */
/* Log-n space */
/* Do not use, since it is buggy and has not been updated lately */  
void AlignedReadQuickSort(AlignedRead *a,
		int32_t sortOrder,
		int32_t showPercentComplete)
{
	int32_t i;
	for(i=0;i<a->numEnds;i++) {
		AlignedEndQuickSort(&a->ends[i],
				sortOrder,
				showPercentComplete);
	}
}

/* TODO */
/* O(n) space, but really double */
void AlignedReadMergeSort(AlignedRead *a,
		int32_t sortOrder,
		int32_t showPercentComplete)
{
	double percentComplete;
	int32_t i;
	/* Sort the first entry */
	for(i=0;i<a->numEnds;i++) {
		percentComplete=0.0;
		AlignedEndMergeSort(&a->ends[i],
				0,
				showPercentComplete);
	}
}

/* TODO */
void AlignedReadReallocate(AlignedRead *a,
		int32_t numEnds)
{
	char *FnName = "AlignedReadReallocate";
	int32_t i;

	/* we have to free if we are reducing the number of entries */
	if(numEnds < a->numEnds) {
		for(i=numEnds;i<a->numEnds;i++) {
			AlignedEndFree(&a->ends[i]);
		}
	}
	a->numEnds = numEnds;

	/* Allocate memory for the entries */ 
	a->ends = realloc(a->ends, sizeof(AlignedEnd)*a->numEnds);
	if(a->numEnds > 0 && NULL==a->ends) {
		if(NULL == a->ends) {
			PrintError(FnName,
					"a->ends",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
}

/* TODO */
void AlignedReadAllocate(AlignedRead *a,
		char *readName,
		int32_t numEnds,
		int32_t space)
{
	char *FnName = "AlignedReadAllocate";
	int32_t i;

	a->numEnds = numEnds;
	a->space = space;
	a->readNameLength = (int)strlen(readName);
	a->readName = malloc(sizeof(char)*(a->readNameLength+1));
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

	/* Allocate memory for the entries */ 
	a->ends = malloc(sizeof(AlignedEnd)*a->numEnds);
	if(0 < a->numEnds && a->ends == NULL) {

		PrintError(FnName,
				"a->ends",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Initialize */
	for(i=0;i<a->numEnds;i++) {
		AlignedEndInitialize(&a->ends[i]);
	}
}

/* TODO */
void AlignedReadFree(AlignedRead *a)
{
	int32_t i;
	for(i=0;i<a->numEnds;i++) {
		AlignedEndFree(&a->ends[i]);
	}
	free(a->readName);
	AlignedReadInitialize(a);
}

/* TODO */
void AlignedReadInitialize(AlignedRead *a) 
{
	a->readNameLength=0;
	a->readName=NULL;
	a->numEnds=0;
	a->ends=NULL;
	a->space=NTSpace;
}

void AlignedReadCopy(AlignedRead *dest, AlignedRead *src) 
{
	int32_t i;

	/* Free and Allocate destination */
	AlignedReadFree(dest);
	AlignedReadAllocate(dest,
			src->readName,
			src->numEnds,
			src->space);
	/* Copy over */
	for(i=0;i<src->numEnds;i++) {
		AlignedEndCopy(&dest->ends[i], &src->ends[i]);
	}
}

/* O(n) space, but really double */
/* Should be a list of pointers to individual AlignedRead */
void AlignedReadMergeSortAll(AlignedRead **a,
		int64_t low, 
		int64_t high)
{
	int64_t mid = (low + high)/2;

	if(low >= high) {
		return;
	}

	/* Split */
	AlignedReadMergeSortAll(a,
			low,
			mid);
	AlignedReadMergeSortAll(a,
			mid+1,
			high);

	AlignedReadMergeAll(a,
			low,
			mid,
			high);
}

void AlignedReadMergeAll(AlignedRead **a,
		int64_t low, 
		int64_t mid,
		int64_t high)
{
	char *FnName="AlignedReadMergeAll";
	int64_t startLower = low;
	int64_t endLower = mid;
	int64_t startUpper = mid+1;
	int64_t endUpper = high;
	int64_t ctr, i;
	AlignedRead **tmp=NULL;

	/* Merge */
	tmp = malloc(sizeof(AlignedRead*)*(high-low+1));
	if(NULL == tmp) {
		PrintError(FnName,
				"tmp",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Merge the two lists */
	ctr=0;
	while( (startLower <= endLower) && (startUpper <= endUpper)) {
		if(AlignedReadCompareAll(a[startLower], a[startUpper]) <= 0) {
			tmp[ctr] = a[startLower];
			startLower++;
		}
		else {
			tmp[ctr] = a[startUpper];
			startUpper++;
		}
		ctr++;
	}
	while(startLower <= endLower) {
		tmp[ctr] = a[startLower];
		startLower++;
		ctr++;
	}
	while(startUpper <= endUpper) {
		tmp[ctr] = a[startUpper];
		startUpper++;
		ctr++;
	}
	/* Copy back */
	for(i=low, ctr=0;
			i<=high;
			i++, ctr++) {
		a[i] = tmp[ctr];
	}

	free(tmp);
	tmp=NULL;
}

/* TODO */
int32_t AlignedReadCompareAll(AlignedRead *one, AlignedRead *two)
{
	char *FnName="AlignedReadCompareAll";
	/* Compare by chr/pos */ 
	int32_t cmp=0;
	int32_t numLeft, i;
	int32_t minIndexOne, minIndexTwo;
	AlignedEnd **oneA=NULL;
	AlignedEnd **twoA=NULL;

	if(1 == one->numEnds &&
			1 == two->numEnds) {
		assert(1 == one->ends[0].numEntries);
		assert(1 == two->ends[0].numEntries);

		return AlignedEndCompare(&one->ends[0], &two->ends[0], AlignedEntrySortByContigPos);
	}
	else {
		assert(one->numEnds == two->numEnds);
		oneA = malloc(sizeof(AlignedEnd*)*one->numEnds);
		if(NULL == oneA) {
			PrintError(FnName,
					"oneA",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		for(i=0;i<one->numEnds;i++) {
			assert(1 == one->ends[i].numEntries);
			oneA[i] = &one->ends[i];
		}
		twoA = malloc(sizeof(AlignedEnd*)*two->numEnds);
		if(NULL == twoA) {
			PrintError(FnName,
					"twoA",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		for(i=0;i<two->numEnds;i++) {
			assert(1 == two->ends[i].numEntries);
			twoA[i] = &two->ends[i];
		}

		numLeft = one->numEnds;
		while(0 < numLeft) {
			/* Get min on one */
			minIndexOne=0;
			for(i=1;i<numLeft;i++) {
				if(0 < AlignedEndCompare(oneA[i], oneA[minIndexOne], AlignedEntrySortByContigPos)) {
					minIndexOne = i;
				}
			}
			/* Get min on two */
			minIndexTwo=0;
			for(i=1;i<numLeft;i++) {
				if(0 < AlignedEndCompare(twoA[i], twoA[minIndexTwo], AlignedEntrySortByContigPos)) {
					minIndexTwo = i;
				}
			}
			/* Compare */
			cmp = AlignedEndCompare(oneA[minIndexOne], twoA[minIndexTwo], AlignedEntrySortByContigPos);
			if(cmp != 0) {
				/* Exit out of the loop */
				numLeft=0;
			}
			else {
				numLeft--;
			}
			/* Remove */
			if(numLeft != minIndexOne) {
				oneA[minIndexOne] = oneA[numLeft];
				oneA[numLeft]=NULL;
			}
			if(numLeft != minIndexTwo) {
				twoA[minIndexTwo] = twoA[numLeft];
				twoA[numLeft]=NULL;
			}
		}

		free(oneA);
		free(twoA);

		return cmp;
	}
}
