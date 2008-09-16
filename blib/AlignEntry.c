#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "BError.h"
#include "BLib.h"
#include "AlignEntry.h"

/* TODO */
int AlignEntryPrint(AlignEntry *aEntry,
		FILE *outputFP)
{
	assert(NULL != aEntry->read);
	assert(NULL != aEntry->reference);
	assert(NULL != aEntry->colorError);

	/* Print the read name, alignment length, chromosome, position, strand, score */
	if(fprintf(outputFP, "%d\t%d\t%d\t%c\t%lf\n",
				aEntry->length,
				aEntry->chromosome,
				aEntry->position,
				aEntry->strand,
				aEntry->score) < 0) {
		return EOF;
	}

	/* Print the reference and read alignment */
	if(fprintf(outputFP, "%s\n%s\n%s\n",
				aEntry->reference,
				aEntry->read,
				aEntry->colorError) < 0) {
		return EOF;
	}

	return 1;
}

/* TODO */
int AlignEntryRead(AlignEntry *aEntry,
		FILE *inputFP)
{
	char *FnName = "AlignEntryRead";

	if(aEntry->read == NULL) {
		aEntry->read = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == aEntry->read) {
			PrintError(FnName,
					"aEntry->read",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	if(aEntry->reference == NULL) {
		aEntry->reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == aEntry->reference) {
			PrintError(FnName,
					"aEntry->reference",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	if(aEntry->colorError == NULL) {
		aEntry->colorError = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == aEntry->colorError) {
			PrintError(FnName,
					"aEntry->colorError",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}

	/* Read the read name, alignment length, chromosome, position, strand, score */
	if(fscanf(inputFP, "%d %d %d %c %lf\n",
				&aEntry->length,
				&aEntry->chromosome,
				&aEntry->position,
				&aEntry->strand,
				&aEntry->score)==EOF) {
		return EOF;
	}

	/* Read the reference and read alignment */
	if(fscanf(inputFP, "%s %s %s", 
				aEntry->reference,
				aEntry->read,
				aEntry->colorError)==EOF) {
		return EOF;
	}

	/*
	assert(((int)strlen(aEntry->read)) == aEntry->length);
	assert(strlen(aEntry->reference) == aEntry->length);
	assert((int)strlen(aEntry->colorError) == aEntry->length);
	*/

	return 1;
}

/* TODO */
int AlignEntryRemoveDuplicates(AlignEntry **a,
		int length,
		int sortOrder)
{
	int i, prevIndex;

	if(length > 0) {
		/* Sort the array */
		AlignEntryQuickSort(a, 0, length-1, sortOrder, 0, NULL, 0);
		/*
		AlignEntryMergeSort(a, 0, length-1, sortOrder, 0, NULL, 0);
		*/

		/* Check sort */
		/*
		   for(i=1;i<length;i++) {
		   assert(AlignEntryCompareAtIndex((*a), i-1, (*a), i, sortOrder)<=0);
		   }
		   */

		/* Remove duplicates */
		prevIndex=0;
		for(i=1;i<length;i++) {
			if(AlignEntryCompareAtIndex((*a), prevIndex, (*a), i, sortOrder)==0) {
				/* Do nothing */
			}
			else {
				/* Increment prevIndex */
				prevIndex++;
				/* Copy to prevIndex (incremented) */
				AlignEntryCopyAtIndex((*a), i, (*a), prevIndex);
			}
		}

		/* Free duplicates */
		for(i=prevIndex+1;i<length;i++) {
			AlignEntryFree(&((*a)[i]));
		}
		/* Update length */
		length = prevIndex+1;
		/* Reallocate based on new length */
		(*a) = realloc((*a), sizeof(AlignEntry)*length);
		if(NULL == (*a)) {
			PrintError("AlignEntryRemoveDuplicates",
					"(*a)",
					"Could not reallocate Align Entries while removing duplicates",
					Exit,
					ReallocMemory);
		}
	}
	return length;
}

/* TODO */
/* Log-n space */
/* Do not use, since it is buggy and has not been updated lately */  
void AlignEntryQuickSort(AlignEntry **a,
		int low,
		int high,
		int sortOrder,
		int showPercentComplete,
		double *curPercent,
		int total)
{
	char *FnName = "AlignEntryQuickSort";
	int i;
	int pivot=-1;
	AlignEntry *temp=NULL;

	if(low < high) {
		/* Allocate memory for the temp used for swapping */
		temp=malloc(sizeof(AlignEntry));
		if(NULL == temp) {
			PrintError(FnName,
					"temp",
					"Could not allocate temp",
					Exit,
					MallocMemory);
		}
		AlignEntryInitialize(temp);

		pivot = AlignEntryGetPivot((*a),
				sortOrder,
				low,
				high);
		if(showPercentComplete == 1 && VERBOSE >= 0) {
			assert(NULL!=curPercent);
			if((*curPercent) < 100.0*((double)low)/total) {
				while((*curPercent) < 100.0*((double)low)/total) {
					(*curPercent) += SORT_ROTATE_INC;
				}
				PrintPercentCompleteShort((*curPercent));
			}
		}

		AlignEntryCopyAtIndex((*a), pivot, temp, 0);
		AlignEntryCopyAtIndex((*a), high, (*a), pivot);
		AlignEntryCopyAtIndex(temp, 0, (*a), high);

		pivot = low;

		for(i=low;i<high;i++) {
			if(AlignEntryCompareAtIndex((*a), i, (*a), high, sortOrder) <= 0) {
				if(i!=pivot) {
					AlignEntryCopyAtIndex((*a), i, temp, 0);
					AlignEntryCopyAtIndex((*a), pivot, (*a), i);
					AlignEntryCopyAtIndex(temp, 0, (*a), pivot);
				}
				pivot++;
			}
		}
		AlignEntryCopyAtIndex((*a), pivot, temp, 0);
		AlignEntryCopyAtIndex((*a), high, (*a), pivot);
		AlignEntryCopyAtIndex(temp, 0, (*a), high);

		/* Free temp before the recursive call, otherwise we have a worst
		 * case of O(n) space (NOT IN PLACE) 
		 * */
		AlignEntryFree(temp);
		free(temp);
		temp=NULL;

		AlignEntryQuickSort(a, low, pivot-1, sortOrder, showPercentComplete, curPercent, total);
		if(showPercentComplete == 1 && VERBOSE >= 0) {
			assert(NULL!=curPercent);
			if((*curPercent) < 100.0*((double)pivot)/total) {
				while((*curPercent) < 100.0*((double)pivot)/total) {
					(*curPercent) += SORT_ROTATE_INC;
				}
				PrintPercentCompleteShort((*curPercent));
			}
		}
		AlignEntryQuickSort(a, pivot+1, high, sortOrder, showPercentComplete, curPercent, total);
		if(showPercentComplete == 1 && VERBOSE >= 0) {
			assert(NULL!=curPercent);
			if((*curPercent) < 100.0*((double)high)/total) {
				while((*curPercent) < 100.0*((double)high)/total) {
					(*curPercent) += SORT_ROTATE_INC;
				}
				PrintPercentCompleteShort((*curPercent));
			}
		}
	}
}

/* TODO */
/* O(n) space, but really double */
void AlignEntryMergeSort(AlignEntry **a,
		int low,
		int high,
		int sortOrder,
		int showPercentComplete,
		double *curPercent,
		int total)
{
	char *FnName = "AlignEntryMergeSort";
	int i, ctr;
	int mid = (low + high)/2;
	int startLower =  low;
	int endLower = mid;
	int startUpper = mid + 1;
	int endUpper = high;
	AlignEntry *tempEntries=NULL;

	if(low >= high) {
		if(VERBOSE >= 0 &&
				showPercentComplete == 1) { 
			assert(NULL!=curPercent);
			if((*curPercent) < 100.0*((double)low)/total) {
				while((*curPercent) < 100.0*((double)low)/total) {
					(*curPercent) += SORT_ROTATE_INC;
				}
				PrintPercentCompleteShort((*curPercent));
			}
		}
		return;
	}


	/* Partition the list into two lists and sort them recursively */
	AlignEntryMergeSort(a,
			low,
			mid,
			sortOrder,
			showPercentComplete,
			curPercent,
			total);
	AlignEntryMergeSort(a,
			mid+1,
			high,
			sortOrder,
			showPercentComplete,
			curPercent,
			total);

	/* Allocate pointers */
	tempEntries = malloc(sizeof(AlignEntry)*(high-low+1));
	if(NULL == tempEntries) {
		PrintError(FnName,
				"tempEntries",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Initialize */
	for(i=0;i<high-low+1;i++) {
		AlignEntryInitialize(&tempEntries[i]);
	}

	/* Merge the two lists */
	ctr=0;
	while( (startLower <= endLower) && (startUpper <= endUpper)) {
		if(AlignEntryCompareAtIndex((*a), startLower, (*a), startUpper, sortOrder) <= 0) {
			AlignEntryCopyAtIndex((*a), startLower, tempEntries, ctr);
			startLower++;
		}
		else {
			AlignEntryCopyAtIndex((*a), startUpper, tempEntries, ctr);
			startUpper++;
		}
		ctr++;
	}
	while(startLower <= endLower) {
		AlignEntryCopyAtIndex((*a), startLower, tempEntries, ctr);
		startLower++;
		ctr++;
	}
	while(startUpper <= endUpper) {
		AlignEntryCopyAtIndex((*a), startUpper, tempEntries, ctr);
		startUpper++;
		ctr++;
	}
	/* Copy back */
	for(i=low, ctr=0;
			i<=high;
			i++, ctr++) {
		AlignEntryCopyAtIndex(tempEntries, ctr, (*a), i);
	}

	/* Free memory */
	for(i=0;i<high-low+1;i++) {
		AlignEntryFree(&tempEntries[i]);
	}
	free(tempEntries);
	tempEntries=NULL;

	/* Test sort */
	/* 
	   for(i=low+1;i<=high;i++) {
	   assert(AlignEntryCompareAtIndex((*a), i-1, (*a), i, sortOrder) <= 0);
	   }
	   */
}

/* TODO */
int AlignEntryCompareAtIndex(AlignEntry *a, int indexA, AlignEntry *b, int indexB, int sortOrder)
{
	int cmp[5];
	int i;
	int top;

	if(sortOrder == AlignEntrySortByAll) {

		/* Old 
		   cmp[0] = strcmp(a[indexA].read, b[indexB].read);
		   cmp[1] = strcmp(a[indexA].reference, b[indexB].reference);
		   cmp[2] = (a[indexA].chromosome <= b[indexB].chromosome)?((a[indexA].chromosome<b[indexB].chromosome)?-1:0):1;
		   cmp[3] = (a[indexA].position <= b[indexB].position)?((a[indexA].position<b[indexB].position)?-1:0):1;
		   cmp[4] = (a[indexA].strand <= b[indexB].strand)?((a[indexA].strand<b[indexB].strand)?-1:0):1;
		   */

		/* If there are multiple alignments to the same starting chr/pos/strand with the same score,
		 * this will pick ensure that we will only pick one of them.
		 * */
		cmp[0] = (a[indexA].chromosome <= b[indexB].chromosome)?((a[indexA].chromosome<b[indexB].chromosome)?-1:0):1;
		cmp[1] = (a[indexA].position <= b[indexB].position)?((a[indexA].position<b[indexB].position)?-1:0):1;
		cmp[2] = (a[indexA].strand <= b[indexB].strand)?((a[indexA].strand<b[indexB].strand)?-1:0):1;
		cmp[3] = (a[indexA].score <= b[indexB].score)?((a[indexA].score<b[indexB].score)?-1:0):1;

		top = 4;
	}
	else {
		assert(sortOrder == AlignEntrySortByChrPos);
		cmp[0] = (a[indexA].chromosome <= b[indexB].chromosome)?((a[indexA].chromosome<b[indexB].chromosome)?-1:0):1;
		cmp[1] = (a[indexA].position <= b[indexB].position)?((a[indexA].position<b[indexB].position)?-1:0):1;

		top = 2;

		/*
		   fprintf(stderr, "a[%d].chromosome=%d\ta[%d].position=%d\n",
		   indexA,
		   a[indexA].chromosome,
		   indexA,
		   a[indexA].position);
		   fprintf(stderr, "a[%d].chromosome=%d\ta[%d].position=%d\n",
		   indexB,
		   a[indexB].chromosome,
		   indexB,
		   a[indexB].position);
		   fprintf(stderr, "cmp[%d]=%d\ncmp[%d]=%d\n",
		   0,
		   cmp[0],
		   1,
		   cmp[1]);
		   */
	}

	/* ingenious */
	for(i=0;i<top;i++) {
		if(cmp[i] != 0) {
			return cmp[i];
		}
	}

	return 0;
}

/* TODO */
void AlignEntryCopyAtIndex(AlignEntry *src, int srcIndex, AlignEntry *dest, int destIndex)
{
	if(dest != src || srcIndex != destIndex) {
		AlignEntryCopy(&(src[srcIndex]), &(dest[destIndex]));
	}
}

/* TODO */
void AlignEntryCopy(AlignEntry *src, AlignEntry *dest)
{
	char *FnName = "AlignEntryCopy";
	if(src != dest) {
		if(NULL == dest->read) {
			dest->read = malloc(sizeof(char)*SEQUENCE_LENGTH);
			if(NULL == dest->read) {
				PrintError(FnName,
						"dest->read",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}
		assert(src->read != NULL);
		strcpy(dest->read, src->read);
		if(NULL == dest->reference) {
			dest->reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
			if(NULL == dest->reference) {
				PrintError(FnName,
						"dest->reference",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}
		assert(src->reference!= NULL);
		strcpy(dest->reference, src->reference);
		if(NULL == dest->colorError) {
			dest->colorError = malloc(sizeof(char)*SEQUENCE_LENGTH);
			if(NULL == dest->colorError) {
				PrintError(FnName,
						"dest->colorError",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}
		assert(src->colorError!= NULL);
		strcpy(dest->colorError, src->colorError);
		dest->length = src->length;
		dest->chromosome = src->chromosome;
		dest->position = src->position;
		dest->strand = src->strand;
		dest->score = src->score;
	}
}

void AlignEntryFree(AlignEntry *aEntry)
{
	free(aEntry->read);
	free(aEntry->reference);
	free(aEntry->colorError);
	AlignEntryInitialize(aEntry);
}

void AlignEntryInitialize(AlignEntry *aEntry) 
{
	aEntry->read=NULL;
	aEntry->reference=NULL;
	aEntry->colorError=NULL;
	aEntry->length=0;
	aEntry->chromosome=0;
	aEntry->position=0;
	aEntry->strand=0;
	aEntry->score=0.0;
}

/* TODO */
/* Debugging function */
void AlignEntryCheckReference(AlignEntry *aEntry, RGBinary *rg)
{
	char *FnName = "AlignEntryCheckReference";
	int i;
	int curPos;
	char rgBase;
	char reference[SEQUENCE_LENGTH]="\0";

	if(aEntry->strand == REVERSE) {
		GetReverseComplimentAnyCase(aEntry->reference, reference, aEntry->length);
	}
	else {
		strcpy(reference, aEntry->reference);
	}

	for(i=0, curPos = aEntry->position;i<aEntry->length;i++) {
		if(reference[i] != GAP) {
			rgBase = RGBinaryGetBase(rg, aEntry->chromosome, curPos);	
			if(rgBase != reference[i]) {
				fprintf(stderr, "\n[%d]\t[%d]\n[%c]\t[%c]\n[%s]\n",
						curPos,
						i,
						reference[i],
						rgBase,
						reference);
				AlignEntryPrint(aEntry, stderr);
				PrintError(FnName,
						NULL,
						"Reference in the align entry does not match the reference genome",
						Exit,
						OutOfRange);
			}
			curPos++;
		}
	}
}

int AlignEntryGetPivot(AlignEntry *a,
		int sortOrder,
		int low,
		int high) 
{
	int cmp[3];
	int pivot = (low + high)/2;
	cmp[0] = AlignEntryCompareAtIndex(a, low, a, pivot, sortOrder); 
	cmp[1] = AlignEntryCompareAtIndex(a, low, a, high, sortOrder); 
	cmp[2] = AlignEntryCompareAtIndex(a, pivot, a, high, sortOrder); 

	if(cmp[0] <= 0) {
		/* low <= pivot */
		if(cmp[1] >= 0) {
			/* high <= low */
			/* so high <= low <= pivot */
			pivot = low;
		}
		else {
			/* low < high */
			if(cmp[2] <= 0) {
				/* pivot <= high */
				/* so low <= pivot <= high */
				/* choose pivot */
			}
			else {
				/* high < pivot */
				/* so low < high < pivot */
				pivot = high;
			}
		}
	}
	else {
		/* pivot < low */
		if(cmp[1] <= 0) {
			/* low <= high */
			/* so pivot < low <= high */
			pivot = low;
		}
		else {
			/* high < low */
			if(cmp[2] <= 0) {
				/* pivot <= high */
				/* so pivot <= high < low */
				pivot = high;
			}
			else {
				/* high < pivot */
				/* so high < pivot < low */
				/* choose pivot */
			}
		}
	}
	return pivot;
}
