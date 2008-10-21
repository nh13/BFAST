#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "BError.h"
#include "BLib.h"
#include "AlignEntry.h"

/* TODO */
int AlignEntryPrint(AlignEntry *a,
		FILE *outputFP,
		int space,
		int binaryOutput)
{
	assert(NULL != a->read);
	assert(NULL != a->reference);
	assert(space == NTSpace ||
			(space == ColorSpace && NULL != a->colorError));

	/* Print contig name */
	if(BStringPrint(&a->contigName, outputFP, binaryOutput) < 0) {
		return EOF;
	}

	if(binaryOutput == TextOutput) {
		if(fprintf(outputFP, "%u\t%u\t%c\t%lf\t%u\t%u\n",
					a->contig,
					a->position,
					a->strand,
					a->score,
					a->referenceLength,
					a->length) < 0) {
			return -1;
		}
	}
	else {
		if(fwrite(&a->contig, sizeof(uint32_t), 1, outputFP) != 1 ||
				fwrite(&a->position, sizeof(uint32_t), 1, outputFP) != 1 ||
				fwrite(&a->strand, sizeof(char), 1, outputFP) != 1 ||
				fwrite(&a->score, sizeof(double), 1, outputFP) != 1 ||
				fwrite(&a->referenceLength, sizeof(uint32_t), 1, outputFP) != 1 ||
				fwrite(&a->length, sizeof(uint32_t), 1, outputFP) != 1) {
			return -1;
		}
	}

	/* Print the reference and read alignment */
	if(BStringPrint(&a->reference, outputFP, binaryOutput) < 0 ||
			BStringPrint(&a->read, outputFP, binaryOutput) < 0) {
		return -1;
	}

	/* Print the color errors if necessary */
	if(space == ColorSpace) {
		if(BStringPrint(&a->colorError, outputFP, binaryOutput)<0) {
			return -1;
		}
	}

	return 1;
}

/* TODO */
int AlignEntryRead(AlignEntry *a,
		FILE *inputFP,
		int space,
		int binaryInput)
{
	char *FnName = "AlignEntryRead";

	/* Read in contig name */
	if(BStringRead(&a->contigName, inputFP, binaryInput)==EOF) {
		return EOF;
	}

	if(binaryInput == TextOutput) {

		if(fscanf(inputFP, "%u\t%u\t%c\t%lf\t%u\t%u\n",
					&a->contig,
					&a->position,
					&a->strand,
					&a->score,
					&a->referenceLength,
					&a->length) < 0) {
			return EOF;
		}

	}
	else {
		if(fread(&a->contig, sizeof(uint32_t), 1, inputFP) != 1 ||
				fread(&a->position, sizeof(uint32_t), 1, inputFP) != 1 ||
				fread(&a->strand, sizeof(char), 1, inputFP) != 1 ||
				fread(&a->score, sizeof(double), 1, inputFP) != 1 ||
				fread(&a->referenceLength, sizeof(uint32_t), 1, inputFP) != 1 ||
				fread(&a->length, sizeof(uint32_t), 1, inputFP) != 1) {
			return EOF;
		}
	}
	
	/* Read the reference and read alignment */
	if(BStringRead(&a->reference, inputFP, binaryInput) < 0 ||
			BStringRead(&a->read, inputFP, binaryInput) < 0) {
		return EOF;
	}

		/* Read the color errors if necessary */
	if(space == ColorSpace) {
		if(BStringRead(&a->colorError, inputFP, binaryInput)<0) {
			return EOF;
		}
	}

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

		/* If there are multiple alignments to the same starting chr/pos/strand with the same score,
		 * this will pick ensure that we will only pick one of them.
		 * */
		cmp[0] = (a[indexA].contig <= b[indexB].contig)?((a[indexA].contig<b[indexB].contig)?-1:0):1;
		cmp[1] = (a[indexA].position <= b[indexB].position)?((a[indexA].position<b[indexB].position)?-1:0):1;
		cmp[2] = (a[indexA].strand <= b[indexB].strand)?((a[indexA].strand<b[indexB].strand)?-1:0):1;
		cmp[3] = (a[indexA].score <= b[indexB].score)?((a[indexA].score<b[indexB].score)?-1:0):1;

		top = 4;
	}
	else {
		assert(sortOrder == AlignEntrySortByContigPos);
		cmp[0] = (a[indexA].contig <= b[indexB].contig)?((a[indexA].contig<b[indexB].contig)?-1:0):1;
		cmp[1] = (a[indexA].position <= b[indexB].position)?((a[indexA].position<b[indexB].position)?-1:0):1;

		top = 2;

		/*
		   fprintf(stderr, "a[%d].contig=%d\ta[%d].position=%d\n",
		   indexA,
		   a[indexA].contig,
		   indexA,
		   a[indexA].position);
		   fprintf(stderr, "a[%d].contig=%d\ta[%d].position=%d\n",
		   indexB,
		   a[indexB].contig,
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
		/* Contig name length */
		BStringCopy(&dest->contigName, &src->contigName);
		BStringCopy(&dest->read, &src->read);
		BStringCopy(&dest->reference, &src->reference);
		/* Color error, if necessary */
		if(src->colorError != NULL) {
			BStringCopy(&dest->colorError, &src->colorError);
		}
		/* Metadata */
		dest->referenceLength = src->referenceLength;
		dest->length = src->length;
		dest->contig = src->contig;
		dest->position = src->position;
		dest->strand = src->strand;
		dest->score = src->score;
	}
}

void AlignEntryFree(AlignEntry *a)
{
	BStringFree(a->contigName);
	BStringFree(a->read);
	BStringFree(a->reference);
	BStringFree(a->colorError);
	AlignEntryInitialize(a);
}

void AlignEntryInitialize(AlignEntry *a) 
{
	BStringInitialize(&a->contigName);
	a->contig=0;
	a->position=0;
	a->strand=0;
	a->score=0.0;
	a->referenceLength=0;
	a->length=0;
	BStringInitialize(&a->read);
	BStringInitialize(&a->reference);
	BStringInitialize(&a->colorError);
}

/* TODO */
/* Debugging function */
void AlignEntryCheckReference(AlignEntry *a, RGBinary *rg, int space)
{
	char *FnName = "AlignEntryCheckReference";
	int i;
	int curPos;
	char rgBase;

	if(a->strand == REVERSE) {
		GetReverseComplimentAnyCase(&a->reference);
	}

	for(i=0, curPos = a->position;i<a->length;i++) {
		if(a->reference.string[i] != GAP) {
			rgBase = RGBinaryGetBase(rg, a->contig, curPos);	
			if(rgBase != a->reference.string[i]) {
				fprintf(stderr, "\n[%d]\t[%d]\n[%c]\t[%c]\n[%s]\n",
						curPos,
						i,
						a->reference.string[i],
						rgBase,
						a->reference.string);
				AlignEntryPrint(a, stderr, space, 0);
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
