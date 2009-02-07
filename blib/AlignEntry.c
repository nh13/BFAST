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

	/*
	   assert(((int)strlen(a->contigName)) == a->contigNameLength);
	   assert(((int)strlen(a->read)) == a->length);
	   assert(strlen(a->reference) == a->length);
	   assert((int)strlen(a->colorError) == a->length);
	   */

	if(binaryOutput == TextOutput) {

		if(fprintf(outputFP, "%s\t%u\t%u\t%c\t%lf\t%u\t%u\n",
					a->contigName,
					a->contig,
					a->position,
					a->strand,
					a->score,
					a->referenceLength,
					a->length) < 0) {
			return EOF;
		}

		/* Print the reference and read alignment */
		if(fprintf(outputFP, "%s\n%s\n",
					a->reference,
					a->read) < 0) {
			return EOF;
		}

		/* Print the color errors if necessary */
		if(space == ColorSpace) {
			if(fprintf(outputFP, "%s\n",
						a->colorError) < 0) {
				return EOF;
			}
		}
	}
	else {
		if(fwrite(&a->contigNameLength, sizeof(int32_t), 1, outputFP) != 1 ||
				fwrite(a->contigName, sizeof(char), a->contigNameLength, outputFP) != a->contigNameLength ||
				fwrite(&a->contig, sizeof(uint32_t), 1, outputFP) != 1 ||
				fwrite(&a->position, sizeof(uint32_t), 1, outputFP) != 1 ||
				fwrite(&a->strand, sizeof(char), 1, outputFP) != 1 ||
				fwrite(&a->score, sizeof(double), 1, outputFP) != 1 ||
				fwrite(&a->referenceLength, sizeof(uint32_t), 1, outputFP) != 1 ||
				fwrite(&a->length, sizeof(uint32_t), 1, outputFP) != 1 ||
				fwrite(a->read, sizeof(char), a->length, outputFP) != a->length ||
				fwrite(a->reference, sizeof(char), a->length, outputFP) != a->length) {
			return EOF;
		}
		if(ColorSpace==space) {
			if(fwrite(a->colorError, sizeof(char), a->length, outputFP) != a->length) {
				return EOF;
			}
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
	char tempContigName[MAX_CONTIG_NAME_LENGTH]="\0";

	/* Allocate memory for the alignment */
	if(a->read == NULL) {
		a->read = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == a->read) {
			PrintError(FnName,
					"a->read",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	if(a->reference == NULL) {
		a->reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == a->reference) {
			PrintError(FnName,
					"a->reference",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	if(space == ColorSpace) {
		if(a->colorError == NULL) {
			a->colorError = malloc(sizeof(char)*SEQUENCE_LENGTH);
			if(NULL == a->colorError) {
				PrintError(FnName,
						"a->colorError",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}
	}

	if(binaryInput == TextInput) {

		if(fscanf(inputFP, "%s\t%u\t%u\t%c\t%lf\t%u\t%u\n",
					tempContigName,
					&a->contig,
					&a->position,
					&a->strand,
					&a->score,
					&a->referenceLength,
					&a->length) < 0) {
			return EOF;
		}

		/* Copy over contig name */
		a->contigNameLength = (int)strlen(tempContigName);
		a->contigName = malloc(sizeof(char)*(a->contigNameLength+1));
		if(NULL==a->contigName) {
			PrintError(FnName,
					"a->contigName",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		strcpy(a->contigName, tempContigName);

		/* Read the reference and read alignment */
		if(fscanf(inputFP, "%s %s", 
					a->reference,
					a->read)==EOF) {
			return EOF;
		}

		/* Read the color errors if necessary */
		if(space == ColorSpace) {
			if(fscanf(inputFP, "%s",
						a->colorError)==EOF) {
				return EOF;
			}
		}
	}
	else {
		if(fread(&a->contigNameLength, sizeof(int32_t), 1, inputFP) != 1) {
			return EOF;
		}
		/* Copy over contig name */
		a->contigName = malloc(sizeof(char)*(a->contigNameLength+1));
		if(NULL==a->contigName) {
			PrintError(FnName,
					"a->contigName",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		if(fread(a->contigName, sizeof(char), a->contigNameLength, inputFP) != a->contigNameLength ||
				fread(&a->contig, sizeof(uint32_t), 1, inputFP) != 1 ||
				fread(&a->position, sizeof(uint32_t), 1, inputFP) != 1 ||
				fread(&a->strand, sizeof(char), 1, inputFP) != 1 ||
				fread(&a->score, sizeof(double), 1, inputFP) != 1 ||
				fread(&a->referenceLength, sizeof(uint32_t), 1, inputFP) != 1 ||
				fread(&a->length, sizeof(uint32_t), 1, inputFP) != 1 ||
				fread(a->read, sizeof(char), a->length, inputFP) != a->length ||
				fread(a->reference, sizeof(char), a->length, inputFP) != a->length) {
			return EOF;
		}
		/* Add the null terminator to strings */
		a->contigName[a->contigNameLength]='\0';
		a->read[a->length]='\0';
		a->reference[a->length]='\0';
		if(ColorSpace==space) {
			if(fread(a->colorError, sizeof(char), a->length, inputFP) != a->length) {
				return EOF;
			}
			a->colorError[a->length]='\0';
		}
	}

	/* Reallocate to conserve memory */
	assert(a->length > 0);
	a->read = realloc(a->read, sizeof(char)*(a->length+1));
	if(NULL == a->read) {
		PrintError(FnName,
				"a->read",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}
	/* Reference */
	a->reference = realloc(a->reference, sizeof(char)*(a->length+1));
	if(NULL == a->reference) {
		PrintError(FnName,
				"a->reference",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}
	/* Color error, if necessary */
	if(space == ColorSpace) {
		a->colorError = realloc(a->colorError, sizeof(char)*(a->length+1));
		if(NULL == a->colorError) {
			PrintError(FnName,
					"a->colorError",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
	}
	else {
		assert(NULL == a->colorError);
	}

	/*
	   assert(((int)strlen(a->contigName)) == a->contigNameLength);
	   assert(((int)strlen(a->read)) == a->length);
	   assert(strlen(a->reference) == a->length);
	   assert((int)strlen(a->colorError) == a->length);
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
		   cmp[2] = (a[indexA].contig <= b[indexB].contig)?((a[indexA].contig<b[indexB].contig)?-1:0):1;
		   cmp[3] = (a[indexA].position <= b[indexB].position)?((a[indexA].position<b[indexB].position)?-1:0):1;
		   cmp[4] = (a[indexA].strand <= b[indexB].strand)?((a[indexA].strand<b[indexB].strand)?-1:0):1;
		   */

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
		assert(src->contigNameLength > 0);
		dest->contigNameLength = src->contigNameLength;
		/* Contig name */
		dest->contigName = realloc(dest->contigName, sizeof(char)*(dest->contigNameLength+1));
		if(NULL == dest->contigName) {
			PrintError(FnName,
					"dest->contigName",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		assert(src->contigName!= NULL);
		strcpy(dest->contigName, src->contigName);
		/* Read */
		assert(src->length > 0);
		dest->read = realloc(dest->read, sizeof(char)*(src->length+1));
		if(NULL == dest->read) {
			PrintError(FnName,
					"dest->read",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		assert(src->read != NULL);
		strcpy(dest->read, src->read);
		/* Reference */
		dest->reference = realloc(dest->reference, sizeof(char)*(src->length+1));
		if(NULL == dest->reference) {
			PrintError(FnName,
					"dest->reference",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		assert(src->reference!= NULL);
		strcpy(dest->reference, src->reference);
		/* Color error, if necessary */
		if(src->colorError != NULL) {
			assert(src->length > 0);
			dest->colorError = realloc(dest->colorError, sizeof(char)*(src->length+1));
			if(NULL == dest->colorError) {
				PrintError(FnName,
						"dest->colorError",
						"Could not reallocate memory",
						Exit,
						ReallocMemory);
			}
			assert(src->colorError!= NULL);
			strcpy(dest->colorError, src->colorError);
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
	free(a->contigName);
	free(a->read);
	free(a->reference);
	free(a->colorError);
	AlignEntryInitialize(a);
}

void AlignEntryInitialize(AlignEntry *a) 
{
	a->contigNameLength=0;
	a->contigName=NULL;
	a->contig=0;
	a->position=0;
	a->strand=0;
	a->score=0.0;
	a->referenceLength=0;
	a->length=0;
	a->read=NULL;
	a->reference=NULL;
	a->colorError=NULL;
}

/* TODO */
/* Debugging function */
void AlignEntryCheckReference(AlignEntry *a, RGBinary *rg, int space)
{
	char *FnName = "AlignEntryCheckReference";
	int i;
	int curPos;
	char rgBase;
	char reference[SEQUENCE_LENGTH]="\0";

	if(a->strand == REVERSE) {
		GetReverseComplimentAnyCase(a->reference, reference, a->length);
	}
	else {
		strcpy(reference, a->reference);
	}

	for(i=0, curPos = a->position;i<a->length;i++) {
		if(reference[i] != GAP) {
			rgBase = RGBinaryGetBase(rg, a->contig, curPos);	
			if(rgBase != reference[i]) {
				fprintf(stderr, "\n[%d]\t[%d]\n[%c]\t[%c]\n[%s]\n",
						curPos,
						i,
						reference[i],
						rgBase,
						reference);
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

int64_t AlignEntryGetSize(AlignEntry *a)
{
	int64_t size = 0;

	size += sizeof(AlignEntry);

	if(0 < a->contigNameLength) {
		size += sizeof(char)*(a->contigNameLength+1);
	}
	if(0 < a->length) {
		size += 2*sizeof(char)*(a->length+1); /* Read and reference */
		if(NULL != a->colorError) {
			size += sizeof(char)*(a->length+1); /* Color error */
		}
	}

	return size;
}
