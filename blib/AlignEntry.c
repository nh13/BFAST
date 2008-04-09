#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "BError.h"
#include "AlignEntry.h"

/* TODO */
void AlignEntryPrint(AlignEntry *aEntry,
		FILE *outputFP)
{
	/* Print the read name, alignment length, chromosome, position, strand, score */
	fprintf(outputFP, "%s\t%d\t%d\t%d\t%c\t%lf\n",
			aEntry->readName,
			aEntry->length,
			aEntry->chromosome,
			aEntry->position,
			aEntry->strand,
			aEntry->score);

	/* Print the reference and read alignment */
	fprintf(outputFP, "%s\n%s\n",
			aEntry->reference,
			aEntry->read);
}

/* TODO */
int AlignEntryRead(AlignEntry *aEntry,
		FILE *inputFP)
{
	/* Read the read name, alignment length, chromosome, position, strand, score */
	if(fscanf(inputFP, "%s %d %d %d %c %lf\n",
				aEntry->readName,
				&aEntry->length,
				&aEntry->chromosome,
				&aEntry->position,
				&aEntry->strand,
				&aEntry->score)==EOF) {
		return EOF;
	}

	/* Read the reference and read alignment */
	if(fscanf(inputFP, "%s", aEntry->reference)==EOF) {
		PrintError("AlignEntryRead", 
				NULL, 
				"Could not read reference alignent",
				Exit,
				EndOfFile);
	}
	if(fscanf(inputFP, "%s", aEntry->read)==EOF) {
		PrintError("AlignEntryRead", 
				NULL, 
				"Could not read 'read' alignent",
				Exit,
				EndOfFile);
	}

	return 1;
}

/* TODO */
int AlignEntryRemoveDuplicates(AlignEntry **a,
		int length)
{
	int i, prevIndex;

	if(length > 0) {
		/* Sort the array */
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Sorting... %d\n", length);
		}
		AlignEntryQuickSort(a, 0, length-1);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Sorted\n");
		}

		/* Remove duplicates */
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Removing duplicates\n");
		}
		prevIndex=0;
		for(i=1;i<length;i++) {
			if(AlignEntryCompareAtIndex((*a), prevIndex, (*a), i)==0) {
				/* Free memory */
				assert((*a)[i].read!=NULL);
				free((*a)[i].read);
				assert((*a)[i].reference!=NULL);
				free((*a)[i].reference);
			}
			else {
				/* Increment prevIndex */
				prevIndex++;
				/* Copy to prevIndex (incremented) */
				AlignEntryCopyAtIndex((*a), i, (*a), prevIndex);
			}
		}

		/* Reallocate pair */
		length = prevIndex+1;
		(*a) = realloc((*a), sizeof(AlignEntry)*length);
		if(NULL == (*a)) {
			PrintError("AlignEntryRemoveDuplicates",
					"(*a)",
					"Could not reallocate Align Entries while removing duplicates",
					Exit,
					ReallocMemory);
		}

		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Duplicates removed\n");
		}
	}
	return length;
}

/* TODO */
void AlignEntryQuickSort(AlignEntry **a,
		int low,
		int high)
{
	int i;
	int pivot=-1;
	AlignEntry *temp;

	if(low < high) {
		/* Allocate memory for the temp used for swapping */
		temp=malloc(sizeof(AlignEntry));
		if(NULL == temp) {
			PrintError("AlignEntryQuickSort",
					"temp",
					"Could not allocate temp",
					Exit,
					MallocMemory);
		}
		
		pivot = (low+high)/2;

		AlignEntryCopyAtIndex((*a), pivot, temp, 0);
		AlignEntryCopyAtIndex((*a), high, (*a), pivot);
		AlignEntryCopyAtIndex(temp, 0, (*a), high);

		pivot = low;

		for(i=low;i<high;i++) {
			if(AlignEntryCompareAtIndex((*a), i, (*a), high) <= 0) {
				AlignEntryCopyAtIndex((*a), i, temp, 0);
				AlignEntryCopyAtIndex((*a), pivot, (*a), i);
				AlignEntryCopyAtIndex(temp, 0, (*a), pivot);
				pivot++;
			}
		}
		AlignEntryCopyAtIndex((*a), pivot, temp, 0);
		AlignEntryCopyAtIndex((*a), high, (*a), pivot);
		AlignEntryCopyAtIndex(temp, 0, (*a), high);

		/* Free temp before the recursive call, otherwise we have a worst
		 * case of O(n) space (NOT IN PLACE) 
		 * */
		free(temp);

		AlignEntryQuickSort(a, low, pivot-1);
		AlignEntryQuickSort(a, pivot+1, high);
	}
}

void AlignEntryCopyAtIndex(AlignEntry *src, int srcIndex, AlignEntry *dest, int destIndex)
{
	strcpy(dest[destIndex].readName, src[srcIndex].readName);
	dest[destIndex].read = src[srcIndex].read;
	dest[destIndex].reference = src[srcIndex].reference;
	dest[destIndex].length = src[srcIndex].length;
	dest[destIndex].chromosome = src[srcIndex].chromosome;
	dest[destIndex].position = src[srcIndex].position;
	dest[destIndex].strand = src[srcIndex].strand;
	dest[destIndex].score = src[srcIndex].score;

}

int AlignEntryCompareAtIndex(AlignEntry *a, int indexA, AlignEntry *b, int indexB)
{
	int cmp[8];
	int result;
	int i;

	cmp[0] = strcmp(a[indexA].readName, b[indexB].readName);
	cmp[1] = strcmp(a[indexA].read, b[indexB].read);
	cmp[2] = strcmp(a[indexA].reference, b[indexB].reference);
	cmp[3] = (a[indexA].length <= b[indexB].length)?((a[indexA].length<b[indexB].length)?-1:0):1;
	cmp[4] = (a[indexA].chromosome <= b[indexB].chromosome)?((a[indexA].chromosome<b[indexB].chromosome)?-1:0):1;
	cmp[5] = (a[indexA].position <= b[indexB].position)?((a[indexA].position<b[indexB].position)?-1:0):1;
	cmp[6] = (a[indexA].strand <= b[indexB].strand)?((a[indexA].strand<b[indexB].strand)?-1:0):1;
	/* Do not compare score */
	/*
	cmp[7] = (a[indexA].score <= b[indexB].score)?((a[indexA].score<b[indexB].score)?-1:0):1;
	*/

	/* ingenious */
	result = 0;
	for(i=0;i<7;i++) {
		result += pow(10, 8-i-1)*cmp[i];
	}

	if(result < 0) {
		return -1;
	}
	else if(result == 0) {
		return 0;
	}
	else {
		return 1;
	}
}

