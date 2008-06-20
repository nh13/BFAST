#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "BError.h"
#include "BLib.h"
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
	int i;

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

	/* Update reference length */
	aEntry->referenceLength = aEntry->length;
	for(i=0;i<aEntry->length;i++) {
		if(GAP == aEntry->reference[i]) {
			aEntry->referenceLength--;
		}
	}
	assert(aEntry->referenceLength >= 0 && aEntry->referenceLength <= aEntry->length);

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
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Sorting... %d\n", length);
		}
		AlignEntryQuickSort(a, 0, length-1, sortOrder, 0, NULL, 0);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Sorted\n");
		}

		/* Remove duplicates */
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Removing duplicates\n");
		}
		prevIndex=0;
		for(i=1;i<length;i++) {
			if(AlignEntryCompareAtIndex((*a), prevIndex, (*a), i, sortOrder)==0) {
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
/* Log-n space */
void AlignEntryQuickSort(AlignEntry **a,
		int low,
		int high,
		int sortOrder,
		int showPercentComplete,
		double *curPercent,
		int total)
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
				fprintf(stderr, "\r%3.2lf percent complete", 100.0*((double)low)/total);
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
		free(temp);
		temp=NULL;

		AlignEntryQuickSort(a, low, pivot-1, sortOrder, showPercentComplete, curPercent, total);
		if(showPercentComplete == 1 && VERBOSE >= 0) {
			assert(NULL!=curPercent);
			if((*curPercent) < 100.0*((double)pivot)/total) {
				while((*curPercent) < 100.0*((double)pivot)/total) {
					(*curPercent) += SORT_ROTATE_INC;
				}
				fprintf(stderr, "\r%3.2lf percent complete", 100.0*((double)pivot)/total);
			}
		}
		AlignEntryQuickSort(a, pivot+1, high, sortOrder, showPercentComplete, curPercent, total);
		if(showPercentComplete == 1 && VERBOSE >= 0) {
			assert(NULL!=curPercent);
			if((*curPercent) < 100.0*((double)high)/total) {
				while((*curPercent) < 100.0*((double)high)/total) {
					(*curPercent) += SORT_ROTATE_INC;
				}
				fprintf(stderr, "\r%3.2lf percent complete", 100.0*((double)high)/total);
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
	AlignEntry *tempEntries;

	if(low >= high) {
		if(VERBOSE >= 0 &&
				showPercentComplete == 1) { 
			assert(NULL!=curPercent);
			if((*curPercent) < 100.0*((double)low)/total) {
				while((*curPercent) < 100.0*((double)low)/total) {
					(*curPercent) += SORT_ROTATE_INC;
				}
				fprintf(stderr, "\r%3.2lf percent complete", 100.0*((double)low)/total);
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
	for(i=0;i<high-low+1;i++) {
		AlignEntryInitialize(&tempEntries[i]);
		tempEntries[i].read = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == tempEntries[i].read) {
			PrintError(FnName,
					"tempEntries[i].read",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		tempEntries[i].reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == tempEntries[i].reference) {
			PrintError(FnName,
					"tempEntries[i].reference",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
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
}

/* TODO */
int AlignEntryCompareAtIndex(AlignEntry *a, int indexA, AlignEntry *b, int indexB, int sortOrder)
{
	int cmp[6];
	int i;
	int top;

	if(sortOrder == AlignEntrySortByAll) {

		cmp[0] = strcmp(a[indexA].readName, b[indexB].readName);
		cmp[1] = strcmp(a[indexA].read, b[indexB].read);
		cmp[2] = strcmp(a[indexA].reference, b[indexB].reference);
		cmp[3] = (a[indexA].chromosome <= b[indexB].chromosome)?((a[indexA].chromosome<b[indexB].chromosome)?-1:0):1;
		cmp[4] = (a[indexA].position <= b[indexB].position)?((a[indexA].position<b[indexB].position)?-1:0):1;
		cmp[5] = (a[indexA].strand <= b[indexB].strand)?((a[indexA].strand<b[indexB].strand)?-1:0):1;

		top = 6;
	}
	else {
		assert(sortOrder == AlignEntrySortByChrPos);
		cmp[0] = (a[indexA].chromosome <= b[indexB].chromosome)?((a[indexA].chromosome<b[indexB].chromosome)?-1:0):1;
		cmp[1] = (a[indexA].position <= b[indexB].position)?((a[indexA].position<b[indexB].position)?-1:0):1;

		top = 2;
	}

	/* ingenious */
	for(i=0;i<top;i++) {
		if(cmp[i] != 0) {
			return cmp[0];
		}
	}

	return 0;
}

/* TODO */
int AlignEntryGetOneRead(AlignEntry **entries, FILE *fp)
{
	char *FnName="AlignEntryGetOne";
	int cont = 1;
	int numEntries = 0;
	fpos_t position;
	AlignEntry temp;

	/* Initialize temp */
	temp.read = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL==temp.read) {
		PrintError(FnName,
				"temp.read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	temp.reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL==temp.reference) {
		PrintError(FnName,
				"temp.reference",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Try to read in the first entry */
	if(EOF != AlignEntryRead(&temp, fp)) {

		/* Initialize the first entry */
		numEntries++;
		(*entries) = malloc(sizeof(AlignEntry)*numEntries);
		if(NULL == (*entries)) {
			PrintError(FnName,
					"(*entries)",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		(*entries)[numEntries-1].read = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL==(*entries)[numEntries-1].read) {
			PrintError(FnName,
					"(*entries)[numEntries-1].read",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		(*entries)[numEntries-1].reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL==(*entries)[numEntries-1].reference) {
			PrintError(FnName,
					"(*entries)[numEntries-1].reference",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Copy over from temp */
		AlignEntryCopy(&temp, &(*entries)[numEntries-1]);

		/* Save current position */
		fgetpos(fp, &position);

		while(cont == 1 && EOF != AlignEntryRead(&temp, fp)) {
			/* Check that the read name matches the previous */
			if(strcmp(temp.readName, (*entries)[numEntries-1].readName)==0) {
				/* Save current position */
				fgetpos(fp, &position);
				/* Copy over */
				numEntries++;
				(*entries) = realloc((*entries), sizeof(AlignEntry)*numEntries);
				if(NULL == (*entries)) {
					PrintError(FnName,
							"(*entries)",
							"Could not reallocate memory",
							Exit,
							ReallocMemory);
				}
				(*entries)[numEntries-1].read = malloc(sizeof(char)*SEQUENCE_LENGTH);
				if(NULL==(*entries)[numEntries-1].read) {
					PrintError(FnName,
							"(*entries)[numEntries-1].read",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}
				(*entries)[numEntries-1].reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
				if(NULL==(*entries)[numEntries-1].reference) {
					PrintError(FnName,
							"(*entries)[numEntries-1].reference",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}
				/* Copy over from temp */
				AlignEntryCopy(&temp, &(*entries)[numEntries-1]);
			}
			else {
				/* Exit the loop */
				cont = 0;
				/* Reset position in the file */
				fsetpos(fp, &position);
			}
		}
	}

	/* Free memory */
	assert(temp.read!=NULL);
	free(temp.read);
	assert(temp.reference!=NULL);
	free(temp.reference);

	return numEntries;
}

/* TODO */
int AlignEntryGetAll(AlignEntry **entries, FILE *fp)
{
	char *FnName="AlignEntryGetAll";
	int numEntries = 0;
	AlignEntry temp;

	/* Initialize temp */
	temp.read = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL==temp.read) {
		PrintError(FnName,
				"temp.read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	temp.reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL==temp.reference) {
		PrintError(FnName,
				"temp.reference",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Try to read in the first entry */
	numEntries=0;
	while(EOF != AlignEntryRead(&temp, fp)) {

		/* Initialize the first entry */
		numEntries++;
		(*entries) = realloc((*entries), sizeof(AlignEntry)*numEntries);
		if(NULL == (*entries)) {
			PrintError(FnName,
					"(*entries)",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		(*entries)[numEntries-1].read = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL==(*entries)[numEntries-1].read) {
			PrintError(FnName,
					"(*entries)[numEntries-1].read",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		(*entries)[numEntries-1].reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL==(*entries)[numEntries-1].reference) {
			PrintError(FnName,
					"(*entries)[numEntries-1].reference",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Copy over from temp */
		AlignEntryCopy(&temp, &(*entries)[numEntries-1]);
	}

	/* Free memory */
	assert(NULL!=temp.read);
	free(temp.read);
	assert(NULL!=temp.reference);
	free(temp.reference);

	return numEntries;
}

/* TODO */
void AlignEntryCopyAtIndex(AlignEntry *src, int srcIndex, AlignEntry *dest, int destIndex)
{
	if(dest != src || srcIndex != destIndex) {
		strcpy(dest[destIndex].readName, src[srcIndex].readName);
		strcpy(dest[destIndex].reference, src[srcIndex].reference);
		strcpy(dest[destIndex].read, src[srcIndex].read);
		dest[destIndex].length = src[srcIndex].length;
		dest[destIndex].referenceLength = src[srcIndex].referenceLength;
		dest[destIndex].chromosome = src[srcIndex].chromosome;
		dest[destIndex].position = src[srcIndex].position;
		dest[destIndex].strand = src[srcIndex].strand;
		dest[destIndex].score = src[srcIndex].score;
	}

}

/* TODO */
void AlignEntryCopy(AlignEntry *src, AlignEntry *dst)
{
	if(src != dst) {
		strcpy(dst->readName, src->readName);
		strcpy(dst->read, src->read);
		strcpy(dst->reference, src->reference);
		dst->length = src->length;
		dst->referenceLength = src->referenceLength;
		assert(dst->referenceLength >= 0 && dst->referenceLength <= dst->length);
		dst->chromosome = src->chromosome;
		dst->position = src->position;
		dst->strand = src->strand;
		dst->score = src->score;
	}
}

void AlignEntryFree(AlignEntry *aEntry)
{
	assert(aEntry->read!=NULL);
	free(aEntry->read);
	assert(aEntry->reference!=NULL);
	free(aEntry->reference);
	AlignEntryInitialize(aEntry);
}

void AlignEntryInitialize(AlignEntry *aEntry) 
{
	aEntry->readName[0] = '\0';
	aEntry->read=NULL;
	aEntry->reference=NULL;
	aEntry->length=0;
	aEntry->referenceLength=0;
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
