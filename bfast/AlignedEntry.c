#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <zlib.h>
#include "BError.h"
#include "BLib.h"
#include "AlignedEntry.h"

/* TODO */
int32_t AlignedEntryPrint(AlignedEntry *a,
		gzFile outputFP,
		int32_t space)
{
	assert(NULL != a->read);
	assert(NULL != a->reference);
	assert(space == NTSpace ||
			(space == ColorSpace && NULL != a->colorError));
	
	if(gzwrite64(outputFP, &a->contigNameLength, sizeof(int32_t))!=sizeof(int32_t)||
			gzwrite64(outputFP, a->contigName, sizeof(char)*a->contigNameLength)!=sizeof(char)*a->contigNameLength||
			gzwrite64(outputFP, &a->contig, sizeof(uint32_t))!=sizeof(uint32_t)||
			gzwrite64(outputFP, &a->position, sizeof(uint32_t))!=sizeof(uint32_t)||
			gzwrite64(outputFP, &a->strand, sizeof(char))!=sizeof(char)||
			gzwrite64(outputFP, &a->score, sizeof(double))!=sizeof(double)||
			gzwrite64(outputFP, &a->mappingQuality, sizeof(int32_t))!=sizeof(int32_t)||
			gzwrite64(outputFP, &a->referenceLength, sizeof(uint32_t))!=sizeof(uint32_t)||
			gzwrite64(outputFP, &a->length, sizeof(uint32_t))!=sizeof(uint32_t)||
			gzwrite64(outputFP, a->read, sizeof(char)*a->length)!=sizeof(char)*a->length||
			gzwrite64(outputFP, a->reference, sizeof(char)*a->length)!=sizeof(char)*a->length) {
		return EOF;
	}
	if(ColorSpace==space) {
		if(gzwrite64(outputFP, a->colorError, sizeof(char)*a->length)!=sizeof(char)*a->length) {
			return EOF;
		}
	}
	
	return 1;
}

/* TODO */
int32_t AlignedEntryPrintText(AlignedEntry *a,
		FILE *outputFP,
		int32_t space)
{
	assert(NULL != a->read);
	assert(NULL != a->reference);
	assert(space == NTSpace ||
			(space == ColorSpace && NULL != a->colorError));

	if(fprintf(outputFP, "%s\t%u\t%u\t%c\t%lf\t%d\t%u\t%u\n",
				a->contigName,
				a->contig,
				a->position,
				a->strand,
				a->score,
				a->mappingQuality,
				a->referenceLength,
				a->length) < 0) {
		return EOF;
	}

	/* Print32_t the reference and read alignment */
	if(fprintf(outputFP, "%s\n%s\n",
				a->reference,
				a->read) < 0) {
		return EOF;
	}

	/* Print32_t the color errors if necessary */
	if(space == ColorSpace) {
		if(fprintf(outputFP, "%s\n",
					a->colorError) < 0) {
			return EOF;
		}
	}

	return 1;
}

/* TODO */
int32_t AlignedEntryRead(AlignedEntry *a,
		gzFile inputFP,
		int32_t space)
{
	char *FnName = "AlignedEntryRead";

	assert(NULL != a);

	/* Allocate memory for the alignment */
	if(a->read == NULL) {
		a->read = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == a->read) {
			PrintError(FnName, "a->read", "Could not allocate memory", Exit, MallocMemory);
		}
	}
	if(a->reference == NULL) {
		a->reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == a->reference) {
			PrintError(FnName, "a->reference", "Could not allocate memory", Exit, MallocMemory);
		}
	}
	if(space == ColorSpace) {
		if(a->colorError == NULL) {
			a->colorError = malloc(sizeof(char)*SEQUENCE_LENGTH);
			if(NULL == a->colorError) {
				PrintError(FnName, "a->colorError", "Could not allocate memory", Exit, MallocMemory);
			}
		}
	}

	if(gzread64(inputFP, &a->contigNameLength, sizeof(int32_t))!=sizeof(int32_t)) {
		return EOF;
	}
	/* Copy over contig name */
	a->contigName = malloc(sizeof(char)*(a->contigNameLength+1));
	if(NULL==a->contigName) {
		PrintError(FnName, "a->contigName", "Could not allocate memory", Exit, MallocMemory);
	}
	if(gzread64(inputFP, a->contigName, sizeof(char)*a->contigNameLength)!=sizeof(char)*a->contigNameLength||
			gzread64(inputFP, &a->contig, sizeof(uint32_t))!=sizeof(uint32_t)||
			gzread64(inputFP, &a->position, sizeof(uint32_t))!=sizeof(uint32_t)||
			gzread64(inputFP, &a->strand, sizeof(char))!=sizeof(char)||
			gzread64(inputFP, &a->score, sizeof(double))!=sizeof(double)||
			gzread64(inputFP, &a->mappingQuality, sizeof(int32_t))!=sizeof(int32_t)||
			gzread64(inputFP, &a->referenceLength, sizeof(uint32_t))!=sizeof(uint32_t)||
			gzread64(inputFP, &a->length, sizeof(uint32_t))!=sizeof(uint32_t)||
			gzread64(inputFP, a->read, sizeof(char)*a->length)!=sizeof(char)*a->length||
			gzread64(inputFP, a->reference, sizeof(char)*a->length)!=sizeof(char)*a->length) {
		return EOF;
	}
	/* Add the null terminator to strings */
	a->contigName[a->contigNameLength]='\0';
	a->read[a->length]='\0';
	a->reference[a->length]='\0';
	if(ColorSpace==space) {
		if(gzread64(inputFP, a->colorError, sizeof(char)*a->length)!=sizeof(char)*a->length) {
			return EOF;
		}
		a->colorError[a->length]='\0';
	}

	/* Reallocate to conserve memory */
	assert(a->length > 0);
	a->read = realloc(a->read, sizeof(char)*(a->length+1));
	if(NULL == a->read) {
		PrintError(FnName, "a->read", "Could not reallocate memory", Exit, ReallocMemory);
	}
	/* Reference */
	a->reference = realloc(a->reference, sizeof(char)*(a->length+1));
	if(NULL == a->reference) {
		PrintError(FnName, "a->reference", "Could not reallocate memory", Exit, ReallocMemory);
	}
	/* Color error, if necessary */
	if(space == ColorSpace) {
		a->colorError = realloc(a->colorError, sizeof(char)*(a->length+1));
		if(NULL == a->colorError) {
			PrintError(FnName, "a->colorError", "Could not reallocate memory", Exit, ReallocMemory);
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
int32_t AlignedEntryReadText(AlignedEntry *a,
		FILE *inputFP,
		int32_t space)
{
	char *FnName = "AlignedEntryReadText";
	char tempContigName[MAX_CONTIG_NAME_LENGTH]="\0";

	assert(NULL != a);

	/* Allocate memory for the alignment */
	if(a->read == NULL) {
		a->read = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == a->read) {
			PrintError(FnName, "a->read", "Could not allocate memory", Exit, MallocMemory);
		}
	}
	if(a->reference == NULL) {
		a->reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == a->reference) {
			PrintError(FnName, "a->reference", "Could not allocate memory", Exit, MallocMemory);
		}
	}
	if(space == ColorSpace) {
		if(a->colorError == NULL) {
			a->colorError = malloc(sizeof(char)*SEQUENCE_LENGTH);
			if(NULL == a->colorError) {
				PrintError(FnName, "a->colorError", "Could not allocate memory", Exit, MallocMemory);
			}
		}
	}

	if(fscanf(inputFP, "%s\t%u\t%u\t%c\t%lf\t%d\t%u\t%u\n",
				tempContigName,
				&a->contig,
				&a->position,
				&a->strand,
				&a->score,
				&a->mappingQuality,
				&a->referenceLength,
				&a->length) < 0) {
		return EOF;
	}

	/* Copy over contig name */
	a->contigNameLength = (int)strlen(tempContigName);
	a->contigName = malloc(sizeof(char)*(a->contigNameLength+1));
	if(NULL==a->contigName) {
		PrintError(FnName, "a->contigName", "Could not allocate memory", Exit, MallocMemory);
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

	/* Reallocate to conserve memory */
	assert(a->length > 0);
	a->read = realloc(a->read, sizeof(char)*(a->length+1));
	if(NULL == a->read) {
		PrintError(FnName, "a->read", "Could not reallocate memory", Exit, ReallocMemory);
	}
	/* Reference */
	a->reference = realloc(a->reference, sizeof(char)*(a->length+1));
	if(NULL == a->reference) {
		PrintError(FnName, "a->reference", "Could not reallocate memory", Exit, ReallocMemory);
	}
	/* Color error, if necessary */
	if(space == ColorSpace) {
		a->colorError = realloc(a->colorError, sizeof(char)*(a->length+1));
		if(NULL == a->colorError) {
			PrintError(FnName, "a->colorError", "Could not reallocate memory", Exit, ReallocMemory);
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
/* Log-n space */
void AlignedEntryQuickSort(AlignedEntry **a,
		int32_t low,
		int32_t high,
		int32_t sortOrder,
		int32_t showPercentComplete,
		double *curPercent,
		int32_t total)
{
	char *FnName = "AlignedEntryQuickSort";
	int32_t i;
	int32_t pivot=-1;
	AlignedEntry *temp=NULL;

	if(low < high) {

		if(high - low + 1 <= ALIGNEDENTRY_SHELL_SORT_MAX) {
			AlignedEntryShellSort(a, low, high, sortOrder);
			return;
		}

		/* Allocate memory for the temp used for swapping */
		temp=malloc(sizeof(AlignedEntry));
		if(NULL == temp) {
			PrintError(FnName, "temp", "Could not allocate temp", Exit, MallocMemory);
		}
		AlignedEntryInitialize(temp);

		pivot = AlignedEntryGetPivot((*a),
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

		AlignedEntryCopyAtIndex(temp, 0, (*a), pivot);
		AlignedEntryCopyAtIndex((*a), pivot, (*a), high);
		AlignedEntryCopyAtIndex((*a), high, temp, 0);

		pivot = low;

		for(i=low;i<high;i++) {
			if(AlignedEntryCompareAtIndex((*a), i, (*a), high, sortOrder) <= 0) {
				if(i!=pivot) {
					AlignedEntryCopyAtIndex(temp, 0, (*a), i);
					AlignedEntryCopyAtIndex((*a), i, (*a), pivot);
					AlignedEntryCopyAtIndex((*a), pivot, temp, 0);
				}
				pivot++;
			}
		}
		AlignedEntryCopyAtIndex(temp, 0, (*a), pivot);
		AlignedEntryCopyAtIndex((*a), pivot, (*a), high);
		AlignedEntryCopyAtIndex((*a), high, temp, 0);

		/* Free temp before the recursive call, otherwise we have a worst
		 * case of O(n) space (NOT IN PLACE) 
		 * */
		AlignedEntryFree(temp);
		free(temp);
		temp=NULL;

		AlignedEntryQuickSort(a, low, pivot-1, sortOrder, showPercentComplete, curPercent, total);
		if(showPercentComplete == 1 && VERBOSE >= 0) {
			assert(NULL!=curPercent);
			if((*curPercent) < 100.0*((double)pivot)/total) {
				while((*curPercent) < 100.0*((double)pivot)/total) {
					(*curPercent) += SORT_ROTATE_INC;
				}
				PrintPercentCompleteShort((*curPercent));
			}
		}
		AlignedEntryQuickSort(a, pivot+1, high, sortOrder, showPercentComplete, curPercent, total);
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

void AlignedEntryShellSort(AlignedEntry **a,
		int32_t low,
		int32_t high,
		int32_t sortOrder)
{
	char *FnName = "AlignedEntryShellSort";
	int32_t i, j, inc;
	AlignedEntry *temp=NULL;

	inc = ROUND((high - low + 1) / 2);

	/* Allocate memory for the temp used for swapping */
	temp=malloc(sizeof(AlignedEntry));
	if(NULL == temp) {
		PrintError(FnName, "temp", "Could not allocate temp", Exit, MallocMemory);
	}
	AlignedEntryInitialize(temp);

	while(0 < inc) {
		for(i=inc + low;i<=high;i++) {
			AlignedEntryCopyAtIndex(temp, 0, (*a), i);
			j = i;
			while(inc + low <= j && AlignedEntryCompareAtIndex(temp, 0, (*a), j - inc, sortOrder) < 0) {
				AlignedEntryCopyAtIndex((*a), j, (*a), j - inc);
				j -= inc;
			}
			AlignedEntryCopyAtIndex((*a), j, temp, 0);
		}
		inc = ROUND(inc / SHELL_SORT_GAP_DIVIDE_BY);
	}
	AlignedEntryFree(temp);
	free(temp);
	temp=NULL;

}

/* TODO */
int32_t AlignedEntryCompareAtIndex(AlignedEntry *a, int32_t indexA, AlignedEntry *b, int32_t indexB, int32_t sortOrder)
{
	return AlignedEntryCompare(&(a[indexA]), &(b[indexB]), sortOrder);
}

/* TODO */
int32_t AlignedEntryCompare(AlignedEntry *a, AlignedEntry *b, int32_t sortOrder)
{
	int32_t cmp[5];
	int32_t i;
	int32_t top;

	if(sortOrder == AlignedEntrySortByAll) {

		/* If there are multiple alignments to the same starting chr/pos/strand with the same score,
		 * this will pick ensure that we will only pick one of them.
		 * */
		cmp[0] = (a->contig <= b->contig)?((a->contig<b->contig)?-1:0):1;
		cmp[1] = (a->position <= b->position)?((a->position<b->position)?-1:0):1;
		cmp[2] = (a->strand <= b->strand)?((a->strand<b->strand)?-1:0):1;
		cmp[3] = (a->score <= b->score)?((a->score<b->score)?-1:0):1;

		top = 4;
	}
	else {
		assert(sortOrder == AlignedEntrySortByContigPos);
		cmp[0] = (a->contig <= b->contig)?((a->contig<b->contig)?-1:0):1;
		cmp[1] = (a->position <= b->position)?((a->position<b->position)?-1:0):1;

		top = 2;
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
void AlignedEntryCopyAtIndex(AlignedEntry *dest, int32_t destIndex, AlignedEntry *src, int32_t srcIndex)
{
	if(dest != src || srcIndex != destIndex) {
		AlignedEntryCopy(&(dest[destIndex]), &(src[srcIndex]));
	}
}

/* TODO */
void AlignedEntryCopy(AlignedEntry *dest, AlignedEntry *src)
{
	char *FnName = "AlignedEntryCopy";
	if(src != dest) {
		/* Contig name length */
		assert(src->contigNameLength > 0);
		dest->contigNameLength = src->contigNameLength;
		/* Contig name */
		dest->contigName = realloc(dest->contigName, sizeof(char)*(dest->contigNameLength+1));
		if(NULL == dest->contigName) {
			PrintError(FnName, "dest->contigName", "Could not reallocate memory", Exit, ReallocMemory);
		}
		assert(src->contigName != NULL);
		strcpy(dest->contigName, src->contigName);
		/* Read */
		assert(src->length > 0);
		dest->read = realloc(dest->read, sizeof(char)*(src->length+1));
		if(NULL == dest->read) {
			PrintError(FnName, "dest->read", "Could not reallocate memory", Exit, ReallocMemory);
		}
		assert(src->read != NULL);
		strcpy(dest->read, src->read);
		/* Reference */
		dest->reference = realloc(dest->reference, sizeof(char)*(src->length+1));
		if(NULL == dest->reference) {
			PrintError(FnName, "dest->reference", "Could not reallocate memory", Exit, ReallocMemory);
		}
		assert(src->reference!= NULL);
		strcpy(dest->reference, src->reference);
		/* Color error, if necessary */
		if(src->colorError != NULL) {
			assert(src->length > 0);
			dest->colorError = realloc(dest->colorError, sizeof(char)*(src->length+1));
			if(NULL == dest->colorError) {
				PrintError(FnName, "dest->colorError", "Could not reallocate memory", Exit, ReallocMemory);
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
		dest->mappingQuality = src->mappingQuality;
	}
}

void AlignedEntryFree(AlignedEntry *a)
{
	free(a->contigName);
	free(a->read);
	free(a->reference);
	free(a->colorError);
	AlignedEntryInitialize(a);
}

void AlignedEntryInitialize(AlignedEntry *a) 
{
	a->contigNameLength=0;
	a->contigName=NULL;
	a->contig=0;
	a->position=0;
	a->strand=0;
	a->score=0.0;
	a->mappingQuality=0;
	a->referenceLength=0;
	a->length=0;
	a->read=NULL;
	a->reference=NULL;
	a->colorError=NULL;
}

/* TODO */
/* Debugging function */
void AlignedEntryCheckReference(AlignedEntry *a, RGBinary *rg, int32_t space)
{
	char *FnName = "AlignedEntryCheckReference";
	int32_t i;
	int32_t curPos;
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
				AlignedEntryPrintText(a, stderr, space);
				PrintError(FnName, NULL, "Reference in the align entry does not match the reference genome", Exit, OutOfRange);
			}
			curPos++;
		}
	}
}

int32_t AlignedEntryGetPivot(AlignedEntry *a,
		int32_t sortOrder,
		int32_t low,
		int32_t high) 
{
	int32_t cmp[3];
	int32_t pivot = (low + high)/2;
	cmp[0] = AlignedEntryCompareAtIndex(a, low, a, pivot, sortOrder); 
	cmp[1] = AlignedEntryCompareAtIndex(a, low, a, high, sortOrder); 
	cmp[2] = AlignedEntryCompareAtIndex(a, pivot, a, high, sortOrder); 

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

void AlignedEntryUpdateAlignment(AlignedEntry *a,
		uint32_t position,
		double score,
		int32_t referenceLength,
		int32_t length,
		char *read,
		char *reference,
		char *colorError) 
{
	char *FnName="AlignedEntryAllocate";

	a->position = position;
	a->score = score;
	a->referenceLength = referenceLength;
	a->length = length;
	a->read = strdup(read);
	if(NULL == a->read) {
		PrintError(FnName, "a->read", "Could not allocate memory", Exit, MallocMemory);
	}
	a->reference = strdup(reference);
	if(NULL == a->reference) {
		PrintError(FnName, "a->reference", "Could not allocate memory", Exit, MallocMemory);
	}
	if(NULL != colorError) {
		a->colorError = strdup(colorError);
		if(NULL == a->colorError) {
			PrintError(FnName, "a->colorError", "Could not allocate memory", Exit, MallocMemory);
		}
	}
}
