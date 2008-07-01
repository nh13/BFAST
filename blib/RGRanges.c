#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <string.h>
#include "BLibDefinitions.h"
#include "BError.h"
#include "RGMatch.h"
#include "RGRanges.h"

/* TODO */
/* Note: this exploits the fact that the ranges in the indexes are non-overlapping.
 * Thus does not take the union of all ranges otherwise.
 * */
void RGRangesRemoveDuplicates(RGRanges *r)
{
	int32_t i;
	int32_t prevIndex=0;

	if(r->numEntries > 0) {
		/* Quick sort the data structure */
		RGRangesQuickSort(r, 0, r->numEntries-1);

		/* Remove duplicates */
		prevIndex=0;
		for(i=1;i<r->numEntries;i++) {
			if(RGRangesCompareAtIndex(r, prevIndex, r, i)==0) {
				/* ignore */
			}
			else {
				prevIndex++;
				/* Copy to prevIndex (incremented) */
				RGRangesCopyAtIndex(r, i, r, prevIndex);
			}
		}

		/* Reallocate pair */
		/* does not make sense if there are no entries */
		RGRangesReallocate(r, prevIndex+1);

	}
}

/* TODO */
void RGRangesCopyToRGMatch(RGRanges *r,
		RGIndex *index,
		RGMatch *m)
{
	int32_t i, k;
	int64_t j;

	if(r->numEntries > 0) {

		/* Allocate memory for the matches */
		RGMatchAllocate(m, r->numEntries);

		/* Copy over for each range */
		for(i=0, k=0;i<r->numEntries;i++) {
			/* Copy over for the given range */
			for(j=r->startIndex[i];j<=r->endIndex[i];j++) {
				m->positions[k] = index->positions[j];
				m->chromosomes[k] = index->chromosomes[j];
				m->strand[k] = r->strand[i];
				k++;
			}
		}
	}
}

/* TODO */
void RGRangesQuickSort(RGRanges *r, int32_t low, int32_t high)
{
	int32_t i;
	int32_t pivot=-1;
	RGRanges *temp;

	if(low < high) {
		/* Allocate memory for the temp used for swapping */
		temp=malloc(sizeof(RGRanges));
		RGRangesInitialize(temp);
		if(NULL == temp) {
			PrintError("RGRangesQuickSort",
					"temp",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		RGRangesAllocate(temp, 1);

		pivot = (low+high)/2;

		RGRangesCopyAtIndex(r, pivot, temp, 0);
		RGRangesCopyAtIndex(r, high, r, pivot);
		RGRangesCopyAtIndex(temp, 0, r, high);

		pivot = low;

		for(i=low;i<high;i++) {
			if(RGRangesCompareAtIndex(r, i, r, high) <= 0) {
				if(i!=pivot) {
					RGRangesCopyAtIndex(r, i, temp, 0);
					RGRangesCopyAtIndex(r, pivot, r, i);
					RGRangesCopyAtIndex(temp, 0, r, pivot);
				}
				pivot++;
			}
		}
		RGRangesCopyAtIndex(r, pivot, temp, 0);
		RGRangesCopyAtIndex(r, high, r, pivot);
		RGRangesCopyAtIndex(temp, 0, r, high);

		/* Free temp before the recursive call, otherwise we have a worst
		 * case of O(n) space (NOT IN PLACE) 
		 * */
		RGRangesFree(temp);
		free(temp);
		temp=NULL;

		RGRangesQuickSort(r, low, pivot-1);
		RGRangesQuickSort(r, pivot+1, high);
	}
}

/* TODO */
int32_t RGRangesCompareAtIndex(RGRanges *rOne, int32_t indexOne, RGRanges *rTwo, int32_t indexTwo) 
{
	assert(indexOne >= 0 && indexOne < rOne->numEntries);
	assert(indexTwo >= 0 && indexTwo < rTwo->numEntries);
	if(rOne->startIndex[indexOne] < rTwo->startIndex[indexTwo] ||
			(rOne->startIndex[indexOne] == rTwo->startIndex[indexTwo] && rOne->endIndex[indexOne] < rTwo->endIndex[indexTwo]) ||
			(rOne->startIndex[indexOne] == rTwo->startIndex[indexTwo] && rOne->endIndex[indexOne] == rTwo->endIndex[indexTwo] && rOne->strand[indexOne] < rTwo->strand[indexTwo])) {
		return -1;
	}
	else if(rOne->startIndex[indexOne] ==  rTwo->startIndex[indexTwo] && rOne->endIndex[indexOne] == rTwo->endIndex[indexTwo] && rOne->strand[indexOne] == rTwo->strand[indexTwo]) {
		return 0;
	}
	else {
		return 1;
	}
}

/* TODO */
void RGRangesAppend(RGRanges *src, RGRanges *dest)
{
	int32_t i, start;

	assert(src != dest);
	start = dest->numEntries;
	RGRangesReallocate(dest, dest->numEntries + src->numEntries);
	assert(dest->numEntries == start + src->numEntries);
	assert(start <= dest->numEntries);
	for(i=start;i<dest->numEntries;i++) {
		RGRangesCopyAtIndex(src, i-start, dest, i);
	}
}

/* TODO */
void RGRangesCopyAtIndex(RGRanges *src, int32_t srcIndex, RGRanges *dest, int32_t destIndex)
{
	assert(srcIndex >= 0 && srcIndex < src->numEntries);
	assert(destIndex >= 0 && destIndex < dest->numEntries);

	dest->startIndex[destIndex] = src->startIndex[srcIndex];
	dest->endIndex[destIndex] = src->endIndex[srcIndex];
	dest->strand[destIndex] = src->strand[srcIndex];
}

void RGRangesAllocate(RGRanges *r, int32_t numEntries)
{
	assert(r->numEntries==0);
	r->numEntries = numEntries;
	assert(r->startIndex==NULL);
	r->startIndex = malloc(sizeof(int64_t)*numEntries); 
	if(NULL == r->startIndex) {
		PrintError("RGRangesAllocate",
				"r->startIndex",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(r->endIndex==NULL);
	r->endIndex = malloc(sizeof(int64_t)*numEntries); 
	if(NULL == r->endIndex) {
		PrintError("RGRangesAllocate",
				"r->endIndex",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(r->strand==NULL);
	r->strand = malloc(sizeof(int8_t)*numEntries); 
	if(NULL == r->strand) {
		PrintError("RGRangesAllocate",
				"r->strand",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
}

void RGRangesReallocate(RGRanges *r, int32_t numEntries)
{
	if(numEntries > 0) {
		r->numEntries = numEntries;
		r->startIndex = realloc(r->startIndex, sizeof(int64_t)*numEntries); 
		if(numEntries > 0 && NULL == r->startIndex) {
			PrintError("RGRangesReallocate",
					"r->startIndex",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		r->endIndex = realloc(r->endIndex, sizeof(int64_t)*numEntries); 
		if(numEntries > 0 && NULL == r->endIndex) {
			PrintError("RGRangesReallocate",
					"r->endIndex",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		r->strand = realloc(r->strand, sizeof(int8_t)*numEntries); 
		if(numEntries > 0 && NULL == r->strand) {
			PrintError("RGRangesReallocate",
					"r->strand",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
	}
	else {
		RGRangesFree(r);
	}
}

void RGRangesFree(RGRanges *r) 
{
	if(r->numEntries>0) {
		free(r->startIndex);
		free(r->endIndex);
		free(r->strand);
	}
	RGRangesInitialize(r);
}

void RGRangesInitialize(RGRanges *r)
{
	r->numEntries=0;
	r->startIndex=NULL;
	r->endIndex=NULL;
	r->strand=NULL;
}
