#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <string.h>
#include "BLibDefinitions.h"
#include "BError.h"
#include "BMatch.h"
#include "BRanges.h"

/* TODO */
/* Note: this exploits the fact that the ranges in the indexes are non-overlapping.
 * Thus does not take the union of all ranges otherwise.
 * */
void BRangesRemoveDuplicates(BRanges *r)
{
	int32_t i;
	int32_t prevIndex=0;

	if(r->numEntries > 0) {
		/* Quick sort the data structure */
		BRangesQuickSort(r, 0, r->numEntries-1);

		/* Remove duplicates */
		prevIndex=0;
		for(i=1;i<r->numEntries;i++) {
			if(BRangesCompareAtIndex(r, prevIndex, r, i)==0) {
				/* ignore */
			}
			else {
				prevIndex++;
				/* Copy to prevIndex (incremented) */
				BRangesCopyAtIndex(r, i, r, prevIndex);
			}
		}

		/* Reallocate pair */
		/* does not make sense if there are no entries */
		BRangesReallocate(r, prevIndex+1);
	}
}

/* TODO */
void BRangesCopyToBMatch(BRanges *r,
		BIndex *index,
		BMatch *m)
{
	int64_t i;
	int64_t j;
	int64_t counter;

	assert(m->numEntries > 0);
	assert(r->numEntries > 0);

	if(r->numEntries > 0) {
		/* Copy over for each range */
		counter =0;
		for(i=0;i<r->numEntries;i++) {
			/* Copy over for the given range */
			for(j=r->startIndex[i];j<=r->endIndex[i];j++) {
				assert(j>=0 && j<index->length);
				assert(counter >= 0 && counter < m->numEntries);
				/* Get contig number */ 
				if(index->contigType == Contig_8) {
					m->contigs[counter] = index->contigs_8[j];
				}
				else {
					m->contigs[counter] = index->contigs_32[j];
				}
				/* Adjust position with the offset */
				m->positions[counter] = (uint32_t)(index->positions[j] - r->offset[i]);
				m->strand[counter] = r->strand[i];
				counter++;
			}
		}
	}
}

/* TODO */
void BRangesQuickSort(BRanges *r, int32_t low, int32_t high)
{
	int32_t i;
	int32_t pivot=-1;
	BRanges *temp;

	if(low < high) {
		/* Allocate memory for the temp used for swapping */
		temp=malloc(sizeof(BRanges));
		BRangesInitialize(temp);
		if(NULL == temp) {
			PrintError("BRangesQuickSort",
					"temp",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		BRangesAllocate(temp, 1);

		pivot = (low+high)/2;

		BRangesCopyAtIndex(r, pivot, temp, 0);
		BRangesCopyAtIndex(r, high, r, pivot);
		BRangesCopyAtIndex(temp, 0, r, high);

		pivot = low;

		for(i=low;i<high;i++) {
			if(BRangesCompareAtIndex(r, i, r, high) <= 0) {
				if(i!=pivot) {
					BRangesCopyAtIndex(r, i, temp, 0);
					BRangesCopyAtIndex(r, pivot, r, i);
					BRangesCopyAtIndex(temp, 0, r, pivot);
				}
				pivot++;
			}
		}
		BRangesCopyAtIndex(r, pivot, temp, 0);
		BRangesCopyAtIndex(r, high, r, pivot);
		BRangesCopyAtIndex(temp, 0, r, high);

		/* Free temp before the recursive call, otherwise we have a worst
		 * case of O(n) space (NOT IN PLACE) 
		 * */
		BRangesFree(temp);
		free(temp);
		temp=NULL;

		BRangesQuickSort(r, low, pivot-1);
		BRangesQuickSort(r, pivot+1, high);
	}
}

/* TODO */
int32_t BRangesCompareAtIndex(BRanges *rOne, int32_t indexOne, BRanges *rTwo, int32_t indexTwo) 
{
	int i;
	int cmp[4] = {0,0,0,0};
	int top = 4;

	assert(indexOne >= 0 && indexOne < rOne->numEntries);
	assert(indexTwo >= 0 && indexTwo < rTwo->numEntries);

	cmp[0] = (rOne->startIndex[indexOne]<=rTwo->startIndex[indexTwo])?((rOne->startIndex[indexOne]<rTwo->startIndex[indexTwo])?-1:0):1;
	cmp[1] = (rOne->endIndex[indexOne]<=rTwo->endIndex[indexTwo])?((rOne->endIndex[indexOne]<rTwo->endIndex[indexTwo])?-1:0):1;
	cmp[2] = (rOne->strand[indexOne]<=rTwo->startIndex[indexTwo])?((rOne->startIndex[indexOne]<rTwo->startIndex[indexTwo])?-1:0):1;
	cmp[3] = (rOne->offset[indexOne]<=rTwo->startIndex[indexTwo])?((rOne->startIndex[indexOne]<rTwo->startIndex[indexTwo])?-1:0):1;

	for(i=0;i<top;i++) {
		if(cmp[i] < 0) {
			return cmp[i];
		}
		else if(cmp[i] > 0) {
			return cmp[i];
		}
	}

	return 0;
}

/* TODO */
void BRangesAppend(BRanges *src, BRanges *dest)
{
	int32_t i, start;

	assert(src != dest);
	start = dest->numEntries;
	BRangesReallocate(dest, dest->numEntries + src->numEntries);
	assert(dest->numEntries == start + src->numEntries);
	assert(start <= dest->numEntries);
	for(i=start;i<dest->numEntries;i++) {
		BRangesCopyAtIndex(src, i-start, dest, i);
	}
}

/* TODO */
void BRangesCopyAtIndex(BRanges *src, int32_t srcIndex, BRanges *dest, int32_t destIndex)
{
	assert(srcIndex >= 0 && srcIndex < src->numEntries);
	assert(destIndex >= 0 && destIndex < dest->numEntries);

	if(src != dest || srcIndex != destIndex) {
		dest->startIndex[destIndex] = src->startIndex[srcIndex];
		dest->endIndex[destIndex] = src->endIndex[srcIndex];
		dest->strand[destIndex] = src->strand[srcIndex];
		dest->offset[destIndex] = src->offset[srcIndex];
	}
}

void BRangesAllocate(BRanges *r, int32_t numEntries)
{
	assert(r->numEntries==0);
	r->numEntries = numEntries;
	assert(r->startIndex==NULL);
	r->startIndex = malloc(sizeof(int64_t)*numEntries); 
	if(NULL == r->startIndex) {
		PrintError("BRangesAllocate",
				"r->startIndex",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(r->endIndex==NULL);
	r->endIndex = malloc(sizeof(int64_t)*numEntries); 
	if(NULL == r->endIndex) {
		PrintError("BRangesAllocate",
				"r->endIndex",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(r->strand==NULL);
	r->strand = malloc(sizeof(int8_t)*numEntries); 
	if(NULL == r->strand) {
		PrintError("BRangesAllocate",
				"r->strand",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(r->offset==NULL);
	r->offset = malloc(sizeof(int32_t)*numEntries); 
	if(NULL == r->offset) {
		PrintError("BRangesAllocate",
				"r->offset",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
}

void BRangesReallocate(BRanges *r, int32_t numEntries)
{
	if(numEntries > 0) {
		r->numEntries = numEntries;
		r->startIndex = realloc(r->startIndex, sizeof(int64_t)*numEntries); 
		if(numEntries > 0 && NULL == r->startIndex) {
			PrintError("BRangesReallocate",
					"r->startIndex",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		r->endIndex = realloc(r->endIndex, sizeof(int64_t)*numEntries); 
		if(numEntries > 0 && NULL == r->endIndex) {
			PrintError("BRangesReallocate",
					"r->endIndex",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		r->strand = realloc(r->strand, sizeof(int8_t)*numEntries); 
		if(numEntries > 0 && NULL == r->strand) {
			PrintError("BRangesReallocate",
					"r->strand",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		r->offset = realloc(r->offset, sizeof(int32_t)*numEntries); 
		if(numEntries > 0 && NULL == r->offset) {
			PrintError("BRangesReallocate",
					"r->offset",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
	}
	else {
		BRangesFree(r);
	}
}

void BRangesFree(BRanges *r) 
{
	if(r->numEntries>0) {
		free(r->startIndex);
		free(r->endIndex);
		free(r->strand);
		free(r->offset);
	}
	BRangesInitialize(r);
}

void BRangesInitialize(BRanges *r)
{
	r->numEntries=0;
	r->startIndex=NULL;
	r->endIndex=NULL;
	r->strand=NULL;
	r->offset=NULL;
}
