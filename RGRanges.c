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
int RGRangesRemoveDuplicates(RGRanges *r)
{
	int32_t i;
	int32_t prevIndex=0;

	int32_t numEntries = 0;

	if(r->numEntries > 0) {
		/* Quick sort the data structure */
		RGRangesQuickSort(r, 0, r->numEntries-1);

		/* Remove duplicates */
		prevIndex=0;
		numEntries = r->endIndex[0] - r->startIndex[0] + 1;
		for(i=1;i<r->numEntries;i++) {
			if(RGRangesCompareAtIndex(r, prevIndex, r, i)==0) {
				/* ignore */
			}
			else {
				prevIndex++;
				/* Copy to prevIndex (incremented) */
				RGRangesCopyAtIndex(r, prevIndex, r, i);
				numEntries += r->endIndex[i] - r->startIndex[i] + 1; 
			}
		}

		/* Reallocate pair */
		/* does not make sense if there are no entries */
		RGRangesReallocate(r, prevIndex+1);
	}
	return numEntries;
}

/* TODO */
void RGRangesCopyToRGMatch(RGRanges *r,
		RGIndex *index,
		RGMatch *m)
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
				if(FORWARD == r->strand[i]) {
					m->positions[counter] = index->positions[j] - r->offset[i];
				}
				else {
					m->positions[counter] = index->positions[j] + r->offset[i] - m->readLength;
				}
				/* Default position to one if the end of the read matched
				 * the beginning of the contig
				 * */
				if(m->positions[counter] <= 0) {
					m->positions[counter] = 1;
				}
				m->strands[counter] = r->strand[i];
				counter++;
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

		RGRangesCopyAtIndex(temp, 0, r, pivot);
		RGRangesCopyAtIndex(r, pivot, r, high);
		RGRangesCopyAtIndex(r, high, temp, 0);

		pivot = low;

		for(i=low;i<high;i++) {
			if(RGRangesCompareAtIndex(r, i, r, high) <= 0) {
				if(i!=pivot) {
					RGRangesCopyAtIndex(temp, 0, r, i);
					RGRangesCopyAtIndex(r, i, r, pivot);
					RGRangesCopyAtIndex(r, pivot, temp, 0);
				}
				pivot++;
			}
		}
		RGRangesCopyAtIndex(temp, 0, r, pivot);
		RGRangesCopyAtIndex(r, pivot, r, high);
		RGRangesCopyAtIndex(r, high, temp, 0);

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
	int i;
	int cmp[4] = {0,0,0,0};
	int top = 4;

	assert(indexOne >= 0 && indexOne < rOne->numEntries);
	assert(indexTwo >= 0 && indexTwo < rTwo->numEntries);

	cmp[0] = (rOne->startIndex[indexOne]<=rTwo->startIndex[indexTwo])?((rOne->startIndex[indexOne]<rTwo->startIndex[indexTwo])?-1:0):1;
	cmp[1] = (rOne->endIndex[indexOne]<=rTwo->endIndex[indexTwo])?((rOne->endIndex[indexOne]<rTwo->endIndex[indexTwo])?-1:0):1;
	cmp[2] = (rOne->strand[indexOne]<=rTwo->strand[indexTwo])?((rOne->strand[indexOne]<rTwo->strand[indexTwo])?-1:0):1;
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
void RGRangesAppend(RGRanges *dest, RGRanges *src)
{
	int32_t i, start;

	assert(src != dest);
	start = dest->numEntries;
	RGRangesReallocate(dest, dest->numEntries + src->numEntries);
	assert(dest->numEntries == start + src->numEntries);
	assert(start <= dest->numEntries);
	for(i=start;i<dest->numEntries;i++) {
		RGRangesCopyAtIndex(dest, i, src, i-start);
	}
}

/* TODO */
void RGRangesCopyAtIndex(RGRanges *dest, int32_t destIndex, RGRanges *src, int32_t srcIndex)
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
	r->strand = malloc(sizeof(char)*numEntries); 
	if(NULL == r->strand) {
		PrintError("RGRangesAllocate",
				"r->strand",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(r->offset==NULL);
	r->offset = malloc(sizeof(int32_t)*numEntries); 
	if(NULL == r->offset) {
		PrintError("RGRangesAllocate",
				"r->offset",
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
		r->strand = realloc(r->strand, sizeof(char)*numEntries); 
		if(numEntries > 0 && NULL == r->strand) {
			PrintError("RGRangesReallocate",
					"r->strand",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		r->offset = realloc(r->offset, sizeof(int32_t)*numEntries); 
		if(numEntries > 0 && NULL == r->offset) {
			PrintError("RGRangesReallocate",
					"r->offset",
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
		free(r->offset);
	}
	RGRangesInitialize(r);
}

void RGRangesInitialize(RGRanges *r)
{
	r->numEntries=0;
	r->startIndex=NULL;
	r->endIndex=NULL;
	r->strand=NULL;
	r->offset=NULL;
}
