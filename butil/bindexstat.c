#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BLib.h"
#include "../blib/BError.h"
#include "../blib/RGIndex.h"
#include "bindexstat.h"

#define NUM_MISMATCHES_START 0
#define NUM_MISMATCHES_END 1

int main(int argc, char *argv[]) 
{
	FILE *fp=NULL;
	char indexFileName[MAX_FILENAME_LENGTH]="\0";
	char rgFileName[MAX_FILENAME_LENGTH]="\0";
	int algorithm;
	int numMismatchesStart = NUM_MISMATCHES_START;
	int numMismatchesEnd = NUM_MISMATCHES_END;

	if(argc == 4) {
		RGBinary rg;
		RGIndex index;

		strcpy(rgFileName, argv[1]);
		strcpy(indexFileName, argv[2]);
		algorithm = atoi(argv[3]);
		assert(0 <= algorithm && algorithm <= 2);


		/* Read in the rg binary file */
		RGBinaryReadBinary(&rg, rgFileName);

		/* Read the index */
		fprintf(stderr, "Reading in index from %s.\n",
				indexFileName);
		if(!(fp=fopen(indexFileName, "rb"))) {
			PrintError("bfixhash",
					indexFileName,
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}
		RGIndexRead(fp, &index, 1);
		fclose(fp);

		/* Get the desired stats */
		switch(algorithm) {
			case 0:
				/* One */
				PrintSummary(&index, &rg);
				break;
			case 1:
				/* The other */
				PrintHistogram(&index, 
						&rg, 
						numMismatchesStart,
						numMismatchesEnd,
						stderr);
				break;
			case 2:
				/* Both */
				PrintSummary(&index, &rg);
				PrintHistogram(&index, 
						&rg, 
						numMismatchesStart,
						numMismatchesEnd,
						stderr);
				break;
			default:
				PrintError("main",
						NULL,
						"Control reached a place not intended.",
						Exit,
						OutOfRange);
				break;
		}


		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Cleaning up.\n");
		/* Delete the index */
		RGIndexDelete(&index);
		/* Delete the rg */
		RGBinaryDelete(&rg);
	}
	else {
		fprintf(stdout, "Usage: bindexstat [OPTIONS]\n\t\t<reference genome file name>\n\t\t<index file name>\n\t\t<algorithm: 0-Summary 1-Histogram 2-Both>\n");
	}

	return 0;
}

void PrintSummary(RGIndex *index, RGBinary *rg)
{
	int64_t start, end;
	int64_t numEntries;
	int64_t max, min;
	long double sum;
	long double mean, variance;

	/* Get the mean */
	fprintf(stderr, "%s", BREAK_LINE);
	fprintf(stderr, "Getting the mean. Out of %lld, currently on:\n0",
			(long long int)index->length);
	numEntries = 0;
	start = 0;
	for(end=1;end < index->length;end++) {
		if(end%RGINDEX_ROTATE_NUM==0) {
			fprintf(stderr, "\r%lld", (long long int)end);
		}
		int cmp = RGIndexCompareAt(index, rg, start, end, 0);
		assert(cmp <= 0);
		if(cmp < 0) {
			numEntries++;
			/* Update */
			start = end;
		}
	}
	fprintf(stderr, "\r%lld\n", (long long int)end);
	/* Times two because we have both forward and reverse strands */
	mean = (index->length*2.0)/numEntries;

	/* Get the variance, max, and min */
	fprintf(stderr, "Getting the variance. Out of %lld, currently on:\n0",
			(long long int)index->length);
	min = UINT_MAX;
	max = -1;
	sum = 0.0;
	start = 0;
	for(end=1;end < index->length;end++) {
		if(end%RGINDEX_ROTATE_NUM==0) {
			fprintf(stderr, "\r%lld", (long long int)end);
		}
		int cmp = RGIndexCompareAt(index, rg, start, end, 0);
		assert(cmp <= 0);
		if(cmp < 0) {
			int val = ( (end-1) - start + 1);
			/* Get max */
			if(val > max) {
				max = val;
			}
			/* Get min */
			if(val < min) {
				min = val;
			}
			/* Update variance sum */
			sum += (val - mean)*(val - mean);
			/* Update */
			start = end;
		}
	}
	fprintf(stderr, "\r%lld\n", (long long int)end);
	variance = sum/numEntries;

	fprintf(stderr, "The mean was: %Lf\nThe variance was: %Lf\nThe max was: %lld\nThe min was: %lld\n",
			mean,
			variance,
			(long long int)max,
			(long long int)min);
}

void PrintHistogram(RGIndex *index, 
		RGBinary *rg,
		int numMismatchesStart,
		int numMismatchesEnd,
		FILE *fp)
{
	char *FnName = "PrintHistogram";
	int i;
	int64_t j;
	int* counts=NULL;
	int maxCount = 0;
	int64_t curIndex;
	int curChr, curPos, curNum;
	char read[SEQUENCE_LENGTH]="\0";
	char reverseRead[SEQUENCE_LENGTH]="\0";
	int readLength = index->totalLength;
	int returnLength, returnPosition;
	int64_t startIndex=-1;
	int64_t endIndex=-1;
	int64_t counter;
	int64_t numDifferent = 0;

	for(i=numMismatchesStart;i<=numMismatchesEnd;i++) {
		fprintf(fp, "Out of %lld, currently on %d mismatches:\n0",
				(long long int)index->length,
				i);
		curIndex = 0;
		counter = 0;
		while(curIndex < index->length) {
			if(counter >= RGINDEX_ROTATE_NUM) {
				fprintf(fp, "\r%lld", (long long int)curIndex);
				counter -= RGINDEX_ROTATE_NUM;
			}
			curNum = 0;
			/* Get the current chromosome and position from which
			 * to draw the read. */
			curChr = index->chromosomes[curIndex];
			curPos = index->positions[curIndex];
			/* Get the read */
			RGBinaryGetSequence(rg,
					curChr,
					curPos,
					FORWARD,
					0,
					read,
					readLength,
					&returnLength,
					&returnPosition);
			assert(returnLength == readLength);
			assert(returnPosition == curPos);
			/* Compute the number of places with this read */
			curNum += GetNumberOfMatches(index,
					rg,
					read,
					readLength,
					startIndex, 
					&endIndex,
					curIndex);
			/* Do the same thing but for the reverse compliment */
			GetReverseComplimentAnyCase(read,
					reverseRead,
					readLength);
			curNum += GetNumberOfMatches(index,
					rg,
					reverseRead,
					readLength,
					-1,
					NULL,
					-1);
			/* Add to our list.  We may have to reallocate this array */
			if(curNum > maxCount) {
				j = maxCount; /* Save previous length */
				/* Reallocate */
				maxCount = curNum;
				counts = realloc(counts, sizeof(int)*maxCount);
				if(NULL == counts) {
					PrintError(FnName,
							"counts",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}
				/* Initialize from j to maxCount */
				while(j<curNum) {
					counts[j] = 0;
					j++;
				}
			}
			/* Add the number of places to the appropriate counts */
			assert(curNum >= 0 && curNum <= maxCount);
			counts[curNum-1] += curNum;
			/* Update curIndex */
			counter += curNum;
			curIndex = endIndex + 1;
			/* Update stats */
			numDifferent++;
		}
		fprintf(fp, "\n");
		fprintf(fp, "Number of unique places was: %lld\nThe mean number of CALs was: %lld/%lld=%lf\n",
				(long long int)numDifferent,
				(long long int)numDifferent,
				(long long int)2*index->length,
				(double)(index->length*2.0)/numDifferent);
		fprintf(fp, "Found counts for %d mismatches ranging from %d to %d.\n",
				i,
				1,
				maxCount);
		for(j=0;j<maxCount;j++) {
			fprintf(fp, "%lld\t%d\n",
					(long long int)j,
					counts[j]);
		}
	}
}

int64_t GetNumberOfMatches(RGIndex *index,
		RGBinary *rg,
		char *read,
		int readLength,
		int64_t startIndex,
		int64_t *endIndex,
		int64_t curIndex)
{
	int64_t tmpEndIndex=-1;
	int64_t foundIndex=0;
	uint32_t hashIndex=0;

	/* Get the hash index for the read */
	hashIndex = RGIndexGetHashIndexFromRead(index, rg, read, readLength, 0);
	assert(hashIndex >= 0 && hashIndex < index->hashLength);
	assert(index->starts[hashIndex] != UINT_MAX);
	assert(index->ends[hashIndex] != UINT_MAX);
	assert(index->starts[hashIndex] >=0 && index->starts[hashIndex] < index->length);
	assert(index->ends[hashIndex] >=0 && index->ends[hashIndex] < index->length);

	/* Search the index using the bounds from the hash */
	foundIndex=RGIndexGetIndex(index,
			rg,
			index->starts[hashIndex],
			index->ends[hashIndex],
			read,
			&startIndex,
			&tmpEndIndex);
	assert(foundIndex==1);
	assert(startIndex < 0 || curIndex < 0 || startIndex == curIndex);
	if(endIndex != NULL) {
		(*endIndex) = tmpEndIndex;
	}
	/* Return the number of places with this read */
	return (tmpEndIndex - startIndex + 1);
}
