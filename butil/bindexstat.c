#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BLib.h"
#include "../blib/BError.h"
#include "../blib/RGIndex.h"
#include "../blib/RGMatch.h"
#include "../blib/RGReads.h"
#include "bindexstat.h"

#define BINDEXSTAT_ROTATE_NUM 100000
#define NUM_MISMATCHES_START 0
#define NUM_MISMATCHES_END 4

typedef struct {
	int **counts;
	int *maxCount;
} Counts;

int main(int argc, char *argv[]) 
{
	FILE *fp=NULL;
	char indexFileName[MAX_FILENAME_LENGTH]="\0";
	char rgFileName[MAX_FILENAME_LENGTH]="\0";
	int algorithm;
	int numMismatchesStart = NUM_MISMATCHES_START;
	int numMismatchesEnd = NUM_MISMATCHES_END;

	if(argc == 5) {
		RGBinary rg;
		RGIndex index;

		strcpy(rgFileName, argv[1]);
		strcpy(indexFileName, argv[2]);
		algorithm = atoi(argv[3]);
		numMismatchesEnd = atoi(argv[4]);
		assert(0 <= algorithm && algorithm <= 2);
		assert(numMismatchesEnd >= numMismatchesStart);


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
				fprintf(stderr, "%s", BREAK_LINE);
				PrintSummary(&index, &rg);
				fprintf(stderr, "%s", BREAK_LINE);
				break;
			case 1:
				/* The other */
				fprintf(stderr, "%s", BREAK_LINE);
				PrintHistogram(&index, 
						&rg, 
						numMismatchesStart,
						numMismatchesEnd,
						stdout);
				fprintf(stderr, "%s", BREAK_LINE);
				break;
			case 2:
				/* Both */
				fprintf(stderr, "%s", BREAK_LINE);
				PrintSummary(&index, &rg);
				fprintf(stderr, "%s", BREAK_LINE);
				PrintHistogram(&index, 
						&rg, 
						numMismatchesStart,
						numMismatchesEnd,
						stdout);
				fprintf(stderr, "%s", BREAK_LINE);
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
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully!\n");
		fprintf(stderr, "%s", BREAK_LINE);
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

	fprintf(stderr, "The mean was: %Lf\nThe variance was: %Lf\nThe max was: %lld\nThe min was: %lld\nThe number of unique reads was: %lld\n",
			mean,
			variance,
			(long long int)max,
			(long long int)min,
			(long long int)numEntries);
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
	int64_t curIndex, nextIndex;
	int curChr, curPos, curNum;
	char read[SEQUENCE_LENGTH]="\0";
	int readLength = index->totalLength;
	int returnLength, returnPosition;
	int64_t counter;
	int64_t numDifferent = 0;
	int64_t numReadsNoMismatches = 0;
	Counts c; /* Used to store the histogram data */
	RGMatch m; /* Used to store the matches returned */
	int offsets[1] = {0}; /* Only use one offset */
	int *curCounts = NULL;

	/* Not implemented for numMismatchesStart > 0 */
	assert(numMismatchesStart == 0);

	/* Allocate memory to hold histogram data */
	c.counts = malloc(sizeof(int*)*(numMismatchesEnd - numMismatchesStart + 1));
	if(NULL == c.counts) {
		PrintError(FnName,
				"c.counts",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	c.maxCount = malloc(sizeof(int)*(numMismatchesEnd - numMismatchesStart + 1));
	if(NULL == c.maxCount) {
		PrintError(FnName,
				"c.maxCount",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	curCounts = malloc(sizeof(int)*(numMismatchesEnd - numMismatchesStart + 1));
	if(NULL == curCounts) {
		PrintError(FnName,
				"curCounts",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Initialize counts */
	for(i=0;i<(numMismatchesEnd - numMismatchesStart + 1);i++) {
		c.counts[i] = NULL;
		c.maxCount[i] = 0;
	}

	/* Go through every possible read in the genome using the index */
	fprintf(stderr, "Out of %lld, currently on\n0",
			(long long int)index->length);
	curIndex = 0;
	nextIndex = 0;
	numDifferent = 0;
	for(curIndex=0, nextIndex=0, counter=0, numDifferent=0;
			curIndex < index->length;
			curIndex = nextIndex) {
		if(counter >= BINDEXSTAT_ROTATE_NUM) {
			fprintf(stderr, "\r%lld", (long long int)curIndex);
			counter -= BINDEXSTAT_ROTATE_NUM;
		}
		/* Initialize */
		curIndex = nextIndex;
		numReadsNoMismatches = 0;
		/* Try each mismatch */
		for(i=numMismatchesStart;i<=numMismatchesEnd;i++) {
			/* Initialize variables */
			RGMatchInitialize(&m);
			curCounts[i] = 0;

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

			/* Get the matches */
			RGReadsFindMatches(index,
					rg,
					&m,
					read,
					readLength,
					offsets,
					1, /* One offset */
					i, /* The number of mismatches */
					0,
					0,
					0,
					0,
					INT_MAX);
			/* Update the counts */
			curNum = m.numEntries;
			assert(curNum > 0);
			/*
			fprintf(stderr, "curIndex=%lld\ti=%d\tcurNum=%d\n",
					(long long int)curIndex,
					i,
					curNum);
			*/
			/* Update the value of numReadsNoMismatches if necessary */
			if(i==0) {
				numReadsNoMismatches = curNum; /* This will be the basis for update c.counts */
			}

			/* curNum is now the x-axis, or the second index in the c.counts array */
			/* The value we should add is actually numReadsNoMismatches sinc we wish to skip over these,
			 * i.e. don't search again since they have the same answer */

			/* Add to our list.  We may have to reallocate this array */
			if(curNum > c.maxCount[i]) {
				j = c.maxCount[i]; /* Save previous length */
				/* Reallocate */
				c.maxCount[i] = curNum;
				assert(c.maxCount[i] > 0);
				c.counts[i] = realloc(c.counts[i], sizeof(int)*c.maxCount[i]);
				if(NULL == c.counts[i]) {
					PrintError(FnName,
							"counts",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}
				/* Initialize from j to maxCount */
				while(j<curNum) {
					c.counts[i][j] = 0;
					j++;
				}
			}
			/* Add the number of places to the appropriate counts */
			assert(curNum >= 0 && curNum <= c.maxCount[i]);
			c.counts[i][curNum-1] += numReadsNoMismatches; /* Add the number of reads for which we will skip over */
			/* Update curIndex */
			if(i==0) {
				/* Add the range since we will be skipping over them */
				nextIndex += curNum;
				counter += curNum;
				/* Update stats */
				numDifferent++;
			}
			/* Free variables */
			RGMatchFree(&m);
		}
	}
	fprintf(stderr, "\n");

	/* Print results */
	fprintf(fp, "# Number of unique places was: %lld\nThe mean number of CALs was: %lld/%lld=%lf\n",
			(long long int)numDifferent,
			(long long int)numDifferent,
			(long long int)2*index->length, /* Times two for both strands */
			(double)(index->length*2.0)/numDifferent); /* Times two for both strands */
	for(i=numMismatchesStart;i<=numMismatchesEnd;i++) {
		fprintf(fp, "# Found counts for %d mismatches ranging from %d to %d.\n",
				i,
				1,
				c.maxCount[i]);
		for(j=0;j<c.maxCount[i];j++) {
			fprintf(fp, "%lld\t%d\n",
					(long long int)j,
					c.counts[i][j]);
		}
	}

	/* Free memory */
	for(i=0;i<(numMismatchesEnd - numMismatchesStart + 1);i++) {
		free(c.counts[i]);
		c.counts[i] = NULL;
	}
	free(c.counts);
	c.counts = NULL;
	free(c.maxCount);
	c.maxCount = NULL;
	free(curCounts);
	curCounts = NULL;
}

int64_t GetNumberOfMatches(RGIndex *index,
		RGBinary *rg,
		char *read,
		int readLength,
		int64_t *endIndex,
		int64_t curIndex)
{
	char *FnName = "GetNumberOfMatches";
	int64_t tmpEndIndex=-1;
	int64_t foundIndex=0;
	int64_t startIndex=0;
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
	if(1!=foundIndex) {
		fprintf(stderr, "\nread:[%s]\n", read);
		fprintf(stderr, "\nstarts:[%lld]\tends:[%lld]\n",
				(long long int)index->starts[hashIndex],
				(long long int)index->ends[hashIndex]);
		PrintError(FnName,
				"foundIndex",
				"Did not find the read",
				Exit,
				OutOfRange);
	}
	/* Reverse compliment will have startIndex == -1 */
	assert(startIndex == -1 || startIndex == curIndex);
	assert(startIndex <= tmpEndIndex);
	if(endIndex != NULL) {
		(*endIndex) = tmpEndIndex;
	}
	/* Return the number of places with this read */
	return (tmpEndIndex - startIndex + 1);
}
