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
	char histogramFileName[MAX_FILENAME_LENGTH]="\0";
	char rgFileName[MAX_FILENAME_LENGTH]="\0";
	int algorithm;
	int numMismatchesStart = NUM_MISMATCHES_START;
	int numMismatchesEnd = NUM_MISMATCHES_END;

	if(argc >= 4 && argc <= 5) {
		RGBinary rg;
		RGIndex index;

		strcpy(rgFileName, argv[1]);
		strcpy(indexFileName, argv[2]);
		algorithm = atoi(argv[3]);
		if(argc == 5) {
			numMismatchesEnd = atoi(argv[4]);
		}
		assert(0 <= algorithm && algorithm <= 2);
		assert(numMismatchesEnd >= numMismatchesStart);

		/* Create the histogram file name */
		strcpy(histogramFileName, indexFileName);
		strcat(histogramFileName, ".hist");

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
						histogramFileName);
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
						histogramFileName);
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
		fprintf(stdout, "Usage: bindexstat [OPTIONS]\n\t\t<reference genome file name>\n\t\t<index file name>\n\t\t<algorithm: 0-Summary 1-Histogram 2-Both>\n\t\t<number of mismatches for histogram>\n");
	}

	return 0;
}

void PrintSummary(RGIndex *index, RGBinary *rg)
{
	int64_t max, min;
	long double sum;
	long double mean, variance;
	int64_t curIndex, nextIndex;
	int64_t counter;
	int64_t numDifferent = 0;
	int numForward;
	RGMatch m; /* Used to store the matches returned */

	/* Get the mean */
	fprintf(stderr, "%s", BREAK_LINE);
	fprintf(stderr, "Getting the mean. Out of %lld, currently on:\n0",
			(long long int)index->length);
	for(curIndex=0, nextIndex=0, counter=0, numDifferent=0;
			curIndex < index->length;
			curIndex = nextIndex) {
		if(counter >= BINDEXSTAT_ROTATE_NUM) {
			fprintf(stderr, "\r%lld", (long long int)curIndex);
			counter -= BINDEXSTAT_ROTATE_NUM;
		}

		/* Initialize variables */
		RGMatchInitialize(&m);

		/* Get the matches for the chr/pos */
		numForward = GetMatchesFromChrPos(index,
				rg,
				index->chromosomes[curIndex],
				index->positions[curIndex],
				0,
				&m);

		/* We add two since a there is symmetry between the forward and reverse strands */
		numDifferent += 2; 
		/* Update the counter and the next index.  We only want to skip over the forward
		 * matches. 
		 * */
		counter += numForward;
		nextIndex = curIndex + numForward;

		/* Free matches */
		RGMatchFree(&m);
	}
	fprintf(stderr, "\r%lld\n", (long long int)curIndex);

	/* Multiply by two in the numerator since it is only the length of the forward strand and
	 * we considered both the forward and reverse strands */
	mean = ((double)(2.0*index->length))/numDifferent;

	/* Get the variance, max, and min */
	fprintf(stderr, "Getting the variance. Out of %lld, currently on:\n0",
			(long long int)index->length);
	min = UINT_MAX;
	max = -1;
	sum = 0.0;
	for(curIndex=0, nextIndex=0, counter=0;
			curIndex < index->length;
			curIndex = nextIndex) {
		if(counter >= BINDEXSTAT_ROTATE_NUM) {
			fprintf(stderr, "\r%lld", (long long int)curIndex);
			counter -= BINDEXSTAT_ROTATE_NUM;
		}

		/* Initialize variables */
		RGMatchInitialize(&m);

		/* Get the matches for the chr/pos */
		numForward = GetMatchesFromChrPos(index,
				rg,
				index->chromosomes[curIndex],
				index->positions[curIndex],
				0,
				&m);

		/* Update the counter and the next index.  We only want to skip over the forward
		 * matches. 
		 * */
		counter += numForward;
		nextIndex = curIndex + numForward;

		/* Get max */
		if(m.numEntries > max) {
			max = m.numEntries;
		}
		/* Get min */
		if(m.numEntries < min) {
			min = m.numEntries;
		}
		/* Update variance sum */
		/* Do this twice, since the forward and reverse strand 
		 * will be symmetric in their matches. 
		 * */
		sum += 2*(m.numEntries - mean)*(m.numEntries - mean);

		/* Free matches */
		RGMatchFree(&m);
	}
	fprintf(stderr, "\r%lld\n", (long long int)curIndex);
	variance = sum/numDifferent;

	fprintf(stderr, "The mean was: %Lf\nThe variance was: %Lf\nThe max was: %lld\nThe min was: %lld\nThe number of unique reads was: %lld out of %lld\n",
			mean,
			variance,
			(long long int)max,
			(long long int)min,
			(long long int)numDifferent,
			(long long int)2*index->length);
}

void PrintHistogram(RGIndex *index, 
		RGBinary *rg,
		int numMismatchesStart,
		int numMismatchesEnd,
		char *histogramFileName)
{
	char *FnName = "PrintHistogram";
	int i;
	int64_t j;
	int64_t curIndex, nextIndex;
	int32_t curNum;
	int64_t counter;
	int64_t numDifferent = 0;
	int64_t numReadsNoMismatches = 0;
	int numForward;
	Counts c; /* Used to store the histogram data */
	RGMatch m; /* Used to store the matches returned */
	int *curCounts = NULL;
	FILE *fp;

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

			/* Get the matches for the chr/pos */
			numForward = GetMatchesFromChrPos(index,
					rg,
					index->chromosomes[curIndex],
					index->positions[curIndex],
					i,
					&m);

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
			/* The value we should add is actually numReadsNoMismatches */

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
			/* Add the number of original reads, this should equal curNum for no mismatches, but will be less when
			 * we consider mismatches. 
			 * Add it twice since the forward and reverse strands will be symmetric. */
			c.counts[i][curNum-1] += 2*numReadsNoMismatches; 
			/* Update curIndex */
			if(i==0) {
				/* Add the range since we will be skipping over them */
				nextIndex += numForward;
				counter += numForward;
				/* Update stats */
				numDifferent+=2; /* Two for both strands */
			}
			/* Free variables */
			RGMatchFree(&m);
		}
	}
	fprintf(stderr, "\r%lld\n", (long long int)curIndex);

	/* Print results */
	if(!(fp = fopen(histogramFileName, "w"))) {
		PrintError(FnName,
				histogramFileName,
				"Could not open file for writing",
				Exit,
				OpenFileError);
	}
	fprintf(fp, "# Number of unique places was: %lld\nThe mean number of CALs was: %lld/%lld=%lf\n",
			(long long int)numDifferent,
			(long long int)numDifferent,
			(long long int)2*index->length, /* Times two for both strands */
			((double)index->length*2.0)/numDifferent); /* Times two for both strands */
	for(i=numMismatchesStart;i<=numMismatchesEnd;i++) {
		fprintf(fp, "# Found counts for %d mismatches ranging from %d to %d.\n",
				i,
				1,
				c.maxCount[i]);
		for(j=0;j<c.maxCount[i];j++) {
			fprintf(fp, "%lld\t%d\n",
					(long long int)j+1,
					c.counts[i][j]);
		}
	}
	fclose(fp);

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

/* Get the matches for the chr/pos */
int GetMatchesFromChrPos(RGIndex *index,
		RGBinary *rg,
		uint32_t curChr,
		uint32_t curPos,
		int numMismatches,
		RGMatch *m)
{
	char *FnName = "GetMatchesFromChrPos";
	char read[SEQUENCE_LENGTH]="\0";
	int readLength = index->totalLength;
	int returnLength, returnPosition;
	int i;
	int numForward;
	RGReads reads;

	/* Initialiez reads */
	RGReadsInitialize(&reads);

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

	/* First generate the perfect match */
	RGReadsGeneratePerfectMatch(read,
			readLength,
			FORWARD,
			0,
			index->numTiles,
			index->tileLengths,
			index->gaps,
			index->totalLength,
			&reads);

	if(numMismatches > 0) {
		/* Generate reads with the necessary mismatches */
		RGReadsGenerateMismatches(read,
				readLength,
				FORWARD,
				0,
				index->numTiles,
				index->tileLengths,
				index->gaps,
				index->totalLength,
				numMismatches,
				&reads);
	}

	for(i=0;i<reads.numReads;i++) {
		/* Get the matches for the read */
		RGIndexGetMatches(index,
				rg,
				reads.reads[i],
				reads.readLength[i],
				reads.strand[i],
				reads.offset[i],
				m,
				INT_MAX);
	}

	/* Free memory */
	RGReadsFree(&reads);

	/* Error check */
	if(m->numEntries <= 0 ||
			m->maxReached == 1) {
		RGMatchPrint(stderr,
				"Name",
				"Sequence",
				NULL,
				m,
				NULL,
				0,
				0);
		PrintError(FnName,
				"m",
				"Returned zero matches or maximum matches was reached",
				Exit,
				OutOfRange);
	}

	/* Return the number of FORWARD strand matches so that we can skip over */
	numForward = 0;
	for(i=0;i<m->numEntries;i++) {
		if(m->strand[i] == FORWARD) {
			numForward++;
		}
	}
	assert(numForward>0);
	return numForward;
}
