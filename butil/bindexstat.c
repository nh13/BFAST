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

int main(int argc, char *argv[]) 
{
	FILE *fp=NULL;
	char indexFileName[MAX_FILENAME_LENGTH]="\0";
	char histogramFileName[MAX_FILENAME_LENGTH]="\0";
	char rgFileName[MAX_FILENAME_LENGTH]="\0";
	int numMismatchesStart = NUM_MISMATCHES_START;
	int numMismatchesEnd = NUM_MISMATCHES_END;

	if(argc == 4) {
		RGBinary rg;
		RGIndex index;

		strcpy(rgFileName, argv[1]);
		strcpy(indexFileName, argv[2]);
		numMismatchesEnd = atoi(argv[4]);
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

		fprintf(stderr, "%s", BREAK_LINE);
		PrintHistogram(&index, 
				&rg, 
				numMismatchesStart,
				numMismatchesEnd,
				histogramFileName);
		fprintf(stderr, "%s", BREAK_LINE);

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
		fprintf(stderr, "Usage: bindexstat [OPTIONS]\n");
		fprintf(stderr, "\t\t<reference genome file name>\n");
		fprintf(stderr, "\t\t<index file name>\n");
		fprintf(stderr, "\t\t<number of mismatches for histogram>\n");
	}

	return 0;
}

void PrintHistogram(RGIndex *index, 
		RGBinary *rg,
		int numMismatchesStart,
		int numMismatchesEnd,
		char *histogramFileName)
{
	char *FnName = "PrintHistogram";
	char tmpFileName[2048]="\0";
	int i;
	int64_t j;
	int64_t curIndex, nextIndex;
	int64_t counter;
	int64_t numDifferent = 0;
	int64_t numReadsNoMismatches = 0;
	int numForward, numReverse;
	Counts c; /* Used to store the histogram data */
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
	/* Go through the index */
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

			/* Get the matches for the chr/pos */
			GetMatchesFromChrPos(index,
					rg,
					index->chromosomes[curIndex],
					index->positions[curIndex],
					i,
					&numForward, 
					&numReverse);

			/* Update the value of numReadsNoMismatches and numDifferent
			 * if we have the rsults for no mismatches */
			if(i==0) {
				/* This will be the basis for update c.counts */
				numReadsNoMismatches = numForward + numReverse;
				numDifferent++;
				/* If the reverse compliment does not match the + strand then 
				 * count it as unique. */
				if(numReverse == 0) {
					numDifferent++;
				}
				/* Add the range since we will be skipping over them */
				nextIndex += numForward;
				counter += numForward;
			}

			/* Add to our list.  We may have to reallocate this array */
			if(numForward + numReverse > c.maxCount[i]) {
				j = c.maxCount[i]+1; /* This will determine where we begin initialization after reallocation */
				/* Reallocate */
				c.maxCount[i] = numForward + numReverse;
				assert(c.maxCount[i] > 0);
				c.counts[i] = realloc(c.counts[i], sizeof(int)*(c.maxCount[i]+1));
				if(NULL == c.counts[i]) {
					PrintError(FnName,
							"counts",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}
				/* Initialize from j to maxCount */
				while(j<=c.maxCount[i]) {
					c.counts[i][j] = 0;
					j++;
				}
			}
			/* Add twice the number of reads when we searched with no mismatches, since these are the original seeds.
			 * Add it twice since the forward and reverse strands will be symmetric. */
			assert(numReadsNoMismatches > 0);
			c.counts[i][numForward+numReverse] += 2*numReadsNoMismatches;
		}
	}
	fprintf(stderr, "\r%lld\n", (long long int)index->length);

	/* Print results */
	for(i=numMismatchesStart;i<=numMismatchesEnd;i++) {
		/* Create file name */
		sprintf(tmpFileName, "%s.%d",
				histogramFileName,
				i);
		if(!(fp = fopen(tmpFileName, "w"))) {
			PrintError(FnName,
					tmpFileName,
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}
		fprintf(fp, "# Number of unique places was: %lld\n# The mean number of CALs was: %lld/%lld=%lf\n",
				(long long int)numDifferent,
				(long long int)numDifferent,
				(long long int)2*index->length, /* Times two for both strands */
				((double)index->length*2.0)/numDifferent); /* Times two for both strands */
		fprintf(fp, "# Found counts for %d mismatches ranging from %d to %d.\n",
				i,
				1,
				c.maxCount[i]);
		for(j=0;j<=c.maxCount[i];j++) {
			fprintf(fp, "%lld\t%d\n",
					(long long int)j,
					c.counts[i][j]);
		}
		fclose(fp);
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
}

/* Get the matches for the chr/pos */
void GetMatchesFromChrPos(RGIndex *index,
		RGBinary *rg,
		uint32_t curChr,
		uint32_t curPos,
		int numMismatches,
		int *numForward,
		int *numReverse)
{
	char *FnName = "GetMatchesFromChrPos";
	char read[SEQUENCE_LENGTH]="\0";
	char reverseRead[SEQUENCE_LENGTH]="\0";
	int readLength = index->totalLength;
	int returnLength, returnPosition;
	int i;
	RGReads reads;
	RGRanges ranges;

	/* Initialiez reads */
	RGReadsInitialize(&reads);
	RGRangesInitialize(&ranges);

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

	/* First generate the perfect match for the forward and
	 * reverse strand */
	GetReverseComplimentAnyCase(read,
			reverseRead,
			readLength);
	RGReadsGeneratePerfectMatch(read,
			readLength,
			FORWARD,
			0,
			index->numTiles,
			index->tileLengths,
			index->gaps,
			index->totalLength,
			&reads);
	RGReadsGeneratePerfectMatch(reverseRead,
			readLength,
			REVERSE,
			0,
			index->numTiles,
			index->tileLengths,
			index->gaps,
			index->totalLength,
			&reads);

	if(numMismatches > 0) {
		/* Generate reads with the necessary mismatches for 
		 * both the forward and reverse strands */
		RGReadsGenerateMismatches(reverseRead,
				readLength,
				REVERSE,
				0,
				index->numTiles,
				index->tileLengths,
				index->gaps,
				index->totalLength,
				numMismatches,
				&reads);
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
		RGIndexGetRanges(index,
				rg,
				reads.reads[i],
				reads.readLength[i],
				reads.strand[i],
				reads.offset[i],
				&ranges);
	}

	/* Remove duplicates */
	RGRangesRemoveDuplicates(&ranges);

	/* Error check */
	if(ranges.numEntries <= 0) {
		PrintError(FnName,
				"ranges",
				"Returned zero ranges",
				Exit,
				OutOfRange);
	}

	/* Return the number of FORWARD strand matches so that we can skip over */
	(*numForward) = (*numReverse) = 0;
	for(i=0;i<ranges.numEntries;i++) {
		switch(ranges.strand[i]) {
			case FORWARD:
				break;
				(*numForward) += ranges.endIndex[i] - ranges.startIndex[i] + 1;
			case REVERSE:
				(*numReverse) += ranges.endIndex[i] - ranges.startIndex[i] + 1;
				break;
			default:
				PrintError(FnName,
						"m->strand[i]",
						"Could not understand strand",
						Exit,
						OutOfRange);
				break;
		}
	}
	assert(numForward>0);

	/* Free memory */
	RGReadsFree(&reads);
	RGRangesFree(&ranges);
}
