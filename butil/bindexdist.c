#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <pthread.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BLib.h"
#include "../blib/BError.h"
#include "../blib/RGIndex.h"
#include "../blib/RGRanges.h"
#include "../blib/RGMatch.h"
#include "../blib/RGReads.h"
#include "bindexdist.h"

#define Name "bindexdist"
#define BINDEXDIST_ROTATE_NUM 1000000

int main(int argc, char *argv[]) 
{
	FILE *fp=NULL;
	char indexFileName[MAX_FILENAME_LENGTH]="\0";
	char distributionFileName[MAX_FILENAME_LENGTH]="\0";
	char rgFileName[MAX_FILENAME_LENGTH]="\0";

	if(argc == 3) {
		RGBinary rg;
		RGIndex index;

		strcpy(rgFileName, argv[1]);
		strcpy(indexFileName, argv[2]);

		/* Create the distribution file name */
		strcpy(distributionFileName, indexFileName);
		strcat(distributionFileName, ".dist");

		/* Read in the rg binary file */
		RGBinaryReadBinary(&rg, rgFileName);

		/* Read the index */
		fprintf(stderr, "Reading in index from %s.\n",
				indexFileName);
		if(!(fp=fopen(indexFileName, "rb"))) {
			PrintError(Name,
					indexFileName,
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}
		RGIndexRead(fp, &index, 1);
		fclose(fp);

		fprintf(stderr, "%s", BREAK_LINE);
		PrintDistribution(&index, 
				&rg, 
				distributionFileName);
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
		fprintf(stderr, "Usage: bindexdist [OPTIONS]\n");
		fprintf(stderr, "\t\t<reference genome file name>\n");
		fprintf(stderr, "\t\t<index file name>\n");
	}

	return 0;
}

void PrintDistribution(RGIndex *index, 
		RGBinary *rg,
		char *distributionFileName)
{
	char *FnName = "PrintDistribution";
	FILE *fp;
	int64_t startIndex = 0;
	int64_t endIndex = index->length-1;
	int64_t curIndex=0, nextIndex=0;
	int64_t counter=0;
	int64_t numDifferent = 0;
	int64_t numForward, numReverse;
	char read[SEQUENCE_LENGTH] = "\0";

	/* Open the file */
	if(!(fp = fopen(distributionFileName, "w"))) {
		PrintError(FnName,
				distributionFileName,
				"Could not open file for writing",
				Exit,
				OpenFileError);
	}

	/* Print header */
	fprintf(fp, "Read\tTotal\tForward\tReverse\n");

	/* Go through every possible read in the genome using the index */
	for(curIndex=startIndex, nextIndex=startIndex, counter=0, numDifferent=0;
			curIndex <= endIndex;
			curIndex = nextIndex) {
		if(counter >= BINDEXDIST_ROTATE_NUM) {
			fprintf(stderr, "\r%10lld", 
					(long long int)(curIndex-startIndex));
			counter -= BINDEXDIST_ROTATE_NUM;
		}
		/* Get the matches for the chr/pos */
		GetMatchesFromChrPos(index,
				rg,
				index->chromosomes[curIndex],
				index->positions[curIndex],
				&numForward, 
				&numReverse,
				read);
		assert(numForward > 0);

		nextIndex += numForward;
		counter += numForward;

		/* Print the number of matches to file */
		fprintf(fp, "%s\t%lld\t%lld\t%lld\n",
				read,
				(long long int)(numForward+numReverse),
				(long long int)numForward,
				(long long int)numReverse);
	}
	fprintf(stderr, "\r%10lld\n", 
			(long long int)(curIndex-startIndex));

	/* Close the file */
	fclose(fp);
}

/* Get the matches for the chr/pos */
void GetMatchesFromChrPos(RGIndex *index,
		RGBinary *rg,
		uint32_t curChr,
		uint32_t curPos,
		int64_t *numForward,
		int64_t *numReverse, 
		char *read)
{
	char *FnName = "GetMatchesFromChrPos";
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
				(*numForward) += ranges.endIndex[i] - ranges.startIndex[i] + 1;
				break;
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
	assert((*numForward)>0);

	/* Free memory */
	RGReadsFree(&reads);
	RGRangesFree(&ranges);
}
