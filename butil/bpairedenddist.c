#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <pthread.h>
#include <zlib.h>

#include "../BLibDefinitions.h"
#include "../BLib.h"
#include "../BError.h"
#include "../RGMatch.h"
#include "../RGMatches.h"
#include "AlignedRead.h"
#include "bpairedenddist.h"

#define Name "bpairedenddist"
#define ROTATE_NUM 100000

/* Prints the distribution of the distance between paired-end reads
 * using reads that have both ends matching only one location on 
 * the same strand.
 * */

int main(int argc, char *argv[]) 
{
	gzFile fpIn=NULL;
	FILE *fpOut=NULL;
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	char outputID[MAX_FILENAME_LENGTH]="\0";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	int32_t i;
	Bins b;

	if(6 <= argc) {
		strcpy(outputID, argv[1]);
		b.minDistance = atoi(argv[2]);
		b.maxDistance = atoi(argv[3]);
		b.binSize = atoi(argv[4]);

		assert(b.minDistance <= b.maxDistance);

		/* Allocate memory */
		b.numCounts = (b.maxDistance - b.minDistance + 1) % b.binSize;
		b.numCounts += (int32_t)(b.maxDistance - b.minDistance + 1)/b.binSize;
		b.counts = calloc(b.numCounts, sizeof(int32_t));
		if(NULL == b.counts) {
			PrintError(Name,
					"b.counts",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}

		for(i=5;i<argc;i++) {
			strcpy(inputFileName, argv[i]);

			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Reading in from %s.\n",
					inputFileName);
			if(!(fpIn=gzopen(inputFileName, "rb"))) {
				PrintError(Name,
						inputFileName,
						"Could not open file for reading",
						Exit,
						OpenFileError);
			}
			fprintf(stderr, "%s", BREAK_LINE);

			if(NULL!=strstr(inputFileName, BFAST_MATCHES_FILE_EXTENSION)) {
				PrintDistributionFromBMF(fpIn,
						&b);
			}
			else if(NULL!=strstr(inputFileName, BFAST_ALIGNED_FILE_EXTENSION)) {
				PrintDistributionFromBAF(fpIn,
						&b);
			}
			else {
				PrintError(Name,
						"input file",
						"Could not recognize input file extension",
						Warn,
						OutOfRange);
			}
			fprintf(stderr, "%s", BREAK_LINE);

			/* Close the file */
			gzclose(fpIn);
		}

		/* Output */
		sprintf(outputFileName, "%s.paired.end.distribution.%s.txt",
				PROGRAM_NAME,
				outputID);
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Writing to %s.\n",
				outputFileName);
		if(!(fpOut=fopen(outputFileName, "wb"))) {
			PrintError(Name,
					outputFileName,
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}
		BinsPrint(&b, fpOut);
		fclose(fpOut);

		/* Free memory */
		free(b.counts);

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully!\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<output id>\n");
		fprintf(stderr, "\t<minimum distance>\n");
		fprintf(stderr, "\t<maximum distance>\n");
		fprintf(stderr, "\t<bin size>\n");
		fprintf(stderr, "\t<bfast matches, aligned, or reported file name(s)>\n");
	}

	return 0;
}

void PrintDistributionFromBMF(gzFile fpIn,
		Bins *b)
{
	char *FnName = "PrintDistributionFromBMF";
	RGMatches m;
	int64_t posOne, posTwo, difference;
	int64_t counter=0, numFound=0;
	int32_t warnUnpaired=0;

	/* Initialize */
	RGMatchesInitialize(&m);

	fprintf(stderr, "Currently on:\n0");
	while(EOF != RGMatchesRead(fpIn,
				&m)) {
		if(0==counter%ROTATE_NUM) {
			fprintf(stderr, "\r%lld",
					(long long int)counter);
		}
		if(2 != m.numEnds && 0 == warnUnpaired) {
			warnUnpaired=1;
			PrintError(FnName,
					"m.numEnds",
					"Data includes non paired end",
					Warn,
					OutOfRange);
		}
		else {
			/* Only use found sequences on the same contig and strand */
			if(1 == m.ends[0].numEntries &&
					1 == m.ends[1].numEntries &&
					m.ends[0].contigs[0] == m.ends[1].contigs[0] &&
					m.ends[0].strands[0] == m.ends[1].strands[0]) {
				/* Simple way to avoid overflow */
				posOne = m.ends[0].positions[0];
				posTwo = m.ends[1].positions[0];
				difference = posTwo - posOne;
				/* Print */
				/*
				   fprintf(fpOut, "%lld\t%lld\t%lld\n",
				   (long long int)posOne,
				   (long long int)posTwo,
				   (long long int)difference);
				   */
				if(1==BinsInsert(b, difference)) {
					numFound++;
				}
			}
		}
		counter++;

		/* Free matches */
		RGMatchesFree(&m);
	}
	fprintf(stderr, "\r%lld\n",
			(long long int)counter);

	fprintf(stderr, "number found was %lld out of %lld.\n",
			(long long int)numFound,
			(long long int)counter);
}

void PrintDistributionFromBAF(gzFile fpIn,
		Bins *b)
{
	char *FnName = "PrintDistributionFromBAF";
	AlignedRead a;
	int64_t posOne, posTwo, difference;
	int64_t counter=0, numFound=0;
	int32_t warnUnpaired=0;

	/* Initialize */
	AlignedReadInitialize(&a);

	fprintf(stderr, "Currently on:\n0");
	while(EOF != AlignedReadRead(&a,
				fpIn)) {
		if(0==counter%ROTATE_NUM) {
			fprintf(stderr, "\r%lld",
					(long long int)counter);
		}
		if(2 != a.numEnds && 0 == warnUnpaired) {
			warnUnpaired = 1;
			PrintError(FnName,
					"a.numEnds",
					"Data includes non paired end",
					Warn,
					OutOfRange);
		}
		else {
			/* Only use found sequences on the same contig and strand */
			if(1 == a.ends[0].numEntries &&
					1 == a.ends[1].numEntries &&
					a.ends[0].entries[0].contig == a.ends[1].entries[0].contig &&
					a.ends[0].entries[0].strand == a.ends[1].entries[0].strand) {
				/* Simple way to avoid overflow */
				posOne = a.ends[0].entries[0].position;
				posTwo = a.ends[1].entries[0].position;
				difference = posTwo - posOne;
				if(1==BinsInsert(b, difference)) {
					numFound++;
				}
			}
		}
		counter++;

		/* Free memory */
		AlignedReadFree(&a);
	}
	fprintf(stderr, "\r%lld\n",
			(long long int)counter);

	fprintf(stderr, "number found was %lld out of %lld.\n",
			(long long int)numFound,
			(long long int)counter);

}

int BinsInsert(Bins *b,
		int32_t difference)
{
	int32_t index;

	if(b->minDistance <= difference &&
			difference <= b->maxDistance) {
		index = (difference - b->minDistance); 
		index -= (difference % b->binSize);
		index /= b->binSize;
		assert(0 <= index &&
				index < b->numCounts);
		b->counts[index]++;
		return 1;
	}
	return 0;
}

void BinsPrint(Bins *b,
		FILE *fpOut)
{
	int32_t i;

	for(i=0;i<b->numCounts;i++) {
		fprintf(fpOut, "%10d\t%10d\t%10d\n",
				b->minDistance + i*b->binSize,
				MIN(b->minDistance + (i+1)*b->binSize-1, b->maxDistance),
				b->counts[i]);
	}
}
