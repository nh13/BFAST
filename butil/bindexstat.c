#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/RGIndex.h"
#include "bindexstat.h"

int main(int argc, char *argv[]) 
{
	FILE *fp=NULL;
	char indexFileName[MAX_FILENAME_LENGTH]="\0";
	char rgFileName[MAX_FILENAME_LENGTH]="\0";

	if(argc == 3) {
		RGBinary rg;
		RGIndex index;

		strcpy(rgFileName, argv[1]);
		strcpy(indexFileName, argv[2]);

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
		PrintMeanAndVarianceOfCAL(&index, &rg);

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Cleaning up.\n");
		/* Delete the index */
		RGIndexDelete(&index);
		/* Delete the rg */
		RGBinaryDelete(&rg);
	}
	else {
		fprintf(stdout, "Please give a reference genome file name then an index file name.  Terminating!\n");
	}

	return 0;
}

void PrintMeanAndVarianceOfCAL(RGIndex *index, RGBinary *rg)
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
		int cmp = RGIndexCompareAt(index, rg, start, end);
		assert(cmp <= 0);
		if(cmp < 0) {
			numEntries++;
			/* Update */
			start = end;
		}
	}
	mean = (index->length*1.0)/numEntries;

	/* Get the variance, max, and min */
	fprintf(stderr, "Getting the mean. Out of %lld, currently on:\n0",
			(long long int)index->length);
	min = UINT_MAX;
	max = -1;
	sum = 0.0;
	start = 0;
	for(end=1;end < index->length;end++) {
		if(end%RGINDEX_ROTATE_NUM==0) {
			fprintf(stderr, "\r%lld", (long long int)end);
		}
		int cmp = RGIndexCompareAt(index, rg, start, end);
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
	variance = sum/numEntries;

	fprintf(stderr, "The mean was: %Lf\nThe variance was: %Lf\nThe max was: %lld\nThe min was: %lld\n",
			mean,
			variance,
			(long long int)max,
			(long long int)min);
}

