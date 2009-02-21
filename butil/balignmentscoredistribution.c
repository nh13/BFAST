#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "../blib/AlignEntries.h"
#include "../blib/AlignEntry.h"
#include "../blib/BError.h"
#include "balignmentscoredistribution.h"

#define Name "balignmentscoredistribution"
#define MIN_SCORE -2500
#define MAX_SCORE 2500
#define ROTATE_NUM 100000

int main(int argc, char *argv[])
{
	FILE *fpIn=NULL;
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	AlignEntries a;
	int64_t counter;
	double from, by, to;
	Dist dist;

	if(5 == argc) {
		strcpy(inputFileName, argv[1]);
		from=atof(argv[2]);
		by=atof(argv[3]);
		to=atof(argv[4]);

		/* Open the input file */
		if(!(fpIn=fopen(inputFileName, "rb"))) {
			PrintError(Name,
					inputFileName,
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}

		DistInitialize(&dist);
		AlignEntriesInitialize(&a);
		counter = 0;
		fprintf(stderr, "Currently on:\n0");
		/* Read in each match */
		while(EOF != AlignEntriesRead(&a, fpIn, PairedEndDoesNotMatter, SpaceDoesNotMatter, BinaryInput)) {
			if(counter%ROTATE_NUM==0) {
				fprintf(stderr, "\r%lld",
						(long long int)counter);
			}
			counter++;
			/* Add */
			DistAdd(&dist, a.entriesOne, a.numEntriesOne, from, by, to);
			if(PairedEnd == a.pairedEnd) {
				DistAdd(&dist, a.entriesTwo, a.numEntriesTwo, from, by, to);
			}
			/* Free */
			AlignEntriesFree(&a);
		}
		fprintf(stderr, "\r%lld\n",
				(long long int)counter);

		/* Close the input file */
		fclose(fpIn);

		/* Print */
		DistPrint(&dist, stdout);

		/* Free */
		DistFree(&dist);

		fprintf(stderr, "Terminating successfully!\n");
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast aligned file>\n");
		fprintf(stderr, "\t<bins: from>\n");
		fprintf(stderr, "\t<bins: by>\n");
		fprintf(stderr, "\t<bins: to>\n");
	}
	return 0;
}

void DistInitialize(Dist *dist)
{
	dist->max_cals = 0;
	dist->cals = NULL;
}

int32_t DistAdd(Dist *dist, 
		AlignEntry *a,
		int32_t numEntries,
		double from,
		double by,
		double to)
{
	char *FnName="DistAdd";
	int32_t i, prev;

	if(numEntries <= 0) {
		return 0;
	}

	/* Check if we need to reallocate */
	if(dist->max_cals < numEntries) {
		prev = dist->max_cals;
		dist->cals = realloc(dist->cals, sizeof(CAL)*numEntries);
		if(NULL == dist->cals) {
			PrintError(FnName,
					"dist->cals",
					"Could not allocate memory",
					Exit,
					OutOfRange);
		}
		for(i=prev;i<numEntries;i++) {
			CALInitialize(&dist->cals[i], i+1, from, by, to);
		}
		dist->max_cals = numEntries;
	}
	/* Copy over */
	return CALAdd(&dist->cals[numEntries-1], a);
}

void DistPrint(Dist *dist,
		FILE *fpOut)
{
	int32_t i;

	for(i=0;i<dist->max_cals;i++) {
		CALPrint(&dist->cals[i], fpOut);
	}
}

void DistFree(Dist *dist)
{
	int32_t i;

	for(i=0;i<dist->max_cals;i++) {
		CALFree(&dist->cals[i]);
	}
	DistInitialize(dist);
}

void CALInitialize(CAL *c, int32_t num_cal, double from, double by, double to) 
{
	char *FnName="CALInitialize";
	int i;

	if(num_cal < 0) {
		c->num_cal = 0;
		c->from = 0.0;
		c->by = 0.0;
		c->to = 0.0;
		c->length = 0;
		c->hist = NULL;
	}
	else {
		c->num_cal = num_cal;
		c->from = from;
		c->by = by;
		c->to = to;
		c->length = (int32_t)((to - from)/by)+ 1;
		c->hist = malloc(sizeof(int32_t)*c->length);
		if(NULL == c->hist) {
			PrintError(FnName,
					"c->hist",
					"Could not allocate memory",
					Exit,
					OutOfRange);
		}
		for(i=0;i<c->length;i++) {
			c->hist[i] = 0;
		}
	}
}

int32_t CALAdd(CAL *c, 
		AlignEntry *a)
{
	int32_t i, index;
	int32_t numOmitted=0;

	for(i=0;i<c->num_cal;i++) {
		/* Check if we should reallocate */
		/* Add in */
		index = (int32_t)((a[i].score - c->from)/c->by);
		if(0 <= index && index < c->length) {
			c->hist[index]++;
			numOmitted++;
		}
	}
	return numOmitted;
}

void CALPrint(CAL *c, FILE *fpOut) 
{
	int32_t i;

	fprintf(fpOut, "%d\t%lf\t%lf\t%lf", 
			c->num_cal,
			c->from,
			c->by,
			c->to);
	for(i=0;i<c->length;i++) { 
		fprintf(fpOut, "\t%d", c->hist[i]);
	}
	fprintf(fpOut, "\n");
}

void CALFree(CAL *c)
{
	free(c->hist);
	CALInitialize(c, -1, 0.0, 0.0, 0.0);
}
