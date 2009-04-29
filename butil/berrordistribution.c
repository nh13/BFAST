#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <limits.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "../blib/AlignedRead.h"
#include "../blib/AlignedEntry.h"
#include "berrordistribution.h"

#define Name "berrordistribution"
#define ROTATE_NUM 100000

/* Counts the number of errors at each position in the read
 * for both nucleotide and color space, as well as counting the 
 * number of reads containing a given number of errors.
 */

int main(int argc, char *argv[])
{
	int32_t i;
	int32_t countGaps=0;
	int32_t trimEndGaps=0;
	char outputID[MAX_FILENAME_LENGTH]="\0";
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	Errors e;

	if(5 <= argc) {

		/* Get arguments */
		countGaps=atoi(argv[1]);
		trimEndGaps=atoi(argv[2]);
		strcpy(outputID, argv[3]);

		fprintf(stderr, "%s", BREAK_LINE);
		ErrorsInitialize(&e);
		for(i=4;i<argc;i++) {
			strcpy(inputFileName, argv[i]);
			/* Run the program */
			ErrorDistribution(inputFileName,
					countGaps,
					trimEndGaps,
					&e);
		}
		ErrorDistributionPrint(outputID, countGaps, &e);
		fprintf(stderr, "%s%s%s", 
				BREAK_LINE, 
				"Terminating successfully!\n",
				BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<count gaps as errors>\n");
		fprintf(stderr, "\t<trim gaps at the end of the read>\n");
		fprintf(stderr, "\t<output id>\n");
		fprintf(stderr, "\t<BFAST reported files>\n");
	}
	return 0;
}

void ErrorDistribution(char *inputFileName,
		int32_t countGaps,
		int32_t trimEndGaps,
		Errors *e)
{
	char *FnName="ErrorDistribution";
	int count;
	FILE *FPin=NULL;
	AlignedRead a;

	fprintf(stderr, "Reading from %s.\n%s",
			inputFileName,
			BREAK_LINE);

	/* Open the input file */
	if(!(FPin = fopen(inputFileName, "rb"))) {
		PrintError(FnName,
				inputFileName,
				"Could not open file for reading",
				Exit,
				OpenFileError);
	}

	/* Go through the input file */
	count=0;
	AlignedReadInitialize(&a);
	fprintf(stderr, "Currently on:\n0");
	while(EOF != AlignedReadRead(&a, FPin, BinaryInput)) {
		if(e->numEnds < a.numEnds) {
			e->numEnds = a.numEnds;
		}
		if(count % ROTATE_NUM == 0) {
			fprintf(stderr, "\r%d", count);
		}
		/* Update */
		ErrorsUpdate(e, &a, countGaps, trimEndGaps);
		/* Free memory */
		AlignedReadFree(&a);
		count++;
	}
	fprintf(stderr, "\r%d\n", count);

	/* Close files */
	fclose(FPin);
}

void ErrorDistributionPrint(char *outputID,
		int32_t countGaps,
		Errors *e)
{
	char *FnName="ErrorDistributionPrintAndClose";
	int i;
	FILE *FPsout[6]={NULL,NULL,NULL,NULL,NULL,NULL};
	char outputFileNames[6][MAX_FILENAME_LENGTH]={"\0","\0","\0","\0","\0","\0"};

	/* Create output file names */
	sprintf(outputFileNames[0], "%s.error.distribution.nt.by.position.%s.txt",
			PROGRAM_NAME,
			outputID);
	sprintf(outputFileNames[1], "%s.error.distribution.color.by.position.%s.txt",
			PROGRAM_NAME,
			outputID);
	sprintf(outputFileNames[2], "%s.error.distribution.gaps.by.position.%s.txt",
			PROGRAM_NAME,
			outputID);
	sprintf(outputFileNames[3], "%s.error.distribution.nt.across.reads.%s.txt",
			PROGRAM_NAME,
			outputID);
	sprintf(outputFileNames[4], "%s.error.distribution.color.across.reads.%s.txt",
			PROGRAM_NAME,
			outputID);
	sprintf(outputFileNames[5], "%s.error.distribution.gaps.across.reads.%s.txt",
			PROGRAM_NAME,
			outputID);

	/* Open the output file */
	for(i=0;i<3;i++) {
		if(0 == countGaps || i != 2) {
			if(!(FPsout[i] = fopen(outputFileNames[i], "wb"))) {
				PrintError(FnName,
						outputFileNames[i],
						"Could not open file for writing",
						Exit,
						OpenFileError);
			}
			if(!(FPsout[i+3] = fopen(outputFileNames[i+3], "wb"))) {
				PrintError(FnName,
						outputFileNames[i+3],
						"Could not open file for writing",
						Exit,
						OpenFileError);
			}
		}
	}

	/* Print output */
	ErrorsPrint(e, FPsout, countGaps);

	/* Close files */
	for(i=0;i<3;i++) {
		if(0 == countGaps || i != 2) {
			fclose(FPsout[i]);
			fclose(FPsout[i+3]);
		}
	}

	/* Free memory */
	ErrorsFree(e);
}

void ErrorsPrint(Errors *e, FILE **fps, int countGaps)
{
	/* Print each out to a separate file */
	assert(fps!=NULL);
	int i;

	for(i=0;i<3;i++) {
		if(0 == countGaps || i != 2) {
			fprintf(fps[i], "# The number of mapped reads was %d\n", e->numReads);
			fprintf(fps[i+3], "# The number of mapped reads was %d\n", e->numReads);
		}
	}

	/* Print out */
	CountPrint(&e->by[0], fps[0]);
	CountPrint(&e->by[1], fps[1]);
	CountPrint(&e->across[0], fps[3]);
	CountPrint(&e->across[1], fps[4]);
	if(0 == countGaps) {
		CountPrint(&e->by[2], fps[2]);
		CountPrint(&e->across[2], fps[5]);
	}
}

void ErrorsUpdate(Errors *e, AlignedRead *a, int32_t countGaps, int32_t trimEndGaps)
{
	int32_t i;

	for(i=0;i<a->numEnds;i++) {
		if(1 == a->ends[i].numEntries) {
			ErrorsUpdateHelper(e, &a->ends[i].entries[0], countGaps, trimEndGaps, a->space, i+1);
		}
	}

	e->numReads++;
}

void ErrorsUpdateHelper(Errors *e, AlignedEntry *a, int countGaps, int trimEndGaps, int space, int which) 
{
	int i, ctr;
	int numNT, numColor, numGap;
	assert(0 < which);
	assert(ColorSpace == space || NTSpace == space);
	int32_t length;

	/* Initialize */
	numNT = numColor = numGap = 0;

	length = a->length;
	if(1 == trimEndGaps) {
		while(0 < length && GAP == a->reference[length-1]) {
			length--;
		}
	}

	/* Go through each position in the alignment */
	for(i=0,ctr=0;i<length;i++) {

		/* Skip gaps in the read */
		if(GAP != a->read[i]) {

			/* nt space */
			/* Check if it is a gap and we are not to count gaps as errors */
			if(0 == countGaps && a->reference[i] == GAP) {
				CountUpdate(&e->by[2],
						CountOnly,
						ctr,
						which);
				numGap++;
			}
			else if(ToLower(a->read[i]) !=  ToLower(a->reference[i])) {
				CountUpdate(&e->by[0],
						CountOnly,
						ctr,
						which);
				numNT++;
			}
			CountUpdate(&e->by[2],
					CountTotal,
					ctr,
					which);
			CountUpdate(&e->by[0],
					CountTotal,
					ctr,
					which);

			/* color space */
			if(ColorSpace == space) {
				if(GAP != a->colorError[i]) {
					CountUpdate(&e->by[1],
							CountBoth,
							ctr,
							which);
					numColor++;
				}
				else {
					CountUpdate(&e->by[1],
							CountTotal,
							ctr,
							which);
				}

			}
			ctr++;
		}
	}

	/* Update across */

	/* nt space */
	CountUpdate(&e->across[0],
			CountBoth,
			numNT,
			which);
	/* color space */
	if(ColorSpace == space) {
		CountUpdate(&e->across[1],
				CountBoth,
				numColor,
				which);
	}
	/* gaps */
	if(0 == countGaps) {
		CountUpdate(&e->across[2],
				CountBoth,
				numGap,
				which);
	}
}

void ErrorsInitialize(Errors *e) 
{
	e->numReads=0;
	CountInitialize(&e->by[0]);
	CountInitialize(&e->by[1]);
	CountInitialize(&e->by[2]);
	CountInitialize(&e->across[0]);
	CountInitialize(&e->across[1]);
	CountInitialize(&e->across[2]);
}

void ErrorsFree(Errors *e)
{
	CountFree(&e->by[0]);
	CountFree(&e->by[1]);
	CountFree(&e->by[2]);
	CountFree(&e->across[0]);
	CountFree(&e->across[1]);
	CountFree(&e->across[2]);
}

void CountPrint(Count *c, FILE *fp)
{
	int i, j, maxLength;
	
	maxLength = INT_MIN;
	for(i=0;i<c->numEnds;i++) {
		if(maxLength < c->lengths[i]) {
			maxLength = c->lengths[i];
		}
	}

	for(i=0;i<maxLength;i++) {
		/* Index - 0-based */
		fprintf(fp, "%d", i);
		/* Counts and totals */
		for(j=0;j<c->numEnds;j++) {
			if(i < c->lengths[j]) {
				fprintf(fp, "\t%d\t%d", c->counts[j][i], c->totals[j][i]);
			}
			else {
				fprintf(fp, "\t%d\t%d", 0, 0);
			}
		}
		fprintf(fp, "\n");
	}
}

void CountUpdate(Count *c,
		int type,
		int index,
		int which)
{
	char *FnName="CountUpdate";
	assert(0<=index);
	int prev, i;

	assert(0 < which);

	/* Reallocate if necessary based on the number of ends */
	if(c->numEnds < which) {
		prev = c->numEnds;
		c->numEnds = which;
		c->lengths = realloc(c->lengths, sizeof(int)*c->numEnds);
		if(NULL==c->lengths) {
			PrintError(FnName,
					"c->lengths",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		c->counts = realloc(c->counts, sizeof(int*)*c->numEnds);
		if(NULL==c->counts) {
			PrintError(FnName,
					"c->counts",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		c->totals = realloc(c->totals, sizeof(int*)*c->numEnds);
		if(NULL==c->totals) {
			PrintError(FnName,
					"c->totals",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		/* Initialize */
		for(i=prev;i<c->numEnds;i++) {
			c->lengths[i] = 0;
			c->counts[i] = c->totals[i] = NULL;
		}
	}
	/* Reallocate if necessary based on the length of the end */
	if(c->lengths[which-1] <= index) {
		prev = c->lengths[which-1];
		c->lengths[which-1] = index+1;

		c->counts[which-1] = realloc(c->counts[which-1], sizeof(int)*(c->lengths[which-1]));
		if(NULL==c->counts[which-1]) {
			PrintError(FnName,
					"c->counts[which-1]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		for(i=prev;i<c->lengths[which-1];i++) {
			c->counts[which-1][i] = 0;
		}
		c->totals[which-1] = realloc(c->totals[which-1], sizeof(int)*(c->lengths[which-1]));
		if(NULL==c->totals[which-1]) {
			PrintError(FnName,
					"c->totals[which-1]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		for(i=prev;i<c->lengths[which-1];i++) {
			c->totals[which-1][i] = 0;
		}
	}

	/* Increment */
	switch(type) {
		case CountOnly:
			c->counts[which-1][index]++;
			break;
		case CountTotal:
			c->totals[which-1][index]++;
			break;
		case CountBoth:
			c->counts[which-1][index]++;
			c->totals[which-1][index]++;
			break;
		default:
			PrintError(FnName,
					"type",
					"Could not understand type",
					Exit,
					OutOfRange);
	}
}

void CountInitialize(Count *c)
{
	c->numEnds=0;
	c->lengths=NULL;
	c->counts=c->totals=NULL;
}

void CountFree(Count *c)
{
	int32_t i;
	for(i=0;i<c->numEnds;i++) {
		free(c->counts[i]);
		free(c->totals[i]);
	}
	free(c->counts);
	free(c->totals);
	free(c->lengths);
	CountInitialize(c);
}
