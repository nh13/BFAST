#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "../blib/AlignEntries.h"
#include "../blib/AlignEntry.h"
#include "bgeterrordistribution.h"

#define Name "bgeterrordistribution"
#define ROTATE_NUM 100000

/* Counts the number of errors at each position in the read
 * for both nucleotide and color space, as well as counting the 
 * number of reads containing a given number of errors.
 */

int main(int argc, char *argv[])
{
	if(5 == argc) {
		char inputFileName[MAX_FILENAME_LENGTH]="\0";
		int32_t countGaps=0;
		int32_t trimEndGaps=0;
		char outputID[MAX_FILENAME_LENGTH]="\0";

		/* Get arguments */
		strcpy(inputFileName, argv[1]);
		countGaps=atoi(argv[2]);
		trimEndGaps=atoi(argv[3]);
		strcpy(outputID, argv[4]);

		/* Run the program */
		GetErrorDistribution(inputFileName,
				countGaps,
				trimEndGaps,
				outputID);

		fprintf(stderr, "%s%s%s", 
				BREAK_LINE, 
				"Terminating successfully!\n",
				BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<BFAST reported file>\n");
		fprintf(stderr, "\t<count gaps as errors>\n");
		fprintf(stderr, "\t<trim gaps at the end of the read>\n");
		fprintf(stderr, "\t<output id>\n");
	}
	return 0;
}

void GetErrorDistribution(char *inputFileName,
		int32_t countGaps,
		int32_t trimEndGaps,
		char *outputID)
{
	char *FnName="GetErrorDistribution";
	int i, count;
	FILE *FPin=NULL;
	FILE *FPsout[6]={NULL,NULL,NULL,NULL,NULL,NULL};
	char outputFileNames[6][MAX_FILENAME_LENGTH]={"\0","\0","\0","\0","\0","\0"};
	AlignEntries a;
	Errors e;
	int pairedEnd = SingleEnd;

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

	/* Open the input file */
	if(!(FPin = fopen(inputFileName, "rb"))) {
		PrintError(FnName,
				inputFileName,
				"Could not open file for reading",
				Exit,
				OpenFileError);
	}

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

	/* Go through the input file */
	count=0;
	AlignEntriesInitialize(&a);
	ErrorsInitialize(&e);
	fprintf(stderr, "Currently on:\n0");
	while(EOF != AlignEntriesRead(&a, FPin, PairedEndDoesNotMatter, SpaceDoesNotMatter, BinaryInput)) {
		pairedEnd = a.pairedEnd;
		if(count % ROTATE_NUM == 0) {
			fprintf(stderr, "\r%d", count);
		}
		/* Update */
		ErrorsUpdate(&e, &a, countGaps, trimEndGaps);
		/* Free memory */
		AlignEntriesFree(&a);
		count++;
	}
	fprintf(stderr, "\r%d\n", count);

	/* Print output */
	ErrorsPrint(&e, FPsout, countGaps, pairedEnd);

	/* Close files */
	fclose(FPin);
	for(i=0;i<3;i++) {
		if(0 == countGaps || i != 2) {
			fclose(FPsout[i]);
			fclose(FPsout[i+3]);
		}
	}

	/* Free memory */
	ErrorsFree(&e);
}

void ErrorsPrint(Errors *e, FILE **fps, int countGaps, int pairedEnd)
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
	CountPrint(&e->by[0], fps[0], pairedEnd);
	CountPrint(&e->by[1], fps[1], pairedEnd);
	CountPrint(&e->across[0], fps[3], pairedEnd);
	CountPrint(&e->across[1], fps[4], pairedEnd);
	if(0 == countGaps) {
		CountPrint(&e->by[2], fps[2], pairedEnd);
		CountPrint(&e->across[2], fps[5], pairedEnd);
	}
}

void ErrorsUpdate(Errors *e, AlignEntries *a, int32_t countGaps, int32_t trimEndGaps)
{
	assert(a->numEntriesOne == 1 &&
			(a->pairedEnd == SingleEnd || a->numEntriesTwo == 1));

	ErrorsUpdateHelper(e, &a->entriesOne[0], countGaps, trimEndGaps, a->space, a->pairedEnd, 0);
	if(PairedEnd == a->pairedEnd) {
		ErrorsUpdateHelper(e, &a->entriesTwo[0], countGaps, trimEndGaps, a->space, a->pairedEnd, 1);
	}

	e->numReads++;
}

void ErrorsUpdateHelper(Errors *e, AlignEntry *a, int countGaps, int trimEndGaps, int space, int pairedEnd, int which) 
{
	int i, ctr;
	int numNT, numColor, numGap;
	assert(which == 0 || which == 1);
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
						pairedEnd,
						which);
				numGap++;
			}
			else if(ToLower(a->read[i]) !=  ToLower(a->reference[i])) {
				CountUpdate(&e->by[0],
						CountOnly,
						ctr,
						pairedEnd,
						which);
				numNT++;
			}
			CountUpdate(&e->by[2],
					CountTotal,
					ctr,
					pairedEnd,
					which);
			CountUpdate(&e->by[0],
					CountTotal,
					ctr,
					pairedEnd,
					which);

			/* color space */
			if(ColorSpace == space) {
				if('1' == a->colorError[i]) {
					CountUpdate(&e->by[1],
							CountBoth,
							ctr,
							pairedEnd,
							which);
					numColor++;
				}
				else {
					CountUpdate(&e->by[1],
							CountTotal,
							ctr,
							pairedEnd,
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
			pairedEnd,
			which);
	/* color space */
	if(ColorSpace == space) {
		CountUpdate(&e->across[1],
				CountBoth,
				numColor,
				pairedEnd,
				which);
	}
	/* gaps */
	if(0 == countGaps) {
		CountUpdate(&e->across[2],
				CountBoth,
				numGap,
				pairedEnd,
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

void CountPrint(Count *c, FILE *fp, int pairedEnd)
{
	int i;

	if(SingleEnd == pairedEnd) {
		for(i=0;i<c->length;i++) {
			fprintf(fp, "%d\t%d\t%d\n", i, c->countOne[i], c->totalOne[i]);
		}
	}
	else {
		for(i=0;i<c->length;i++) {
			fprintf(fp, "%d\t%d\t%d\t%d\t%d\n", 
					i, 
					c->countOne[i],
					c->totalOne[i],
					c->countTwo[i],
					c->totalTwo[i]);
		}
	}
}

void CountUpdate(Count *c,
		int type,
		int index,
		int pairedEnd,
		int which)
{
	char *FnName="CountUpdate";
	assert(0<=index);
	int prev, i;

	/* Reallocate if necessary */
	if(index >= c->length) {
		prev = c->length;
		c->length = index+1;

		c->countOne = realloc(c->countOne, sizeof(int)*(c->length));
		if(NULL==c->countOne) {
			PrintError(FnName,
					"c->countOne",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		for(i=prev;i<c->length;i++) {
			c->countOne[i] = 0;
		}
		c->totalOne = realloc(c->totalOne, sizeof(int)*(c->length));
		if(NULL==c->totalOne) {
			PrintError(FnName,
					"c->totalOne",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		for(i=prev;i<c->length;i++) {
			c->totalOne[i] = 0;
		}

		if(PairedEnd == pairedEnd) {
			c->countTwo = realloc(c->countTwo, sizeof(int)*(c->length));
			if(NULL==c->countTwo) {
				PrintError(FnName,
						"c->countTwo",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
			for(i=prev;i<c->length;i++) {
				c->countTwo[i] = 0;
			}
			c->totalTwo = realloc(c->totalTwo, sizeof(int)*(c->length));
			if(NULL==c->totalTwo) {
				PrintError(FnName,
						"c->totalTwo",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
			for(i=prev;i<c->length;i++) {
				c->totalTwo[i] = 0;
			}
		}
	}

	/* Increment */
	if(0 == which) {
		switch(type) {
			case CountOnly:
				c->countOne[index]++;
				break;
			case CountTotal:
				c->totalOne[index]++;
				break;
			case CountBoth:
				c->countOne[index]++;
				c->totalOne[index]++;
				break;
			default:
				PrintError(FnName,
						"type",
						"Could not understand type",
						Exit,
						OutOfRange);
		}
	}
	else {
		switch(type) {
			case CountOnly:
				c->countTwo[index]++;
				break;
			case CountTotal:
				c->totalTwo[index]++;
				break;
			case CountBoth:
				c->countTwo[index]++;
				c->totalTwo[index]++;
				break;
			default:
				PrintError(FnName,
						"type",
						"Could not understand type",
						Exit,
						OutOfRange);
		}
	}
}

void CountInitialize(Count *c)
{
	c->length=0;
	c->countOne=c->totalOne=c->countTwo=c->totalTwo=NULL;
}

void CountFree(Count *c)
{
	free(c->countOne);
	free(c->totalOne);
	free(c->countTwo);
	free(c->totalTwo);
	CountInitialize(c);
}
