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

/* Counts the number of errors at each position in the read
 * for both nucleotide and color space, as well as the 
 * number of reads containing a given number of errors.
 */

int main(int argc, char *argv[])
{
	if(3 == argc) {
		char inputFileName[MAX_FILENAME_LENGTH]="\0";
		char outputID[MAX_FILENAME_LENGTH]="\0";

		/* Get arguments */
		strcpy(inputFileName, argv[1]);
		strcpy(outputID, argv[2]);

		/* Run the program */
		GetErrorDistribution(inputFileName,
				outputID);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<BFAST reported file>\n");
		fprintf(stderr, "\t<output id>\n");
	}
	return 0;
}

void GetErrorDistribution(char *inputFileName,
		char *outputID)
{
	char *FnName="GetErrorDistribution";
	int i;
	FILE *FPin=NULL;
	FILE *FPsout[4]={NULL,NULL,NULL,NULL};
	char outputFileNames[4][MAX_FILENAME_LENGTH]={"\0","\0","\0","\0"};
	AlignEntries a;
	Errors e;

	/* Create output file names */
	sprintf(outputFileNames[0], "%s.error.distribution.nt.by.position.%s.txt",
			PROGRAM_NAME,
			outputID);
	sprintf(outputFileNames[1], "%s.error.distribution.color.by.position.%s.txt",
			PROGRAM_NAME,
			outputID);
	sprintf(outputFileNames[2], "%s.error.distribution.nt.across.reads.%s.txt",
			PROGRAM_NAME,
			outputID);
	sprintf(outputFileNames[3], "%s.error.distribution.color.across.reads.%s.txt",
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
	for(i=0;i<4;i++) {
		if(!(FPsout[i] = fopen(outputFileNames[i], "wb"))) {
			PrintError(FnName,
					outputFileNames[i],
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}
	}

	/* Go through the input file */
	AlignEntriesInitialize(&a);
	ErrorsInitialize(&e);
	while(EOF != AlignEntriesRead(&a, FPin, PairedEndDoesNotMatter, SpaceDoesNotMatter, BinaryInput)) {
		/* Update */
		ErrorsUpdate(&e, &a);
		/* Free memory */
		AlignEntriesFree(&a);
	}

	/* Print output */
	ErrorsPrint(&e, FPsout);

	/* Close files */
	fclose(FPin);
	for(i=0;i<4;i++) {
		fclose(FPsout[i]);
	}
	
	/* Free memory */
	ErrorsFree(&e);
}

void ErrorsPrint(Errors *e, FILE **fps)
{
	/* Print each out to a separate file */
	assert(fps!=NULL);
	int i;

	for(i=0;i<4;i++) {
		CountPrint(&e->by[i%2], fps[i]);
	}
}

void ErrorsUpdate(Errors *e, AlignEntries *a)
{
	assert(a->numEntriesOne == 1 &&
			(a->pairedEnd == SingleEnd || a->numEntriesTwo == 1));

	ErrorsUpdateHelper(e, &a->entriesOne[0], a->space, 0);
	ErrorsUpdateHelper(e, &a->entriesTwo[0], a->space, 1);

	e->numReads++;
}

void ErrorsUpdateHelper(Errors *e, AlignEntry *a, int space, int which) 
{
	int i;
	int numNT, numColor;
	assert(which == 0 || which == 1);
	assert(ColorSpace == space || NTSpace == space);

	/* Initialize */
	numNT = numColor = 0;

	/* Go through each position in the alignment */
	for(i=0;i<a->length;i++) {

		/* Skip gaps in the read */
		if(GAP != a->read[i]) {
		
			/* nt space */
			if(ToLower(a->read[i]) !=  ToLower(a->reference[i])) {
				CountUpdate(&e->by[0],
						i,
						which);
				numNT++;
			}
		
			/* color space */
			if(ColorSpace == space &&
					1 == a->colorError[i]) {
				CountUpdate(&e->by[1],
						i,
						which);
				numColor++;
			}
		}
	}

	/* Update across */

	/* nt space */
	CountUpdate(&e->across[0],
			numNT,
			which);
	/* color space */
	if(PairedEnd == space) {
		CountUpdate(&e->across[1],
				numColor,
				which);
	}
}

void ErrorsInitialize(Errors *e) 
{
	CountInitialize(&e->by[0]);
	CountInitialize(&e->by[1]);
	CountInitialize(&e->across[0]);
	CountInitialize(&e->across[1]);
}

void ErrorsFree(Errors *e)
{
	CountFree(&e->by[0]);
	CountFree(&e->by[1]);
	CountFree(&e->across[0]);
	CountFree(&e->across[1]);
}

void CountPrint(Count *c, FILE *fp)
{
	int i;

	if(SingleEnd == c->pairedEnd) {
		for(i=0;i<c->length;i++) {
			fprintf(fp, "%d\t%d\n", i, c->countOne[i]);
		}
	}
	else {
		for(i=0;i<c->length;i++) {
			fprintf(fp, "%d\t%d\t%d\n", 
					i, 
					c->countOne[i],
					c->countTwo[i]);
		}
	}
}

void CountUpdate(Count *c,
		int index,
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

		if(PairedEnd == c->pairedEnd) {
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
		}
	}

	/* Increment */
	if(SingleEnd == which) {
		c->countOne[index]++;
	}
	else {
		c->countTwo[index]++;
	}
}

void CountInitialize(Count *c)
{
	c->pairedEnd=SingleEnd;
	c->length=0;
	c->countOne=c->countTwo=NULL;
}

void CountFree(Count *c)
{
	free(c->countOne);
	free(c->countTwo);
	CountInitialize(c);
}
