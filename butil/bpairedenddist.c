#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <pthread.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BLib.h"
#include "../blib/BError.h"
#include "../blib/RGMatch.h"
#include "../blib/RGMatches.h"
#include "../blib/AlignEntries.h"
#include "bpairedenddist.h"

#define Name "bpairedenddist"
#define BINDEXDIST_ROTATE_NUM 1000000

/* Prints the distribution of the distance between paired-end reads
 * using reads that have both ends matching only one location on 
 * the same strand.
 * */

int main(int argc, char *argv[]) 
{
	FILE *fpIn=NULL;
	FILE *fpOut=NULL;
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	char outputID[MAX_FILENAME_LENGTH]="\0";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";

	if(argc == 3) {
		strcpy(inputFileName, argv[1]);
		strcpy(outputID, argv[2]);
		sprintf(outputFileName, "%s.paired.end.distribution.%s.txt",
				PROGRAM_NAME,
				outputID);

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Reading in from %s.\n",
				inputFileName);
		if(!(fpIn=fopen(inputFileName, "rb"))) {
			PrintError(Name,
					inputFileName,
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}
		fprintf(stderr, "%s", BREAK_LINE);

		fprintf(stderr, "Writing to %s.\n",
				outputFileName);
		if(!(fpOut=fopen(outputFileName, "rb"))) {
			PrintError(Name,
					outputFileName,
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}
		fprintf(stderr, "%s", BREAK_LINE);

		if(NULL!=strstr(inputFileName, BFAST_MATCHES_FILE_EXTENSION)) {
			PrintDistributionFromBMF(fpIn,
					fpOut);
		}
		else if(NULL!=strstr(inputFileName, BFAST_ALIGNED_FILE_EXTENSION)) {
			PrintDistributionFromBAF(fpIn,
					fpOut);
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
		fclose(fpIn);
		fclose(fpOut);

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully!\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t\t<bfast matches, aligned, or reported file name>\n");
		fprintf(stderr, "\t\t<output id>\n");
	}

	return 0;
}

void PrintDistributionFromBMF(FILE *fpIn,
		FILE *fpOut)
{
	char *FnName = "PrintDistributionFromBMF";
	RGMatches m;
	int64_t posOne, posTwo, difference;
	int64_t counter=0, numUnique=0;

	/* Initialize */
	RGMatchesInitialize(&m);

	while(EOF != RGMatchesRead(fpIn,
				&m,
				PairedEndDoesNotMatter,
				BinaryInput)) {
		if(m.pairedEnd != PairedEnd) {
			PrintError(FnName,
					"m.pairedEnd",
					"Data was not paired end",
					Exit,
					OutOfRange);
		}
		/* Only use unique sequences on the same contig and strand */
		if(1 == m.matchOne.numEntries &&
				1 == m.matchTwo.numEntries &&
				m.matchOne.contigs[0] == m.matchTwo.contigs[0] &&
				m.matchOne.strand[0] == m.matchTwo.strand[0]) {
			/* Simple way to avoid overflow */
			posOne = m.matchOne.positions[0];
			posTwo = m.matchTwo.positions[0];
			difference = posOne - posTwo;
			/* Print */
			/*
			   fprintf(fpOut, "%lld\t%lld\t%lld\n",
			   (long long int)posOne,
			   (long long int)posTwo,
			   (long long int)difference);
			   */
			fprintf(fpOut, "%lld\n",
					(long long int)difference);
			numUnique++;
		}
		counter++;

		/* Free matches */
		RGMatchesFree(&m);
	}

	fprintf(stderr, "number unique was %lld out of %lld.\n",
			(long long int)numUnique,
			(long long int)counter);

}

void PrintDistributionFromBAF(FILE *fpIn,
		FILE *fpOut)
{
	char *FnName = "PrintDistributionFromBAF";
	AlignEntries a;
	int64_t posOne, posTwo, difference;
	int64_t counter=0, numUnique=0;

	/* Initialize */
	AlignEntriesInitialize(&a);

	while(EOF != AlignEntriesRead(&a,
				fpIn,
				PairedEndDoesNotMatter,
				SpaceDoesNotMatter,
				BinaryInput)) {
		if(a.pairedEnd != PairedEnd) {
			PrintError(FnName,
					"a.pairedEnd",
					"Data was not paired end",
					Exit,
					OutOfRange);
		}
		/* Only use unique sequences on the same contig and strand */
		if(1 == a.numEntriesOne &&
				1 == a.numEntriesTwo &&
				a.entriesOne[0].contig == a.entriesTwo[0].contig &&
				a.entriesOne[0].strand == a.entriesTwo[0].strand) {
			/* Simple way to avoid overflow */
			posOne = a.entriesOne[0].position;
			posTwo = a.entriesTwo[0].position;
			difference = posOne - posTwo;
			fprintf(fpOut, "%lld\n",
					(long long int)difference);
			numUnique++;
		}
		counter++;

		/* Free memory */
		AlignEntriesFree(&a);
	}

	fprintf(stderr, "number unique was %lld out of %lld.\n",
			(long long int)numUnique,
			(long long int)counter);

}
