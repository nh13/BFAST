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
#include "../blib/RGMatches.h"
#include "../blib/RGReads.h"
#include "bpairedenddist.h"

#define Name "bpairedenddist"
#define BINDEXDIST_ROTATE_NUM 1000000

/* Prints the distribution of the distance between paired-end reads
 * using reads that have both ends matching only one location on 
 * the same strand.
 * */

int main(int argc, char *argv[]) 
{
	FILE *fp=NULL;
	char matchesFileName[MAX_FILENAME_LENGTH]="\0";
	int binaryInput;

	if(argc == 3) {
		strcpy(matchesFileName, argv[1]);
		binaryInput = atoi(argv[2]);
		assert(binaryInput == 0 || binaryInput == 1);

		/* Read the index */
		fprintf(stderr, "Reading in bfast matches file from %s.\n",
				matchesFileName);
		if(!(fp=fopen(matchesFileName, "rb"))) {
			PrintError(Name,
					matchesFileName,
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}

		fprintf(stderr, "%s", BREAK_LINE);
		PrintDistribution(fp,
				binaryInput,
				stdout);
		fprintf(stderr, "%s", BREAK_LINE);

		/* Close the file */
		fclose(fp);

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully!\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t\t<bfast matches file name>\n");
		fprintf(stderr, "\t\t<binary Input: 0-text 1-binary>\n");
	}

	return 0;
}

void PrintDistribution(FILE *fpIn,
		int binaryInput,
		FILE *fpOut)
{
	RGMatches m;
	int64_t posOne, posTwo, difference;
	int64_t counter=0, numUnique=0;

	/* Initialize */
	RGMatchesInitialize(&m);

	while(EOF != RGMatchesRead(fpIn,
				&m,
				1,
				binaryInput)) {
		/* Only use unique sequences on the same chromosome and strand */
			if(1 == m.matchOne.numEntries &&
					1 == m.matchTwo.numEntries &&
					m.matchOne.chromosomes[0] == m.matchTwo.chromosomes[0] &&
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
