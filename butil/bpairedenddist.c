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

int main(int argc, char *argv[]) 
{
	FILE *fp=NULL;
	char matchesFileName[MAX_FILENAME_LENGTH]="\0";

	if(argc == 2) {
		strcpy(matchesFileName, argv[1]);

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
				stderr);
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
	}

	return 0;
}

void PrintDistribution(FILE *fpIn,
		FILE *fpOut)
{
	RGMatches m;

	/* Initialize */
	RGMatchesInitialize(&m);

	while(EOF != RGMatchesRead(fpIn,
				&m,
				1,
				1)) {
		/* Only use unique sequences on the same chromosome and strand */
		if(1 == m.matchOne.numEntries &&
				1== m.matchTwo.numEntries &&
				m.matchOne.chromosomes[0] == m.matchTwo.chromosomes[0] &&
				m.matchOne.strand[0] == m.matchTwo.strand[0]) {
			/* Print */
			fprintf(stderr, "%d",
					(m.matchOne.positions[0] - m.matchTwo.positions[0]));
		}

		/* Free matches */
		RGMatchesFree(&m);
	}
}
