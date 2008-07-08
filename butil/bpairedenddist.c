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
	RGMatch m1, m2;
	char name[SEQUENCE_NAME_LENGTH]="\0";
	char read1[SEQUENCE_LENGTH]="\0";
	char read2[SEQUENCE_LENGTH]="\0";

	/* Initialize */
	RGMatchInitialize(&m1);
	RGMatchInitialize(&m2);

	/* Wow, a do-while loop */
	while(EOF != RGMatchRead(fpIn,
				name,
				read1, 
				read2,
				&m1,
				&m2,
				1,
				1)) {
		/* Only use unique sequences on the same chromosome and strand */
		if(1 == m1.numEntries &&
				1== m2.numEntries &&
				m1.chromosomes[0] == m2.chromosomes[0] &&
				m1.strand[0] == m2.strand[0]) {
			/* Print */
			fprintf(stderr, "%d",
					(m1.positions[0] - m2.positions[0]));
		}

		/* Free matches */
		RGMatchFree(&m1);
		RGMatchFree(&m2);
	}
}
