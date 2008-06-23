#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "Definitions.h"
#include "../blib/BError.h"
#include "../blib/BLibDefinitions.h"
#include "ReadInputFiles.h"

/* TODO */
int ReadScoringMatrix(char *scoringMatrixFileName, ScoringMatrix *sm)
{
	int i, j;
	FILE *fp;

	/* Open the scoring matrix file */
	if((fp=fopen(scoringMatrixFileName, "r"))==0) {
		PrintError("ReadScoringMatrix",
				scoringMatrixFileName,
				"Could not open scoringMatrixFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Read in the gap open penalty */
	if(fscanf(fp, "%lf", &sm->gapOpenPenalty)==EOF) {
		PrintError("ReadScoringMatrix",
				scoringMatrixFileName,
				"Could not read in the gap open penalty",
				Exit,
				OutOfRange);
	}

	/* Read in the gap close penalty */
	if(fscanf(fp, "%lf", &sm->gapExtensionPenalty)==EOF) {
		PrintError("ReadScoringMatrix",
				scoringMatrixFileName,
				"Could not read in the gap extension penalty",
				Exit,
				OutOfRange);
	}

	/* Assume the key is acgt */
	assert(ALPHABET_SIZE==4);
	/* Allocate memory for the key */
	sm->key = (char*)malloc(sizeof(char)*(ALPHABET_SIZE+1));
	if(NULL == sm->key) {
		PrintError("ReadScoringMatrix",
				"sm->key",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Read in the key */
	/* Assume the key is acgt */
	sm->key[0] = 'a';
	sm->key[1] = 'c';
	sm->key[2] = 'g';
	sm->key[3] = 't';
	sm->key[4] = 'n';

	/* Allocate memory for the scores */
	sm->scores = (double**)malloc(sizeof(double*)*(ALPHABET_SIZE+1));
	if(NULL==sm->scores) {
		PrintError("ReadScoringMatrix",
			"sm->scores",
			"Could not allocate memory",
			Exit,
			MallocMemory);
	}
	for(i=0;i<ALPHABET_SIZE+1;i++) {
		sm->scores[i] = (double*)malloc(sizeof(double)*(ALPHABET_SIZE+1));
		if(NULL==sm->scores[i]) {
			PrintError("ReadScoringMatrix",
					"sm->scores[i]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	/* Read in the score matrix */
	for(i=0;i<ALPHABET_SIZE+1;i++) { /* Read row */
		for(j=0;j<ALPHABET_SIZE+1;j++) { /* Read column */
			if(fscanf(fp, "%lf", &sm->scores[i][j])==EOF) {
				PrintError("ReadScoringMatrix",
						scoringMatrixFileName,
						"Could not read in the scoring matrix",
						Exit,
						OutOfRange);
			}
		}
	}

	/* Close the file */
	fclose(fp);

	return 1;
}
