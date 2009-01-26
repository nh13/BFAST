#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "Definitions.h"
#include "../blib/BError.h"
#include "../blib/BLibDefinitions.h"
#include "ScoringMatrix.h"

/* TODO */
int ScoringMatrixRead(char *scoringMatrixFileName, 
		ScoringMatrix *sm,
		int space)
{
	char *FnName="ScoringMatrixRead";
	int i, j;
	int32_t tempInt;
	FILE *fp;

	/* Open the scoring matrix file */
	if((fp=fopen(scoringMatrixFileName, "r"))==0) {
		PrintError(FnName,
				scoringMatrixFileName,
				"Could not open scoringMatrixFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Read in the gap open penalty */
	if(fscanf(fp, "%d", &tempInt)==EOF) {
		PrintError(FnName,
				scoringMatrixFileName,
				"Could not read in the gap open penalty",
				Exit,
				OutOfRange);
	}
	sm->gapOpenPenalty = tempInt;

	/* Read in the gap close penalty */
	if(fscanf(fp, "%d", &tempInt)==EOF) {
		PrintError(FnName,
				scoringMatrixFileName,
				"Could not read in the gap extension penalty",
				Exit,
				OutOfRange);
	}
	sm->gapExtensionPenalty = tempInt;

	/* Assume the NT key is acgt */
	assert(ALPHABET_SIZE==4);
	/* Allocate memory for the NT key */
	sm->NTKeys = (char*)malloc(sizeof(char)*(ALPHABET_SIZE+1));
	if(NULL == sm->NTKeys) {
		PrintError(FnName,
				"sm->NTKeys",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Read in the NT key */
	/* Assume the NT key is acgt */
	sm->NTKeys[0] = 'a';
	sm->NTKeys[1] = 'c';
	sm->NTKeys[2] = 'g';
	sm->NTKeys[3] = 't';
	sm->NTKeys[4] = 'n';

	/* Allocate memory for the scores */
	sm->NTScores = (int32_t**)malloc(sizeof(int32_t*)*(ALPHABET_SIZE+1));
	if(NULL==sm->NTScores) {
		PrintError(FnName,
				"sm->NTScores",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	for(i=0;i<ALPHABET_SIZE+1;i++) {
		sm->NTScores[i] = (int32_t*)malloc(sizeof(int32_t)*(ALPHABET_SIZE+1));
		if(NULL==sm->NTScores[i]) {
			PrintError(FnName,
					"sm->NTScores[i]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	/* Read in the score matrix */
	sm->maxNTScore = INT_MIN;
	for(i=0;i<ALPHABET_SIZE+1;i++) { /* Read row */
		for(j=0;j<ALPHABET_SIZE+1;j++) { /* Read column */
			if(fscanf(fp, "%d", &tempInt)==EOF) {
				PrintError(FnName,
						scoringMatrixFileName,
						"Could not read in the scoring matrix (NTScores)",
						Exit,
						OutOfRange);
			}
			sm->NTScores[i][j]= tempInt;
			if(sm->maxNTScore < sm->NTScores[i][j]) {
				sm->maxNTScore = sm->NTScores[i][j];
			}
		}
	}

	if(space == 1) {
		/* Assume the color key is 0123 */
		assert(ALPHABET_SIZE==4);
		/* Allocate memory for the color key */
		sm->ColorKeys = malloc(sizeof(int)*(ALPHABET_SIZE+1));
		if(NULL == sm->ColorKeys) {
			PrintError(FnName,
					"sm->ColorKeys",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Read in the color key */
		/* Assume the color key is acgt */
		sm->ColorKeys[0] = 0;
		sm->ColorKeys[1] = 1;
		sm->ColorKeys[2] = 2;
		sm->ColorKeys[3] = 3;
		sm->ColorKeys[4] = 4;

		/* Allocate memory for the scores */
		sm->ColorScores = (int32_t**)malloc(sizeof(int32_t*)*(ALPHABET_SIZE+1));
		if(NULL==sm->ColorScores) {
			PrintError(FnName,
					"sm->ColorScores",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		for(i=0;i<ALPHABET_SIZE+1;i++) {
			sm->ColorScores[i] = (int32_t*)malloc(sizeof(int32_t)*(ALPHABET_SIZE+1));
			if(NULL==sm->ColorScores[i]) {
				PrintError(FnName,
						"sm->ColorScores[i]",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}
		/* Read in the score matrix */
		sm->maxColorScore = INT_MIN;
		for(i=0;i<ALPHABET_SIZE+1;i++) { /* Read row */
			for(j=0;j<ALPHABET_SIZE+1;j++) { /* Read column */
				if(fscanf(fp, "%d", &tempInt)==EOF) {
					PrintError(FnName,
							scoringMatrixFileName,
							"Could not read in the scoring matrix (ColorScores)",
							Exit,
							OutOfRange);
				}
				sm->ColorScores[i][j] = tempInt;
				if(sm->maxColorScore < sm->ColorScores[i][j]) {
					sm->maxColorScore = sm->ColorScores[i][j];
				}
			}
		}
	}
	else {
		sm->ColorKeys=NULL;
		sm->ColorScores=NULL;
	}
	
	ScoringMatrixCheck(sm, space);

	/* Close the file */
	fclose(fp);

	return 1;
}

/* TODO */
void ScoringMatrixInitialize(ScoringMatrix *sm)
{
	sm->gapOpenPenalty=0;
	sm->gapExtensionPenalty=0;
	sm->NTKeys=NULL;
	sm->NTScores=NULL;
	sm->maxNTScore=0;
	sm->ColorKeys=NULL;
	sm->ColorScores=NULL;
	sm->maxColorScore=0;
}

/* TODO */
void ScoringMatrixFree(ScoringMatrix *sm)
{
	int i;
	free(sm->NTKeys);
	free(sm->ColorKeys);
	for(i=0;i<ALPHABET_SIZE+1;i++) {
		free(sm->NTScores[i]);
		if(NULL != sm->ColorScores) {
			free(sm->ColorScores[i]);
		}
	}
	free(sm->NTScores);
	free(sm->ColorScores);
	ScoringMatrixInitialize(sm);
}

/* TODO */
int32_t ScoringMatrixGetNTScore(char a,
		char b,
		ScoringMatrix *sm)
{
	int indexA=-1;
	int indexB=-1;

	/* Get index for a */
	switch(a) {
		case 'A':
		case 'a':
			indexA=0;
			break;
		case 'C':
		case 'c':
			indexA=1;
			break;
		case 'G':
		case 'g':
			indexA=2;
			break;
		case 'T':
		case 't':
			indexA=3;
			break;
		case 'N':
		case 'n':
			indexA=4;
			break;
		default:
			fprintf(stderr, "\n[%c]\n", a);
			PrintError("GetNTScore",
					NULL,
					"Could not understand key",
					Exit,
					OutOfRange);
	}
	/* Get index for b */
	switch(b) {
		case 'A':
		case 'a':
			indexB=0;
			break;
		case 'C':
		case 'c':
			indexB=1;
			break;
		case 'G':
		case 'g':
			indexB=2;
			break;
		case 'T':
		case 't':
			indexB=3;
			break;
		case 'N':
		case 'n':
			indexB=4;
			break;
		default:
			fprintf(stderr, "b key:[%c]\n", b);
			PrintError("GetNTScore",
					NULL,
					"Could not understand 'b' key",
					Exit,
					OutOfRange);
			break;
	}

	return sm->NTScores[indexA][indexB];
}

/* TODO */
int32_t ScoringMatrixGetColorScore(uint8_t a,
		uint8_t b,
		ScoringMatrix *sm)
{
	char *FnName="ScoringMatrixGetColorScore";
	int indexA=(int)a;
	int indexB=(int)b;
	if(indexA < 0 || indexA > 4) {
		fprintf(stderr, "indexA=%d\n", indexA);
		PrintError(FnName,
				NULL,
				"indexA out of range",
				Exit,
				OutOfRange);
	}
	if(indexB < 0 || indexB > 4) {
		fprintf(stderr, "indexB=%d\n", indexB);
		PrintError(FnName,
				NULL,
				"indexB out of range",
				Exit,
				OutOfRange);
	}

	return sm->ColorScores[indexA][indexB];
}

/* TODO */
/* For color space only */
int32_t ScoringMatrixCheck(ScoringMatrix *sm,
		int space) {
	char *FnName="ScoringMatrixCheck";
	int i, j;

	if(0 < sm->gapOpenPenalty) {
		PrintError(FnName,
				"sm->gapOpenPenalty",
				"Must be less than or equal to zero",
				Exit,
				OutOfRange);
	}
	if(0 < sm->gapExtensionPenalty) {
		PrintError(FnName,
				"sm->gapExtensionPenalty",
				"Must be less than or equal to zero",
				Exit,
				OutOfRange);
	}

	for(i=0;i<ALPHABET_SIZE+1;i++) {
		for(j=0;j<ALPHABET_SIZE+1;j++) {
			if(i==j) {
				if(sm->NTScores[i][i] < 0) {
					PrintError(FnName,
							"sm->NTScores[i][i]",
							"Must be greater than or equal to zero",
							Exit,
							OutOfRange);
				}
				if(ColorSpace == space && sm->ColorScores[i][i] < 0) {
					PrintError(FnName,
							"sm->ColorScores[i][i]",
							"Must be greater than or equal to zero",
							Exit,
							OutOfRange);
				}
			}
			else {
				if(0 < sm->NTScores[i][j]) {
					PrintError(FnName,
							"sm->NTScores[i][j]",
							"Must be less than or equal to zero",
							Exit,
							OutOfRange);
				}
				if(ColorSpace == space && 0 < sm->ColorScores[i][j]) {
					PrintError(FnName,
							"sm->ColorScores[i][j]",
							"Must be less than or equal to zero",
							Exit,
							OutOfRange);
				}
			}
		}
	}
	return 1;
}
