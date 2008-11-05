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
		int colorSpace)
{
	char *FnName="ScoringMatrixRead";
	int i, j;
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
	if(fscanf(fp, "%lf", &sm->gapOpenPenalty)==EOF) {
		PrintError(FnName,
				scoringMatrixFileName,
				"Could not read in the gap open penalty",
				Exit,
				OutOfRange);
	}

	/* Read in the gap close penalty */
	if(fscanf(fp, "%lf", &sm->gapExtensionPenalty)==EOF) {
		PrintError(FnName,
				scoringMatrixFileName,
				"Could not read in the gap extension penalty",
				Exit,
				OutOfRange);
	}

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
	sm->NTScores = (double**)malloc(sizeof(double*)*(ALPHABET_SIZE+1));
	if(NULL==sm->NTScores) {
		PrintError(FnName,
				"sm->NTScores",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	for(i=0;i<ALPHABET_SIZE+1;i++) {
		sm->NTScores[i] = (double*)malloc(sizeof(double)*(ALPHABET_SIZE+1));
		if(NULL==sm->NTScores[i]) {
			PrintError(FnName,
					"sm->NTScores[i]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}
	/* Read in the score matrix */
	for(i=0;i<ALPHABET_SIZE+1;i++) { /* Read row */
		for(j=0;j<ALPHABET_SIZE+1;j++) { /* Read column */
			if(fscanf(fp, "%lf", &sm->NTScores[i][j])==EOF) {
				PrintError(FnName,
						scoringMatrixFileName,
						"Could not read in the scoring matrix (NTScores)",
						Exit,
						OutOfRange);
			}
		}
	}

	if(colorSpace == 1) {
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
		sm->ColorScores = (double**)malloc(sizeof(double*)*(ALPHABET_SIZE+1));
		if(NULL==sm->ColorScores) {
			PrintError(FnName,
					"sm->ColorScores",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		for(i=0;i<ALPHABET_SIZE+1;i++) {
			sm->ColorScores[i] = (double*)malloc(sizeof(double)*(ALPHABET_SIZE+1));
			if(NULL==sm->ColorScores[i]) {
				PrintError(FnName,
						"sm->ColorScores[i]",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}
		/* Read in the score matrix */
		for(i=0;i<ALPHABET_SIZE+1;i++) { /* Read row */
			for(j=0;j<ALPHABET_SIZE+1;j++) { /* Read column */
				if(fscanf(fp, "%lf", &sm->ColorScores[i][j])==EOF) {
					PrintError(FnName,
							scoringMatrixFileName,
							"Could not read in the scoring matrix (ColorScores)",
							Exit,
							OutOfRange);
				}
			}
		}
	}
	else {
		sm->ColorKeys=NULL;
		sm->ColorScores=NULL;
	}

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
	sm->ColorKeys=NULL;
	sm->ColorScores=NULL;
}

/* TODO */
void ScoringMatrixFree(ScoringMatrix *sm)
{
	free(sm->NTKeys);
	free(sm->NTScores);
	free(sm->ColorKeys);
	free(sm->ColorScores);
	ScoringMatrixInitialize(sm);
}

/* TODO */
double ScoringMatrixGetNTScore(char a,
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
double ScoringMatrixGetColorScore(uint8_t a,
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
