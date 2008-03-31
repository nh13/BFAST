#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "ReadInputFiles.h"
#include "Align.h"

/* TODO */
/* Note: we could do this in log-space if necessary */
/* Note: we should be aware that we are adding a lot of negative
 * numbers, and therefore we might an overflow problem if we keep
 * adding negative numbers (long sequences)
 * */
/* Return the offset of the alignment (gaps at the beginning of
 * the read are OK) in relation to the start of the reference
 * sequence 
 * */
int AlignmentGetScore(char *read,
		int readLength,
		char *reference,
		int referenceLength,
		ScoringMatrix *sm,
		AlignOutput *ao)
{
	int i, j;
	double curMismatchScore;
	double hScore, vScore, dScore, maxScore;
	int endRow, endCol;
	int numRows=readLength+1;
	int numCols=referenceLength+1;
	MatrixEntry **Entries; /* store the dynamic programming array */

	/* Allocate memory - remember to include an extra row and column */
	Entries = (MatrixEntry**)malloc(sizeof(MatrixEntry*)*numRows);
	for(i=0;i<numRows;i++) {
		Entries[i] = (MatrixEntry*)malloc(sizeof(MatrixEntry)*numCols);
	}

	/* Initialize:
	 * 
	 * Since we are aligning a read to the reference, we want a local 
	 * alignment within the reference that uses the entire read.  
	 * */
	/* Scoring Matrix: intialize the first columns for all rows except for the first */
	for(i=1;i<numRows;i++) {
		Entries[i][0].hScore = NEGATIVE_INFINITY;
		Entries[i][0].vScore = NEGATIVE_INFINITY;
		Entries[i][0].dScore = NEGATIVE_INFINITY; 
		Entries[i][0].prevRow = -1;
		Entries[i][0].prevCol = -1;
	}
	/* Scoring Matrix: initialize columns for the first row */
	for(j=0;j<numCols;j++) {
			Entries[0][j].hScore = NEGATIVE_INFINITY;
			Entries[0][j].vScore = NEGATIVE_INFINITY;
			Entries[0][j].dScore = 0;
			Entries[0][j].prevRow = -1;
			Entries[0][j].prevCol = -1;
	}

	/* Perform alignment - Dynamic programming */
	for(i=1;i<numRows;i++) { /* for each row */
		for(j=1;j<numCols;j++) { /* for each column */
			/* i = row, j = col */

			/* Update horizontal */
			/* Take the max of either extending horizontally, or starting a new gap from the previous diagonal */
			Entries[i][j].hScore = GetMaximumOfTwoDoubles(Entries[i][j-1].hScore + sm->gapExtensionPenalty, Entries[i][j-1].dScore + sm->gapExtensionPenalty + sm->gapOpenPenalty);

			/* Update vertical */
			/* Take the max of either extending vertically, or starting a new gap from the previous diagonal */
			Entries[i][j].vScore = GetMaximumOfTwoDoubles(Entries[i-1][j].vScore + sm->gapExtensionPenalty, Entries[i-1][j].dScore + sm->gapExtensionPenalty + sm->gapOpenPenalty);

			/* Update diagonal */
			/* Get mismatch score */
			curMismatchScore=GetScoreFromMatrix(read[i-1],
					ToLower(reference[j-1]),
					sm);
			hScore = Entries[i][j].hScore;
			vScore = Entries[i][j].vScore;
			dScore = Entries[i-1][j-1].dScore + curMismatchScore;
			/* Get the maximum score of the three cases: horizontal, vertical and diagonal */
			/* Intialize the maximum score to be the horizontal */
			maxScore = hScore;
			/* Initialize the previous cell to be from the previous column */
			Entries[i][j].prevRow = i;
			Entries[i][j].prevCol = j-1;
			/* Check if vertical is greater */
			if(vScore > maxScore) {
				maxScore = vScore;
				maxScore = dScore;
				/* Update the previous cell to be from the previous row */
				Entries[i][j].prevRow = i-1;
				Entries[i][j].prevCol = j;
			}
			/* Check if the diagonal is greater */
			if(dScore > maxScore) {
				maxScore = dScore;
				/* Update the previous cell to be from the previous diagonal */
				Entries[i][j].prevRow = i-1;
				Entries[i][j].prevCol = j-1;
			}
			/* Set the score for the diagonal */
			Entries[i][j].dScore = maxScore;

		}
	}

	/* Get the best alignment.  We can find the best score in the last row and then
	 * trace back.  We choose the best score from the last row since we want to 
	 * align the read completely and only locally to the reference. */
	/* Initialize the best score to be from the first column and last row */ 
	endRow = numRows-1;
	endCol = 1;
	maxScore = Entries[endRow][endCol].dScore; 
	for(i=endCol+1;i<numCols;i++) { /* For each column */
		/* Check to see if the current cell has a better score */ 
		if(Entries[endRow][i].dScore > maxScore) {
			endCol = i;
			maxScore = Entries[endRow][endCol].dScore; 
		}
	}

	/* Store results */
	/* TODO */
	
	curCol=endCol;
	curRow=endRow;
	while(curCol != -1 && curRow != -1) {
		/* Update the row and column */
		tempCol = curCol;
		curCol = Entries[curRow][curCol];
		curRow = Entries[curRow][tempCol];
	}

	/* Free memory */
	for(i=0;i<numRows;i++) {
		free(Entries[i]);
	}
	free(Entries);

	return 0.0;
}

/* TODO */
double GetScoreFromMatrix(char a, 
		char b, 
		ScoringMatrix *sm)
{
	int indexA;
	int indexB;

	/* Get index for a */
	switch(a) {
		case 'a':
			indexA=0;
			break;
		case 'c':
			indexA=1;
			break;
		case 'g':
			indexA=2;
			break;
		case 't':
			indexA=3;
			break;
		default:
			fprintf(stderr, "Error.  GetScoreFromMatrix could not understand [%c].  Terminating!\n",
					a);
			exit(1);
			break;
	}

	/* Get index for b */
	switch(b) {
		case 'a':
			indexB=0;
			break;
		case 'c':
			indexB=1;
			break;
		case 'g':
			indexB=2;
			break;
		case 't':
			indexB=3;
			break;
		default:
			fprintf(stderr, "Error.  GetScoreFromMatrix could not understand [%c].  Terminating!\n",
					a);
			exit(1);
			break;
	}

	return sm->scores[indexA][indexB];
}

/* TODO */
double GetMaximumOfTwoDoubles(double a, double b)
{
	return ((a<b)?a:b); 
}
