#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/AlignEntry.h"
#include "ScoringMatrix.h"
#include "Align.h"
#include "AlignNTSpace.h"

int AlignNTSpace(char *read,
		int readLength,
		char *reference,
		int referenceLength,
		ScoringMatrix *sm,
		AlignEntry *a,
		int type)
{
	char *FnName="AlignNTSpace";

	switch(type) {
		case MismatchesOnly:
			return AlignNTSpaceMismatchesOnly(read,
					readLength,
					reference,
					referenceLength,
					sm,
					a);
			break;
		case FullAlignment:
			return AlignNTSpaceFull(read,
					readLength,
					reference,
					referenceLength,
					sm,
					a);
			break;
		default:
			PrintError(FnName,
					"type",
					"Could not understand alignment type",
					Exit,
					OutOfRange);
	}
	return -1;
}

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
int AlignNTSpaceFull(char *read,
		int readLength,
		char *reference,
		int referenceLength,
		ScoringMatrix *sm,
		AlignEntry *a)
{
	char *FnName = "AlignNTSpaceFull";
	/* read goes on the rows, reference on the columns */
	AlignMatrix **matrix=NULL;
	int offset = 0;
	int i, j, k, l;

	/* Allocate memory for the matrix */
	matrix = malloc(sizeof(AlignMatrix*)*(readLength+1));
	if(NULL==matrix) {
		PrintError(FnName,
				"matrix",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	for(i=0;i<readLength+1;i++) {
		matrix[i] = malloc(sizeof(AlignMatrix)*(referenceLength+1));
		if(NULL==matrix[i]) {
			PrintError(FnName,
					"matrix[i]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}

	/* Initialize "the matrix" */
	/* Row i (i>0) column 0 should be negative infinity since we want to
	 * align the full read */
	for(i=1;i<readLength+1;i++) {
		for(k=0;k<ALIGNMATRIXCELL_NUM_SUB_CELLS;k++) {
			matrix[i][0].score[k] = NEGATIVE_INFINITY;
			matrix[i][0].from[k] = Start;
			matrix[i][0].length[k] = 0;
			matrix[i][0].colorError[k] = '0';
		}
	}
	/* Row 0 column j should be zero since we want to find the best
	 * local alignment within the reference */
	for(j=0;j<referenceLength+1;j++) {
		for(k=0;k<ALIGNMATRIXCELL_NUM_SUB_CELLS;k++) {
			if(k==0) {
				matrix[0][j].score[k] = 0.0;
			}
			else {
				matrix[0][j].score[k] = NEGATIVE_INFINITY;
			}
			matrix[0][j].from[k] = Start;
			matrix[0][j].length[k] = 0;
			matrix[0][j].colorError[k] = '0';
		}
	}

	/* Fill in the matrix according to the recursive rules */
	for(i=0;i<readLength;i++) { /* read/rows */
		for(j=0;j<referenceLength;j++) { /* reference/columns */
			/* No color errors */
			matrix[i+1][j+1].colorError[0] = '0';
			matrix[i+1][j+1].colorError[1] = '0';
			matrix[i+1][j+1].colorError[2] = '0';

			/* Deletion relative to reference across a column */
			/* Insertion relative to reference is down a row */
			/* Match/Mismatch is a diagonal */

			/* Update deletion */
			/* Deletion extension */
			matrix[i+1][j+1].score[1] = matrix[i+1][j].score[1] + sm->gapExtensionPenalty; 
			matrix[i+1][j+1].length[1] = matrix[i+1][j].length[1] + 1;
			matrix[i+1][j+1].from[1] = DeletionExt;
			/* Check if starting a new deletion is better */
			if(matrix[i+1][j].score[0] + sm->gapOpenPenalty > matrix[i+1][j+1].score[1]) {
				matrix[i+1][j+1].score[1] = matrix[i+1][j].score[0] + sm->gapOpenPenalty;
				matrix[i+1][j+1].length[1] = matrix[i+1][j].length[0] + 1;
				matrix[i+1][j+1].from[1] = DeletionA;
			}

			/* Update insertion */
			/* Insertion extension */
			matrix[i+1][j+1].score[2] = matrix[i][j+1].score[2] + sm->gapExtensionPenalty; 
			matrix[i+1][j+1].length[2] = matrix[i][j+1].length[2] + 1;
			matrix[i+1][j+1].from[2] = InsertionExt;
			/* Check if starting a new insertion is better */
			if(matrix[i][j+1].score[0] + sm->gapOpenPenalty > matrix[i+1][j+1].score[2]) {
				matrix[i+1][j+1].score[2] = matrix[i][j+1].score[0] + sm->gapOpenPenalty;
				matrix[i+1][j+1].length[2] = matrix[i][j+1].length[0] + 1;
				matrix[i+1][j+1].from[2] = InsertionA;
			}

			/* Update diagonal */
			/* Get mismatch score */
			matrix[i+1][j+1].score[0] = matrix[i][j].score[0] + ScoringMatrixGetNTScore(read[i], reference[j], sm);
			matrix[i+1][j+1].length[0] = matrix[i][j].length[0] + 1;
			matrix[i+1][j+1].from[0] = DiagA;
			/* Get the maximum score of the three cases: horizontal, vertical and diagonal */
			if(matrix[i+1][j+1].score[1] > matrix[i+1][j+1].score[0]) {
				matrix[i+1][j+1].score[0] = matrix[i+1][j+1].score[1];
				matrix[i+1][j+1].length[0] = matrix[i+1][j+1].length[1];
				matrix[i+1][j+1].from[0] = matrix[i+1][j+1].from[1];
			}
			if(matrix[i+1][j+1].score[2] > matrix[i+1][j+1].score[0]) {
				matrix[i+1][j+1].score[0] = matrix[i+1][j+1].score[2];
				matrix[i+1][j+1].length[0] = matrix[i+1][j+1].length[2];
				matrix[i+1][j+1].from[0] = matrix[i+1][j+1].from[2];
			}

			/* Update the rest to negative infinity */
			for(l=3;l<ALIGNMATRIXCELL_NUM_SUB_CELLS;l++) {
				matrix[i+1][j+1].score[l] = NEGATIVE_INFINITY;
				matrix[i+1][j+1].length[l] = -1;
				matrix[i+1][j+1].from[l] = Start;
			}
		}
	}

	offset = FillAlignEntryFromMatrix(a,
			matrix,
			read,
			readLength,
			reference,
			referenceLength,
			0,
			0);

	/* Free the matrix, free your mind */
	for(i=0;i<readLength+1;i++) {
		free(matrix[i]);
		matrix[i]=NULL;
	}
	free(matrix);
	matrix=NULL;

	/* The return is the number of gaps at the beginning of the reference */
	return offset;
}

int AlignNTSpaceMismatchesOnly(char *read,
		int readLength,
		char *reference,
		int referenceLength,
		ScoringMatrix *sm,
		AlignEntry *a)
{
	char *FnName = "AlignNTSpaceMismatchesOnly";
	/* Read goes on the second row, reference on the first */
	int i, j;
	int32_t maxScore = NEGATIVE_INFINITY;
	int offset=-1;
	int32_t curScore = 0.0;
	
	assert(readLength <= referenceLength);
	
	for(i=0;i<referenceLength-readLength+1;i++) { /* Starting position */
		curScore = 0.0;
		for(j=0;j<readLength;j++) { /* Position in the alignment */
			curScore += ScoringMatrixGetNTScore(read[j], reference[i+j], sm);
		}
		if(maxScore < curScore) {
			maxScore = curScore;
			offset = i;
		}
	}

	/* Copy over */
	a->referenceLength = readLength;
	a->length = readLength;
	a->score = maxScore;
	/* Allocate memory */
	assert(NULL==a->read);
	a->read = malloc(sizeof(char)*(a->length+1));
	if(NULL==a->read) {
		PrintError(FnName,
				"a->read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(NULL==a->reference);
	a->reference = malloc(sizeof(char)*(a->length+1));
	if(NULL==a->reference) {
		PrintError(FnName,
				"a->reference",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(NULL==a->colorError);
	/*
	   a->colorError = malloc(sizeof(char)*SEQUENCE_LENGTH);
	   if(NULL==a->colorError) {
	   PrintError(FnName,
	   "a->colorError",
	   "Could not allocate memory",
	   Exit,
	   MallocMemory);
	   }
	   */
	/* Copy over */
	for(i=0;i<a->length;i++) {
		a->read[i] = read[i];
		a->reference[i] = reference[i+offset];
	}
	a->read[a->length] = '\0';
	a->reference[a->length] = '\0';
	
	return offset;
}
