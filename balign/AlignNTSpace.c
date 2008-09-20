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
int AlignNTSpace(char *read,
		int readLength,
		char *reference,
		int referenceLength,
		ScoringMatrix *sm,
		AlignEntry *aEntry)
{
	char *FnName = "AlignNTSpace";
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
			matrix[0][j].score[k] = 0.0;
			matrix[0][j].from[k] = Start;
			matrix[0][j].length[k] = 0;
			matrix[i][0].colorError[k] = '0';
		}
	}

	/* Fill in the matrix according to the recursive rules */
	for(i=0;i<readLength;i++) { /* read/rows */
		for(j=0;j<referenceLength;j++) { /* reference/columns */
			/* No color errors */
			matrix[i+1][i+1].colorError[0] = '0';
			matrix[i+1][i+1].colorError[1] = '0';
			matrix[i+1][i+1].colorError[2] = '0';

			/* Deletion is down a row, insertion is across a column, and
			 * match/mismatch is a diagonal */

			/* Update deletion */
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
				matrix[i+1][j+1].from[0] = DeletionA;
			}
			if(matrix[i+1][j+1].score[2] > matrix[i+1][j+1].score[0]) {
				matrix[i+1][j+1].score[0] = matrix[i+1][j+1].score[2];
				matrix[i+1][j+1].length[0] = matrix[i+1][j+1].length[2];
				matrix[i+1][j+1].from[0] = InsertionA;
			}

			/* Update the rest to negative infinity */
			for(l=3;l<ALIGNMATRIXCELL_NUM_SUB_CELLS;l++) {
				matrix[i+1][j+1].score[l] = NEGATIVE_INFINITY;
				matrix[i+1][j+1].length[l] = -1;
				matrix[i+1][j+1].from[l] = Start;
			}
			/* HERE */
			/*
			   for(l=0;l<ALIGNMATRIXCELL_NUM_SUB_CELLS;l++) {
			   fprintf(stderr, "(i,j,l,score)=(%d,%d,%d,%lf)\n",
			   i+1,
			   j+1,
			   l,
			   matrix[i+1][j+1].score[l]);
			   }
			   */
		}
	}

	/* HERE */
	/*
	   fprintf(stderr, "read=%s\nreference=%s\n",
	   read,
	   reference);
	   */

	offset = FillAlignEntryFromMatrix(aEntry,
			matrix,
			read,
			readLength,
			reference,
			referenceLength,
			0,
			0);

	/* HERE */
	/*
	   fprintf(stderr, "(%d,%d)\n",
	   aEntry->length,
	   (int)strlen(aEntry->read));
	   fprintf(stderr, "aEntry->reference=[%s]\naEntry->read=[%s]\n%lf\n%d\n\n", aEntry->reference, aEntry->read, aEntry->score, offset);
	   exit(1);
	   */

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
