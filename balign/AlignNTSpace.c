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
	AlignMatrixNT **matrix=NULL;
	int offset = 0;
	int i, j;

	/* Allocate memory for the matrix */
	matrix = malloc(sizeof(AlignMatrixNT*)*(readLength+1));
	if(NULL==matrix) {
		PrintError(FnName,
				"matrix",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	for(i=0;i<readLength+1;i++) {
		matrix[i] = malloc(sizeof(AlignMatrixNT)*(referenceLength+1));
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
		matrix[i][0].h.score = matrix[i][0].s.score = matrix[i][0].v.score = NEGATIVE_INFINITY;
		matrix[i][0].h.from = matrix[i][0].s.from = matrix[i][0].v.from = StartNT;
		matrix[i][0].h.length = matrix[i][0].s.length = matrix[i][0].v.length = 0;
	}
	/* Row 0 column j should be zero since we want to find the best
	 * local alignment within the reference */
	for(j=0;j<referenceLength+1;j++) {
		matrix[0][j].s.score = 0;
		matrix[0][j].h.score = matrix[0][j].v.score = NEGATIVE_INFINITY;
		matrix[0][j].h.from = matrix[0][j].s.from = matrix[0][j].v.from = StartNT;
		matrix[0][j].h.length = matrix[0][j].s.length = matrix[0][j].v.length = 0;
	}

	/* Fill in the matrix according to the recursive rules */
	for(i=0;i<readLength;i++) { /* read/rows */
		for(j=0;j<referenceLength;j++) { /* reference/columns */
			/* Deletion relative to reference across a column */
			/* Insertion relative to reference is down a row */
			/* Match/Mismatch is a diagonal */

			/* Update deletion */
			/* Deletion extension */
			matrix[i+1][j+1].h.score = matrix[i+1][j].h.score + sm->gapExtensionPenalty; 
			matrix[i+1][j+1].h.length = matrix[i+1][j].h.length + 1;
			matrix[i+1][j+1].h.from = DeletionExtension;
			/* Check if starting a new deletion is better */
			if(matrix[i+1][j+1].h.score < matrix[i+1][j].s.score + sm->gapOpenPenalty) {
				matrix[i+1][j+1].h.score = matrix[i+1][j].s.score + sm->gapOpenPenalty;
				matrix[i+1][j+1].h.length = matrix[i+1][j].s.length + 1;
				matrix[i+1][j+1].h.from = DeletionStart;
			}

			/* Update insertion */
			/* Insertion extension */
			matrix[i+1][j+1].v.score = matrix[i][j+1].v.score + sm->gapExtensionPenalty; 
			matrix[i+1][j+1].v.length = matrix[i][j+1].v.length + 1;
			matrix[i+1][j+1].v.from = InsertionExtension;
			/* Check if starting a new insertion is better */
			if(matrix[i+1][j+1].v.score < matrix[i][j+1].s.score + sm->gapOpenPenalty) {
				matrix[i+1][j+1].v.score = matrix[i][j+1].s.score + sm->gapOpenPenalty;
				matrix[i+1][j+1].v.length = matrix[i][j+1].s.length + 1;
				matrix[i+1][j+1].v.from = InsertionStart;
			}

			/* Update diagonal */
			/* Get mismatch score */
			matrix[i+1][j+1].s.score = matrix[i][j].s.score + ScoringMatrixGetNTScore(read[i], reference[j], sm);
			matrix[i+1][j+1].s.length = matrix[i][j].s.length + 1;
			matrix[i+1][j+1].s.from = Match;
			/* Get the maximum score of the three cases: horizontal, vertical and diagonal */
			if(matrix[i+1][j+1].s.score < matrix[i+1][j+1].h.score) {
				matrix[i+1][j+1].s.score = matrix[i+1][j+1].h.score;
				matrix[i+1][j+1].s.length = matrix[i+1][j+1].h.length;
				matrix[i+1][j+1].s.from = matrix[i+1][j+1].h.from;
			}
			if(matrix[i+1][j+1].s.score < matrix[i+1][j+1].v.score) {
				matrix[i+1][j+1].s.score = matrix[i+1][j+1].v.score;
				matrix[i+1][j+1].s.length = matrix[i+1][j+1].v.length;
				matrix[i+1][j+1].s.from = matrix[i+1][j+1].v.from;
			}
		}
	}

	offset = FillAlignEntryFromMatrixNTSpace(a,
			matrix,
			read,
			readLength,
			reference,
			referenceLength,
			NTSpace,
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

int FillAlignEntryFromMatrixNTSpace(AlignEntry *a,
		AlignMatrixNT **matrix,
		char *read,
		int readLength,
		char *reference,
		int referenceLength,
		int scoringType,
		int debug)
{
	char *FnName="FillAlignEntryFromMatrixNTSpace";
	int curRow, curCol, startRow, startCol;
	char curReadBase;
	int nextRow, nextCol;
	char nextReadBase;
	int curFrom;
	double maxScore;
	int i;
	int offset;

	curReadBase = nextReadBase = 'X';
	nextRow = nextCol = -1;

	/* Get the best alignment.  We can find the best score in the last row and then
	 * trace back.  We choose the best score from the last row since we want to 
	 * align the read completely and only locally to the reference. */
	startRow=-1;
	startCol=-1;
	maxScore = NEGATIVE_INFINITY;
	for(i=0;i<referenceLength+1;i++) {
		/* Check only the first cell */
		if(maxScore < matrix[readLength][i].s.score) {
			maxScore = matrix[readLength][i].s.score;
			startRow = readLength;
			startCol = i;
		}
	}
	assert(startRow >= 0 && startCol >= 0);

	/* Initialize variables for the loop */
	curRow=startRow;
	curCol=startCol;
	curFrom = Match;
	a->referenceLength=0;
	i=matrix[curRow][curCol].s.length-1; /* Get the length of the alignment */
	a->length=matrix[curRow][curCol].s.length; /* Copy over the length */
	a->score = maxScore; /* Copy over score */

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

	/* Now trace back the alignment using the "from" member in the matrix */
	while(curRow > 0 && curCol > 0) {

		/* Where did the current cell come from */
		switch(curFrom) {
			case DeletionStart:
				curFrom = matrix[curRow][curCol].s.from;
				assert(curFrom == Match || curFrom == InsertionExtension);
				break;
			case DeletionExtension:
				curFrom = matrix[curRow][curCol].h.from;
				assert(curFrom == DeletionStart || curFrom == DeletionExtension);
				break;
			case Match:
				curFrom = matrix[curRow][curCol].s.from;
				break;
			case InsertionStart:
				curFrom = matrix[curRow][curCol].s.from;
				assert(curFrom == Match || curFrom == DeletionExtension);
				break;
			case InsertionExtension:
				curFrom = matrix[curRow][curCol].v.from;
				assert(curFrom == InsertionStart || curFrom == InsertionExtension);
				break;
			default:
				PrintError(FnName,
						"curFrom",
						"Could not recognize curFrom",
						Exit,
						OutOfRange);
		}

		assert(i>=0);

		/* Update alignment */
		switch(curFrom) {
			case DeletionStart:
			case DeletionExtension:
				a->read[i] = GAP;
				a->reference[i] = reference[curCol-1];
				a->referenceLength++;
				nextRow = curRow;
				nextCol = curCol-1;
				break;
			case Match:
				a->read[i] = read[curRow-1];
				a->reference[i] = reference[curCol-1];
				a->referenceLength++;
				nextRow = curRow-1;
				nextCol = curCol-1;
				break;
			case InsertionStart:
			case InsertionExtension:
				a->read[i] = read[curRow-1];
				a->reference[i] = GAP;
				nextRow = curRow-1;
				nextCol = curCol;
				break;
			default:
				fprintf(stderr, "curFrom=%d\n", curFrom);
				PrintError(FnName,
						"curFrom",
						"Could not understand curFrom",
						Exit,
						OutOfRange);
		}

		assert(a->read[i] != GAP || a->read[i] != a->reference[i]);

		/* Update for next loop iteration */
		curRow = nextRow;
		curCol = nextCol;
		i--;

		assert(i!=-1 || a->read[0] != GAP);
	} /* End Loop */
	assert(-1==i);
	assert(a->length >= a->referenceLength);

	a->read[a->length]='\0';
	a->reference[a->length]='\0';

	offset = curCol;

	return offset;
}
