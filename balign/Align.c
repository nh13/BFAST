#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "ReadInputFiles.h"
#include "../blib/BError.h"
#include "../blib/AlignEntry.h"
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
		AlignEntry *aEntry)
{
	int i, j;
	double curMismatchScore;
	double hScore, vScore, dScore, maxScore;
	int endRow, endCol;
	int numRows=readLength+1;
	int numCols=referenceLength+1;
	int curCol=-1;
	int curRow=-1;
	int prevRow, prevCol;
	int offset=0;
	MatrixEntry **Entries; /* store the dynamic programming array */

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "In AlignmentGetScore\nread:[%s] length [%d,%d]\nreference:[%s] length [%d,%d]\n",
				read,
				readLength,
				(int)strlen(read),
				reference,
				referenceLength,
				(int)strlen(reference));
		/*
		   for(i=0;i<ALPHABET_SIZE;i++) {
		   for(j=0;j<ALPHABET_SIZE;j++) {
		   fprintf(stderr, "%lf\t", sm->scores[i][j]);
		   }
		   fprintf(stderr, "\n");
		   }
		   */
	}

	/* Allocate memory - remember to include an extra row and column */
	Entries = (MatrixEntry**)malloc(sizeof(MatrixEntry*)*numRows);
	if(NULL == Entries) {
		PrintError("AlignmentGetScore",
				"Entries",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	for(i=0;i<numRows;i++) {
		Entries[i] = (MatrixEntry*)malloc(sizeof(MatrixEntry)*numCols);
	if(NULL == Entries[i]) {
		PrintError("AlignmentGetScore",
				"Entries[i]",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
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
				/* Update the previous cell to be from the previous row */
				Entries[i][j].prevRow = i-1;
				Entries[i][j].prevCol = j;
			}
			/* Check if the diagonal is greater */
			if(dScore >= maxScore) { /* greater or equal so diagonal always takes precedence! */
				maxScore = dScore;
				/* Update the previous cell to be from the previous diagonal */
				Entries[i][j].prevRow = i-1;
				Entries[i][j].prevCol = j-1;
			}
			/* Set the score for the diagonal */
			Entries[i][j].dScore = maxScore;

		}
	}

	if(VERBOSE >= DEBUG) {
		/*
		for(i=0;i<numRows;i++) { 
			for(j=0;j<numCols;j++) { 
				fprintf(stderr, "(row,col)=[%d,%d]\t(v,d,h)=[%.2lf,%.2lf,%.2lf]\t(prevRow,prevCol)=[%d,%d]\n",
						i,
						j,
						Entries[i][j].vScore,
						Entries[i][j].dScore,
						Entries[i][j].hScore,
						Entries[i][j].prevRow,
						Entries[i][j].prevCol);
			}
		}
		*/
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
	aEntry->score = maxScore;

	/* First allocate the maximum length of the alignment, we can update later */
	aEntry->read = (char*)malloc(sizeof(char)*(endRow+endCol+1));
	if(NULL==aEntry->read) {
		PrintError("AlignmentGetScore",
				"aEntry->read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	aEntry->reference = (char*)malloc(sizeof(char)*(endRow+endCol+1));
	if(NULL==aEntry->reference) {
		PrintError("AlignmentGetScore",
				"aEntry->reference",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	aEntry->length = 0; /* initialize length */
	prevRow=endRow;
	prevCol=endCol;
	while(prevRow != -1 && prevCol != -1) {
		assert(prevRow < numRows && prevCol < numCols);
		curRow = Entries[prevRow][prevCol].prevRow;
		curCol = Entries[prevRow][prevCol].prevCol;
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "(curRow,curCol)=[%d,%d]\t(%c,%c)\n",
					curRow,
					curCol,
					read[curRow],
					reference[curCol]);
		}
		assert(curRow < readLength || curCol < referenceLength);
		if(curRow >= readLength) {
			assert(curRow != -1 && curCol != -1);
			aEntry->read[aEntry->length] = '-';
			aEntry->reference[aEntry->length] = reference[curCol];
			/* Update the offset */
			offset=curCol;
			/* Update the length */
			aEntry->length++;
		}
		else if(curCol >= referenceLength) {
			assert(curRow != -1 && curCol != -1);
			aEntry->read[aEntry->length] = read[curRow];
			aEntry->reference[aEntry->length] = '-';
			/* Update the offset */
			offset=curCol;
			/* Update the length */
			aEntry->length++;
		}
		else if(curRow != -1 && curCol != -1 && curRow < readLength && curCol < referenceLength) {

			/* We have three cases:
			 * 1. Vertical move: curRow == prevRow && curCol == prevCol-1
			 * 2. Horizontal move: curRow == prevRow-1 && curCol == prevCol
			 * 3. Diagonal move: curRow == prevRow-1 && curCol == prevCol-1
			 * */

			if(curRow == prevRow && curCol == prevCol-1) {
				/* Vertical - gap in reference */
				aEntry->read[aEntry->length] = '-';
				aEntry->reference[aEntry->length] = reference[curCol];
			}
			else if(curRow == prevRow-1 && curCol == prevCol) {
				/* Horizontal - gap in read */
				aEntry->read[aEntry->length] = read[curCol];
				aEntry->reference[aEntry->length] = '-';
			}
			else if(curRow == prevRow-1 && curCol == prevCol-1) {
				/* Diagonal */
				aEntry->read[aEntry->length] = read[curRow]; 
				aEntry->reference[aEntry->length] = reference[curCol];
			}
			else {
				PrintError("AlignmentGetScore",
						NULL,
						"Could not recover alignment path",
						Exit,
						OutOfRange);
			}

		/* Update the offset */
			offset=curCol;
			/* Update the length */
			aEntry->length++;
		}

		/* Update the row and column */
		prevRow=curRow;
		prevCol=curCol;
	}
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "\n");
	}
	/* Reallocate memory for the read and reference */
	aEntry->read = (char*)realloc(aEntry->read, sizeof(char)*(aEntry->length+1));
	if(NULL==aEntry->read) {
		PrintError("AlignmentGetScore",
				"aEntry->read",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}
	aEntry->read[aEntry->length]='\0'; /* null terminator */
	aEntry->reference = (char*)realloc(aEntry->reference, sizeof(char)*(aEntry->length+1));
	if(NULL==aEntry->reference) {
		PrintError("AlignmentGetScore",
				"aEntry->reference",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}
	aEntry->reference[aEntry->length]='\0'; /* null terminator */
	/* Reverse read and reference since we got them from the end of the path */
	ReverseSequence(aEntry->read, aEntry->length);
	ReverseSequence(aEntry->reference, aEntry->length);

	/* Free memory */
	for(i=0;i<numRows;i++) {
		free(Entries[i]);
	}
	free(Entries);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "aligned read:%s\naligned reference:%s\n",
				aEntry->read,
				aEntry->reference);
		fprintf(stderr, "Exiting AlignmentGetScore returning offset %d\n", offset);
	}
	/* The return is the number of gaps at the beginning of the reference */
	return offset;
}

/* TODO */
double GetScoreFromMatrix(char a, 
		char b, 
		ScoringMatrix *sm)
{
	int indexA=-1;
	int indexB=-1;

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
		case 'n':
			indexA=4;
			break;
		default:
			PrintError("GetScoreFromMatrix",
					NULL,
					"Could not understand 'a' key",
					Exit,
					OutOfRange);
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
		case 'n':
			indexB=4;
			break;
		default:
			PrintError("GetScoreFromMatrix",
					NULL,
					"Could not understand 'b' key",
					Exit,
					OutOfRange);
			break;
	}

	return sm->scores[indexA][indexB];
}

/* TODO */
double GetMaximumOfTwoDoubles(double a, double b)
{
	return ((a>b)?a:b); 
}

/* TODO */
void ReverseSequence(char *a, int length) 
{
	int i, j;
	char c;

	/* In-place */
	for(i=0,j=length-1;i<j;i++,j--) {
		c = a[i];
		a[i] = a[j];
		a[j] = c;
	}
}
