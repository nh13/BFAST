#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "ReadInputFiles.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/AlignEntry.h"
#include "Align.h"

#define AlignmentGetScore_DEBUG_ON 0

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
	char *FnName = "AlignmentGetScore";
	int i, j;
	double curMismatchScore;
	double maxScore;
	int endRow, endCol;
	int numRows=readLength+1;
	int numCols=referenceLength+1;
	int curCol=-1;
	int curRow=-1;
	int startRow, startCol;
	int type=-1;
	int offset=0;
	MatrixEntry **Entries; /* store the dynamic programming array */
	int *path=NULL;
	int pathLength=0;

	/* Allocate memory - remember to include an extra row and column */
	Entries = (MatrixEntry**)malloc(sizeof(MatrixEntry*)*numRows);
	if(NULL == Entries) {
		PrintError(FnName,
				"Entries",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	for(i=0;i<numRows;i++) {
		Entries[i] = (MatrixEntry*)malloc(sizeof(MatrixEntry)*numCols);
		if(NULL == Entries[i]) {
			PrintError(FnName,
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
		Entries[i][0].dType = NO_TYPE;
	}
	/* Scoring Matrix: initialize columns for the first row */
	for(j=0;j<numCols;j++) {
		Entries[0][j].hScore = NEGATIVE_INFINITY;
		Entries[0][j].vScore = NEGATIVE_INFINITY;
		Entries[0][j].dScore = 0;
		Entries[0][j].prevRow = -1;
		Entries[0][j].prevCol = -1;
		Entries[0][j].dType = NO_TYPE;
	}

	/* Perform alignment - Dynamic programming */
	for(i=1;i<numRows;i++) { /* for each row */
		for(j=1;j<numCols;j++) { /* for each column */
			/* i = row, j = col */

			/* Update horizontal */
			/* Take the max of either extending horizontally, or starting a new gap from the previous diagonal */
			GetGapScore(Entries[i][j-1].hScore + sm->gapExtensionPenalty,
					Entries[i][j-1].dScore + sm->gapOpenPenalty,
					&Entries[i][j].hScore,
					&Entries[i][j].hType);

			if(i==7 && j==5) {
			}

			/* Update vertical */
			/* Take the max of either extending vertically, or starting a new gap from the previous diagonal */
			GetGapScore(Entries[i-1][j].vScore + sm->gapExtensionPenalty,
					Entries[i-1][j].dScore + sm->gapOpenPenalty,
					&Entries[i][j].vScore,
					&Entries[i][j].vType);

			/* Update diagonal */
			/* Get mismatch score */
			curMismatchScore=GetScoreFromMatrix(read[i-1],
					ToLower(reference[j-1]),
					sm);
			/* Get the maximum score of the three cases: horizontal, vertical and diagonal */
			Entries[i][j].dScore = Entries[i-1][j-1].dScore + curMismatchScore;
			Entries[i][j].dType = DIAGONAL;
			Entries[i][j].prevRow = i-1;
			Entries[i][j].prevCol = j-1;
			if(Entries[i][j].hScore > Entries[i][j].dScore) {
				Entries[i][j].dScore = Entries[i][j].hScore;
				Entries[i][j].dType = HORIZONTAL;
				Entries[i][j].prevRow = i;
				Entries[i][j].prevCol = j-1;
			}
			if(Entries[i][j].vScore > Entries[i][j].dScore) {
				Entries[i][j].dScore = Entries[i][j].vScore;
				Entries[i][j].dType = VERTICAL;
				Entries[i][j].prevRow = i-1;
				Entries[i][j].prevCol = j;
			}
		}
	}

	/* Debugging */
	if(AlignmentGetScore_DEBUG_ON) {
		for(i=0;i<numRows;i++) { 
			for(j=0;j<numCols;j++) { 
				fprintf(stderr, "(row,col)=[%d,%d]\t(v,d,h)=[%.2lf,%.2lf,%.2lf]\t(prevRow,prevCol)=[%d,%d]\t(vType, dType, hType)=[%d,%d,%d]",
						i,
						j,
						Entries[i][j].vScore,
						Entries[i][j].dScore,
						Entries[i][j].hScore,
						Entries[i][j].prevRow,
						Entries[i][j].prevCol,
						Entries[i][j].vType,
						Entries[i][j].dType,
						Entries[i][j].hType
					   );
				if(j>0 && i>0) {
					fprintf(stderr, "\t(ref,read)=[%c,%c]",
							reference[j-1],
							read[i-1]);
				}
				fprintf(stderr, "\n");
			}
		}
		fprintf(stderr, "[%lf,%lf]\n", sm->gapOpenPenalty, sm->gapExtensionPenalty);
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
	aEntry->read = malloc(sizeof(char)*(endRow+endCol+1));
	if(NULL==aEntry->read) {
		PrintError(FnName,
				"aEntry->read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	aEntry->reference = malloc(sizeof(char)*(endRow+endCol+1));
	if(NULL==aEntry->reference) {
		PrintError(FnName,
				"aEntry->reference",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Allocate memory for the path */
	pathLength = 0;
	path = malloc(sizeof(int)*(endRow + endCol));
	if(NULL == path) {
		PrintError(FnName,
				"path",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	if(AlignmentGetScore_DEBUG_ON == 1) {
		fprintf(stderr, "(read,ref)=[%s,%s]\n",
				read,
				reference);
	}
	curRow=endRow;
	curCol=endCol;
	type = Entries[curRow][curCol].dType;
	curRow--;
	curCol--;
	i=0;
	while(type != NO_TYPE && curRow != -1 && curCol != -1) {
		if(AlignmentGetScore_DEBUG_ON == 1) {
			fprintf(stderr, "i:%d\t(curRow,curCol)=[%d,%d]\ttype:%d\t", 
					i,
					curRow,
					curCol,
					type);
			i++;
		}

		/* What was the previous type ? */
		switch(type) {
			case DIAGONAL:
				path[pathLength] = DIAGONAL;
				pathLength++;
				/* Update row and column */
				switch(Entries[curRow][curCol].dType) {
					case HORIZONTAL:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "Diagonal[horizontal]");
						}
						curCol--;
						break;
					case VERTICAL:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "Diagonal[vertical]");
						}
						curRow--;
						break;
					case DIAGONAL:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "Diagonal[diagonal]");
						}
						curRow--;
						curCol--;
						break;
					case NO_TYPE:
						type = NO_TYPE;
						break;
					default:
						PrintError(FnName,
								NULL,
								"Could not recover alignment path: dType",
								Exit,
								OutOfRange);
						break;
				}
				type = Entries[curRow][curCol].dType;
				break;
			case HORIZONTAL:
				path[pathLength] = HORIZONTAL;
				pathLength++;
				/* Update row and column */
				switch(Entries[curRow][curCol].hType) {
					case GAP_OPEN:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "Horizontal[open]");
						}
						curCol--;
						type = DIAGONAL;
						break;
					case GAP_EXTENSION:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "Horizontal[extension]");
						}
						curCol--;
						type = HORIZONTAL;
						break;
					case NO_TYPE:
						type = NO_TYPE;
						break;
					default:
						PrintError(FnName,
								NULL,
								"Could not recover alignment path: hType",
								Exit,
								OutOfRange);
						break;
				}
				break;
			case VERTICAL:
				path[pathLength] = VERTICAL;
				pathLength++;
				/* Update row and column */
				switch(Entries[curRow][curCol].vType) {
					case GAP_OPEN:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "Vertical[open]");
						}
						curRow--;
						type = DIAGONAL;
						break;
					case GAP_EXTENSION:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "Vertical[extension]");
						}
						curRow--;
						type = VERTICAL;
						break;
					case NO_TYPE:
						type = NO_TYPE;
						break;
					default:
						PrintError(FnName,
								NULL,
								"Could not recover alignment path: vType",
								Exit,
								OutOfRange);
						break;
				}
				break;
			default:
				PrintError(FnName,
						NULL,
						"Could not recover alignment path: type",
						Exit,
						OutOfRange);
		}
		if(type != NO_TYPE) {
			/* Update the offset */
			offset=curCol;
		}
		if(AlignmentGetScore_DEBUG_ON == 1) {
			fprintf(stderr, "\n");
		}
	}
	path[pathLength] = -1;
	startRow = curRow;
	startCol = curCol;


	if(AlignmentGetScore_DEBUG_ON == 1) {
		fprintf(stderr, "(read,ref)=[%s,%s]\n",
				read,
				reference);
		fprintf(stderr, "(endRow,endCol)=[%d,%d\n]",
				endRow,
				endCol);
		fprintf(stderr, "(curRow,curCol)=[%d,%d]\tpathLength:%d\n",
				curRow,
				curCol,
				pathLength);
	}

	aEntry->length=0;
	curRow = startRow;
	curCol = startCol;
	for(i=pathLength-1;i>=0 && curRow < readLength;i--) {
		if(AlignmentGetScore_DEBUG_ON == 1) {
			fprintf(stderr, "type:%d\t(curRow,curCol)=[%d,%d]\t", 
					path[i],
					curRow,
					curCol);
		}
		if(!(curRow<readLength) || !(curCol < referenceLength)) {
			int nCols = 0;
			int nRows = 0;
			fprintf(stderr, "\ni:%d\n", i);
			fprintf(stderr, "start[%d,%d]\tend[%d,%d]\tcur[%d,%d]\n",
					startRow,
					startCol,
					endRow,
					endCol,
					curRow,
					curCol);
			for(i=pathLength-1;i>=0;i--) {
				fprintf(stderr, "path[%d]:%d\n",
						i,
						path[i]);
				switch(path[i]) {
					case HORIZONTAL:
						nCols++;
						break;
					case VERTICAL:
						nRows++;
						break;
					case DIAGONAL:
						nRows++;
						nCols++;
						break;
				}
			}
			fprintf(stderr, "nRows:%d\tnCols:%d\n",
					nRows,
					nCols);
			fprintf(stderr, "read:%s\nreference:%s\n",
					read,
					reference);
		}
		assert(curRow < readLength);
		assert(curCol < referenceLength);
		switch(path[i]) {
			case DIAGONAL:
				if(AlignmentGetScore_DEBUG_ON == 1) {
					fprintf(stderr, "[%c,%c]\n",
							read[curRow],
							reference[curCol]);
				}
				aEntry->read[aEntry->length] = read[curRow];
				aEntry->reference[aEntry->length] = reference[curCol];
				curRow++;
				curCol++;
				break;
			case HORIZONTAL:
				if(AlignmentGetScore_DEBUG_ON == 1) {
					fprintf(stderr, "[%c,%c]\n",
							GAP,
							reference[curCol]);
				}
				aEntry->read[aEntry->length] = GAP;
				aEntry->reference[aEntry->length] = reference[curCol];
				curCol++;
				break;
			case VERTICAL:
				if(AlignmentGetScore_DEBUG_ON == 1) {
					fprintf(stderr, "[%c,%c]\n",
							read[curRow],
							GAP);
				}
				aEntry->read[aEntry->length] = read[curRow];
				aEntry->reference[aEntry->length] = GAP;
				curRow++;
				break;
			default:
				PrintError(FnName,
						NULL,
						"Could not navigate path",
						Exit,
						OutOfRange);
		}
		aEntry->length++;
	}

	aEntry->reference[aEntry->length]='\0'; /* null terminator */
	aEntry->read[aEntry->length]='\0';

	if(AlignmentGetScore_DEBUG_ON == 1) {
		fprintf(stderr, "length:%d\n%s\n%s\n",
				aEntry->length,
				aEntry->read,
				aEntry->reference);
	}
	assert(strlen(aEntry->read) == strlen(aEntry->reference));

	/* Free memory */
	free(path);
	for(i=0;i<numRows;i++) {
		free(Entries[i]);
	}
	free(Entries);

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
			fprintf(stderr, "b key:[%c]\n", b);
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

void GetGapScore(double gapScore,
		double diagScore,
		double *curGapScore,
		int *curType)
{
	if(gapScore > diagScore) {
		(*curGapScore) = gapScore;
		(*curType) = GAP_EXTENSION;
	}
	else {
		(*curGapScore) = diagScore;
		(*curType) = GAP_OPEN;
	}
}
