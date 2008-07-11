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
	int prevType=-1;
	int offset=0;
	MatrixEntry **Entries; /* store the dynamic programming array */
	int *rows=NULL;
	int *cols=NULL;
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
		Entries[i][0].vType = NO_TYPE;
		Entries[i][0].hType = NO_TYPE;
	}
	/* Scoring Matrix: initialize columns for the first row */
	for(j=0;j<numCols;j++) {
		Entries[0][j].hScore = NEGATIVE_INFINITY;
		Entries[0][j].vScore = NEGATIVE_INFINITY;
		Entries[0][j].dScore = 0;
		Entries[0][j].prevRow = -1;
		Entries[0][j].prevCol = -1;
		Entries[0][j].dType = NO_TYPE;
		Entries[0][j].vType = NO_TYPE;
		Entries[0][j].hType = NO_TYPE;
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
	assert(NULL==aEntry->read);
	aEntry->read = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL==aEntry->read) {
		PrintError(FnName,
				"aEntry->read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(NULL==aEntry->reference);
	aEntry->reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL==aEntry->reference) {
		PrintError(FnName,
				"aEntry->reference",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Allocate memory for the path */
	pathLength = 0;
	cols = malloc(sizeof(int)*(endRow + endCol));
	if(NULL == cols) {
		PrintError(FnName,
				"cols",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	rows = malloc(sizeof(int)*(endRow + endCol));
	if(NULL == rows) {
		PrintError(FnName,
				"rows",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	if(AlignmentGetScore_DEBUG_ON == 1) {
		fprintf(stderr, "(read,ref)=[%s,%s]\n",
				read,
				reference);
	}

	prevType = DIAGONAL;
	curRow=endRow;
	curCol=endCol;
	i=0;
	while(curRow != -1 && curCol != -1 && prevType != NO_TYPE) {
		if(AlignmentGetScore_DEBUG_ON == 1) {
			fprintf(stderr, "i:%d\t(curRow,curCol)=[%d,%d]\tprevType:%d\t", 
					i,
					curRow,
					curCol,
					prevType);
			i++;
		}

		assert(pathLength < endRow + endCol);
		rows[pathLength] = curRow;
		cols[pathLength] = curCol;
		pathLength++;

		/* We can be in a gap (horizontal or vertical) or going along the diagonal */
		switch(prevType) {
			case DIAGONAL:
				/* We were going along a diagonal, check next */
				switch(Entries[curRow][curCol].dType) {
					case DIAGONAL:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "From Diagonal to Diagonal");
						}
						curRow--;
						curCol--;
						break;
					case HORIZONTAL:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "From Diagonal to Horizontal");
						}
						/* Get prevType */
						switch(Entries[curRow][curCol].hType) {
							case GAP_OPEN:
								prevType = DIAGONAL;
								break;
							case GAP_EXTENSION:
								prevType = HORIZONTAL;
								break;
							default:
								PrintError(FnName,
										NULL,
										"Could not recover alignment path: hType",
										Exit,
										OutOfRange);
								break;
						}
						curCol--;
						break;
					case VERTICAL:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "From Diagonal to Vertical");
						}
						switch(Entries[curRow][curCol].vType) {
							case GAP_OPEN:
								prevType = DIAGONAL;
								break;
							case GAP_EXTENSION:
								prevType = VERTICAL;
								break;
							default:
								PrintError(FnName,
										NULL,
										"Could not recover alignment path: vType",
										Exit,
										OutOfRange);
								break;
						}
						curRow--;
						break;
					case NO_TYPE:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "From Diagonal to No Type");
						}
						prevType = NO_TYPE;
						break;
					default:
						PrintError(FnName,
								NULL,
								"Could not recover alignment path: dType",
								Exit,
								OutOfRange);
						break;
				}
				break;
			case HORIZONTAL:
				/* Update previous type.  If we have a gap extension then we 
				 * previous type is Horizontal, if we have a gap open, then we 
				 * are coming from a diagonal */
				switch(Entries[curRow][curCol].hType) {
					case GAP_OPEN:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "From Horizontal to Diagonal");
						}
						prevType = DIAGONAL;
						break;
					case GAP_EXTENSION:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "From Horizontal to Horizontal");
						}
						prevType = HORIZONTAL;
						break;
					case NO_TYPE:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "From Horizontal to No Type");
						}
						prevType = NO_TYPE;
						break;
					default:
						PrintError(FnName,
								NULL,
								"Could not recover alignment path: hType",
								Exit,
								OutOfRange);
						break;
				}
				curCol--;
				break;
			case VERTICAL:
				/* Update row and column */
				/* Update previous type.  If we have a gap extension then we 
				 * previous type is Vetical, if we have a gap open, then we 
				 * are coming from a diagonal */
				switch(Entries[curRow][curCol].vType) {
					case GAP_OPEN:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "From Vertical to Diagonal");
						}
						prevType = DIAGONAL;
						break;
					case GAP_EXTENSION:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "From Vertical to Vertical");
						}
						prevType = VERTICAL;
						break;
					case NO_TYPE:
						if(AlignmentGetScore_DEBUG_ON == 1) {
							fprintf(stderr, "From Vertical to No Type");
						}
						prevType = NO_TYPE;
						break;
					default:
						fprintf(stderr, "[%d]\n",
								Entries[curRow][curCol].vType);
						PrintError(FnName,
								NULL,
								"Could not recover alignment path: vType",
								Exit,
								OutOfRange);
						break;
				}
				curRow--;
				break;
			default:
				PrintError(FnName,
						NULL,
						"Could not recover alignment path: prevType",
						Exit,
						OutOfRange);
		}
		if(AlignmentGetScore_DEBUG_ON == 1) {
			fprintf(stderr, "\n");
		}
	}

	assert(pathLength-1 < endRow + endCol);
	startRow = rows[pathLength-1];
	startCol = cols[pathLength-1];

	if(AlignmentGetScore_DEBUG_ON == 1) {
		fprintf(stderr, "(read,ref)=[%s,%s]\n",
				read,
				reference);
		fprintf(stderr, "(startRow,startCol)=[%d,%d]\n",
				startRow,
				startCol);
		fprintf(stderr, "(endRow,endCol)=[%d,%d]\n",
				endRow,
				endCol);
		fprintf(stderr, "(curRow,curCol)=[%d,%d]\tpathLength:%d\n",
				curRow,
				curCol,
				pathLength);
	}
	if(!(endRow - startRow == readLength)) {
		fprintf(stderr, "read:%s\nreference:%s\n",
				read,
				reference);
	}
	assert(endRow - startRow == readLength);

	/* Update offset */
	offset = startCol;

	aEntry->length=0;
	for(i=pathLength-2;i>=0;i--) {
		if(AlignmentGetScore_DEBUG_ON == 1) {
			fprintf(stderr, "(row,col)=[%d,%d]\t",
					rows[i],
					cols[i]);
		}
		if(rows[i] == (rows[i+1] + 1) &&
				cols[i] == (cols[i+1] + 1)) {
			if(AlignmentGetScore_DEBUG_ON == 1) {
				fprintf(stderr, "Diagonal");
			}
			aEntry->read[aEntry->length] = read[rows[i] - 1];;
			aEntry->reference[aEntry->length] = reference[cols[i] - 1];
		}
		else if(rows[i] == (rows[i+1] + 1) &&
				cols[i] == cols[i+1]) {
			if(AlignmentGetScore_DEBUG_ON == 1) {
				fprintf(stderr, "Vertical");
			}
				aEntry->read[aEntry->length] = read[rows[i] - 1];;
				aEntry->reference[aEntry->length] = GAP;
		}
		else if(rows[i] == rows[i+1] &&
				cols[i] == (cols[i+1] + 1)) {
			if(AlignmentGetScore_DEBUG_ON == 1) {
				fprintf(stderr, "Horizontal");
			}
				aEntry->read[aEntry->length] = GAP;
				aEntry->reference[aEntry->length] = reference[cols[i] - 1];
		}
		else {
			PrintError(FnName,
					NULL,
					"Could not understand rows and columns",
					Exit,
					OutOfRange);
		}
		aEntry->length++;
	}

	aEntry->reference[aEntry->length]='\0'; /* null terminator */
	aEntry->read[aEntry->length]='\0';
	assert(offset >= 0 && offset <= referenceLength);

	if(AlignmentGetScore_DEBUG_ON == 1) {
		fprintf(stderr, "length:%d\n%s\n%s\n",
				aEntry->length,
				aEntry->read,
				aEntry->reference);
		exit(1);
	}

	/* Free memory */
	free(rows);
	rows=NULL;
	free(cols);
	cols=NULL;
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
			PrintError("GetScoreFromMatrix",
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
