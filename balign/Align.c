#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "../blib/BLib.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/AlignEntry.h"
#include "AlignNTSpace.h"
#include "AlignColorSpace.h"
#include "Align.h"

/* TODO */
/* Assumes the read and reference are both forward strand */
/* If we are in color space, the read should also be in 
 * color space */
int Align(char *read,
		int readLength,
		char *reference,
		int referenceLength,
		ScoringMatrix *sm,
		AlignEntry *a,
		char strand,
		int space)
{
	char *FnName="Align";
	int returnValue=-1;
	char reverseRead[SEQUENCE_LENGTH]="\0";
	char tmpString[SEQUENCE_LENGTH]="\0";
	char *reverseReference=NULL;

	switch(space) {
		case NTSpace:
			/* NT Space */
			switch(strand) {
				case FORWARD:
					/* Matches the forward strand */
					/* Align */
					returnValue = AlignNTSpace(read,
							readLength,
							reference,
							referenceLength,
							sm,
							a);
					break;
				case REVERSE:
					/* Reverse the read to match the forward strand  */
					GetReverseComplimentAnyCase(read, reverseRead, readLength);
					/* Align */
					returnValue = AlignNTSpace(reverseRead,
							readLength,
							reference,
							referenceLength,
							sm,
							a);
					/* We must reverse the alignment to match the REVERSE stand */
					GetReverseComplimentAnyCase(a->read, tmpString, a->length);
					strcpy(a->read, tmpString);
					GetReverseComplimentAnyCase(a->reference, tmpString, a->length);
					strcpy(a->reference, tmpString);
					break;
				default:
					PrintError(FnName,
							NULL,
							"Could not understand strand",
							Exit,
							OutOfRange);
					break;

			}
			break;
		case ColorSpace:
			/* Color Space */
			switch(strand) {
				case FORWARD:
					/* Matches the forward strand */
					/* Align */
					returnValue = AlignColorSpace(read,
							readLength,
							reference,
							referenceLength,
							sm,
							a,
							FORWARD);
					break;
				case REVERSE:
					/* Matches the reverse strand */
					/* Reverse compliment the reference */
					reverseReference = malloc(sizeof(char)*(referenceLength+1));
					if(NULL==reverseReference) {
						PrintError(FnName,
								"reverseReference",
								"Could not allocate memory",
								Exit,
								MallocMemory);
					}
					GetReverseComplimentAnyCase(reference,
							reverseReference,
							referenceLength);
					/* Align */
					returnValue = AlignColorSpace(read,
							readLength,
							reverseReference,
							referenceLength,
							sm,
							a,
							REVERSE);
					/* No need to reverse alignment, since we reversed the reference
					 * to be the reverse strand */

					free(reverseReference);
					reverseReference=NULL;
					break;
				default:
					PrintError(FnName,
							NULL,
							"Could not understand strand",
							Exit,
							OutOfRange);
					break;

			}
			break;
		default:
			PrintError(FnName,
					"space",
					"Could not understand space",
					Exit,
					OutOfRange);
	}
	return returnValue;
}

int FillAlignEntryFromMatrix(AlignEntry *a,
		AlignMatrix **matrix,
		char *read,
		int readLength,
		char *reference,
		int referenceLength,
		int space,
		int debug)
{
	char *FnName="FillAlignEntryFromMatrix";
	int curRow, curCol, curCell, startRow, startCol, startCell; 
	char curReadBase;
	int nextRow, nextCol, nextCell;
	char nextReadBase;
	int curFrom;
	double maxScore;
	double maxScoreNT=NEGATIVE_INFINITY;
	int i, j;
	int offset;

	curReadBase = nextReadBase = 'X';
	nextRow = nextCol = nextCell = -1;

	/* First allocate the maximum length of the alignment, we can update later if necessay */
	assert(NULL==a->read);
	a->read = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL==a->read) {
		PrintError(FnName,
				"a->read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(NULL==a->reference);
	a->reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL==a->reference) {
		PrintError(FnName,
				"a->reference",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(NULL==a->colorError);
	a->colorError = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL==a->colorError) {
		PrintError(FnName,
				"a->colorError",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Get the best alignment.  We can find the best score in the last row and then
	 * trace back.  We choose the best score from the last row since we want to 
	 * align the read completely and only locally to the reference. */
	startRow=-1;
	startCol=-1;
	startCell=-1;
	maxScore = NEGATIVE_INFINITY;
	maxScoreNT = NEGATIVE_INFINITY;
	for(i=0;i<referenceLength+1;i++) {
		if(space == ColorSpace) {
			for(j=0;j<ALIGNMATRIXCELL_NUM_SUB_CELLS;j++) {
				/* Cannot end with a deletion from the read */
				if(j != 4 && matrix[readLength][i].score[j] > maxScore) {
					maxScore = matrix[readLength][i].score[j];
					maxScoreNT = matrix[readLength][i].scoreNT[j];
					startRow = readLength;
					startCol = i;
					startCell = j;
				}
			}
		}
		else {
			assert(space == NTSpace);
			/* Check only the first cell */
			if(matrix[readLength][i].score[0] > maxScore) {
				maxScore = matrix[readLength][i].score[0];
				maxScoreNT = matrix[readLength][i].scoreNT[0];
				startRow = readLength;
				startCol = i;
				startCell = 0;
			}
		}
	}
	assert(startRow >= 0 && startCol >= 0 && startCell >= 0);
	/* Cannot end with a deletion from the read */
	assert(startCell != 4);

	/* Initialize variables for the loop */
	curRow=startRow;
	curCol=startCol;
	curCell=startCell;

	if(space == ColorSpace) {
	/* Color space */
		/*Initialize the current read base */
		switch(curCell) {
			case 0:
				curReadBase = 'A';
				break;
			case 1:
				curReadBase = 'C';
				break;
			case 2:
				curReadBase = 'G';
				break;
			case 3:
				curReadBase = 'T';
				break;
			case 5:
				curReadBase = matrix[curRow][curCol].prevInsertionBase;
				break;
			default:
				PrintError(FnName,
						"curCell",
						"Could not initialize curCell",
						Exit,
						OutOfRange);
		}
		assert(curReadBase != 'X');
		a->referenceLength=0;
		i=matrix[curRow][curCol].length[curCell]-1; /* Get the length of the alignment */
		a->length=matrix[curRow][curCol].length[curCell]; /* Copy over the length */
		a->score = maxScoreNT; /* Copy over score */
		/* Now trace back the alignment using the "from" member in the matrix */
		while(curRow > 0 && curCol > 0) {
			/* Where did the current cell come from */
			curFrom = matrix[curRow][curCol].from[curCell];

			assert(i>=0);

			/* Get if there was a color error */
			a->colorError[i] = matrix[curRow][curCol].colorError[curCell];

			/* Update alignment and next row/col */
			switch(curFrom) {
				case DiagA:
				case DiagC:
				case DiagG:
				case DiagT:
				case DeletionEnd:
				case InsertionEnd:
					a->read[i] = curReadBase;
					a->reference[i] = reference[curCol-1];
					a->referenceLength++;
					nextRow = curRow-1;
					nextCol = curCol-1;
					break;
				case DeletionA:
				case DeletionC:
				case DeletionG:
				case DeletionT:
				case DeletionExt:
					a->read[i] = GAP;
					a->reference[i] = reference[curCol-1];
					a->referenceLength++;
					nextRow = curRow;
					nextCol = curCol-1;
					break;
				case InsertionA:
				case InsertionC:
				case InsertionG:
				case InsertionT:
				case InsertionExt:
					assert(curReadBase != GAP);
					a->read[i] = curReadBase;
					a->reference[i] = GAP;
					nextRow = curRow-1;
					nextCol = curCol;
					break;
				default:
					fprintf(stderr, "curFrom=%d\n", curFrom);
					fprintf(stderr, "InsertionExt=%d\n", InsertionExt);
					PrintError(FnName,
							"curFrom",
							"Could not understand curFrom",
							Exit,
							OutOfRange);
			}

			/* Update previous base (relevant for color errors) and the
			 * next cell */
			switch(curFrom) {
				case DiagA:
				case DeletionA:
				case InsertionA:
					nextReadBase = 'A';
					nextCell = 0;
					break;
				case DiagC:
				case DeletionC:
				case InsertionC:
					nextReadBase = 'C';
					nextCell = 1;
					break;
				case DiagG:
				case DeletionG:
				case InsertionG:
					nextReadBase = 'G';
					nextCell = 2;
					break;
				case DiagT:
				case DeletionT:
				case InsertionT:
					nextReadBase = 'T';
					nextCell = 3;
					break;
				case DeletionExt:
				case DeletionEnd:
					nextReadBase = GAP;
					nextCell = 4;
					break;
				case InsertionExt:
				case InsertionEnd:
					nextReadBase = matrix[curRow][curCol].prevInsertionBase;
					assert(nextReadBase != GAP);
					nextCell = 5;
					break;
				default:
					fprintf(stderr, "curFrom=%d\n", curFrom);
					PrintError(FnName,
							"curFrom",
							"Could not understand curFrom",
							Exit,
							OutOfRange);
			}

			/* Update next row and column */

			assert(a->read[i] != GAP || a->read[i] != a->reference[i]); 

			/* Update for next loop iteration */
			curReadBase = nextReadBase;
			curRow = nextRow;
			curCol = nextCol;
			curCell = nextCell;
			i--;
		} /* End loop */
	}
	else { /* NT space */
		assert(space == NTSpace);
		a->referenceLength=0;
		i=matrix[curRow][curCol].length[curCell]-1; /* Get the length of the alignment */
		a->length=matrix[curRow][curCol].length[curCell]; /* Copy over the length */
		a->score = maxScore; /* Copy over score */
		/* Now trace back the alignment using the "from" member in the matrix */
		while(curRow > 0 && curCol > 0) {

			/* Where did the current cell come from */
			curFrom = matrix[curRow][curCol].from[curCell];

			assert(i>=0);

			/* Get if there was a color error (should not be possible
			 * in nt space) */
			a->colorError[i] = matrix[curRow][curCol].colorError[curCell];
			assert(a->colorError[i] == '0');

			/* Update alignment */
			switch(curFrom) {
				case DiagA:
					a->read[i] = read[curRow-1];
					a->reference[i] = reference[curCol-1];
					a->referenceLength++;
					break;
				case DeletionA:
				case DeletionExt:
					a->read[i] = GAP;
					a->reference[i] = reference[curCol-1];
					a->referenceLength++;
					break;
				case InsertionA:
				case InsertionExt:
					a->read[i] = read[curRow-1];
					a->reference[i] = GAP;
					break;
				default:
					fprintf(stderr, "curFrom=%d\n", curFrom);
					fprintf(stderr, "InsertionExt=%d\n", InsertionExt);
					PrintError(FnName,
							"curFrom",
							"Could not understand curFrom",
							Exit,
							OutOfRange);
			}

			/* Next row and column */
			switch(curFrom) {
				case DiagA:
					nextRow = curRow-1;
					nextCol = curCol-1;
					break;
				case DeletionA:
				case DeletionExt:
					nextRow = curRow;
					nextCol = curCol-1;
					break;
				case InsertionA:
				case InsertionExt:
					nextRow = curRow-1;
					nextCol = curCol;
					break;
				default:
					PrintError(FnName,
							"curFrom",
							"Could not understand curFrom (updating row/col)",
							Exit,
							OutOfRange);
			}

			/* Next cell */
			switch(curFrom) {
				case DiagA:
				case DeletionA:
				case InsertionA:
					nextCell = 0;
					break;
				case DeletionExt:
					nextCell = 1;
					break;
				case InsertionExt:
					nextCell = 2;
					break;
				default:
					PrintError(FnName,
							"curFrom",
							"Could not understand curFrom (updating cell)",
							Exit,
							OutOfRange);
			}

			assert(a->read[i] != GAP || a->read[i] != a->reference[i]); 

			/* Update for next loop iteration */
			curRow = nextRow;
			curCol = nextCol;
			curCell = nextCell;
			i--;

			assert(i!=-1 || a->read[0] != GAP);
		} /* End Loop */
	}
	assert(-1==i);
	assert(a->length >= a->referenceLength);

	offset = curCol;
	a->read[a->length]='\0';
	a->reference[a->length]='\0';
	a->colorError[a->length]='\0';

	return offset;
}
