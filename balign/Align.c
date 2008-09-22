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
		AlignEntry *aEntry,
		char strand,
		int colorSpace)
{
	char *FnName="Align";
	int returnValue=-1;
	char reverseRead[SEQUENCE_LENGTH]="\0";
	char tmpString[SEQUENCE_LENGTH]="\0";
	char *reverseReference=NULL;

	/* HERE 42 */
	/*
	   fprintf(stderr, "HERE 42\nread=%s\nreference=%s\n",
	   read,
	   reference);
	   */

	switch(colorSpace) {
		case 0:
			/* Not color space */
			switch(strand) {
				case FORWARD:
					/* Matches the forward strand */
					/* Align */
					returnValue = AlignNTSpace(read,
							readLength,
							reference,
							referenceLength,
							sm,
							aEntry);
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
							aEntry);
					/* We must reverse the alignment to match the REVERSE stand */
					GetReverseComplimentAnyCase(aEntry->read, tmpString, aEntry->length);
					strcpy(aEntry->read, tmpString);
					GetReverseComplimentAnyCase(aEntry->reference, tmpString, aEntry->length);
					strcpy(aEntry->reference, tmpString);
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
			/* HERE 49 */
			/*
			   fprintf(stderr, "HERE 49\nreference=%s\nread=%s\n",
			   reference,
			   read);
			   */
			/* Color Space */
			switch(strand) {
				case FORWARD:
					/* Matches the forward strand */
					/* Align */
					/* HERE 51 */
					/*
					   fprintf(stderr, "HERE 51\nreference=%s\nread=%s\nreverseRead=%s\n",
					   reference,
					   read,
					   reverseRead);
					   */
					returnValue = AlignColorSpace(read,
							readLength,
							reference,
							referenceLength,
							sm,
							aEntry,
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
					/* HERE 50 */
					/*
					   fprintf(stderr, "HERE 50\nreference=%s\nread=%s\nreverseRead=%s\n",
					   reference,
					   read,
					   reverseRead);
					   */
					/* Align */
					returnValue = AlignColorSpace(read,
							readLength,
							reverseReference,
							referenceLength,
							sm,
							aEntry,
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
			/* HERE 39 */
			/*
			   fprintf(stderr, "HERE 39\n");
			   */
			break;
	}
	return returnValue;
}

int FillAlignEntryFromMatrix(AlignEntry *aEntry,
		AlignMatrix **matrix,
		char *read,
		int readLength,
		char *reference,
		int referenceLength,
		int colorSpace,
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

	/* HERE */
	if(debug == 1) {
		fprintf(stderr, "\n%s\n", FnName);
		fprintf(stderr, "read=%s\n", read);
		fprintf(stderr, "reference=%s\n", reference);
	}

	/* First allocate the maximum length of the alignment, we can update later if necessay */
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
	assert(NULL==aEntry->colorError);
	aEntry->colorError = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL==aEntry->colorError) {
		PrintError(FnName,
				"aEntry->colorError",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* HERE */
	if(debug == 1) {
		for(i=0;i<readLength+1;i++) {
			for(j=0;j<referenceLength+1;j++) {
				int k;
				for(k=0;k<ALIGNMATRIXCELL_NUM_SUB_CELLS;k++) {
					fprintf(stderr, "(row,col,cell,score,length,from,colorError)=(%d,%d,%d,%lf,%d,%d,%c)\n",
							i,
							j,
							k,
							matrix[i][j].score[k],
							matrix[i][j].length[k],
							matrix[i][j].from[k],
							matrix[i][j].colorError[k]
						   );
				}
			}
		}
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
	assert(startRow >= 0 && startCol >= 0 && startCell >= 0);
	/* Cannot end with a deletion from the read */
	assert(startCell != 4);

	/* HERE */
	if(debug == 1) {
		fprintf(stderr, "(startRow,startCol,startCell,score,length)=(%d,%d,%d,%lf,%d)\n",
				startRow,
				startCol,
				startCell,
				matrix[startRow][startCol].score[startCell],
				matrix[startRow][startCol].length[startCell]);
	}

	/* Initialize variables for the loop */
	curRow=startRow;
	curCol=startCol;
	curCell=startCell;

	/* Color space */
	if(colorSpace == 1) {
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
		i=matrix[curRow][curCol].length[curCell]-1; /* Get the length of the alignment */
		aEntry->length=matrix[curRow][curCol].length[curCell]; /* Copy over the length */
		aEntry->score = maxScoreNT; /* Copy over score */
		/* HERE */
		if(debug == 1) {
			fprintf(stderr, "[%d,%lf]\n",
					matrix[curRow][curCol].length[curCell],
					maxScore);
		}
		/* Now trace back the alignment using the "from" member in the matrix */
		while(curRow > 0 && curCol > 0) {
			/* Where did the current cell come from */
			curFrom = matrix[curRow][curCol].from[curCell];

			/* HERE */
			if(debug == 1) {
				fprintf(stderr, "(curRow,curCol,i,curFrom,curCell)=(%d,%d,%d,%d,%d)\n",
						curRow,
						curCol,
						i,
						curFrom,
						curCell);
				fprintf(stderr, "cur.length=%d\n",
						matrix[curRow][curCol].length[curCell]);
			}
			assert(i>=0);

			/* Get if there was a color error */
			aEntry->colorError[i] = matrix[curRow][curCol].colorError[curCell];

			/* Update alignment and next row/col */
			switch(curFrom) {
				case DiagA:
				case DiagC:
				case DiagG:
				case DiagT:
				case DeletionEnd:
				case InsertionEnd:
					aEntry->read[i] = curReadBase;
					aEntry->reference[i] = reference[curCol-1];
					nextRow = curRow-1;
					nextCol = curCol-1;
					break;
				case DeletionA:
				case DeletionC:
				case DeletionG:
				case DeletionT:
				case DeletionExt:
					aEntry->read[i] = GAP;
					aEntry->reference[i] = reference[curCol-1];
					nextRow = curRow;
					nextCol = curCol-1;
					break;
				case InsertionA:
				case InsertionC:
				case InsertionG:
				case InsertionT:
				case InsertionExt:
					assert(curReadBase != GAP);
					aEntry->read[i] = curReadBase;
					aEntry->reference[i] = GAP;
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
			/* HERE */
			if(debug == 1) {
				fprintf(stderr, "[%c][%c](%c,%c)\n",
						aEntry->read[i],
						aEntry->reference[i],
						read[curRow-1],
						reference[curCol-1]);
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

			assert(aEntry->read[i] != GAP || aEntry->read[i] != aEntry->reference[i]); 

			/* Update for next loop iteration */
			curReadBase = nextReadBase;
			curRow = nextRow;
			curCol = nextCol;
			curCell = nextCell;
			i--;

			/* HERE */
			if(debug == 1) {
				fprintf(stderr, "next(row,col,curCell,curReadBase)=(%d,%d,%d,%c)\n",
						nextRow,
						nextCol,
						nextCell,
						curReadBase);
				fprintf(stderr, "%s", BREAK_LINE);
			}
		} /* End loop */
	}
	else { /* NT space */

		/*Initialize the current read base */
		switch(curCell) {
			case 0:
			case 2:
				/* Diagonal or insertion */
				curReadBase = read[curRow-1];
				break;
			default:
				/* Deletion not legal to end with */
				PrintError(FnName,
						"curCell",
						"Could not initialize curCell",
						Exit,
						OutOfRange);
		}
		assert(curReadBase != 'X');
		i=matrix[curRow][curCol].length[curCell]-1; /* Get the length of the alignment */
		aEntry->length=matrix[curRow][curCol].length[curCell]; /* Copy over the length */
		aEntry->score = maxScore; /* Copy over score */
		/* HERE */
		if(debug == 1) {
			fprintf(stderr, "reference=%s\nread=%s\n",
					reference,
					read);
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "[%d,%lf]\n",
					matrix[curRow][curCol].length[curCell],
					maxScore);
		}
		/* Now trace back the alignment using the "from" member in the matrix */
		while(curRow > 0 && curCol > 0) {

			/* Where did the current cell come from */
			curFrom = matrix[curRow][curCol].from[curCell];

			/* HERE */
			if(debug == 1) {
				fprintf(stderr, "(curRow,curCol,i,curFrom,curCell,curReadBase)=(%d,%d,%d,%d,%d,%c)\n",
						curRow,
						curCol,
						i,
						curFrom,
						curCell,
						curReadBase);
				fprintf(stderr, "cur.length=%d\n",
						matrix[curRow][curCol].length[curCell]);
			}
			assert(i>=0);

			/* Get if there was a color error (should not be possible
			 * in nt space) */
			aEntry->colorError[i] = matrix[curRow][curCol].colorError[curCell];
			assert(aEntry->colorError[i] == '0');

			/* Update alignment */
			switch(curFrom) {
				case DiagA:
					aEntry->read[i] = curReadBase;
					aEntry->reference[i] = reference[curCol-1];
					break;
				case DeletionA:
				case DeletionExt:
					aEntry->read[i] = GAP;
					aEntry->reference[i] = reference[curCol-1];
					break;
				case InsertionA:
				case InsertionExt:
					assert(curReadBase != GAP);
					aEntry->read[i] = curReadBase;
					aEntry->reference[i] = GAP;
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
			/* HERE */
			if(debug == 1) {
				fprintf(stderr, "aEntry:[%c][%c]\toriginal:(%c,%c)\n",
						aEntry->read[i],
						aEntry->reference[i],
						read[curRow-1],
						reference[curCol-1]);
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

			/* Update next base */
			switch(curFrom) {
				case DiagA:
				case InsertionA:
				case InsertionExt:
					if(curRow >= 2) {
						nextReadBase = read[curRow-2];
					}
					else {
						nextReadBase = 'X';
					}
					break;
				case DeletionA:
				case DeletionExt:
					nextReadBase = GAP;
					break;
				default:
					PrintError(FnName,
							"curFrom",
							"Could not understand curFrom (updating next read base)",
							Exit,
							OutOfRange);
			}

			assert(aEntry->read[i] != GAP || aEntry->read[i] != aEntry->reference[i]); 

			/* Update for next loop iteration */
			curReadBase = nextReadBase;
			curRow = nextRow;
			curCol = nextCol;
			curCell = nextCell;
			i--;

			/* HERE */
			if(debug == 1) {
				fprintf(stderr, "next(row,col,curCell,curReadBase)=(%d,%d,%d,%c)\n",
						nextRow,
						nextCol,
						nextCell,
						curReadBase);
				fprintf(stderr, "%s", BREAK_LINE);
			}
		} /* End Loop */
	}
	/* HERE */
	if(debug == 1) {
		int tempi=i;
		for(i=0;i<readLength+1;i++) {
			for(j=0;j<referenceLength+1;j++) {
				/*
				   fprintf(stderr, "(i,j,length)=(%d,%d,%d)\n",
				   i,
				   j,
				   matrix[i][j].length[0]);
				   */
				fprintf(stderr, "%d ",
						matrix[i][j].length[0]);
			}
			fprintf(stderr, "\n");
		}
		i=tempi;
	}

	/* HERE */
	if(debug == 1) {
		fprintf(stderr, "i=%d\n", i);
	}
	assert(-1==i);

	/* HERE C1 */
	/*
	   fprintf(stderr, "%s\n%s\n%s\n",
	   aEntry->reference,
	   aEntry->read,
	   aEntry->colorError);
	   fprintf(stderr, "HERE C1\n");
	   exit(1);
	   */
	/* This might be off by one */
	offset = curCol;
	aEntry->read[aEntry->length]='\0';
	aEntry->reference[aEntry->length]='\0';
	aEntry->colorError[aEntry->length]='\0';

	return offset;
}
