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
	int returnValue;
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
		int colorSpace)
{
	char *FnName="FillAlignEntryFromMatrix";
	int curRow, curCol, curCell, startRow, startCol, startCell; 
	int curFrom;
	double maxScore;
	int i, j;
	char curReadBase;
	int offset;

	/* HERE */
	if(ALIGN_DEBUG_ON == 1) {
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
	if(ALIGN_DEBUG_ON == 1) {
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
	for(i=0;i<referenceLength+1;i++) {
		for(j=0;j<ALIGNMATRIXCELL_NUM_SUB_CELLS;j++) {
			if(matrix[readLength][i].score[j] > maxScore) {
				maxScore = matrix[readLength][i].score[j];
				startRow = readLength;
				startCol = i;
				startCell = j;
			}
		}
	}
	assert(startRow >= 0 && startCol >= 0 && startCell >= 0);

	/* HERE */
	if(ALIGN_DEBUG_ON == 1) {
		fprintf(stderr, "(startRow,startCol,startCell,score,length)=(%d,%d,%d,%lf,%d)\n",
				startRow,
				startCol,
				startCell,
				matrix[startRow][startCol].score[startCell],
				matrix[startRow][startCol].length[startCell]);
	}

	curRow=startRow;
	curCol=startCol;
	curCell=startCell;
	curReadBase = 'X';
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
		case 4:
		case 5:
			curReadBase = GAP;
			break;
		default:
			PrintError(FnName,
					"curCell",
					"Could not initialize curCell",
					Exit,
					OutOfRange);
	}
	i=matrix[curRow][curCol].length[curCell]-1;
	aEntry->length=matrix[curRow][curCol].length[curCell];
	aEntry->read[i+1]='\0';
	aEntry->reference[i+1]='\0';
	aEntry->colorError[i+1]='\0';
	aEntry->score = maxScore;
	/* HERE */
	if(ALIGN_DEBUG_ON == 1) {
		fprintf(stderr, "[%d,%lf]\n",
				matrix[curRow][curCol].length[curCell],
				maxScore);
	}
	/* Now trace back the alignment using the "from" member in the matrix */
	while(curRow > 0 && curCol > 0) {
		/* HERE */
		if(ALIGN_DEBUG_ON == 1) {
			fprintf(stderr, "(curRow,curCol,i)=(%d,%d,%d)\n",
					curRow,
					curCol,
					i);
			fprintf(stderr, "cur.length=%d,%d,%d\n",
					matrix[curRow][curCol].length[0],
					matrix[curRow][curCol].length[1],
					matrix[curRow][curCol].length[2]);
		}
		assert(i>=0);

		/* Where did the current cell come from */
		curFrom = matrix[curRow][curCol].from[curCell];
		/* Get if there was a color error */
		aEntry->colorError[i] = matrix[curRow][curCol].colorError[curCell];

		/* Update alignment */
		switch(curFrom) {
			case DiagA:
			case DiagC:
			case DiagG:
			case DiagT:
			case DeletionEnd:
			case InsertionEnd:
				aEntry->read[i] = curReadBase;
				aEntry->reference[i] = reference[curCol-1];
				break;
			case DeletionA:
			case DeletionC:
			case DeletionG:
			case DeletionT:
			case DeletionExt:
				aEntry->read[i] = GAP;
				aEntry->reference[i] = reference[curCol-1];
				break;
			case InsertionA:
			case InsertionC:
			case InsertionG:
			case InsertionT:
			case InsertionExt:
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
		if(ALIGN_DEBUG_ON == 1) {
			fprintf(stderr, "[%c][%c](%c,%c)\n",
					aEntry->read[i],
					aEntry->reference[i],
					read[curRow-1],
					reference[curCol-1]);
		}

		/* Update previous base (relevant for color errors) and the
		 * next cell */
		if(colorSpace==1) {
			switch(curFrom) {
				case DiagA:
				case DeletionA:
				case InsertionA:
					curReadBase = 'A';
					curCell = 0;
					break;
				case DiagC:
				case DeletionC:
				case InsertionC:
					curReadBase = 'C';
					curCell = 1;
					break;
				case DiagG:
				case DeletionG:
				case InsertionG:
					curReadBase = 'G';
					curCell = 2;
					break;
				case DiagT:
				case DeletionT:
				case InsertionT:
					curReadBase = 'T';
					curCell = 3;
					break;
				case DeletionExt:
				case DeletionEnd:
					curReadBase = GAP;
					curCell = 4;
					break;
				case InsertionExt:
				case InsertionEnd:
					curReadBase = GAP;
					curCell = 5;
					break;
				default:
					fprintf(stderr, "curFrom=%d\n", curFrom);
					PrintError(FnName,
							"curFrom",
							"Could not understand curFrom",
							Exit,
							OutOfRange);
			}
		}
		else {
			switch(curFrom) {
				case DiagA:
				case DeletionA:
				case InsertionA:
				case DiagC:
				case DeletionC:
				case InsertionC:
				case DiagG:
				case DeletionG:
				case InsertionG:
				case DiagT:
				case DeletionT:
				case InsertionT:
					if(curRow>1) {
						curReadBase = read[curRow-2];
					}
					else {
						curReadBase = 'X';
					}
					break;
				case DeletionExt:
				case InsertionExt:
					curReadBase = GAP;
					break;
				default:
					PrintError(FnName,
							"curFrom",
							"Could not understand curFrom",
							Exit,
							OutOfRange);
			}
		}

		/* Get next row and column */
		switch(curFrom) {
			case DiagA:
			case DiagC:
			case DiagG:
			case DiagT:
			case DeletionEnd:
			case InsertionEnd:
				curRow--;
				curCol--;
				break;
			case DeletionA:
			case DeletionC:
			case DeletionG:
			case DeletionT:
			case DeletionExt:
				curRow--;
				break;
			case InsertionA:
			case InsertionC:
			case InsertionG:
			case InsertionT:
			case InsertionExt:
				curCol--;
				break;
			case Start:
				curRow=-1;
				curCol=-1;
				break;
			default:
				PrintError(FnName,
						"curFrom",
						"Could not understand curFrom",
						Exit,
						OutOfRange);
		}

		i--;
		/* HERE */
		if(ALIGN_DEBUG_ON == 1) {
			fprintf(stderr, "next(row,col,i)=(%d,%d,%d)\n",
					curRow,
					curCol,
					i);
		}
	}
	/* HERE */
	if(ALIGN_DEBUG_ON == 1) {
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
	if(ALIGN_DEBUG_ON == 1) {
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

	return offset;
}
