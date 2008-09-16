#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "ReadInputFiles.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/AlignEntry.h"
#include "AlignNTSpace.h"
#include "AlignColorSpace.h"
#include "Align.h"

int Align(char *read,
		int readLength,
		char *reference,
		int referenceLength,
		ScoringMatrix *sm,
		AlignEntry *aEntry,
		int colorSpace)
{
	switch(colorSpace) {
		case 0:
			return AlignNTSpace(read,
					readLength,
					reference,
					referenceLength,
					sm,
					aEntry);
			break;
		default:
			return AlignColorSpace(read,
					readLength,
					reference,
					referenceLength,
					sm,
					aEntry);
			break;
	}
}

/* TODO */
double GetNTScore(char a,
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
			PrintError("GetNTScore",
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
			PrintError("GetNTScore",
					NULL,
					"Could not understand 'b' key",
					Exit,
					OutOfRange);
			break;
	}

	return sm->scores[indexA][indexB];
}

/* TODO */
double GetColorScore(uint8_t a, 
		uint8_t b,
		ScoringMatrix *sm) 
{
	/* simple for now */
	return (a==b)?(sm->colorMatch):(sm->colorError);
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
	/* HERE */
	if(ALIGN_DEBUG_ON == 1) {
		/*
		for(i=0;i<readLength+1;i++) {
			for(j=0;j<referenceLength+1;j++) {
				int k;
				for(k=0;k<ALIGNMATRIXCELL_NUM_SUB_CELLS;k++) {
					fprintf(stderr, "(row,col,cell,score,length)=(%d,%d,%d,%lf,%d)\n",
							i,
							j,
							k,
							matrix[i][j].score[k],
							matrix[i][j].length[k]);
				}
			}
		}
		*/
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
	curReadBase = read[curRow-1];
	i=matrix[curRow][curCol].length[curCell]-1;
	aEntry->length=matrix[curRow][curCol].length[curCell];
	aEntry->read[i+1]='\0';
	aEntry->reference[i+1]='\0';
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
			fprintf(stderr, "[%c][%c]\n",
					aEntry->read[i],
					aEntry->reference[i]);
		}

		/* Update previous base (relevant for color errors */
		if(colorSpace==1) {
			switch(curFrom) {
				case DiagA:
				case DeletionA:
				case InsertionA:
					curReadBase = 'A';
					break;
				case DiagC:
				case DeletionC:
				case InsertionC:
					curReadBase = 'C';
					break;
				case DiagG:
				case DeletionG:
				case InsertionG:
					curReadBase = 'G';
					break;
				case DiagT:
				case DeletionT:
				case InsertionT:
					curReadBase = 'T';
					break;
				case DeletionExt:
				case InsertionExt:
				case DeletionEnd:
				case InsertionEnd:
					curReadBase = GAP;
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
	/* This might be off by one */
	offset = curCol;

	return offset;
}
