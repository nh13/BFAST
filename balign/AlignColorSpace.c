#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "../blib/BLib.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/AlignEntry.h"
#include "ScoringMatrix.h"
#include "Align.h"
#include "AlignColorSpace.h"

/* TODO */
int AlignColorSpace(char *read,
		int readLength,
		char *reference,
		int referenceLength,
		ScoringMatrix *sm,
		AlignEntry *aEntry,
		char strand)
{
	/* read goes on the rows, reference on the columns */
	char *FnName = "AlignColorSpace";
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
			matrix[i][0].scoreNT[k] = NEGATIVE_INFINITY;
			matrix[i][0].from[k] = Start;
			matrix[i][0].length[k] = 0;
			matrix[i][0].colorError[k] = '0';
		}
		matrix[i][0].prevDeletionBase = COLOR_SPACE_START_NT;
		matrix[i][0].prevInsertionBase = COLOR_SPACE_START_NT;
	}
	/* Row 0 column j should be zero since we want to find the best
	 * local alignment within the reference */
	for(j=0;j<referenceLength+1;j++) {
		for(k=0;k<ALIGNMATRIXCELL_NUM_SUB_CELLS;k++) {
			/* Assumes both DNA and COLOR_SPACE_START_NT are upper case */
			if(DNA[k] == COLOR_SPACE_START_NT) { 
				/* Starting adaptor NT */
				matrix[0][j].score[k] = 0;
				matrix[0][j].scoreNT[k] = 0;
			}
			else {
				matrix[0][j].score[k] = NEGATIVE_INFINITY;
				matrix[0][j].scoreNT[k] = NEGATIVE_INFINITY;
			}
			matrix[0][j].from[k] = Start;
			matrix[0][j].length[k] = 0;
			matrix[0][j].colorError[k] = '0';
		}
		matrix[0][j].prevDeletionBase = COLOR_SPACE_START_NT;
		matrix[0][j].prevInsertionBase = COLOR_SPACE_START_NT;
	}

	/* Fill in the matrix according to the recursive rules */
	for(i=0;i<readLength;i++) { /* read/rows */
		/* Get the current color */
		uint8_t curColor;
		char curReadBase, prevReadBase;
		/* In color space, the first color is determined by the adapter NT */
		curReadBase = read[i];
		prevReadBase = (i==0)?COLOR_SPACE_START_NT:read[i-1];

		/* Get the current color for the read */
		if(0 == ConvertBaseToColorSpace(prevReadBase, curReadBase, &curColor)) {
			fprintf(stderr, "prevReadBase=%c\tcurReadBase=%c\n",
					prevReadBase,
					curReadBase);
			PrintError(FnName,
					"curColor",
					"Could not convert base to color space",
					Exit,
					OutOfRange);
		}

		for(j=0;j<referenceLength;j++) { /* reference/columns */

			for(k=0;k<ALIGNMATRIXCELL_NUM_SUB_CELLS;k++) { /* To NT */
				char DNA[4] = "ACGT";
				double maxScore = NEGATIVE_INFINITY-1;
				double maxScoreNT = NEGATIVE_INFINITY-1;
				int maxFrom = -1;
				char maxColorError = '0';
				int maxLength = 0;
				char maxPrevDeletionBase='X';
				char maxPrevInsertionBase='X';

				for(l=0;l<ALIGNMATRIXCELL_NUM_SUB_CELLS;l++) { /* From NT */
					double curScore=NEGATIVE_INFINITY;
					double curScoreNT=NEGATIVE_INFINITY;
					int curFrom=-1;
					int curLength=-1;
					char curPrevDeletionBase='X';
					char curPrevInsertionBase='X';
					uint8_t convertedColor='X';
					switch(k) {
						case 0:
						case 1:
						case 2:
						case 3:
							/* Check diagonals from NT */
							/* Previous score */ 
							curScore = matrix[i][j].score[l];
							curScoreNT = matrix[i][j].scoreNT[l];
							curLength = matrix[i][j].length[l] + 1;
							/* Plus score for colors */
							switch(l) {
								case 0:
								case 1:
								case 2:
								case 3:
									/* Diagonals */
									if(0 == ConvertBaseToColorSpace(DNA[l], DNA[k], &convertedColor)) {
										fprintf(stderr, "DNA[l=%d]=%c\tDNA[k=%d]=%c\n",
												l,
												DNA[l],
												k,
												DNA[k]);
										PrintError(FnName,
												"convertedColor",
												"Could not convert base to color space",
												Exit,
												OutOfRange);
									}
									break;
								case 4:
									/* Use previous base for deletions */
									if(0 == ConvertBaseToColorSpace(matrix[i][j].prevDeletionBase, DNA[k], &convertedColor)) {
										fprintf(stderr, "matrix[i=%d][j=%d].prevDeletionBase=%c\tDNA[k=%d]=%c\n",
												i,
												j,
												matrix[i][j].prevDeletionBase,
												k,
												DNA[k]);
										PrintError(FnName,
												"convertedColor",
												"Could not convert base to color space",
												Exit,
												OutOfRange);
									}
									break;
								case 5:
									/* Use the previous base carried along by insertions */
									if(0 == ConvertBaseToColorSpace(matrix[i][j].prevInsertionBase, DNA[k], &convertedColor)) {
										fprintf(stderr, "matrix[i=%d][j=%d].prevInsertionBase=%c\tDNA[k=%d]=%c\n",
												i,
												j,
												matrix[i][j].prevInsertionBase,
												k,
												DNA[k]);
										PrintError(FnName,
												"convertedColor",
												"Could not convert base to color space",
												Exit,
												OutOfRange);
									}
									break;
								default:
									PrintError(FnName,
											"l",
											"Could not understand l",
											Exit,
											OutOfRange);
									break;
							}
							/* Add score for color error, if any */
							curScore += ScoringMatrixGetColorScore(curColor,
									convertedColor,
									sm);
							/* Add score for NT */
							curScore += ScoringMatrixGetNTScore(reference[j], DNA[k], sm);
							curScoreNT += ScoringMatrixGetNTScore(reference[j], DNA[k], sm);

							if(curScore < NEGATIVE_INFINITY/2) {
								curScore = NEGATIVE_INFINITY;
								curScoreNT = NEGATIVE_INFINITY;
							}

							/* Check to see if this is better than the max */
							if(curScore > maxScore) {
								maxScore = curScore;
								maxScoreNT = curScoreNT;
								maxLength = curLength;
								maxColorError = (curColor == convertedColor)?'0':'1';
								switch(l) {
									case 0:
										maxFrom = DiagA;
										break;
									case 1:
										maxFrom = DiagC;
										break;
									case 2:
										maxFrom = DiagG;
										break;
									case 3:
										maxFrom = DiagT;
										break;
									case 4:
										maxFrom = DeletionEnd;
										break;
									case 5:
										maxFrom = InsertionEnd;
										break;
									default:
										/*
										   fprintf(stderr, "l=%d\n", l);
										   */
										PrintError(FnName,
												"l",
												"Could not understand character",
												Exit,
												OutOfRange);
								}
							}
							break;
						case 4:
							/* Deletion - previous column */
							curLength = matrix[i+1][j].length[l] + 1;
							switch(l) {
								/* 0-3 gap open */
								/* 4 gap extension */
								case 0:
									curScore = matrix[i+1][j].score[l] + sm->gapOpenPenalty;
									curScoreNT = matrix[i+1][j].scoreNT[l] + sm->gapOpenPenalty;
									curFrom = DeletionA;
									curPrevDeletionBase = DNA[l]; 
									break;
								case 1:
									curScore = matrix[i+1][j].score[l] + sm->gapOpenPenalty;
									curScoreNT = matrix[i+1][j].scoreNT[l] + sm->gapOpenPenalty;
									curFrom = DeletionC;
									curPrevDeletionBase = DNA[l]; 
									break;
								case 2:
									curScore = matrix[i+1][j].score[l] + sm->gapOpenPenalty;
									curScoreNT = matrix[i+1][j].scoreNT[l] + sm->gapOpenPenalty;
									curFrom = DeletionG;
									curPrevDeletionBase = DNA[l]; 
									break;
								case 3:
									curScore = matrix[i+1][j].score[l] + sm->gapOpenPenalty;
									curScoreNT = matrix[i+1][j].scoreNT[l] + sm->gapOpenPenalty;
									curFrom = DeletionT;
									curPrevDeletionBase = DNA[l]; 
									break;
								case 4:
									curScore = matrix[i+1][j].score[l] + sm->gapExtensionPenalty;
									curScoreNT = matrix[i+1][j].scoreNT[l] + sm->gapExtensionPenalty;
									curFrom = DeletionExt;
									curPrevDeletionBase = matrix[i+1][j].prevDeletionBase;
									break;
								case 5:
									break;
								default:
									/*
									   fprintf(stderr, "l=%d\n", l);
									   */
									PrintError(FnName,
											"l",
											"Could not understand character",
											Exit,
											OutOfRange);
							}
							if(curScore < NEGATIVE_INFINITY/2) {
								curScore = NEGATIVE_INFINITY;
								curScoreNT = NEGATIVE_INFINITY;
							}
							if(curScore > maxScore && l != 5) {
								maxScore = curScore;
								maxScoreNT = curScoreNT;
								maxFrom = curFrom;
								maxColorError = '0';
								maxLength = curLength;
								maxPrevDeletionBase=curPrevDeletionBase;
							}
							break;
						case 5:
							/* Insertion - previous row */
							curLength = matrix[i][j+1].length[l] + 1;
							switch(l) {
								/* 0-3 gap open */
								/* 5 gap extension */
								case 0:
									curScore = matrix[i][j+1].score[l] + sm->gapOpenPenalty;
									curScoreNT = matrix[i][j+1].scoreNT[l] + sm->gapOpenPenalty;
									curFrom = InsertionA;
									curPrevInsertionBase = DNA[l];
									break;
								case 1:
									curScore = matrix[i][j+1].score[l] + sm->gapOpenPenalty;
									curScoreNT = matrix[i][j+1].scoreNT[l] + sm->gapOpenPenalty;
									curFrom = InsertionC;
									curPrevInsertionBase = DNA[l];
									break;
								case 2:
									curScore = matrix[i][j+1].score[l] + sm->gapOpenPenalty;
									curScoreNT = matrix[i][j+1].scoreNT[l] + sm->gapOpenPenalty;
									curFrom = InsertionG;
									curPrevInsertionBase = DNA[l];
									break;
								case 3:
									curScore = matrix[i][j+1].score[l] + sm->gapOpenPenalty;
									curScoreNT = matrix[i][j+1].scoreNT[l] + sm->gapOpenPenalty;
									curFrom = InsertionT;
									curPrevInsertionBase = DNA[l];
									break;
								case 4:
									break;
								case 5:
									curScore = matrix[i][j+1].score[l] + sm->gapExtensionPenalty;
									curScoreNT = matrix[i][j+1].scoreNT[l] + sm->gapExtensionPenalty;
									curFrom = InsertionExt;
									if(0 == ConvertBaseAndColor(matrix[i+1][j].prevInsertionBase, curColor, (uint8_t*)&curPrevInsertionBase)) {
										fprintf(stderr, "matrix[i+1=%d][j=%d].prevInsertionBase=%c\tcurColor=%c\n",
												i,
												j,
												matrix[i+1][j].prevInsertionBase, 
												curColor);
										PrintError(FnName,
												"curPrevInsertionBase",
												"Could not convert base and color",
												Exit,
												OutOfRange);
									}
									break;
								default:
									/*
									   fprintf(stderr, "l=%d\n", l);
									   */
									PrintError(FnName,
											"l",
											"Could not understand character",
											Exit,
											OutOfRange);
							}
							if(curScore < NEGATIVE_INFINITY/2) {
								curScore = NEGATIVE_INFINITY;
								curScoreNT = NEGATIVE_INFINITY;
							}
							if(curScore > maxScore && l != 4) {
								maxScore = curScore;
								maxScoreNT = curScoreNT;
								maxFrom = curFrom;
								maxColorError = '0';
								maxLength = curLength;
								maxPrevInsertionBase = curPrevInsertionBase;
							}
							break;
						default:
							PrintError(FnName,
									"k",
									"Could not understand k",
									Exit,
									OutOfRange);
					}
				}
				/* Update */
				/*
				   if(maxFrom < 0) {
				   fprintf(stderr, "(i,j,k)=(%d,%d,%d)\n",
				   i,
				   j,
				   k);
				   }
				   */
				assert(maxFrom >= 0);
				matrix[i+1][j+1].score[k] = maxScore;
				matrix[i+1][j+1].scoreNT[k] = maxScoreNT;
				matrix[i+1][j+1].from[k] = maxFrom;
				matrix[i+1][j+1].colorError[k] = maxColorError;
				matrix[i+1][j+1].length[k] = maxLength;
				/* Update the previous insertion base */
				if(k==4) {
					matrix[i+1][j+1].prevDeletionBase = maxPrevDeletionBase;
				}
				else if(k==5) {
					matrix[i+1][j+1].prevInsertionBase = maxPrevInsertionBase;
				}
				/* Uncomment this to remove insertions */
				/*
				   if(k>=5) {
				   matrix[i+1][j+1].score[k] = NEGATIVE_INFINITY; 
				   }
				   */
			}
		}
	}

	offset = FillAlignEntryFromMatrix(aEntry,
			matrix,
			read,
			readLength,
			reference,
			referenceLength,
			1,
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
