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

	/* HERE 40 */
	/*
	   fprintf(stderr, "HERE 40\nreference=%s\nread=%s\n",
	   reference,
	   read);
	   */

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
			matrix[i][0].prevInsertionBase = COLOR_SPACE_START_NT;
		}
	}
	/* Row 0 column j should be zero since we want to find the best
	 * local alignment within the reference */
	for(j=0;j<referenceLength+1;j++) {
		for(k=0;k<ALIGNMATRIXCELL_NUM_SUB_CELLS;k++) {
			/* Assumes both DNA and COLOR_SPACE_START_NT are upper case */
			if(DNA[k] == COLOR_SPACE_START_NT) { 
				/* Starting adaptor NT */
				matrix[0][j].score[k] = 0;
			}
			else {
				matrix[0][j].score[k] = NEGATIVE_INFINITY;
			}
			matrix[0][j].from[k] = Start;
			matrix[0][j].length[k] = 0;
			matrix[0][j].colorError[k] = '0';
			matrix[0][j].prevInsertionBase = COLOR_SPACE_START_NT;
		}
	}

	/* Fill in the matrix according to the recursive rules */
	for(i=0;i<readLength;i++) { /* read/rows */
		/* Get the current color */
		uint8_t curColor;
		char curReadBase, prevReadBase;
		/* In color space, the first color is determined by the adapter NT */
		prevReadBase = (i<=0)?COLOR_SPACE_START_NT:read[i-1];
		curReadBase = read[i];

		/* Get the current color for the read */
		curColor = ConvertBaseToColorSpace(prevReadBase, curReadBase);

		for(j=0;j<referenceLength;j++) { /* reference/columns */

			for(k=0;k<ALIGNMATRIXCELL_NUM_SUB_CELLS;k++) { /* To NT */
				char DNA[4] = "ACGT";
				double maxScore = NEGATIVE_INFINITY-1;
				int maxFrom = -1;
				char maxColorError = '0';
				int maxLength = 0;
				char maxPrevInsertionBase;

				for(l=0;l<ALIGNMATRIXCELL_NUM_SUB_CELLS;l++) { /* From NT */
					double curScore=NEGATIVE_INFINITY;
					int curFrom=-1;
					int curLength=-1;
					char curPrevInsertionBase;
					uint8_t convertedColor='X';
					switch(k) {
						case 0:
						case 1:
						case 2:
						case 3:
							/* Check diagonals from NT */
							/* Previous score */ 
							curScore = matrix[i][j].score[l];
							curLength = matrix[i][j].length[l] + 1;
							/* Plus score for colors */
							switch(l) {
								case 0:
								case 1:
								case 2:
								case 3:
									/* Diagonals */
									convertedColor=ConvertBaseToColorSpace(DNA[l], DNA[k]);
									break;
								case 4:
									/* Use previous base for deletions */
									convertedColor=ConvertBaseToColorSpace(prevReadBase, DNA[k]);
									break;
								case 5:
									/* Use the previous base carried along by insertions */
									convertedColor=ConvertBaseToColorSpace(matrix[i][j].prevInsertionBase, DNA[k]);
									break;
								default:
									PrintError(FnName,
											"l",
											"Could not understand l",
											Exit,
											OutOfRange);
									break;
							}
							/* HERE A5 */
							/*
							   if(strcmp(read, "AGCTTTTCATTCTGACTGCAACGGT")==0 &&
							   i==readLength-1 &&
							   i==j) {
							   fprintf(stderr, "HERE A5\n");
							   char stuff[6] = "ACGTDI";
							   fprintf(stderr, "(%d,%d,to=%c,from=%c)=>%d,%d,%c,%c,%lf,%lf,%lf,%lf(from=%c,to=%c)\n",
							   i+1,
							   j+1,
							   stuff[k],
							   stuff[l],
							   (int)curColor,
							   (int)convertedColor,
							   (int)reference[j],
							   (int)stuff[k],
							   curScore,
							   ScoringMatrixGetColorScore(curColor, convertedColor, sm),
							   ScoringMatrixGetNTScore(reference[j], DNA[k], sm),
							   curScore + ScoringMatrixGetColorScore(curColor, convertedColor, sm) + ScoringMatrixGetNTScore(reference[j], DNA[k], sm),
							   stuff[l],
							   stuff[k]
							   );
							   }
							   */
							/* Add score for color error, if any */
							curScore += ScoringMatrixGetColorScore(curColor,
									convertedColor,
									sm);
							/* Add score for NT */
							curScore += ScoringMatrixGetNTScore(reference[j], DNA[k], sm);

							if(curScore < NEGATIVE_INFINITY/2) {
								curScore = NEGATIVE_INFINITY;
							}

							/* Check to see if this is better than the max */
							if(curScore > maxScore) {
								maxScore = curScore;
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
										fprintf(stderr, "l=%d\n", l);
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
									curFrom = DeletionA;
									break;
								case 1:
									curScore = matrix[i+1][j].score[l] + sm->gapOpenPenalty;
									curFrom = DeletionC;
									break;
								case 2:
									curScore = matrix[i+1][j].score[l] + sm->gapOpenPenalty;
									curFrom = DeletionG;
									break;
								case 3:
									curScore = matrix[i+1][j].score[l] + sm->gapOpenPenalty;
									curFrom = DeletionT;
									break;
								case 4:
									curScore = matrix[i+1][j].score[l] + sm->gapExtensionPenalty;
									curFrom = DeletionExt;
									break;
								case 5:
									break;
								default:
									fprintf(stderr, "l=%d\n", l);
									PrintError(FnName,
											"l",
											"Could not understand character",
											Exit,
											OutOfRange);
							}
							if(curScore < NEGATIVE_INFINITY/2) {
								curScore = NEGATIVE_INFINITY;
							}
							if(curScore > maxScore && l != 5) {
								maxScore = curScore;
								maxFrom = curFrom;
								maxColorError = '0';
								maxLength = curLength;
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
									curFrom = InsertionA;
									curPrevInsertionBase = DNA[l];
									break;
								case 1:
									curScore = matrix[i][j+1].score[l] + sm->gapOpenPenalty;
									curFrom = InsertionC;
									curPrevInsertionBase = DNA[l];
									break;
								case 2:
									curScore = matrix[i][j+1].score[l] + sm->gapOpenPenalty;
									curFrom = InsertionG;
									curPrevInsertionBase = DNA[l];
									break;
								case 3:
									curScore = matrix[i][j+1].score[l] + sm->gapOpenPenalty;
									curFrom = InsertionT;
									curPrevInsertionBase = DNA[l];
									break;
								case 4:
									break;
								case 5:
									curScore = matrix[i][j+1].score[l] + sm->gapExtensionPenalty;
									curFrom = InsertionExt;
									curPrevInsertionBase = ConvertBaseAndColor(matrix[i+1][j].prevInsertionBase, 
											curColor);
									break;
								default:
									fprintf(stderr, "l=%d\n", l);
									PrintError(FnName,
											"l",
											"Could not understand character",
											Exit,
											OutOfRange);
							}
							if(curScore < NEGATIVE_INFINITY/2) {
								curScore = NEGATIVE_INFINITY;
							}
							if(curScore > maxScore && l != 4) {
								maxScore = curScore;
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
				if(maxFrom < 0) {
					fprintf(stderr, "(i,j,k)=(%d,%d,%d)\n",
							i,
							j,
							k);
				}
				assert(maxFrom >= 0);
				matrix[i+1][j+1].score[k] = maxScore;
				matrix[i+1][j+1].from[k] = maxFrom;
				matrix[i+1][j+1].colorError[k] = maxColorError;
				matrix[i+1][j+1].length[k] = maxLength;
				/* Update the previous insertion base */
				if(k==5) {
					matrix[i+1][j+1].prevInsertionBase = maxPrevInsertionBase;
				}
				/* Uncomment this to remove insertions */
				/*
				if(k>=5) {
					matrix[i+1][j+1].score[k] = NEGATIVE_INFINITY; 
				}
				*/

				/* HERE A4 */
				/*
				   if(strcmp(read, "AGCTTTTCATTCTGACTGCAACGGT")==0 &&
				   i==readLength-1 &&
				   i==j) {
				   fprintf(stderr, "(row,col,cell,score,length,from,colorError)=(%d,%d,%d,%lf,%d,%d,%c)\n",
				   i+1,
				   j+1,
				   k,
				   matrix[i+1][j+1].score[k],
				   matrix[i+1][j+1].length[k],
				   matrix[i+1][j+1].from[k],
				   matrix[i+1][j+1].colorError[k]
				   );
				   }
				   */
				/*
				   if(i == 18 && j== 21 && k == 5) {
				   fprintf(stderr, "HERE A4\n");
				   exit(1);
				   }
				   */
			}
		}
	}

	/* TODO */
	/* Get results and store them in the align entry */

	/* HERE A6 */
	/*
	   fprintf(stderr, "HERE A6\n");
	   exit(1);
	   */

	offset = FillAlignEntryFromMatrix(aEntry,
			matrix,
			read,
			readLength,
			reference,
			referenceLength,
			1,
			/* HERE */
			1);

	/* HERE E2 */
	int count = 0;
	for(i=0;i<aEntry->length;i++) {
		if(aEntry->read[i] != GAP) {
			count++;
		}
	}
	if(count != strlen(read)) {
		AlignEntryFree(aEntry);
		offset = FillAlignEntryFromMatrix(aEntry,
				matrix,
				read,
				readLength,
				reference,
				referenceLength,
				1,
				1);
		fprintf(stderr, "read=%s\naEntry->read=%s\n",
				read,
				aEntry->read);
	}
	assert(count == strlen(read));

	/* Free the matrix, free your mind */
	for(i=0;i<readLength+1;i++) {
		free(matrix[i]);
		matrix[i]=NULL;
	}
	free(matrix);
	matrix=NULL;

	/* HERE */
	/*
	   fprintf(stderr, "%s\n%s\n%s\n%lf\n",
	   aEntry->reference,
	   aEntry->read,
	   aEntry->colorError,
	   aEntry->score);
	   fprintf(stderr, "Exiting HERE\n");
	   exit(1);
	   */

	/* The return is the number of gaps at the beginning of the reference */
	return offset;
}
