#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include "../blib/BLibDefinitions.h"

#define NEGATIVE_INFINITY INT_MIN/2 /* cannot make this too small, otherwise we will not have numerical stability, i.e. become positive */
#define VERY_NEGATIVE_INFINITY (INT_MIN/2)-1000 /* cannot make this too small, otherwise we will not have numerical stability, i.e. become positive */

/* Algorithm command line options:
 * 0: Dynamic programming 
 * */
#define MIN_ALGORITHM 0
#define MAX_ALGORITHM 1
#define ALIGNMATRIXCELL_NUM_SUB_CELLS 6
#define COLOR_MATCH 0
#define COLOR_ERROR -1

/* Align.c specific definitions */
typedef struct {
	double gapOpenPenalty;
	double gapExtensionPenalty;
	char *NTKeys;
	double **NTScores;
	int *ColorKeys;
	double **ColorScores;
} ScoringMatrix;

/* TODO */
typedef struct {
	/* add three for the insertion, deletion and best score */
	double score[ALIGNMATRIXCELL_NUM_SUB_CELLS]; /* current score */
	double scoreNT[ALIGNMATRIXCELL_NUM_SUB_CELLS]; /* current score */
	int from[ALIGNMATRIXCELL_NUM_SUB_CELLS]; /* previous arc */
	int length[ALIGNMATRIXCELL_NUM_SUB_CELLS]; /* path length */
	char colorError[ALIGNMATRIXCELL_NUM_SUB_CELLS]; /* color error */
	char prevInsertionBase; /* When we create a run of insertions, we must keep track of the bases in the read */
} AlignMatrix;

/* For the "from" in the struct "AlignMatrixCell" */
enum {
	Start, /* 0 */
	/* From diagonals */
	DiagA, /* 1 */
	DiagC, /* 2 */
	DiagG, /* 3 */
	DiagT, /* 4 */
	/* From deletion */
	/* Start of deletion */
	DeletionA, /* 5 */
	DeletionC, /* 6 */
	DeletionG, /* 7 */
	DeletionT, /* 8 */
	/* Extension of deletion */
	DeletionExt, /* 9 */
	/* End of deletion */
	DeletionEnd, /* 10 */
	/* From insertion */
	InsertionA, /* 11 */
	InsertionC, /* 12 */
	InsertionG, /* 13 */
	InsertionT, /* 14 */
	InsertionExt, /* 15 */
	InsertionEnd, /* 16 */
};

#endif
