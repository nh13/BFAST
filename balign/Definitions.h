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

enum {FullAlignment, MismatchesOnly};

/* Align.c specific definitions */
typedef struct {
	int32_t gapOpenPenalty;
	int32_t gapExtensionPenalty;
	char *NTKeys;
	int32_t **NTScores;
	int32_t *ColorKeys;
	int32_t **ColorScores;
} ScoringMatrix;

/* TODO */
typedef struct {
	/* add three for the insertion, deletion and best score */
	int32_t score[ALIGNMATRIXCELL_NUM_SUB_CELLS]; /* current score */
	int32_t scoreNT[ALIGNMATRIXCELL_NUM_SUB_CELLS]; /* current score */
	int32_t from[ALIGNMATRIXCELL_NUM_SUB_CELLS]; /* previous arc */
	int32_t length[ALIGNMATRIXCELL_NUM_SUB_CELLS]; /* path length */
	char colorError[ALIGNMATRIXCELL_NUM_SUB_CELLS]; /* color error */
	char prevDeletionBase; /* When we create a run of deletions, we must keep track of the base used right before the deletion */
	char prevInsertionBase; /* When we create a run of insertions, we must keep track of the base used right before the insertion */
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
