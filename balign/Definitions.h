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

enum {
	HORIZONTAL,		/* 0 */ 
	VERTICAL,		/* 1 */ 
	DIAGONAL,		/* 2 */ 
	GAP_OPEN,		/* 3 */ 
	GAP_EXTENSION, 	/* 4 */
	NO_TYPE 		/* 5 */
};

/* Align.c specific definitions */
typedef struct {
	double gapOpenPenalty;
	double gapExtensionPenalty;
	int alphabetSize; /* = ALPHABET_SIZE */
	char *key;
	double **scores;
	double colorMatch;
	double colorError;
} ScoringMatrix;

/* Structure for the dynamic programming with affine gap penalties */
typedef struct {
	double hScore; /* horizontal score */
	int hType;
	double vScore; /* vertical score */
	int vType;
	double dScore; /* diagonal score */
	int dType;
	int prevRow; /* previous row */
	int prevCol; /* previous column */
} MatrixEntry;

/* TODO */
typedef struct {
	/* add three for the insertion, deletion and best score */
	double score[ALIGNMATRIXCELL_NUM_SUB_CELLS];
	int from[ALIGNMATRIXCELL_NUM_SUB_CELLS];
	int length[ALIGNMATRIXCELL_NUM_SUB_CELLS];
} AlignMatrix;

/* For the "from" in the struct "AlignMatrixCell" */
enum {
	Start,
	/* From diagonals */
	DiagA,
	DiagC,
	DiagG,
	DiagT,
	/* From deletion */
	/* Start of deletion */
	DeletionA,
	DeletionC,
	DeletionG,
	DeletionT,
	/* Extension of deletion */
	DeletionExt,
	/* End of deletion */
	DeletionEnd,
	/* From insertion */
	InsertionA,
	InsertionC,
	InsertionG,
	InsertionT,
	InsertionExt,
	InsertionEnd,
};

#endif
