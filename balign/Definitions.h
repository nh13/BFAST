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
	int32_t score;
	int32_t from;
	int32_t length;
} AlignCellNT;

/* TODO */
typedef struct {
	AlignCellNT h;
	AlignCellNT s;
	AlignCellNT v;
} AlignMatrixNT;

typedef struct {
	/* For A, C, G, T, and N */
	int32_t score[ALPHABET_SIZE+1];
	int32_t scoreNT[ALPHABET_SIZE+1];
	int32_t from[ALPHABET_SIZE+1];
	int32_t length[ALPHABET_SIZE+1];
	char colorError[ALPHABET_SIZE+1];
} AlignCellCS;

/* TODO */
typedef struct {
	/* Deletion */
	AlignCellCS h;
	/* Match/Mismatch  */
	AlignCellCS s;
	/* Insertion */
	AlignCellCS v;
} AlignMatrixCS;

/* For the "from" in the struct "AlignCellNT" */
enum {StartNT, /* 0 */
	DeletionStart, /* 1 */
	DeletionExtension, /* 2 */
	Match, /* 3 */
	InsertionStart, /* 4 */
	InsertionExtension}; /* 5 */

/* For the "from" in the struct "AlignCellCS" */
enum {
	Start, /* 0 */
	DeletionA, /* 1 */
	DeletionC, /* 2 */
	DeletionG, /* 3 */
	DeletionT, /* 4 */
	DeletionN, /* 5 */
	MatchA, /* 6 */
	MatchC, /* 7 */
	MatchG, /* 8 */
	MatchT, /* 9 */
	MatchN, /* 10 */
	InsertionA, /* 11 */
	InsertionC, /* 12 */
	InsertionG, /* 13 */
	InsertionT, /* 14 */
	InsertionN  /* 15 */
};

#endif
