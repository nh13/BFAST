#ifndef ALIGN_H_
#define ALIGN_H_
#include "BLibDefinitions.h"

/* Align.c specific definitions */

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
	InsertionExtension, /* 5 */
	NoFromNT}; /* 6 */

/* For the "from" in the struct "AlignCellCS" */
enum {
	StartCS, /* 0 */
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
	InsertionN, /* 15 */
	NoFromCS /* 16 */
};

int AlignRGMatches(RGMatches*, RGBinary*, AlignedRead*, int32_t, int32_t, ScoringMatrix*, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, AlignMatrixNT***, AlignMatrixCS***, int32_t*, int32_t*);
void AlignRGMatchesOneEnd(RGMatch*, RGBinary*, AlignedEnd*, int32_t, int32_t, ScoringMatrix*, int32_t, int32_t, int32_t, double*, int32_t*, AlignMatrixNT***, AlignMatrixCS***, int32_t*, int32_t*);
int32_t AlignExact(char*, int32_t, char*, int32_t, ScoringMatrix*, AlignedEntry*, char, int32_t, int32_t);
void AlignUngapped(char*, int32_t, char*, int32_t, int32_t, ScoringMatrix*, AlignedEntry*, int32_t, int32_t, uint32_t, char);
void AlignFullWithBound(char*, int32_t, char*, int32_t, ScoringMatrix*, AlignedEntry*, int32_t, char, int32_t, double, AlignMatrixNT***, AlignMatrixCS***);
int32_t AlignRGMatchesKeepBestScore(AlignedEnd*, double);

#endif
