#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#define NEGATIVE_INFINITY INT_MIN/2 /* cannot make this too small, otherwise we will not have numerical stability, i.e. become positive */

/* Algorithm command line options:
 * 0: Dynamic programming 
 * */
#define MIN_ALGORITHM 0
#define MAX_ALGORITHM 1
#define DEFAULT_ALGORITHM 0

/* Align.c specific definitions */
typedef struct {
	double gapOpenPenalty;
	double gapExtensionPenalty;
	int alphabetSize; /* = ALPHABET_SIZE */
	char *key;
	double **scores;
} ScoringMatrix;

/* Structure for the dynamic programming with affine gap penalties */
typedef struct {
	int hScore; /* horizontal score */
	int vScore; /* vertical score */
	int dScore; /* diagonal score */
	int prevRow; /* previous row */
	int prevCol; /* previous column */
} MatrixEntry;

/* TODO */
typedef struct {
	char *read; /* The read */
	char *reference;
	int length; /* The length of the alignment */
	int chromosome;
	int position;
	int strand;
	double score;
} AlignOutput;

#endif
