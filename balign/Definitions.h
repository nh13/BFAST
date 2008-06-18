#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#define NEGATIVE_INFINITY INT_MIN/2 /* cannot make this too small, otherwise we will not have numerical stability, i.e. become positive */

/* Algorithm command line options:
 * 0: Dynamic programming 
 * */
#define MIN_ALGORITHM 0
#define MAX_ALGORITHM 1
#define DEFAULT_ALGORITHM 0

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

#endif
