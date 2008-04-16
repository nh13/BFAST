#ifndef BLIBDEFINITIONS_H_
#define BLIBDEFINITIONS_H_
#define DEFAULT_FILENAME "Default.txt"
#define MAX_FILENAME_LENGTH 2048
#define DEFAULT_OUTPUT_ID "OutputID"
#define DEFAULT_OUTPUT_DIR "\0"
#define BLATTER_TREE_FILE_EXTENSION "btf"
#define BLATTER_INDEX_FILE_EXTENSION "bif"
#define BLATTER_MATCHES_FILE_EXTENSION "bmf"
#define BLATTER_ALIGN_FILE_EXTENSION "baf"
#define DEFAULT_MATCH_LENGTH 11
#define READ_ROTATE_NUM 1000000
#define RGINDEX_ROTATE_NUM 1000000
#define SORT_ROTATE_INC 0.01
#define BREAK_LINE "************************************************************\n"

#define VERBOSE 0
#define DEBUG 10

#define ALPHABET_SIZE 4

#define SEQUENCE_LENGTH 512
#define SEQUENCE_NAME_LENGTH 4028

#define FORWARD '+'
#define REVERSE '-'

/* For RGIndex.c */
enum {KILOBYTES, MEGABYTES, GIGABYTES};

/************************************/
/* 		Data structures 			*/
/************************************/

/* TODO */
typedef struct {
	int *positions;
	unsigned char *chromosomes;
	char *strand;
	int numEntries;
} RGMatch;

/* TODO */
typedef struct {
	int numReads;
	char **reads;
	char *strand;
	int *offset;
} RGReads;

typedef struct {
	unsigned int chromosome;
	unsigned int startPos;
	unsigned int endPos;
	unsigned char *sequence; 
} RGBinaryChr;

typedef struct {
	RGBinaryChr *chromosomes;
	int numChrs;
	int startChr;
	int startPos;
	int endChr;
	int endPos;
} RGBinary;

/* TODO */
typedef struct {
	/* Index storage */
	unsigned int *positions;
	unsigned char *chromosomes;
	unsigned int length;

	/* Index definition */
	unsigned int totalLength;
	unsigned int numTiles;
	unsigned int *tileLengths;
	unsigned int *gaps; /* There should be numTiles - 1 gaps */

	/* Index range */
	unsigned int startChr;
	unsigned int startPos;
	unsigned int endChr;
	unsigned int endPos;
} RGIndex;

/* TODO */
typedef struct {
	int numIndexes;
	int *numTiles;
	int **tileLengths;
	int **gaps;
} RGIndexLayout;

/* TODO */
typedef struct {
	RGIndex *index;
	RGBinary *rg;
	unsigned int low;
	unsigned int high;
	int threadID;
	int showPercentComplete;
} ThreadRGIndexSortData;

#endif
