#ifndef BLIBDEFINITIONS_H_
#define BLIBDEFINITIONS_H_

#include <sys/types.h>
#include <stdint.h>

/* Program defaults */
#define PROGRAM_NAME "bfast"
#define DEFAULT_FILENAME "Default.txt"
#define MAX_FILENAME_LENGTH 2048
#define DEFAULT_OUTPUT_ID "OutputID"
#define DEFAULT_OUTPUT_DIR "\0"
#define BREAK_LINE "************************************************************\n"
#define SEQUENCE_LENGTH 2048
#define SEQUENCE_NAME_LENGTH 4028
/* 0 - quick sort in place
 * 1 - merge sort with tmp file I/O (not fully debugged) */
#define SORT_TYPE 1
#define ONE_GIGABYTE 1073741824

/* File extensions */
#define BFAST_RG_FILE_EXTENSION "brg"
#define BFAST_INDEX_FILE_EXTENSION "bif"
#define BFAST_MATCHES_FILE_EXTENSION "bmf"
#define BFAST_ALIGN_FILE_EXTENSION "baf"
#define BFAST_NOT_ALIGNED_FILE_EXTENSION "bnf"

#define RGMATCH_MERGE_ROTATE_NUM 1000000
#define READ_ROTATE_NUM 1000000
#define RGINDEX_ROTATE_NUM 1000000
#define SORT_ROTATE_INC 0.01
#define RGINDEX_SORT_ROTATE_INC 0.001
#define ALIGN_ROTATE_NUM 100000
#define PARTITION_MATCHES_ROTATE_NUM 100000
#define BFAST_TMP_TEMPLATE ".bfast.tmp.XXXXXX"

/* For printing to stderr */
#define VERBOSE 0
#define DEBUG 10

/* Algorithm defaults */
#define DEFAULT_MATCH_LENGTH 11
#define ALPHABET_SIZE 4
#define FORWARD '+'
#define REVERSE '-'
#define GAP '-'
#define NULL_LETTER 'N'

/* For RGIndex.c */
enum {KILOBYTES, MEGABYTES, GIGABYTES};

/************************************/
/* 		Data structures 			*/
/************************************/

/* TODO */
typedef struct {
	uint8_t *chromosomes;
	uint32_t *positions;
	int8_t *strand;
	int32_t numEntries;
	int32_t maxReached;
} RGMatch;

/* TODO */
typedef struct { 
	int64_t *startIndex;
	int64_t *endIndex;
	int8_t *strand;
	int32_t numEntries;
} RGRanges;

/* TODO */
typedef struct {
	int32_t numReads;
	char **reads;
	int32_t *readLength;
	int8_t *strand;
	int32_t *offset;
} RGReads;

/* TODO */
typedef struct {
	int32_t chromosome;
	int32_t startPos;
	int32_t endPos;
	uint32_t numBytes;
	uint8_t *sequence; 
} RGBinaryChr;

/* TODO */
typedef struct {
	RGBinaryChr *chromosomes;
	int32_t numChrs;
	int32_t startChr;
	int32_t startPos;
	int32_t endChr;
	int32_t endPos;
} RGBinary;

/* TODO */
typedef struct {
	/* Index storage */
	uint32_t *positions;
	uint8_t *chromosomes;
	uint32_t length;

	/* Hash storage */
	uint32_t hashWidth; /* in bases */
	int64_t hashLength; 
	uint32_t *starts;
	uint32_t *ends;

	/* Index definition */
	int32_t totalLength;
	int32_t numTiles;
	int32_t *tileLengths;
	int32_t *gaps; /* There should be numTiles - 1 gaps */
	int32_t repeatMasker;

	/* Index range */
	int32_t startChr;
	int32_t startPos;
	int32_t endChr;
	int32_t endPos;
} RGIndex;

/* TODO */
typedef struct {
	int32_t numIndexes;
	int32_t *hashLengths;
	int32_t *numTiles;
	int32_t **tileLengths;
	int32_t **gaps;
} RGIndexLayout;

/* TODO */
typedef struct {
	RGIndex *index;
	RGBinary *rg;
	int64_t low;
	int64_t high;
	int32_t threadID;
	int32_t showPercentComplete;
	char *tmpDir;
} ThreadRGIndexSortData;

#endif
