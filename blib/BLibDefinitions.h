#ifndef BLIBDEFINITIONS_H_
#define BLIBDEFINITIONS_H_
#define DEFAULT_FILENAME "Default.txt"
#define MAX_FILENAME_LENGTH 2048
#define DEFAULT_OUTPUT_ID "OutputID"
#define DEFAULT_OUTPUT_DIR ""
#define BLATTER_TREE_FILE_EXTENSION "btf"
#define BLATTER_MATCHES_FILE_EXTENSION "bmf"
#define DEFAULT_MATCH_LENGTH 11
#define BREAK_LINE "************************************************************\n"

#define IN_PLACE 0
#define MERGE_SORT_MIN_PERCENT 0.01

#define VERBOSE 11
#define DEBUG 10

#define ALPHABET_SIZE 4

#define SEQUENCE_LENGTH 512
#define SEQUENCE_NAME_LENGTH 4028

/* For SRTree.c */
#define SRT_SEQUENCE_NAME_LENGTH 1024 
#define SRT_SEQUENCE_LENGTH 1024

/* For BError.c  */
enum {Exit, Warn, LastActionType};
enum {
	Dummy,
	OutOfRange, /* e.g. command line args */
	IllegalFileName,   /*  KeepAdding */
	LastErrorType
};

/* For RGTree.c */
enum {RGT_KILOBYTES, RGT_MEGABYTES, RGT_GIGABYTES};

/************************************/
/* 		Data structures 			*/
/************************************/

/* TODO*/
/* We could package this better, i.e. pack
 * numEntries and indexOne and indexTwo.
 * */
typedef struct {
	int indexOne;
	int indexTwo;
	int numEntries;
	int *positions;
	char *chromosomes;
} RGNode;

/* TODO*/
typedef struct {
	RGNode *nodes;
	int numNodes;
	int gap;
	int matchLength;
	int startChr;
	int startPos;
	int endChr;
	int endPos;
} RGTree;

/* TODO */
typedef struct {
	int *positions;
	char *chromosomes;
	char *strand;
	int numEntries;
} RGMatch;

/* TODO */
typedef struct {
	int numPairs;
	int *indexOne;
	int *indexTwo;
	char *strand;
} RGSeqPair;

/* TODO */
typedef struct {
	char *sequence;
	int startPos;
	int endPos;
	int chromosome;
} RGChr;

/* TODO */
typedef struct {
	RGChr *chromosomes;
	int numChrs;
	int startChr;
	int startPos;
	int endChr;
	int endPos;
} RGList;

#endif
