#ifndef BLIBDEFINITIONS_H_
#define BLIBDEFINITIONS_H_
#define DEFAULT_FILENAME "Default.txt"
#define MAX_FILENAME_LENGTH 2048
#define DEFAULT_OUTPUT_ID "OutputID"
#define DEFAULT_OUTPUT_DIR ""
#define BLATTER_TREE_FILE_EXTENSION "btf"
#define BLATTER_INDEX_FILE_EXTENSION "bif"
#define BLATTER_MATCHES_FILE_EXTENSION "bmf"
#define BLATTER_ALIGN_FILE_EXTENSION "baf"
#define DEFAULT_MATCH_LENGTH 11
#define BREAK_LINE "************************************************************\n"

#define VERBOSE 0
#define DEBUG 10

#define ALPHABET_SIZE 4

#define SEQUENCE_LENGTH 512
#define SEQUENCE_NAME_LENGTH 4028

#define FORWARD '+'
#define REVERSE '-'

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
	unsigned char *chromosomes;
} RGTreeNode;

/* TODO*/
typedef struct {
	RGTreeNode *nodes;
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
	unsigned char *chromosomes;
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

/* TODO */
typedef struct {
	int numEntries;
	unsigned char *index; 
	int *positions;
	unsigned char *chromosomes;
} RGIndexNode;

/* TODO */
typedef struct {
	RGIndexNode *nodes;
	int numNodes;
	int matchLength;
	int startChr;
	int startPos;
	int endChr;
	int endPos;
} RGIndex;

typedef struct {
	int chromosome;
	int startPos;
	int endPos;
	char *sequence; /* Store in the bytes via two bits - four nt per char (assuming sizeof(char)==1) */
} RGBinaryChr;

typedef struct {
	RGBinaryChr *chromosomes;
	int numChrs;
	int startChr;
	int startPos;
	int endChr;
	int endPos;
} RGBinary;

#endif
