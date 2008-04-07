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
	unsigned int indexOne;
	unsigned int indexTwo;
	unsigned int numEntries;
	unsigned int *positions;
	unsigned char *chromosomes;
} RGTreeNode;

/* TODO*/
typedef struct {
	RGTreeNode *nodes;
	unsigned int numNodes;
	unsigned int gap;
	unsigned int matchLength;
	unsigned int startChr;
	unsigned int startPos;
	unsigned int endChr;
	unsigned int endPos;
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
	int *offset;
} RGSeqPair;

/* TODO */
typedef struct {
	char *sequence;
	unsigned int startPos;
	unsigned int endPos;
	unsigned int chromosome;
} RGChr;

/* TODO */
typedef struct {
	RGChr *chromosomes;
	unsigned int numChrs;
	unsigned int startChr;
	unsigned int startPos;
	unsigned int endChr;
	unsigned int endPos;
} RGList;

/* TODO */
typedef struct {
	unsigned int numEntries;
	unsigned char *index; 
	unsigned int *positions;
	unsigned char *chromosomes;
} RGIndexNode;

/* TODO */
typedef struct {
	RGIndexNode *nodes;
	unsigned int numNodes;
	unsigned int matchLength;
	unsigned int startChr;
	unsigned int startPos;
	unsigned int endChr;
	unsigned int endPos;
} RGIndex;

typedef struct {
	unsigned int chromosome;
	unsigned int startPos;
	unsigned int endPos;
	unsigned char *sequence; /* Store in the bytes via two bits - four nt per char (assuming sizeof(char)==1) */
} RGBinaryChr;

typedef struct {
	RGBinaryChr *chromosomes;
	unsigned int numChrs;
	unsigned int startChr;
	unsigned int startPos;
	unsigned int endChr;
	unsigned int endPos;
} RGBinary;

#endif
