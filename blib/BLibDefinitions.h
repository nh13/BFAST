#ifndef BLIBDEFINITIONS_H_
#define BLIBDEFINITIONS_H_

#include <sys/types.h>
#include <stdint.h>

/* Obsolete ? HERE */
/*
#define SEQUENCE_LENGTH 2048
#define SEQUENCE_NAME_LENGTH 4028
#define MAX_HEADER_LENGTH 2048
*/

/* Program defaults */
#define PROGRAM_NAME "bfast" /* Could just use PACKAGE_NAME */
#define DEFAULT_FILENAME "Default.txt"
#define MAX_FILENAME_LENGTH 2048
#define DEFAULT_OUTPUT_ID "OutputID"
#define DEFAULT_OUTPUT_DIR "\."
#define BREAK_LINE "************************************************************\n"
#define MAX_MASK_LENGTH 1024
#define MAX_CONTIG_NAME_LENGTH 2048
#define MAX_CONTIG_LOG_10 6
#define MAX_POSITION_LOG_10 10
#define ONE_GIGABYTE (int64_t)1073741824
#define MERGE_MEMORY_LIMIT 12*((int64_t)1073741824) /* In Gigabytes */

/* Testing/Debug */
#define TEST_BINDEX_SORT 0

/* Default output */
enum {TextOutput, BinaryOutput};
enum {TextInput, BinaryInput};
#define BPREPROCESS_DEFAULT_OUTPUT 1 /* 0: text 1: binary */
#define BMATCHES_DEFAULT_OUTPUT 1 /* 0: text 1: binary */
#define BALIGN_DEFAULT_OUTPUT 1 /* 0: text 1: binary */

/* File extensions */
#define BFAST_RG_FILE_EXTENSION "brg"
#define BFAST_INDEX_FILE_EXTENSION "bif"
#define BFAST_MATCHES_FILE_EXTENSION "bmf"
#define BFAST_MATCHES_READS_FILTERED_FILE_EXTENSION "fa"
#define BFAST_ALIGN_FILE_EXTENSION "baf"
#define BFAST_NOT_ALIGNED_FILE_EXTENSION "bnaf"
#define BFAST_MAF_FILE_EXTENSION "maf"

#define BMATCH_MERGE_ROTATE_NUM 100000
#define READ_ROTATE_NUM 1000000
#define BINDEX_ROTATE_NUM 1000000
#define SORT_ROTATE_INC 0.01
#define BINDEX_SORT_ROTATE_INC 0.001
#define ALIGN_ROTATE_NUM 10000
#define PARTITION_MATCHES_ROTATE_NUM 100000
#define ALIGNENTRIES_READ_ROTATE_NUM 10000
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
#define COLOR_SPACE_START_NT 'A'
#define BFAST_ID 'B'+'F'+'A'+'S'+'T'

enum {KILOBYTES, MEGABYTES, GIGABYTES};
enum {Contig_8, Contig_32};
enum {SingleEnd, PairedEnd, PairedEndDoesNotMatter};
enum {NTSpace, ColorSpace, SpaceDoesNotMatter};
enum {AlignEntrySortByAll, AlignEntrySortByContigPos};
enum {IgnoreExons, UseExons};

/************************************/
/* 		Data structures 			*/
/************************************/

/* TODO */
typedef struct {
	BString read;
	int32_t maxReached;
	int32_t numEntries;
	uint32_t *contigs;
	int32_t *positions;
	int8_t *strand;
} BMatch;

/* TODO */
typedef struct {
	int32_t pairedEnd;
	BString readName;
	RGMatch matchOne;
	RGMatch matchTwo;
} BMatches;

/* TODO */
typedef struct { 
	int64_t *startIndex;
	int64_t *endIndex;
	int8_t *strand;
	int32_t *offset;
	int32_t numEntries;
} BRanges;

/* TODO */
typedef struct {
	int32_t numReads;
	BString *reads;
	int8_t *strand;
	int32_t *offset;
} BReads;

/* TODO */
typedef struct {
	/* Storage */
	BString contigName;
	/* Metadata */
	uint8_t *sequence; 
	int32_t sequenceLength;
	uint32_t numBytes;
} BReferenceGenomeContig;

/* TODO */
typedef struct {
	/* Storage type */
	int32_t id;
	/* B storage */
	BReferenceGenomeContig *contigs;
	int32_t numContigs;
	/* Metadata */
	int32_t space;
} BReferenceGenome;

/* TODO */
typedef struct {
	/* Storage type */
	int32_t id;
	/* Index storage */
	uint8_t *contigs_8;
	uint32_t *contigs_32;
	int32_t *positions;
	int64_t length;
	int32_t contigType; 
	/* Index range */
	int32_t startContig;
	int32_t startPos;
	int32_t endContig;
	int32_t endPos;
	/* Index layout */
	int32_t width;
	int32_t keysize;
	int32_t *mask;
	/* Index properties */
	int32_t repeatMasker;
	int32_t space;
	/* Hash storage */
	uint32_t hashWidth; /* in bases */
	int64_t hashLength; 
	uint32_t *starts;
	uint32_t *ends;

} BIndex;

/* TODO */
typedef struct {
	int32_t numIndexes;
	int32_t *hashWidths;
	int32_t **masks;
	int32_t *widths;
	int32_t *keysizes;
} BIndexLayout;

/* TODO */
typedef struct {
	uint32_t startContig;
	uint32_t startPos;
	uint32_t endContig;
	uint32_t endPos;
} BIndexExon;

/* TODO */
typedef struct {
	int numExons;
	BIndexExon *exons;
} BIndexExons;

/* TODO */
typedef struct {
	BString contigName;
	uint32_t contig;
	uint32_t position;
	char strand;
	double score;
	uint32_t referenceLength; /* The length of the reference alignment substracting gaps */
	uint32_t length; /* The length of the alignment */
	BString read;
	BString reference;
	BString colorError;
} AlignEntry;

/* TODO */
typedef struct {
	BString readName;
	int32_t pairedEnd;
	int32_t space;
	int32_t numEntriesOne;
	int32_t numEntriesTwo;
	AlignEntry *entriesOne;
	AlignEntry *entriesTwo;
} AlignEntries;

/* TODO */
typedef struct {
	BIndex *index;
	BReferenceGenome *rg;
	int32_t space;
	int64_t low;
	int64_t high;
	int32_t threadID;
	int32_t showPercentComplete;
	char *tmpDir;
	int64_t mergeMemoryLimit;
} ThreadBIndexSortData;

/* TODO */
typedef struct {
	BIndex *index;
	BReferenceGenome *rg;
	int32_t threadID;
	int64_t low;
	int64_t mid;
	int64_t high;
	int64_t mergeMemoryLimit;
	char *tmpDir;
} ThreadBIndexMergeData;

#endif
