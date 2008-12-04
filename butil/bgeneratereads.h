#ifndef BGENERATEREADS_H_
#define BGENERATEREADS_H_

#include "../blib/RGBinary.h"

typedef struct {
	char *readOne;
	int readOneType[SEQUENCE_LENGTH];
	char *readTwo;
	int readTwoType[SEQUENCE_LENGTH];
	int contig;
	int pos;
	int readNum;
	char strand;
	int whichReadVariants;
	int startIndel;
	/* Do not modify */
	int readLength;
	int pairedEnd;
	int pairedEndLength;
	int indelLength;
} Read;

enum {
	Default,				/* 0 */
	Insertion,				/* 1 */
	SNP,					/* 2 */
	Error,					/* 3 */
	InsertionAndSNP,		/* 4 */
	InsertionAndError,		/* 5 */
	SNPAndError,			/* 6 */
	InsertionSNPAndError	/* 7 */
};

void ReadInitialize(Read*);
void ReadDelete(Read*);
void ReadPrint(Read*,
		FILE*);

void GenerateReads(RGBinary*, int, int, int, int, int, int, int, int, int, int);
void GetRandomRead(RGBinary*, int64_t, Read*);
void GetRandomContigPos(RGBinary*, int64_t, int*, int*, char*);
int ModifyRead(RGBinary*, Read*, int, int, int, int, int, int);
int InsertIndel(RGBinary*, Read*, int, int);
void InsertMismatches(Read*, int, int, int);
void InsertMismatchesHelper(char*, int, int*, int, int, int);
void InsertColorErrors(Read*, int, int);
	
#endif
