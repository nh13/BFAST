#ifndef BGENERATEREADS_H_
#define BGENERATEREADS_H_

#include "../blib/RGBinary.h"

typedef struct {
	char *readOne;
	int readOneType[SEQUENCE_LENGTH];
	char *readTwo;
	int readTwoType[SEQUENCE_LENGTH];
	int chr;
	int pos;
	char strand;
	int whichReadIndel;
	int startIndel;
	/* Do not modify */
	int readLength;
	int pairedEnd;
	int pairedEndLength;
	int indelLength;
} Read;

enum {
	Default,
	Insertion,
	SNP,
	Error,
	InsertionAndSNP,
	InsertionAndError,
	SNPAndError,
	InsertionSNPAndError
};

void ReadInitialize(Read*);
void ReadDelete(Read*);
void ReadPrint(Read*,
		FILE*);

void GenerateReads(RGBinary*, int, int, int, int, int, int, int, int, int, int);
void GetRandomRead(RGBinary*, int64_t, Read*);
void GetRandomChrPos(RGBinary*, int64_t, int*, int*, char*);
int ModifyRead(RGBinary*, Read*, int, int, int, int, int, int);
int InsertIndel(RGBinary*, Read*, int, int);
void InsertSNPs(Read*, int, int);
void InsertColorErrors(Read*, int, int);
	
#endif
