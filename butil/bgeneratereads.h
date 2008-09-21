#ifndef BGENERATEREADS_H_
#define BGENERATEREADS_H_

#include "../blib/RGBinary.h"

typedef struct {
	char *readOne;
	char *readTwo;
	int chr;
	int pos;
	char strand;
	int readLength;
	int pairedEnd;
	int pairedEndLength;
} Read;

void GenerateReads(RGBinary*, int, int, int, int, int, int, int, int, int, int);
void GetRandomRead(RGBinary*, int64_t, Read*);
void ReadInitialize(Read*);
void GetRandomChrPos(RGBinary*, int64_t, int*, int*, char*);
	
#endif
