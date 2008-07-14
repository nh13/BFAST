#ifndef BMAFCONVERT_H_
#define BMAFCONVERT_H_

enum {MAFSortByChrPos};

typedef struct {
	int minPos;
	int maxPos;
	int chromosome;
	int numEntries;
	FILE *FP;
	char *FileName;
} TmpFile;

typedef struct {
	double score;
	uint32_t chromosome;
	uint32_t position;
	int8_t strand;
	int32_t alignmentLength;
	int32_t referenceLength;
	int32_t readLength;
	char read[SEQUENCE_LENGTH];
	char reference[SEQUENCE_LENGTH];
} MAF;

void TmpFileOpen(TmpFile*, char*, int);
void TmpFileClose(TmpFile*);
void TmpFileInitialize(TmpFile*);

void MAFPrint(FILE*, MAF*);
int MAFRead(FILE*, MAF*);
void MAFMergeSort(MAF*, int, int, int);
int MAFCompare(MAF*, MAF*, int);
void MAFCopy(MAF*, MAF*);
void MAFInitialize(MAF*);
void MAFPrintToBedAndWig(MAF*, int, int, int64_t, int64_t, FILE*, FILE*);

int SplitIntoTmpFilesByChr(char*, TmpFile**, char*, int, int);
void SplitMAFAndPrint(FILE*, FILE*, TmpFile*, char*, int);

#endif
