#ifndef BEVALSIM_H_
#define BEVALSIM_H_

typedef struct {
	/* Meta data */
	char strand;
	int chr;
	int pos;
	int space;
	int pairedEnd;
	int pairedEndLength;
	int readLength;
	int whichReadVariants;
	int startIndel;
	int indelLength;
	int numSNPs;
	int numErrors;
	int deletionLength;
	int insertionLength;
	/* Actual data */
	int aChr;
	int aPos;
	int aStrand;
} ReadType;

typedef struct {
	/* actual data */
	int numReads;
	int numCorrectlyAligned[5]; /* 0, 10, 100, 1000, 10000 */
	/* meta data */
	int space;
	int pairedEnd;
	int pairedEndLength;
	int readLength;
	int indelLength;
	int numSNPs;
	int numErrors;
	int deletionLength;
	int insertionLength;
} Stats;

void ReadTypeInitialize(ReadType*);
void ReadTypeCopy(ReadType*,ReadType*);
int ReadTypeCompare(ReadType*,ReadType*);
void ReadTypeReadFromRAF(ReadType*, FILE*);

void StatsInitialize(Stats*, ReadType*);
void StatsPrintHeader(FILE*);
void StatsPrint(Stats*, FILE*);
void StatsAdd(Stats*, ReadType*);

void Evaluate(char*, char*);

#endif
