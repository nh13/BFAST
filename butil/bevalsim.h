#ifndef BEVALSIM_H_
#define BEVALSIM_H_

typedef struct {
	/* Meta data */
	char strand;
	int contig;
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
	int aContig;
	int aPos;
	int aStrand;
} ReadType;

typedef struct {
	/* actual data */
	int numReads;
	int numCorrectlyAligned[5]; /* 0, 10, 100, 1000, 10000 */
	ReadType r;
} Stat;

typedef struct {
	Stat *stats;
	int numStats;
} Stats;

void ReadTypeInitialize(ReadType*);
void ReadTypeCopy(ReadType*, ReadType*);
void ReadTypePrint(ReadType*, FILE*);
int ReadTypeCompare(ReadType*,ReadType*);
void ReadTypeReadFromRAF(ReadType*, FILE*);

void StatInitialize(Stat*, ReadType*);
void StatPrint(Stat*, FILE*);
void StatAdd(Stat*, ReadType*);

void StatsInitialize(Stats*);
void StatsPrintHeader(FILE*);
void StatsPrint(Stats*, FILE*);
void StatsAdd(Stats*, ReadType*);
void StatsDelete(Stats*);

void Evaluate(char*, char*);

#endif
