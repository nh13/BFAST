#ifndef SIMREAD_H_
#define SIMREAD_H_

#define SIM_READ_MAX_GETRANDOM_FAILURES 100000
#define SIM_READ_MAX_MODIFY_FAILURES 100000

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
	int readLength;
	int pairedEnd;
	int pairedEndLength;
	int indelLength;
} SimRead;

enum {
	Default,                /* 0 */
	Insertion,              /* 1 */
	SNP,                    /* 2 */
	Error,                  /* 3 */
	InsertionAndSNP,        /* 4 */
	InsertionAndError,      /* 5 */
	SNPAndError,            /* 6 */
	InsertionSNPAndError    /* 7 */
};

void SimReadInitialize(SimRead*);
void SimReadDelete(SimRead*);
char *SimReadGetName(SimRead*);
void SimReadPrint(SimRead*,
		FILE*);

void SimReadGetRandom(RGBinary*, int64_t, SimRead*, int, int, int, int, int, int, int, int, int);
void SimReadGetRandomContigPos(RGBinary*, int64_t, int*, int*, char*);
int SimReadModify(RGBinary*, SimRead*, int, int, int, int, int, int);
int SimReadInsertIndel(RGBinary*, SimRead*, int, int);
void SimReadInsertMismatches(SimRead*, int, int, int);
void SimReadInsertMismatchesHelper(char*, int, int*, int, int, int);
void SimReadInsertColorErrors(SimRead*, int, int);


#endif
