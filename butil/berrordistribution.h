#ifndef BERRORDISTRIBUTION_H_
#define BERRORDISTRIBUTION_H_

enum {CountOnly, CountTotal, CountBoth};

/* TODO */
typedef struct {
	int length;
	int *countOne;
	int *totalOne;
	int *countTwo;
	int *totalTwo;
} Count;

/* Structure to hold data about the error distributions
 * */
typedef struct {
	int numReads;
	int space;
	int pairedEnd;
	Count by[3]; /* nt errors by position and color errors by position */
	Count across[3]; /* nt errors across reads and color errors across reads */
} Errors;

void ErrorDistribution(char*, int32_t, int32_t, Errors*);
void ErrorDistributionPrint(char*, int32_t, Errors*);

void ErrorsPrint(Errors*, FILE**, int);
void ErrorsUpdate(Errors*, AlignEntries *a, int32_t, int32_t);
void ErrorsUpdateHelper(Errors*, AlignEntry *a, int, int, int, int, int);
void ErrorsInitialize(Errors*);
void ErrorsFree(Errors*);

void CountPrint(Count*, FILE*, int);
void CountUpdate(Count*, int, int, int, int);
void CountInitialize(Count*);
void CountFree(Count*);

#endif
