#ifndef BERRORDISTRIBUTION_H_
#define BERRORDISTRIBUTION_H_

enum {CountOnly, CountTotal, CountBoth};

/* TODO */
typedef struct {
	int numEnds;
	int *lengths;
	int **counts;
	int **totals;
} Count;

/* Structure to hold data about the error distributions
 * */
typedef struct {
	int numReads;
	int space;
	int32_t numEnds;
	Count by[3]; /* nt errors by position, color errors by position, and gaps by position */
	Count across[3]; /* nt errors across reads, color errors across reads, and gaps across reads */
} Errors;

void ErrorDistribution(char*, int32_t, int32_t, Errors*);
void ErrorDistributionPrint(char*, int32_t, Errors*);

void ErrorsPrint(Errors*, FILE**, int);
void ErrorsUpdate(Errors*, AlignedRead *a, int32_t, int32_t);
void ErrorsUpdateHelper(Errors*, AlignedEntry *a, int, int, int, int);
void ErrorsInitialize(Errors*);
void ErrorsFree(Errors*);

void CountPrint(Count*, FILE*); 
void CountUpdate(Count*, int, int, int);
void CountInitialize(Count*);
void CountFree(Count*);

#endif

