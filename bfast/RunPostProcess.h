#ifndef RUNPOSTPROCESS_H_
#define RUNPOSTPROCESS_H_

#include "AlignedRead.h"


/* Paired End Distance Bins */
// This distance of the second end minus the first end
typedef struct {
	int32_t minDistance;
	int32_t maxDistance;
	int32_t bins[MAX_PEDBINS_DISTANCE - MIN_PEDBINS_DISTANCE + 1];
	int32_t numDistances;
	double std;
	double avg;
	int32_t inversionCount;
} PEDBins;

typedef struct {
	RGBinary *rg;
	PEDBins *bins;
	int algorithm;
	int unpaired;
	int avgMismatchQuality;
	int mismatchScore;
	int queueLength;
	int outputFormat;
	char *readGroupString;
	char *outputID;
	int32_t **mappedEndCounts;
	int32_t *mappedEndCountsNumEnds;
	int32_t *numReported, *numUnmapped;

	pthread_mutex_t *inputFP_mutex;
	int32_t *input_threadID;
	gzFile fpIn;

	FILE *fpReported;
	gzFile fpReportedGZ;
	gzFile fpUnmapped;
	pthread_mutex_t *outputFP_mutex;
	int32_t threadID;
} PostProcessThreadData;

void ReadInputFilterAndOutput(RGBinary *rg,
		char *inputFileName,
		int algorithm,
		int space,
		int unpaired,
		int avgMismatchQuality,
		char *scoringMatrixFileName,
		int numThreads,
		int queueLength,
		int outputFormat,
		char *outputID,
		char *readGroup,
		char *unmappedFileName,
		FILE *fpOut);

void *ReadInputFilterAndOutputThread(void*);

int32_t GetPEDBins(char*, int, int, PEDBins*);

int32_t GetAlignedReads(gzFile, AlignedRead*, int32_t, pthread_mutex_t*);

int FilterAlignedRead(AlignedRead *a,
		int algorithm,
		int unpairedInfer,
		int avgMismatchQuality,
		int mismatchScore,
		PEDBins *b);

void PEDBinsInitialize(PEDBins*);
void PEDBinsFree(PEDBins*);
void PEDBinsInsert(PEDBins*, char, char, int32_t);
void PEDBinsPrintStatistics(PEDBins*, FILE*);

#endif
