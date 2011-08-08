#ifndef RUNPOSTPROCESS_H_
#define RUNPOSTPROCESS_H_

#include "AlignedRead.h"
#include "AlignMatrix.h"

/* Paired End Distance Bins */
// This distance of the second end minus the first end
typedef struct {
	int32_t minDistance;
	int32_t maxDistance;
	int32_t bins[MAX_PEDBINS_DISTANCE - MIN_PEDBINS_DISTANCE + 1];
	int32_t numDistances;
	double std;
	double avg;
	double invRatio;
	int32_t inversionCount;
	int32_t doCalc;
} PEDBins;

typedef struct {
	PEDBins *bins;
	RGBinary *rg;
	ScoringMatrix *sm;
	int algorithm;
	int unpaired;
	int reversePaired;
	int avgMismatchQuality;
	int randomBest;
	int matchScore;
	int mismatchScore;
	int minimumMappingQuality;
	int minimumNormalizedScore;
	double pairingStandardDeviation;
	int gappedPairingRescue;
	int queueLength;
	int8_t *foundTypes;
	AlignedRead *alignQueue;
	int32_t *alignQueueThreadIDs;
	int32_t **numEntries;
	int32_t *numEntriesN;
	int32_t numThreads;
	int32_t threadID;
} PostProcessThreadData;

void ReadInputFilterAndOutput(RGBinary *rg,
		char *inputFileName,
		int algorithm,
		int space,
		int unpaired,
		int reversePaired,
		int avgMismatchQuality,
		char *scoringMatrixFileName,
		int randomBest,
		int minimumMappingQuality,
		int minimumNormalizedScore,
		double pairingStandardDeviation,
		int insertSizeSpecified,
		double insertSizeAvg,
		double insertSizeStdDev,
		int gappedPairingRescue,
		int numThreads,
		int queueLength,
		int outputFormat,
		char *outputID,
		char *readGroup,
		FILE *fpOut);

void *ReadInputFilterAndOutputThread(void*);

int32_t GetPEDBins(AlignedRead*, int, PEDBins*);

int32_t GetAlignedReads(gzFile, AlignedRead*, int32_t);

int FilterAlignedRead(AlignedRead *a,
		RGBinary *rg,
		AlignMatrix *matrix,
		ScoringMatrix *sm,
		int algorithm,
		int unpaired,
		int reversePaired,
		int avgMismatchQuality,
		int randomBest,
		int matchScore,
		int mismatchScore,
		int minimumMappingQuality,
		int minimumNormalizedScore,
		double pairingStandardDeviation,
		int gappedPairingRescue,
		PEDBins *b);

void PEDBinsInitialize(PEDBins*, int, double, double);
void PEDBinsFree(PEDBins*);
void PEDBinsInsert(PEDBins*, char, char, int32_t);
void PEDBinsPrintStatistics(PEDBins*, FILE*);

#endif
