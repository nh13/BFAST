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
	int32_t doCalc;
} PEDBins;

typedef struct {
	PEDBins *bins;
	RGBinary *rg;
	ScoringMatrix *sm;
	int algorithm;
        int strandedness;
        int positioning;
	int avgMismatchQuality;
	int randomBest;
	int matchScore;
	int mismatchScore;
	int minimumMappingQuality;
	int minimumNormalizedScore;
	int queueLength;
	int8_t *foundTypes;
	AlignedRead *alignQueue;
	int32_t **numEntries;
	int32_t *numEntriesN;
	int32_t numThreads;
	int32_t threadID;
} PostProcessThreadData;

void ReadInputFilterAndOutput(RGBinary *rg,
		char *inputFileName,
		int algorithm,
		int space,
                int strandedness,
                int positioning,
                int unpaired,
		int avgMismatchQuality,
		char *scoringMatrixFileName,
		int randomBest,
		int minimumMappingQuality,
		int minimumNormalizedScore,
		int insertSizeSpecified,
		double insertSizeAvg,
		double insertSizeStdDev,
		int numThreads,
		int queueLength,
		int outputFormat,
		char *outputID,
		char *readGroup,
                int baseQualityType,
		FILE *fpOut);

void *ReadInputFilterAndOutputThread(void*);

int32_t GetPEDBins(AlignedRead*, int, int, int, PEDBins*);

int32_t GetAlignedReads(gzFile, AlignedRead*, int32_t);

int FilterAlignedRead(AlignedRead *a,
		RGBinary *rg,
		AlignMatrix *matrix,
		ScoringMatrix *sm,
		int algorithm,
                int strandedness,
                int positioning,
		int avgMismatchQuality,
		int randomBest,
		int matchScore,
		int mismatchScore,
		int minimumMappingQuality,
		int minimumNormalizedScore,
		PEDBins *b);

void PEDBinsInitialize(PEDBins*, int, double, double);
void PEDBinsFree(PEDBins*);
void PEDBinsInsert(PEDBins*, int32_t);
void PEDBinsPrintStatistics(PEDBins*, FILE*);

#endif
