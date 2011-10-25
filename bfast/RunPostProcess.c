#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <zlib.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include "BLibDefinitions.h"
#include "BLib.h"
#include "BError.h"
#include "AlignedRead.h"
#include "AlignedEnd.h"
#include "AlignedReadConvert.h"
#include "ScoringMatrix.h"
#include "AlignMatrix.h"
#include "Align.h"
#include "RunPostProcess.h"

#define MAXIMUM_RESCUE_MAPQ 30
	
static inline int getStrandDiff(char strandA, char strandB, int strandedness)
{
  int strandDiff = 0;
  if(strandA == strandB) {
      strandDiff = 0;
  }
  else {
      strandDiff = 1;
  }
  if(1 == strandedness) {
      strandDiff = 1 - strandDiff;
  }
  return strandDiff;
}

//static inline int getTemplateLength(int32_t positionA, int32_t positionB, char strand, int32_t positioning)

static inline int getPositionDiff(int32_t positionA, int32_t positionB, char strandA, char strandB, int32_t positioning, int32_t strandedness)
{
  int positionDiff = positionB - positionA;
  if(0 == positioning || 2 == positioning) { // one before two
      if(FORWARD == strandA) {
          positionDiff = positionB - positionA;
      }
      else {
          positionDiff = positionA - positionB;
      }
  }
  else if(1 == positioning) { // two before one
      if(FORWARD == strandA) {
          positionDiff = positionA - positionB;
      }
      else {
          positionDiff = positionB - positionA;
      }
  }
  /*
  fprintf(stderr, "positionA=%d positionB=%d strandA=%c strandB=%c positionDiff=%d positioning=%d strandedness=%d\n",
          positionA,
          positionB,
          strandA,
          strandB,
          positionDiff,
          positioning,
          strandedness);
  */
  return positionDiff;
}

static inline int isDiscordantPair(AlignedEntry *entryOne,
		AlignedEntry *entryTwo,
                int strandedness,
                int positioning,
		PEDBins *b)
{
        int positionDiff;
        int strandDiff;
        int32_t low, high;

        strandDiff = getStrandDiff(entryOne->strand, entryTwo->strand, strandedness);
        positionDiff = getPositionDiff(entryOne->position, entryTwo->position, entryOne->strand, entryTwo->strand, positioning, strandedness);
        low = b->avg - (INSERT_MAX_STD * b->std); 
        high = b->avg + (INSERT_MAX_STD * b->std); 

        if(0 != strandDiff || positionDiff < 0
           || positionDiff < low
           || high < positionDiff) {
            return 1;
        }
        else {
            return 0;
        }
}

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
		FILE *fpOut)
{
	char *FnName="ReadInputFilterAndOutput";
	gzFile fp=NULL;
	int32_t i, j;
	int32_t numUnmapped=0, numReported=0;
	gzFile fpReportedGZ=NULL;
	FILE *fpReported=NULL;
	int32_t *mappedEndCounts=NULL;
	int32_t mappedEndCountsNumEnds=-1;
	int8_t *foundTypes=NULL;
	char *readGroupString=NULL;
	pthread_t *threads=NULL;
	int errCode;
	void *status=NULL;
	PostProcessThreadData *data=NULL;
	ScoringMatrix sm;
	int32_t matchScore ,mismatchScore, numRead, queueIndex;
	int32_t numReadsProcessed = 0;
	AlignedRead *alignQueue=NULL;
	int32_t alignQueueLength = 0;
	int32_t **numEntries=NULL;
	int32_t *numEntriesN=NULL;
	PEDBins bins;

	srand48(1); // to get the same behavior

	/* Read in scoring matrix */
	ScoringMatrixInitialize(&sm);
	if(NULL != scoringMatrixFileName) {
		ScoringMatrixRead(scoringMatrixFileName, &sm, space);
	}
	/* Calculate mismatch score */
	/* Assumes all match scores are the same and all substitution scores are the same */
	if(space == NTSpace) {
                matchScore = sm.ntMatch;
		mismatchScore = sm.ntMatch - sm.ntMismatch;
	}
	else {
                matchScore = sm.colorMatch + sm.ntMatch;
		mismatchScore = sm.colorMatch - sm.colorMismatch;
	}

	if(NULL != readGroup) {
		readGroupString=ParseReadGroup(readGroup);
	}

	if(NULL == inputFileName) {
		PrintError(FnName, "inputFileName", "Pairing from stdin currently not supported", Exit, OutOfRange);
	}

	/* Open the input file */
	if(NULL == inputFileName) {
		if(!(fp=gzdopen(fileno(stdin), "rb"))) {
			PrintError(FnName, "stdin", "Could not open stdin for reading", Exit, OpenFileError);
		}
	}
	else {
		if(!(fp=gzopen(inputFileName, "rb"))) {
			PrintError(FnName, inputFileName, "Could not open inputFileName for reading", Exit, OpenFileError);
		}
	}

	/* Open output files, if necessary */
	if(BAF == outputFormat) {
		if(!(fpReportedGZ=gzdopen(fileno(fpOut), "wb"))) {
			PrintError(FnName, "stdout", "Could not open stdout for writing", Exit, OpenFileError);
		}
	}
	else {
		if(!(fpReported=fdopen(fileno(fpOut), "wb"))) {
			PrintError(FnName, "stdout", "Could not open stdout for writing", Exit, OpenFileError);
		}
	}

	AlignedReadConvertPrintHeader(fpReported, rg, outputFormat, readGroup);

	/* Allocate memory for threads */
	threads=malloc(sizeof(pthread_t)*numThreads);
	if(NULL==threads) {
		PrintError(FnName, "threads", "Could not allocate memory", Exit, MallocMemory);
	}
	/* Allocate memory to pass data to threads */
	data=malloc(sizeof(PostProcessThreadData)*numThreads);
	if(NULL==data) {
		PrintError(FnName, "data", "Could not allocate memory", Exit, MallocMemory);
	}

	alignQueueLength=queueLength;
	alignQueue=malloc(sizeof(AlignedRead)*alignQueueLength);
	if(NULL == alignQueue) {
		PrintError(FnName, "alignQueue", "Could not allocate memory", Exit, MallocMemory);
	}
	numEntries=malloc(sizeof(int32_t*)*alignQueueLength);
	if(NULL == numEntries) {
		PrintError(FnName, "numEntries", "Could not allocate memory", Exit, MallocMemory);
	}
	numEntriesN=malloc(sizeof(int32_t)*alignQueueLength);
	if(NULL == numEntriesN) {
		PrintError(FnName, "numEntriesN", "Could not allocate memory", Exit, MallocMemory);
	}
	foundTypes=malloc(sizeof(int8_t)*alignQueueLength);
	if(NULL == foundTypes) {
		PrintError(FnName, "foundTypes", "Could not allocate memory", Exit, MallocMemory);
	}

	// Initialize
	for(i=0;i<alignQueueLength;i++) {
		numEntries[i] = NULL;
		numEntriesN[i] = 0;
	}

	/* Go through each read */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Postprocessing...\n");
	}
	numRead = 0;
        PEDBinsInitialize(&bins, insertSizeSpecified, insertSizeAvg, insertSizeStdDev);
	while(0 != (numRead = GetAlignedReads(fp, alignQueue, alignQueueLength))) {

		/* Get the PEDBins if necessary */
                if(0 == unpaired) {
	      	  GetPEDBins(alignQueue, numRead, strandedness, positioning, &bins);
                }

		// Store the original # of entries for SAM output
		for(i=0;i<numRead;i++) {
			numEntries[i] = NULL;
			numEntriesN[i] = 0;
			foundTypes[i] = NoneFound;
		}

		/* Initialize thread data */
		for(i=0;i<numThreads;i++) {
			data[i].bins = (0 == unpaired) ? &bins : NULL;
			data[i].rg = rg;
			data[i].sm = &sm;
			data[i].algorithm = algorithm;
			data[i].strandedness = strandedness;
			data[i].positioning = positioning;
			data[i].avgMismatchQuality = avgMismatchQuality;
			data[i].randomBest = randomBest;
			data[i].matchScore = matchScore;
			data[i].mismatchScore = mismatchScore;
			data[i].minimumMappingQuality = minimumMappingQuality; 
			data[i].minimumNormalizedScore = minimumNormalizedScore; 
			data[i].alignQueue =  alignQueue;
			data[i].queueLength = numRead;
			data[i].foundTypes = foundTypes;
			data[i].numEntriesN = numEntriesN;
			data[i].numEntries = numEntries;
			data[i].threadID = i;
			data[i].numThreads = numThreads;
		}

		/* Open threads */
		for(i=0;i<numThreads;i++) {
			/* Start thread */
			errCode = pthread_create(&threads[i], /* thread struct */
					NULL, /* default thread attributes */
					ReadInputFilterAndOutputThread, /* start routine */
					&data[i]); /* data to routine */
			if(0!=errCode) {
				PrintError(FnName, "pthread_create: errCode", "Could not start thread", Exit, ThreadError);
			}
		}
		/* Wait for threads to return */
		for(i=0;i<numThreads;i++) {
			/* Wait for the given thread to return */
			errCode = pthread_join(threads[i],
					&status);
			/* Check the return code of the thread */
			if(0!=errCode) {
				PrintError(FnName, "pthread_join: errCode", "Thread returned an error", Exit, ThreadError);
			}
		}

		/* Print to Output file */
		for(queueIndex=0;queueIndex<numRead;queueIndex++) {
			int32_t numEnds=0;
			if(NoneFound == foundTypes[queueIndex]) {
				/* Free the alignments for output */
				for(i=0;i<alignQueue[queueIndex].numEnds;i++) {
					for(j=0;j<alignQueue[queueIndex].ends[i].numEntries;j++) {
						AlignedEntryFree(&alignQueue[queueIndex].ends[i].entries[j]);
					}
					alignQueue[queueIndex].ends[i].numEntries=0;
				}
			}
			else {
				numReported++;
			}

			// Get the # of ends
			numEnds = 0;
			for(i=0;i<alignQueue[queueIndex].numEnds;i++) {
				if(0 < alignQueue[queueIndex].ends[i].numEntries) {
					numEnds++;
				}
			}
			if(0 == numEnds) {
				numUnmapped++;
			}
			// Clean up
			if(mappedEndCountsNumEnds < numEnds) {
				// Reallocate
				mappedEndCounts = realloc(mappedEndCounts, sizeof(int32_t)*(1+numEnds));
				if(NULL == mappedEndCounts) {
					PrintError(FnName, "mappedEndCounts", "Could not reallocate memory", Exit, ReallocMemory);
				}
				// Initialize
				for(i=1+mappedEndCountsNumEnds;i<=numEnds;i++) {
					mappedEndCounts[i] = 0;
				}
				mappedEndCountsNumEnds = numEnds;
			}
			mappedEndCounts[numEnds]++;

			// Proper pair ? 
			int properPair = 0;
			if(2 == alignQueue[queueIndex].numEnds && 0 == unpaired) {
                            if(1 == alignQueue[queueIndex].ends[0].numEntries && 1 == alignQueue[queueIndex].ends[1].numEntries) {
                                properPair = 1 - isDiscordantPair(&alignQueue[queueIndex].ends[0].entries[0],
                                                                  &alignQueue[queueIndex].ends[1].entries[0],
                                                                  strandedness,
                                                                  positioning,
                                                                  &bins);
                            }
			}
			AlignedReadConvertPrintOutputFormat(&alignQueue[queueIndex], rg, fpReported, fpReportedGZ, (NULL == outputID) ? "" : outputID, readGroupString, algorithm, numEntries[queueIndex], outputFormat, properPair, baseQualityType, BinaryOutput);

			/* Free memory */
			AlignedReadFree(&alignQueue[queueIndex]);
		}

		// Free
		for(i=0;i<numRead;i++) {
			free(numEntries[i]);
			numEntries[i] = NULL;
			numEntriesN[i] = 0;
		}

		numReadsProcessed += numRead;
		if(VERBOSE >= 0) {
			fprintf(stderr, "Reads processed: %d\n%s", numReadsProcessed, BREAK_LINE);
		}
	}
        /* Free */
        PEDBinsFree(&bins);
	if(0 <= VERBOSE) {
		fprintf(stderr, "Alignment complete.\n");
	}


	/* Close output files, if necessary */
	if(BAF == outputFormat) {
		gzclose(fpReportedGZ);
	}
	else {
		fclose(fpReported);
	}
	/* Close the input file */
	gzclose(fp);

	if(VERBOSE>=0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Found %10lld reads with no ends mapped.\n", 
				(long long int)numUnmapped);
		if(!(mappedEndCountsNumEnds < 1 || numUnmapped == mappedEndCounts[0])) {
			fprintf(stderr, "%d < 1 || %d == %d\n",
					mappedEndCountsNumEnds,
					(int)numUnmapped,
					mappedEndCounts[0]);
		}
		assert(mappedEndCountsNumEnds < 1 || numUnmapped == mappedEndCounts[0]);
		for(i=1;i<=mappedEndCountsNumEnds;i++) {
			if(1 == i) fprintf(stderr, "Found %10d reads with %2d end mapped.\n", mappedEndCounts[i], i);
			else fprintf(stderr, "Found %10d reads with %2d ends mapped.\n", mappedEndCounts[i], i);
		}
		fprintf(stderr, "Found %10lld reads with at least one end mapping.\n",
				(long long int)numReported);
		fprintf(stderr, "%s", BREAK_LINE);
	}
	free(mappedEndCounts);
	free(readGroupString);
	free(threads);
	free(data);
	free(alignQueue);
	free(foundTypes);
	free(numEntries);
	free(numEntriesN);
}

void *ReadInputFilterAndOutputThread(void *arg)
{
	char *FnName="ReadInputFilterAndOutputThread";
	PostProcessThreadData *data = (PostProcessThreadData*)arg;
	PEDBins *bins = data->bins;
	RGBinary *rg = data->rg;
	ScoringMatrix *sm = data->sm;
	int algorithm = data->algorithm;
	int strandedness = data->strandedness;
	int positioning = data->positioning;
	int avgMismatchQuality = data->avgMismatchQuality;
	int randomBest = data->randomBest;
	int matchScore = data->matchScore;
	int mismatchScore = data->mismatchScore;
	int minimumMappingQuality = data->minimumMappingQuality;
	int minimumNormalizedScore = data->minimumNormalizedScore;
	AlignedRead *alignQueue = data->alignQueue;
	int queueLength = data->queueLength;
	int8_t *foundTypes = data->foundTypes;
	int32_t threadID = data->threadID;
	int32_t numThreads = data->numThreads;
	int32_t **numEntries = data->numEntries;
	int32_t *numEntriesN = data->numEntriesN;
	int32_t j;
	int32_t queueIndex=0;
	AlignMatrix matrix;
	AlignMatrixInitialize(&matrix); 

	for(queueIndex=threadID;queueIndex<queueLength;queueIndex+=numThreads) {

                if(numEntriesN[queueIndex] < alignQueue[queueIndex].numEnds) {
                        numEntriesN[queueIndex] = alignQueue[queueIndex].numEnds;
                        numEntries[queueIndex]=realloc(numEntries[queueIndex], sizeof(int32_t)*numEntriesN[queueIndex]);
                        if(NULL == numEntries[queueIndex]) {
                                PrintError(FnName, "numEntries[queueIndex]", "Could not reallocate memory", Exit, ReallocMemory);
                        }
                }
                for(j=0;j<alignQueue[queueIndex].numEnds;j++) {
                        numEntries[queueIndex][j] = alignQueue[queueIndex].ends[j].numEntries;
                }

                /* Filter */
                foundTypes[queueIndex] = FilterAlignedRead(&alignQueue[queueIndex],
                                rg,
                                &matrix,
                                sm,
                                algorithm,
                                strandedness,
                                positioning,
                                avgMismatchQuality,
                                randomBest,
                                matchScore,
                                mismatchScore,
                                minimumMappingQuality,
                                minimumNormalizedScore,
                                bins);
	}

	// Free
	AlignMatrixFree(&matrix);

	return arg;
}

int32_t GetAlignedReads(gzFile fp, AlignedRead *alignQueue, int32_t maxToRead) 
{
	int32_t numRead=0;
	while(numRead < maxToRead) {
		AlignedReadInitialize(&alignQueue[numRead]);
		if(EOF == AlignedReadRead(&alignQueue[numRead], fp)) {
			break;
		}
		numRead++;
	}
	return numRead;
}

static int32_t getPairedScore(AlignedEntry *aOne,
                              AlignedEntry *aTwo,
                              int strandedness,
                              int positioning,
                              int avgMismatchQuality,
                              int matchScore,
                              int mismatchScore,
                              PEDBins *b)
{
  int32_t s = (int)(aOne->score + aTwo->score);
  int32_t strandDiff = 0;
  int32_t positionDiff = 0;
  double numStd = INSERT_MAX_STD;

  strandDiff = getStrandDiff(aOne->strand, aTwo->strand, strandedness);
  positionDiff = getPositionDiff(aOne->position, aTwo->position, aOne->strand, aTwo->strand, positioning, strandedness);

  // strandedness
  if(0 < strandDiff) {
      numStd = INSERT_MAX_STD;
  }
  else { 
      // positioning
      if(aOne->contig == aTwo->contig) {
          numStd = fabs(positionDiff - b->avg) / b->std;
          if(INSERT_MAX_STD < numStd) {
              numStd = INSERT_MAX_STD;
          }
      }
  }
  s -= (int)(mismatchScore * -1.0 * log10( erfc(M_SQRT1_2 * numStd)) + 0.499);

  return s;
}

void DoPairing(AlignedRead *tmpA, 
               int algorithm,
               int positioning, 
               int strandedness,
               int randomBest,
               int avgMismatchQuality,
               int matchScore,
               int mismatchScore,
               PEDBins *b)
{
  char *FnName="DoPairing";
  int32_t i, j;
      int32_t bestScore, bestNum;
      int32_t penultimateScore, penultimateNum;
      int32_t *bestScoreIndexesI = NULL;
      int32_t *bestScoreIndexesJ = NULL;
      int32_t bestScoreIndexesMem = 16;

      bestScoreIndexesI = malloc(sizeof(int32_t) * bestScoreIndexesMem);
      if(NULL == bestScoreIndexesI) {
          PrintError(FnName, "bestScoreIndexesI", "Could not allocate memory", Exit, MallocMemory);
      }
      bestScoreIndexesJ = malloc(sizeof(int32_t) * bestScoreIndexesMem);
      if(NULL == bestScoreIndexesJ) {
          PrintError(FnName, "bestScoreIndexesJ", "Could not allocate memory", Exit, MallocMemory);
      }

      bestScore = INT_MIN;
      bestNum = 0;
      penultimateScore = INT_MIN;
      penultimateNum = 0;

      for(i=0;i<tmpA->ends[0].numEntries;i++) {
          for(j=0;j<tmpA->ends[1].numEntries;j++) {
              int32_t score = getPairedScore(&tmpA->ends[0].entries[i], 
                                             &tmpA->ends[1].entries[j],
                                             positioning,
                                             strandedness,
                                             avgMismatchQuality,
                                             matchScore,
                                             mismatchScore,
                                             b);
              if(bestScore < score) {
                  penultimateScore = bestScore;
                  penultimateNum = bestNum;
                  bestScore = score;
                  bestNum = 1;
                  bestScoreIndexesI[bestNum-1] = i;
                  bestScoreIndexesJ[bestNum-1] = j;
              }
              else if(bestScore == score) {
                  bestNum++;
                  if(bestScoreIndexesMem < bestNum) {
                      bestScoreIndexesMem = bestNum;
                      bestScoreIndexesI = realloc(bestScoreIndexesI, sizeof(int32_t) * bestScoreIndexesMem);
                      if(NULL == bestScoreIndexesI) {
                          PrintError(FnName, "bestScoreIndexesI", "Could not allocate memory", Exit, ReallocMemory);
                      }
                      bestScoreIndexesJ = realloc(bestScoreIndexesJ, sizeof(int32_t) * bestScoreIndexesMem);
                      if(NULL == bestScoreIndexesJ) {
                          PrintError(FnName, "bestScoreIndexesJ", "Could not allocate memory", Exit, ReallocMemory);
                      }
                  }
                  bestScoreIndexesI[bestNum-1] = i;
                  bestScoreIndexesJ[bestNum-1] = j;
              }
              else if(penultimateScore < score) {
                  penultimateScore = score;
                  penultimateNum = 1;
              }
              else if(penultimateScore == score) {
                  penultimateNum++;
              }
          }
      }
      
      if(1 < bestNum) {
          if(1 == randomBest) {
              // pick a random one
              i = (int)(drand48() * bestNum);
              // copy over
              if(i != bestScoreIndexesI[i]) {
                  AlignedEntryCopy(&tmpA->ends[0].entries[0], &tmpA->ends[0].entries[bestScoreIndexesI[i]]);
              }
              if(i != bestScoreIndexesJ[i]) {
                  AlignedEntryCopy(&tmpA->ends[1].entries[0], &tmpA->ends[1].entries[bestScoreIndexesJ[i]]);
              }
              // reallocate
              AlignedEndReallocate(&tmpA->ends[0], 1);
              AlignedEndReallocate(&tmpA->ends[1], 1);
          }
          else {
              AlignedEndReallocate(&tmpA->ends[0], 0);
              AlignedEndReallocate(&tmpA->ends[1], 0);
          }
      }
      else {
          if(0 != bestScoreIndexesI[0]) {
              AlignedEntryCopy(&tmpA->ends[0].entries[0], &tmpA->ends[0].entries[bestScoreIndexesI[0]]);
          }
          if(0 != bestScoreIndexesJ[0]) {
              AlignedEntryCopy(&tmpA->ends[1].entries[0], &tmpA->ends[1].entries[bestScoreIndexesJ[0]]);
          }
          // reallocate
          AlignedEndReallocate(&tmpA->ends[0], 1);
          AlignedEndReallocate(&tmpA->ends[1], 1);
      }
          
      // update mapping quality
      if(1 == bestNum) {
          int32_t mapq = MAXIMUM_MAPPING_QUALITY;
          if(penultimateNum <= 0) {
              penultimateScore = GETMAX(tmpA->ends[i].entries[0].score, tmpA->ends[1].entries[0].score);
              penultimateNum = 1;
          }
          /*
             mapq = (bestScore - penultimateScore) * avgMismatchQuality / mismatchScore;
             if(0 == mapq) {
             mapq = 1;
             }
             else if(MAXIMUM_MAPPING_QUALITY < mapq) {
             mapq = MAXIMUM_MAPPING_QUALITY;
             }
             */
          double sf = 0.2;
          sf *= 250.0 / (matchScore * GETMAX(tmpA->ends[0].readLength, tmpA->ends[1].readLength)); // scale based on the best possible alignment score 
          sf *= (bestNum / (1.0 * penultimateNum)); // scale based on number of sub-optimal mappings
          sf *= (double)(bestScore - penultimateScore + 1); // scale based on distance to the sub-optimal mapping
          mapq = (int32_t)(sf + 0.99999);
          if(mapq > MAXIMUM_MAPPING_QUALITY) mapq = MAXIMUM_MAPPING_QUALITY;
          else if(mapq <= 0) mapq = 1;
          assert(0 < mapq);
          tmpA->ends[0].entries[0].mappingQuality = mapq;
          tmpA->ends[1].entries[0].mappingQuality = mapq;
          assert(0 < tmpA->ends[0].entries[0].mappingQuality);
          assert(0 < tmpA->ends[1].entries[0].mappingQuality);
      }
      
      // free
      free(bestScoreIndexesI);
      free(bestScoreIndexesJ);
}

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
		PEDBins *b)
{
	char *FnName="FilterAlignedRead";
	int foundType;
	int32_t *foundTypes=NULL;
	AlignedRead tmpA;
	int32_t i, j, k, ctr;
	int32_t best, bestIndex, numBest;

	AlignedReadInitialize(&tmpA);

	/* We should only modify "a" if it is going to be reported */ 
	/* Copy in case we do not find anything to report */
	AlignedReadCopy(&tmpA, a);

	foundType=NoneFound;
	foundTypes=malloc(sizeof(int32_t)*tmpA.numEnds);
	if(NULL == foundTypes) {
		PrintError(FnName, "foundTypes", "Could not allocate memory", Exit, MallocMemory);
	}
        
        if(NULL != b && 
           0 <= strandedness && 0 <= positioning && 2 == tmpA.numEnds 
           && BestScore == algorithm 
           && 1 <= tmpA.ends[0].numEntries && 1 <= tmpA.ends[1].numEntries) {
            // Do pairing
            DoPairing(&tmpA,
                      algorithm,
                      positioning, 
                      strandedness,
                      randomBest,
                      avgMismatchQuality,
                      matchScore,
                      mismatchScore,
                      b);

            // check if any were found
            for(i=0;i<tmpA.numEnds;i++) {
                if(1 <= tmpA.ends[i].numEntries) {
                    foundTypes[i]=Found;
                }
                else {
                    foundTypes[i]=NoneFound;
                }
            }
        }
        else {
		/* Pick alignment for each end individually (is this a good idea?) */
		for(i=0;i<tmpA.numEnds;i++) {
			/* Choose each end */
			switch(algorithm) {
				case NoFiltering:
				case AllNotFiltered:
					foundTypes[i] = (0<tmpA.ends[i].numEntries)?Found:NoneFound;
					break;
				case Unique:
					foundTypes[i]=(1==tmpA.ends[i].numEntries)?Found:NoneFound;
					break;
				case BestScore:
				case BestScoreAll:
					best = INT_MIN;
					bestIndex = -1;
					numBest = 0;
					for(j=0;j<tmpA.ends[i].numEntries;j++) {
						if(best < tmpA.ends[i].entries[j].score) {
							best = tmpA.ends[i].entries[j].score;
							bestIndex = j;
							numBest = 1;
						}
						else if(best == tmpA.ends[i].entries[j].score) {
							numBest++;
						}
					}
					// Copy all to the front
					ctr=0;
					for(j=0;j<tmpA.ends[i].numEntries;j++) {
						if(tmpA.ends[i].entries[j].score == best) {
							if(ctr != j) {
								AlignedEntryCopy(&tmpA.ends[i].entries[ctr], 
										&tmpA.ends[i].entries[j]);
							}
							ctr++;
						}
					}
					assert(ctr == numBest);
					AlignedEndReallocate(&tmpA.ends[i], numBest);
					// Random
					if(BestScore == algorithm) {
						if(1 < numBest && 1 == randomBest) {
							int32_t keep = (int)(drand48() * numBest);
							AlignedEntryCopy(&tmpA.ends[i].entries[0], &tmpA.ends[i].entries[keep]);
							AlignedEndReallocate(&tmpA.ends[i], 1);
							tmpA.ends[i].entries[0].mappingQuality = 0; // ambiguous
							numBest = 1;
						}
						if(1 == numBest) {
							foundTypes[i] = Found;
						}
						else {
							foundTypes[i] = NoneFound;
						}
					}
					else if(BestScoreAll == algorithm) {
						if( 1 <= numBest) {
							foundTypes[i] = Found;
						}
						else {
							foundTypes[i] = NoneFound;
						}
					}
					break;
				default:
					PrintError(FnName, "algorithm", "Could not understand algorithm", Exit, OutOfRange);
					break;
			}
			/* Free if not found */
			if(NoneFound == foundTypes[i]) {
				AlignedEndReallocate(&tmpA.ends[i],
						0);
			}
		}
	}

	if(INT_MIN < minimumMappingQuality ||
			INT_MIN < minimumNormalizedScore) {
		for(i=0;i<tmpA.numEnds;i++) {
			for(j=k=0;j<tmpA.ends[i].numEntries;j++) {
				if(minimumMappingQuality <= tmpA.ends[i].entries[j].mappingQuality &&
						minimumNormalizedScore <= tmpA.ends[i].entries[j].score/AlignedEntryGetReadLength(&tmpA.ends[i].entries[j])) {
					// Copy to the front
					if(j != k) {
						AlignedEntryCopy(&tmpA.ends[i].entries[k], &tmpA.ends[i].entries[j]);
					}
					k++;
				}
			}
			// Reallocate
			if(k < tmpA.ends[i].numEntries) {
				AlignedEndReallocate(&tmpA.ends[i], k);
			}
		}
	}

	// Found if one end is found
	foundType=NoneFound;
	for(i=0;NoneFound==foundType && i<tmpA.numEnds;i++) {
		if(Found == foundTypes[i]) {
			foundType=Found;
			break;
		}
	}

	/* copy back */
	if(NoneFound != foundType) {
		AlignedReadFree(a);
		AlignedReadCopy(a, &tmpA);
	}
	AlignedReadFree(&tmpA);
	free(foundTypes);


	return foundType;
}

int32_t GetPEDBins(AlignedRead *alignQueue,
		int queueLength,
                int strandedness,
                int positioning,
		PEDBins *b)
{
	char *FnName="GetPEDBins";
	int64_t counter, foundType;
	int32_t alignQueueLength;
	int32_t numRead, queueIndex;

	alignQueueLength=queueLength;

	/* Go through each read */
	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
		if (1 == b->doCalc) 
			fprintf(stderr, "Estimating paired end distance...\n");
		else
			fprintf(stderr, "Collecting paired end statistics...\n");
	}
	counter = numRead = 0;

	// TODO: make this multi-threaded
	for(queueIndex=0;queueIndex<queueLength;queueIndex++) {
		if(2 == alignQueue[queueIndex].numEnds) { // Only paired end data
			AlignedRead tmpA;

			// Filter base on best scoring
			AlignedReadInitialize(&tmpA);
			AlignedReadCopy(&tmpA, &alignQueue[queueIndex]);

			foundType=FilterAlignedRead(&tmpA,
					NULL,
					NULL,
					NULL,
					BestScore,
                                        -1,
                                        -1,
					INT_MIN,
					0,
                                        INT_MIN,
					INT_MIN,
					INT_MIN,
					INT_MIN,
					NULL);

			if(Found == foundType) {
                                int32_t toInsert = 0, positionDiff = 0;
				assert(2 == tmpA.numEnds);
				/* Must only have one alignment per end and on the same contig.
				 * There is a potential this will be inferred incorrectly under
				 * many scenarios.  Be careful! */
				if(1 == tmpA.ends[0].numEntries &&
						1 == tmpA.ends[1].numEntries &&
						tmpA.ends[0].entries[0].contig == tmpA.ends[1].entries[0].contig) {
                                    // Strands are OK
                                    int32_t strandDiff = getStrandDiff(tmpA.ends[0].entries[0].strand, tmpA.ends[1].entries[0].strand, strandedness);
                                    //fprintf(stderr, "strandDiff=%d\n", strandDiff);
                                    if(0 == strandDiff) {
                                        // Positions are OK
                                        //fprintf(stderr, "positionDiff=%d\n", positionDiff);
                                        positionDiff = getPositionDiff(tmpA.ends[0].entries[0].position, 
                                                                       tmpA.ends[1].entries[0].position, 
                                                                       tmpA.ends[0].entries[0].strand, 
                                                                       tmpA.ends[1].entries[0].strand, 
                                                                       positioning,
                                                                       strandedness);
                                        if(2 == positioning || 0 <= positionDiff) {
                                            toInsert = 1;
                                        }
                                    }
                                }
                                if(1 == toInsert) {
                                    PEDBinsInsert(b, positionDiff); 
				}
			}
			AlignedReadFree(&tmpA);
		}

		/* Increment counter */
		counter++;
	}

	if(1 == b->doCalc && b->numDistances < MIN_PEDBINS_SIZE) {
		fprintf(stderr, "Found only %d distances to infer the insert size distribution\n", b->numDistances);
		PrintError(FnName, "b->numDistances", "Not enough distances to infer insert size distribution", Warn, OutOfRange);
		PEDBinsFree(b);
		return 1;
	}

	if(VERBOSE>=0) {
		// Print Statistics
		PEDBinsPrintStatistics(b, stderr);
	}

	return 0;
}

void PEDBinsInitialize(PEDBins *b, int insertSizeSpecified, double insertSizeAvg, double insertSizeStdDev)
{
	int32_t i;
	b->minDistance = INT_MAX;
	b->maxDistance = INT_MIN;
	for(i=0;i<MAX_PEDBINS_DISTANCE - MIN_PEDBINS_DISTANCE+1;i++) {
		b->bins[i] = 0;
	}
	b->numDistances = 0;
	if (1 == insertSizeSpecified) {
		b->doCalc = 0;
		b->avg = insertSizeAvg;
		b->std = insertSizeStdDev;
	}
	else {
		b->doCalc = 1;
		b->avg = 0;
		b->std = 0;
	}
}

void PEDBinsFree(PEDBins *b) {
	PEDBinsInitialize(b, (1 - b->doCalc), b->avg, b->std);
}

void PEDBinsInsert(PEDBins *b,
		int32_t distance)
{
	if(distance < MIN_PEDBINS_DISTANCE ||
	   MAX_PEDBINS_DISTANCE < distance) {
		return;
	}

	if(distance < b->minDistance) {
		b->minDistance = distance;
		if(0 == b->numDistances) { // First one!
			b->minDistance = b->maxDistance = distance;
		}
	}
	else if(b->maxDistance < distance) {
		b->maxDistance = distance;
		if(0 == b->numDistances) { // First one!
			b->minDistance = distance;
		}
	}

	// Add to bin
	b->bins[distance - b->minDistance]++;
	b->numDistances++;

}

void PEDBinsPrintStatistics(PEDBins *b, FILE *fp)
{
	if (1 == b->doCalc) {
		// Mean, Range, and SD
		int32_t i;

		// Mean
		b->avg = 0.0;
		for(i=0;i<b->maxDistance-b->minDistance+1;i++) {
			b->avg += (b->minDistance + i)*b->bins[i];
		}
		b->avg /= b->numDistances;

		// SD
		b->std = 0.0;
		for(i=0;i<b->maxDistance-b->minDistance+1;i++) {
			b->std += b->bins[i]*((b->minDistance + i) - b->avg)*((b->minDistance + i) - b->avg);
		}
		b->std /= b->numDistances-1;
		b->std = sqrt(b->std);
	}

	if(0<=VERBOSE) {
		if (1 == b->doCalc)
			fprintf(stderr, "Used %d paired end distances to infer the insert size distribution.\n",
				b->numDistances);
		else
			fprintf(stderr, "Collected statistics for %d paired end distances.\n",
				b->numDistances);

		fprintf(stderr, "The paired end distance range was from %d to %d.\n",
			b->minDistance, b->maxDistance);

		if (1 == b->doCalc)
			fprintf(stderr, "The paired end distance mean and standard deviation were %.2lf and %.2lf.\n",
				b->avg, b->std);
		else
			fprintf(stderr, "The paired end distance mean and standard deviation were %.2lf and %.2lf. (Specified.)\n",
				b->avg, b->std);
	}
}
