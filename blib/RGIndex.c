#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <string.h>
#include <pthread.h>
#include "BLibDefinitions.h"
#include "BError.h"
#include "BLib.h"
#include "RGBinary.h"
#include "RGMatch.h"
#include "RGIndex.h"

/* TODO */
void RGIndexCreate(RGIndex *index, 
		RGBinary *rg, 
		RGIndexLayout *rgLayout, 
		int32_t startChr,
		int32_t startPos,
		int32_t endChr,
		int32_t endPos,
		int32_t layoutIndex,
		int32_t numThreads,
		int32_t includeRepeats, 
		int32_t includeNs) 
{

	/* The sort will take care of most of the work.  We just want 
	 * to make sure that we only include sequence that agrees with
	 * includeRepeats and includeNs
	 * */

	char *FnName = "RGIndexCreate";
	int32_t curPos=-1;
	int32_t curStartPos=-1;
	int32_t curEndPos=-1;
	int32_t curChr=-1;
	int32_t insert, i, j, curTilePos;
	int32_t chrIndex = 0;

	/* Initialize the index */
	index->positions=NULL;
	index->chromosomes=NULL;
	index->length=0;
	index->totalLength=0;

	assert(startChr <= endChr);
	assert(startChr < endChr || (startChr == endChr && startPos <= endPos));

	/* Copy over index information from the rg */
	index->startChr = startChr;
	index->startPos = startPos;
	index->endChr = endChr;
	index->endPos = endPos;

	assert(index->startChr > rg->startChr || (index->startChr == rg->startChr && index->startPos >= rg->startPos));
	assert(index->endChr < rg->endChr || (index->endChr == rg->endChr && index->endPos <= rg->endPos));

	/* Copy over index information from the layout */
	index->totalLength = 0;
	index->hashWidth = rgLayout->hashLengths[layoutIndex];
	index->hashLength = pow(4, index->hashWidth);
	index->numTiles = rgLayout->numTiles[layoutIndex];
	/* Allocate memory and copy over tile lengths */
	index->tileLengths = malloc(sizeof(int32_t)*rgLayout->numTiles[layoutIndex]);
	if(NULL == index->tileLengths) {
		PrintError(FnName,
				"index->tileLengths",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	for(i=0;i<rgLayout->numTiles[layoutIndex];i++) {
		index->tileLengths[i] = rgLayout->tileLengths[layoutIndex][i];
		index->totalLength += rgLayout->tileLengths[layoutIndex][i];
	}
	/* Allocate memory and copy over gaps */
	if(rgLayout->numTiles[layoutIndex] > 1) {
		index->gaps = malloc(sizeof(int32_t)*(rgLayout->numTiles[layoutIndex]-1));
		if(NULL == index->gaps) {
			PrintError(FnName,
					"index->gaps",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		for(i=0;i<rgLayout->numTiles[layoutIndex]-1;i++) {
			index->gaps[i] = rgLayout->gaps[layoutIndex][i];
			index->totalLength += rgLayout->gaps[layoutIndex][i];
		}
	}

	assert(index->numTiles > 0);
	assert(index->tileLengths != NULL);

	if(VERBOSE >= 0) {
		fprintf(stderr, "Currently on [chr,pos]:\n");
		fprintf(stderr, "\r[%d,%d]",
				-1,
				-1);
	}
	/* For each chromosome */
	for(curChr=index->startChr;curChr <= index->endChr;curChr++) { 
		/* Update chr index */
		chrIndex = curChr - rg->startChr;
		assert(chrIndex >=0 && chrIndex < rg->numChrs);

		/* Update start and end bounds for this chromosome */
		if(curChr == startChr) {
			curStartPos = startPos;
		}
		else {
			curStartPos = rg->chromosomes[chrIndex].startPos;
		}
		if(curChr == endChr) {
			curEndPos = endPos;
		}
		else {
			curEndPos = rg->chromosomes[chrIndex].endPos;
		}
		/* For each position */
		for(curPos=curStartPos;curPos<=curEndPos;curPos++) {
			if(VERBOSE >= 0) {
				if(curPos%RGINDEX_ROTATE_NUM==0) {
					fprintf(stderr, "\r[%d,%d]",
							curChr,
							curPos);
				}
			}

			/* Forward direction Only */
			insert = 1;
			if(index->totalLength + curPos - 1 <= rg->chromosomes[chrIndex].endPos) {
				curTilePos=curPos;
				for(i=0;i<index->numTiles && insert==1;i++) { /* For each tile */
					for(j=0;insert==1 && j<index->tileLengths[i];j++) { /* For each position in the tile */
						/* Check that we are within bounds */
						assert(curTilePos <= rg->chromosomes[chrIndex].endPos);
						/* Check if there are any repeats in any of the tiles */
						if(0==includeRepeats && 1==RGBinaryIsRepeat(rg, curChr, curTilePos)) {
							insert = 0;
						}
						if(0==includeNs && 1==RGBinaryIsN(rg, curChr, curTilePos)) {
							insert = 0;
						}
						curTilePos++;
					}
					curTilePos--; /* incremented on exit, so decrement */
					if(i<index->numTiles-1) { /* Add gap */
						curTilePos += index->gaps[i];
					}
				}
			}
			else {
				insert = 0;
			}

			/* Insert if desired */
			if(1==insert) {
				/* Insert */
				index->length++;

				/* Allocate memory */
				index->positions = realloc(index->positions, sizeof(uint32_t)*index->length);
				if(NULL == index->positions) {
					PrintError("RGBinaryCreate",
							"index->positions",
							"Could not reallocate memory",
							Exit,
							ReallocMemory);
				}
				index->chromosomes = realloc(index->chromosomes, sizeof(uint8_t)*index->length);
				if(NULL == index->chromosomes) {
					PrintError("RGBinaryCreate",
							"index->chromosomes",
							"Could not reallocate memory",
							Exit,
							ReallocMemory);
				}

				/* Copy over */
				index->positions[index->length-1] = curPos;
				index->chromosomes[index->length-1] = curChr;
			}
		}
	}
	/* Decrement since they were incremented before exiting loops */
	curChr--;
	curPos--;
	if(VERBOSE >= 0) {
		fprintf(stderr, "\r[%d,%d]\n",
				curChr,
				curPos);
	}

	assert(index->length > 0);

	/* Sort the nodes in the index */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Sorting...\n");
	}
	RGIndexSortNodes(index, rg, numThreads);
	if(VERBOSE >= 0) {
		fprintf(stderr, "Sorted.\n");
	}

	/* Create hash table from the index */
	RGIndexCreateHash(index, rg);
}

/* TODO */
void RGIndexCreateHash(RGIndex *index, RGBinary *rg) 
{
	char *FnName = "RGIndexCreateHash";
	uint32_t start, end;
	uint32_t curHash, startHash;
	int64_t i;

	/* Allocate memory for the hash */
	index->starts = malloc(sizeof(uint32_t)*index->hashLength);
	if(NULL==index->starts) {
		PrintError(FnName,
				"index->starts",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	index->ends = malloc(sizeof(uint32_t)*index->hashLength);
	if(NULL==index->ends) {
		PrintError(FnName,
				"index->ends",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* initialize */
	for(i=0;i<index->hashLength;i++) {
		/* Can't use -1, so use UINT_MAX */
		index->starts[i] = UINT_MAX;
		index->ends[i] = UINT_MAX;
	}

	/* Go through index and update the hash */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Creating a hash. Out of %u, currently on:\n0",
				(uint32_t)index->length);
	}
	startHash = RGIndexGetHashIndex(index, rg, 0, 0);
	for(end=1, start=0;end < index->length;end++) {
		if(VERBOSE >= 0 && end%RGINDEX_ROTATE_NUM==0) {
			fprintf(stderr, "\r%u", end);
		}
		curHash = RGIndexGetHashIndex(index, rg, end, 0);
		assert(curHash >= startHash);
		if(curHash == startHash) {
			/* Do nothing */
		}
		else {
			/* Paranoia check */
			assert(startHash < curHash);
			assert(curHash != startHash);
			/* Check that it is within bounds */
			if(startHash < 0 || startHash >= index->hashLength) {
				fprintf(stderr, "%s: %lld\t%lld\n",
						FnName,
						(long long int)startHash,
						(long long int)index->hashLength);
			}
			assert(startHash >= 0 && startHash < index->hashLength);
			/* Check that it has not been already initialized */
			if(index->starts[startHash] != UINT_MAX) {
				fprintf(stderr, "%s: %lld\t%lld\n",
						FnName,
						(long long int)startHash,
						(long long int)index->starts[startHash]);
			}
			assert(index->starts[startHash] == UINT_MAX);
			if(index->ends[startHash] != UINT_MAX) {
				fprintf(stderr, "%s: %lld\t%lld\n",
						FnName,
						(long long int)startHash,
						(long long int)index->ends[startHash]);
			}
			assert(index->ends[startHash] == UINT_MAX);

			/* Store start and end */
			index->starts[startHash] = start;
			index->ends[startHash] = end-1;

			/* Check correctness */
			if(index->starts[startHash] > 0 && index->starts[startHash] != UINT_MAX) {
				assert( RGIndexCompareAt(index, rg, index->starts[startHash]-1, index->starts[startHash], 0) < 0);
			}
			if(index->ends[startHash] < index->length-1 && index->ends[startHash] != UINT_MAX) {
				assert( RGIndexCompareAt(index, rg, index->ends[startHash], index->ends[startHash]+1, 0) < 0);
			}

			/* Update start */
			start = end;
			startHash = curHash;
		}
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\rHash created.\n");
	}
	/* In the boundary case... */
	/* Store start and end */
	index->starts[startHash] = start;
	index->ends[startHash] = end-1;

	/* Check hash creation */
	for(i=0;i<index->hashLength;i++) {
		assert( (index->starts[i] == UINT_MAX && index->ends[i] == UINT_MAX) ||
				(index->starts[i] != UINT_MAX && index->ends[i] != UINT_MAX));
		if(index->starts[i] > 0 && index->starts[i] != UINT_MAX) {
			assert( RGIndexCompareAt(index, rg, index->starts[i]-1, index->starts[i], 0) < 0);
		}
		if(index->ends[i] < index->length-1 && index->ends[i] != UINT_MAX) {
			assert( RGIndexCompareAt(index, rg, index->ends[i], index->ends[i]+1, 0) < 0);
		}
	}
}

/* TODO */
void RGIndexSortNodes(RGIndex *index, RGBinary *rg, int32_t numThreads)
{
	int64_t i;
	ThreadRGIndexSortData *data=NULL;
	ThreadRGIndexSortData tempData;
	pthread_t *threads=NULL;
	int32_t errCode;
	void *status=NULL;
	int64_t *pivots;
	int64_t max, maxIndex;

	/* Only use threads if we want to divide and conquer */
	if(numThreads > 1) {

		/* Allocate memory for the thread arguments */
		data = malloc(sizeof(ThreadRGIndexSortData)*numThreads);
		if(NULL==data) {
			PrintError("RGIndexSortNodes",
					"data",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Allocate memory for the thread point32_ters */
		threads = malloc(sizeof(pthread_t)*numThreads);
		if(NULL==threads) {
			PrintError("RGIndexSortNodes",
					"threads",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}

		/* Should check that the number of threads is a power of 4 */
		assert(IsAPowerOfTwo(numThreads)==1);

		/* Allocate memory for the pivots */
		pivots = malloc(sizeof(int64_t)*(2*numThreads));
		if(NULL == pivots) {
			PrintError("RGIndexSortNodes",
					"pivots",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}

		for(i=0;i<2*numThreads;i++) {
			pivots[i] = -1;
		}

		/* Get the pivots and presort */
		fprintf(stderr, "\rInitializing...");
		RGIndexQuickSortNodesGetPivots(index,
				rg,
				0,
				index->length-1,
				pivots,
				1,
				numThreads);
		/* The last one must be less than index->length */
		pivots[2*numThreads-1]--;

		/* Check pivots */
		for(i=0;i<2*numThreads;i+=2) {
			if(!(pivots[i] >= 0 && pivots[i] < index->length)) {
				fprintf(stderr, "offender\ti:%d\tlength:%d\n",
						(int)i,
						(int)index->length);
				for(i=0;i<2*numThreads;i+=2) {
					fprintf(stderr, "i:%d\t%d\t%d\n",
							(int)i,
							(int)pivots[i],
							(int)pivots[i+1]);
				}
				exit(1);
			}
			assert(pivots[i] >= 0 && pivots[i] < index->length);
			assert(pivots[i+1] >= 0 && pivots[i+1] < index->length);
			assert(pivots[i] <= pivots[i+1]);
			if(i==0) {
				assert(pivots[i] == 0);
			}
			if(i+1==2*numThreads-1) {
				assert(pivots[i+1] == index->length-1);
			}
			if(i>1 && i%2==0) {
				assert(pivots[i] == pivots[i-1] + 1);
			}
			if(i>1) {
				assert(pivots[i] > pivots[i-1]);
			}
		}

		/* Initialize data */
		maxIndex=0;
		max = data[0].high-data[0].low;
		for(i=0;i<numThreads;i++) {
			data[i].index = index;
			data[i].rg = rg;
			data[i].threadID = i;
			data[i].low = pivots[2*i];
			data[i].high = pivots[2*i+1];
			data[i].showPercentComplete = 0;
			assert(data[i].low >= 0 && data[i].high < index->length);
			if(data[i].high - data[i].low >= max) {
				maxIndex = i;
			}
		}
		data[maxIndex].showPercentComplete = 1;

		/* Check that we split correctly */
		for(i=1;i<numThreads;i++) {
			assert(data[i-1].high < data[i].low);
		}

		/* Copy maxIndex to the front so that it is the first that we wait for... */
		tempData.low = data[0].low;
		tempData.high = data[0].high;
		tempData.threadID = data[0].threadID;
		tempData.showPercentComplete = data[0].showPercentComplete;
		data[0].low = data[maxIndex].low;
		data[0].high = data[maxIndex].high;
		data[0].threadID = data[maxIndex].threadID;
		data[0].showPercentComplete = data[maxIndex].showPercentComplete;
		data[maxIndex].low = tempData.low;
		data[maxIndex].high = tempData.high;
		data[maxIndex].threadID = tempData.threadID;
		data[maxIndex].showPercentComplete = tempData.showPercentComplete;

		/* Create threads */
		for(i=0;i<numThreads;i++) {
			/* Start thread */
			errCode = pthread_create(&threads[i], /* thread struct */
					NULL, /* default thread attributes */
					RGIndexQuickSortNodes, /* start routine */
					(void*)(&data[i])); /* data to routine */
			if(0!=errCode) {
				PrintError("RGIndexSortNodes",
						"pthread_create: errCode",
						"Could not start thread",
						Exit,
						ThreadError);
			}
		}

		/* Wait for the threads to finish */
		for(i=0;i<numThreads;i++) {
			/* Wait for the given thread to return */
			errCode = pthread_join(threads[i],
					&status);
			/* Check the return code of the thread */
			if(0!=errCode) {
				PrintError("RGIndexSortNodes",
						"pthread_join: errCode",
						"Thread returned an error",
						Exit,
						ThreadError);
			}
			if(i==maxIndex && VERBOSE >= 0) {
				fprintf(stderr, "\rWaiting for other threads to complete...");
			}
		}

		/* Free memory */
		free(threads);
		free(data);
		free(pivots);
	}
	else {
		if(VERBOSE >= 0) {
			fprintf(stderr, "\r0 percent complete");
		}
		RGIndexQuickSortNodesHelper(index,
				rg,
				0,
				index->length-1,
				1);
		if(VERBOSE >= 0) {
			fprintf(stderr, "\r100.00 percent complete\n");
		}
	}

	/* Test that we sorted correctly */
	/*
	   for(i=1;i<index->length;i++) {
	   assert(RGIndexCompareAt(index, rg, i-1, i, 0) <= 0);
	   }
	   */

}

/* TODO */
void *RGIndexQuickSortNodes(void *arg)
{
	/* thread arguments */
	ThreadRGIndexSortData *data = (ThreadRGIndexSortData*)(arg);

	/* Call helper */
	if(data->showPercentComplete == 1 && VERBOSE >= 0) {
		fprintf(stderr, "\r0 percent complete");
	}
	RGIndexQuickSortNodesHelper(data->index,
			data->rg,
			data->low,
			data->high,
			data->showPercentComplete);
	if(data->showPercentComplete == 1 && VERBOSE >= 0) {
		fprintf(stderr, "\r");
		fprintf(stderr, "thread %3.3lf percent complete", 100.0);
	}

	return arg;
}

/* TODO */
/* Call stack was getting too big, implement non-recursive sort */
void RGIndexQuickSortNodesHelper(RGIndex *index,
		RGBinary *rg,
		int64_t low,
		int64_t high,
		int32_t showPercentComplete)
{
	/* Stack for log n space and non-recursive implementation */
	int64_t *lowStack=NULL;
	int64_t *highStack=NULL;
	int64_t stackLength=0;

	/* Local Variables */
	int64_t i;
	int64_t pivot = 0;
	int64_t total = high-low+1;
	int64_t curLow, curHigh;
	double curPercent = 0.0;

	/* Initialize stack */
	stackLength=1;
	lowStack = malloc(sizeof(int64_t));
	if(NULL==lowStack) {
		PrintError("RGIndexQuickSortNodesHelper",
				"lowStack",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	highStack = malloc(sizeof(int64_t));
	if(NULL==highStack) {
		PrintError("RGIndexQuickSortNodesHelper",
				"highStack",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	lowStack[0] = low;
	highStack[0] = high;

	/* Continue while the stack is not empty */
	while(stackLength > 0) {
		/* Pop off the stack */
		curLow = lowStack[stackLength-1];
		curHigh = highStack[stackLength-1];
		stackLength--;

		/* Reallocate memory */
		lowStack = realloc(lowStack, sizeof(int64_t)*stackLength);
		if(NULL==lowStack && stackLength > 0) {
			PrintError("RGIndexQuickSortNodesHelper",
					"lowStack",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		highStack = realloc(highStack, sizeof(int64_t)*stackLength);
		if(NULL==highStack && stackLength > 0) {
			PrintError("RGIndexQuickSortNodesHelper",
					"highStack",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}

		/* Proceed if we are with range */
		if(curLow < curHigh && curLow >= low && curHigh <= high) {
			/* Choose a new pivot.  We could do this randomly (randomized quick sort)
			 * but lets just choose the median of the front, middle and end 
			 * */
			pivot = RGIndexGetPivot(index, rg, curLow, curHigh);
			assert(pivot >=0 && pivot<index->length);
			assert(curLow >=0 && curLow<index->length);
			assert(curHigh >=0 && curHigh<index->length);

			if(showPercentComplete == 1 && VERBOSE >= 0) {
				if(curPercent < 100.0*((double)(curLow - low))/total) {
					while(curPercent < 100.0*((double)(curLow - low))/total) {
						curPercent += RGINDEX_SORT_ROTATE_INC;
					}
					fprintf(stderr, "\r%3.3lf percent complete", 100.0*((double)(curLow - low))/total);
				}
			}

			/* Swap the node at pivot with the node at curHigh */
			RGIndexSwapAt(index, pivot, curHigh);

			/* Store where the pivot should be */
			pivot = curLow;

			for(i=curLow;i<curHigh;i++) {
				assert(pivot >= 0 && pivot <= curHigh); 
				assert(i>=0 && i <= curHigh);
				if(RGIndexCompareAt(index, rg, i, curHigh, 0) <= 0) {
					/* Swap node at i with node at the new pivot index */
					if(i!=pivot) {
						RGIndexSwapAt(index, pivot, i);
					}
					/* Increment the new pivot index */
					pivot++;
				}
			}

			/* Move pivot element to correct place */
			if(pivot != curHigh) {
				RGIndexSwapAt(index, pivot, curHigh);
			}

			/* Add to the stack */
			stackLength+=2;
			/* Reallocate memory */
			lowStack = realloc(lowStack, sizeof(int64_t)*stackLength);
			if(NULL==lowStack) {
				PrintError("RGIndexQuickSortNodesHelper",
						"lowStack",
						"Could not reallocate memory",
						Exit,
						ReallocMemory);
			}
			highStack = realloc(highStack, sizeof(int64_t)*stackLength);
			if(NULL==highStack) {
				PrintError("RGIndexQuickSortNodesHelper",
						"highStack",
						"Could not reallocate memory",
						Exit,
						ReallocMemory);
			}
			/* Add sub array below */
			lowStack[stackLength-1] = curLow;
			highStack[stackLength-1] = pivot-1;
			/* Add sub array above */
			lowStack[stackLength-2] = pivot+1;
			highStack[stackLength-2] = curHigh;
		}
	}
}

/* TODO */
void RGIndexQuickSortNodesGetPivots(RGIndex *index,
		RGBinary *rg,
		int64_t low,
		int64_t high,
		int64_t *pivots,
		int32_t lowPivot,
		int32_t highPivot)
{
	/* local variables */
	int64_t i;
	int64_t pivot = 0;

	if(low < high ) {
		/* Choose a new pivot.  We could do this randomly (randomized quick sort)
		 * but lets just choose the median of the front, middle and end 
		 * */
		pivot = RGIndexGetPivot(index, rg, low, high);
		assert(pivot >=0 && pivot<index->length);
		assert(low >=0 && low<index->length);
		assert(high >=0 && high<index->length);

		/* Partition the array.
		 * Basically, arrange everything from low to high so that everything that
		 * has value less than or equal to the pivot is on the low of the pivot, and
		 * everthing else (greater than) is on the high side. 
		 * */

		/* Swap the node at pivot with the node at high */
		RGIndexSwapAt(index, pivot, high);

		/* Store where the pivot should be */
		pivot = low;

		for(i=low;i<high;i++) {
			assert(pivot >= 0 && pivot <= high); 
			assert(i>=0 && i <= high);
			if(RGIndexCompareAt(index, rg, i, high, 0) <= 0) {
				/* Swap node at i with node at the new pivot index */
				if(i!=pivot) {
					RGIndexSwapAt(index, pivot, i);
				}
				/* Increment the new pivot index */
				pivot++;
			}
		}

		/* Move pivot element to correct place */
		if(pivot != high) {
			RGIndexSwapAt(index, pivot, high);
		}

		if(lowPivot >= highPivot) {
			assert(pivots!=NULL);
			/* Save pivots if necessary */
			pivots[2*lowPivot-2] = low;
			pivots[2*lowPivot-1] = high+1;
			return;
		}
		else {
			/* Call recursively */

			/* Sort below */
			assert(pivot-1 < high);
			RGIndexQuickSortNodesGetPivots(index, 
					rg, 
					low, 
					pivot-1,
					pivots, 
					lowPivot, 
					(lowPivot+highPivot)/2);
			/* Sort above */
			assert(pivot+1 > low);
			RGIndexQuickSortNodesGetPivots(index, 
					rg, 
					pivot+1, 
					high, 
					pivots, 
					(lowPivot+highPivot)/2 + 1, 
					highPivot);
		}
	}
	else {
		/* Special case when saving pivots */
		assert(pivots!=NULL);
		/* Save pivots if necessary */
		pivots[2*lowPivot-2] = low;
		pivots[2*lowPivot-1] = low;
		return;
	}
}

/* TODO */
/* Do not use this function.  It is too slow. 
 * */
void RGIndexShellSortNodesHelper(RGIndex *index,
		RGBinary *rg,
		int64_t low,
		int64_t high,
		int32_t showPercentComplete,
		double *curPercent,
		int64_t lowTotal,
		int64_t highTotal)
{
	/* local variables */
	int64_t i, j, increment, length;
	uint32_t tempPos;
	uint8_t tempChr;

	length = high - low + 1;

	increment = length/2;

	/*
	   fprintf(stderr, "low:%Ld\thigh:%Ld\tlength:%Ld\n",
	   low,
	   high,
	   length);
	   */

	while(increment > 0) {
		assert(increment < length);
		/* Perform insertion sort with jump size as increment */
		/*
		   if(showPercentComplete==1 && VERBOSE >= 0) {
		   fprintf(stderr, "\rincrement:%Ld\tlength:%Ld\n",
		   increment,
		   length);
		   }
		   */
		for(i=increment+low;i<=high;i+=increment) {
			/*
			   if(showPercentComplete==1 && VERBOSE >= 0 && (i-low-increment)%1000==0) {
			   fprintf(stderr, "\r%Ld",
			   i-low-increment);
			   }
			   */
			j=i;
			while( (j>=increment+low) && RGIndexCompareAt(index, rg, j-increment, j, 0) > 0) {
				/*
				   if(showPercentComplete==1 && VERBOSE >= 0) {
				   fprintf(stderr, "\rincrement:%Ld\ti:%Ld\tj:%Ld\tlow:%Ld\thigh:%Ld\n",
				   increment,
				   i,
				   j,
				   low,
				   high);
				   }
				   */
				/* Swap */
				tempPos = index->positions[j];
				tempChr = index->chromosomes[j];
				index->positions[j] = index->positions[j-increment];
				index->chromosomes[j] = index->chromosomes[j-increment];
				index->positions[j-increment] = tempPos;
				index->chromosomes[j-increment] = tempChr;
				j = j-increment;
			}
		}

		/* Update the increment */
		if(increment == 2) {
			/* Will perform insertion sort */
			increment = 1;
		}
		else { 
			increment = (int64_t) (increment/2.2);
		}
	}
}

/* TODO */
void RGIndexDelete(RGIndex *index)
{
	free(index->positions);
	index->positions = NULL;
	free(index->chromosomes);
	index->chromosomes = NULL;
	index->length = 0;
	index->totalLength = 0;

	index->hashLength=0;
	index->hashWidth=0;
	free(index->starts);
	index->starts=NULL;
	free(index->ends);
	index->ends=NULL;

	/* Only free if we have allocated gaps */
	free(index->gaps);
	index->gaps=NULL;
	index->numTiles=0;
	free(index->tileLengths);
	index->tileLengths=NULL;

	index->startChr=0;
	index->startPos=0;
	index->endChr=0;
	index->endPos=0;
}

/* TODO */
double RGIndexGetSize(RGIndex *index, int32_t outputSize) 
{
	double total=0.0;

	total += sizeof(RGIndex); /* memory used by the index base structure */
	total += sizeof(uint32_t)*index->length;/* memory used by positions */
	total += sizeof(uint8_t)*index->length;/* memory used by positions */
	total += sizeof(uint32_t)*index->hashLength;/* memory used by starts */
	total += sizeof(uint32_t)*index->hashLength;/* memory used by ends */
	total += sizeof(int32_t)*index->numTiles;/* memory used by tileLengths */
	total += sizeof(int32_t)*(index->numTiles-1);/* memory used by gaps */

	switch(outputSize) {
		case KILOBYTES:
			return (total/pow(2, 10));
			break;
		case MEGABYTES:
			return (total/pow(2, 20));
			break;
		case GIGABYTES:
			return (total/pow(2, 30));
			break;
		default:
			return total;
			break;
	}
}

/* TODO */
void RGIndexPrint(FILE *fp, RGIndex *index, int32_t binaryOutput)
{
	int64_t i;

	/* Print header */
	RGIndexPrintHeader(fp, index, binaryOutput);

	if(binaryOutput == 0) {

		/* Print the positions and chromosomes */
		for(i=0;i<index->length;i++) {
			fprintf(fp, "%u\t%u\n", 
					index->positions[i],
					(uint32_t)index->chromosomes[i]);
		}

		/* Print the starts and ends */
		for(i=0;i<index->hashLength;i++) {
			fprintf(fp, "%u\t%u\n",
					index->starts[i],
					index->ends[i]);
		}

		/* Print the tileLengths */
		for(i=0;i<index->numTiles;i++) {
			if(i>0) {
				fprintf(fp, "\t");
			}
			fprintf(fp, "%u",
					index->tileLengths[i]);
		}
		fprintf(fp, "\n");

		/* Print the gaps */
		for(i=0;i<index->numTiles-1;i++) {
			if(i>0) {
				fprintf(fp, "\t");
			}
			fprintf(fp, "%u",
					index->gaps[i]);
		}
		if(index->numTiles>1) {
			fprintf(fp, "\n");
		}
	}
	else {
		/* Print positions */
		if(fwrite(index->positions, sizeof(uint32_t), index->length, fp) != index->length || 
				/* Print chomosomes */
				fwrite(index->chromosomes, sizeof(uint8_t), index->length, fp) != index->length ||
				/* Print the starts */
				fwrite(index->starts, sizeof(uint32_t), index->hashLength, fp) != index->hashLength ||
				/* Print the ends */
				fwrite(index->ends, sizeof(uint32_t), index->hashLength, fp) != index->hashLength || 
				/* Print the tileLengths */
				fwrite(index->tileLengths, sizeof(int32_t), index->numTiles, fp) != index->numTiles ||
				/* Print the gaps */
				fwrite(index->gaps, sizeof(int32_t), index->numTiles-1, fp) != (index->numTiles-1)) {
			PrintError("RGIndexPrint",
					NULL,
					"Could not write index and hash",
					Exit,
					WriteFileError);
		}
	}
}

/* TODO */
void RGIndexRead(FILE *fp, RGIndex *index, int32_t binaryInput)
{
	int64_t i;
	uint32_t tempInt;

	/* Read in the header */
	RGIndexReadHeader(fp, index, binaryInput);

	assert(index->length > 0);

	/* Allocate memory for the positions */
	index->positions = malloc(sizeof(uint32_t)*index->length);
	if(NULL == index->positions) {
		PrintError("RGIndexRead",
				"index->positions",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the chromosomes */
	index->chromosomes = malloc(sizeof(uint8_t)*index->length);
	if(NULL == index->chromosomes) {
		PrintError("RGIndexRead",
				"index->chromosomes",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the starts */
	index->starts = malloc(sizeof(uint32_t)*index->hashLength);
	if(NULL == index->starts) {
		PrintError("RGIndexRead",
				"index->starts",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the ends */
	index->ends = malloc(sizeof(uint32_t)*index->hashLength);
	if(NULL == index->ends) {
		PrintError("RGIndexRead",
				"index->ends",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the tile lengths */
	index->tileLengths = malloc(sizeof(int32_t)*index->numTiles);
	if(NULL == index->tileLengths) {
		PrintError("RGIndexRead",
				"index->tileLengths",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the gaps */
	index->gaps = malloc(sizeof(int32_t)*(index->numTiles-1));
	if(NULL == index->gaps) {
		PrintError("RGIndexRead",
				"index->gaps",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	if(binaryInput == 0) {
		/* Read the positions and chromosomes */
		for(i=0;i<index->length;i++) {
			if(fscanf(fp, "%u\t%u\n",
						&index->positions[i],
						&tempInt)==EOF) {
				PrintError("RGIndexRead",
						NULL,
						"Could not read in chromosome and position",
						Exit,
						EndOfFile);
			}
			index->chromosomes[i] = (uint8_t)tempInt;
		}

		/* Read the positions and chromosomes */
		for(i=0;i<index->hashLength;i++) {
			if(fscanf(fp, "%u\t%u\n",
						&index->starts[i],
						&index->ends[i])==EOF) {
				PrintError("RGIndexRead",
						NULL,
						"Could not read in starts and ends",
						Exit,
						EndOfFile);
			}
		}

		/* Read the tileLengths */
		for(i=0;i<index->numTiles;i++) {
			if(fscanf(fp, "%d",
						&index->tileLengths[i])==EOF) {
				PrintError("RGIndexRead",
						NULL,
						"Could not read in tile length",
						Exit,
						EndOfFile);
			}
		}

		/* Read the gaps */
		for(i=0;i<index->numTiles-1;i++) {
			if(fscanf(fp, "%d",
						&index->gaps[i])==EOF) {
				PrintError("RGIndexRead",
						NULL,
						"Could not read in gap",
						Exit,
						EndOfFile);
			}
		}
	}
	else {
		/* Read in positions */
		if(fread(index->positions, sizeof(uint32_t), index->length, fp)!=index->length) {
			PrintError("RGIndexRead",
					NULL,
					"Could not read in positions",
					Exit,
					ReadFileError);
		}

		/* Read in the chromosomes */
		if(fread(index->chromosomes, sizeof(uint8_t), index->length, fp)!=index->length) {
			PrintError("RGIndexRead",
					NULL,
					"Could not read in chromosomes",
					Exit,
					ReadFileError);
		}

		/* Read in starts */
		if(fread(index->starts, sizeof(uint32_t), index->hashLength, fp)!=index->hashLength) {
			PrintError("RGIndexRead",
					NULL,
					"Could not read in starts",
					Exit,
					ReadFileError);
		}

		/* Read in ends */
		if(fread(index->ends, sizeof(uint32_t), index->hashLength, fp)!=index->hashLength) {
			PrintError("RGIndexRead",
					NULL,
					"Could not read in ends",
					Exit,
					ReadFileError);
		}

		/* Read the tileLengths */
		if(fread(index->tileLengths, sizeof(int32_t), index->numTiles, fp)!=index->numTiles) {
			PrintError("RGIndexRead",
					NULL,
					"Could not read in tile lengths",
					Exit,
					ReadFileError);
		}

		/* Read the gaps */
		if(fread(index->gaps, sizeof(int32_t), index->numTiles-1, fp)!= (index->numTiles-1)) {
			PrintError("RGIndexRead",
					NULL,
					"Could not read in gaps",
					Exit,
					ReadFileError);
		}
	}
}

/* TODO */
void RGIndexPrintHeader(FILE *fp, RGIndex *index, int32_t binaryOutput)
{
	if(binaryOutput == 0) {
		fprintf(fp, "%u\t%u\t%lld\t%d\t%d\t%d\t%d\t%d\t%d\n",
				index->length,
				index->hashWidth,
				(long long int)index->hashLength,
				index->totalLength,
				index->numTiles,
				index->startChr,
				index->startPos,
				index->endChr,
				index->endPos);
	}
	else {
		/* Print Header */
		if(fwrite(&index->length, sizeof(uint32_t), 1, fp) != 1 || 
				fwrite(&index->hashWidth, sizeof(uint32_t), 1, fp) != 1 ||
				fwrite(&index->hashLength, sizeof(int64_t), 1, fp) != 1 ||
				fwrite(&index->totalLength, sizeof(int32_t), 1, fp) != 1 ||
				fwrite(&index->numTiles, sizeof(int32_t), 1, fp) != 1 ||
				fwrite(&index->startChr, sizeof(int32_t), 1, fp) != 1 ||
				fwrite(&index->startPos, sizeof(int32_t), 1, fp) != 1 ||
				fwrite(&index->endChr, sizeof(int32_t), 1, fp) != 1 ||
				fwrite(&index->endPos, sizeof(int32_t), 1, fp) != 1) {
			PrintError("RGIndexPrintHeader",
					NULL,
					"Could not write header",
					Exit,
					WriteFileError);
		}
	}
}

/* TODO */
void RGIndexReadHeader(FILE *fp, RGIndex *index, int32_t binaryInput)
{
	/* Read in header */
	if(binaryInput == 0) {
		if(fscanf(fp, "%u %u %lld %d %d %d %d %d %d",
					&index->length,
					&index->hashWidth,
					(long long int *)&index->hashLength,
					&index->totalLength,
					&index->numTiles,
					&index->startChr,
					&index->startPos,
					&index->endChr,
					&index->endPos)==EOF) {
			PrintError("RGIndexReadHeader",
					NULL,
					"Could not read header",
					Exit,
					EndOfFile);
		}
	}
	else {
		if(fread(&index->length, sizeof(uint32_t), 1, fp)!=1
				|| fread(&index->hashWidth, sizeof(uint32_t), 1, fp)!=1
				|| fread(&index->hashLength, sizeof(int64_t), 1, fp)!=1
				|| fread(&index->totalLength, sizeof(int32_t), 1, fp)!=1
				|| fread(&index->numTiles, sizeof(int32_t), 1, fp)!=1 
				|| fread(&index->startChr, sizeof(int32_t), 1, fp)!=1
				|| fread(&index->startPos, sizeof(int32_t), 1, fp)!=1
				|| fread(&index->endChr, sizeof(int32_t), 1, fp)!=1
				|| fread(&index->endPos, sizeof(int32_t), 1, fp)!=1) {
			PrintError("RGIndexReadHeader",
					NULL,
					"Could not read header",
					Exit,
					ReadFileError);
		}
	}

	/* Error checking */
	assert(index->length > 0);
	assert(index->hashWidth > 0);
	assert(index->hashLength > 0);
	assert(index->totalLength > 0);
	assert(index->numTiles > 0);
	assert(index->startChr > 0);
	assert(index->startPos > 0);
	assert(index->endChr > 0);
	assert(index->endPos > 0);
}

/* TODO */
/* We will append the matches if matches have already been found */
void RGIndexGetMatches(RGIndex *index, RGBinary *rg, char *read, int32_t readLength, int8_t direction, int32_t offset, RGMatch *m, int32_t maxMatches)
{
	int64_t i;
	int64_t startIndex=-1;
	int64_t endIndex=-1;
	int64_t nodeIndex=-1;
	int64_t foundIndex=0;
	int64_t tmpIndex=-1;
	uint32_t hashIndex=0;

	/* Don't search if we have reached the maximum number of matches */
	if(m->maxReached==1) {
		return;
	}

	/* Get the hash index */
	/* The hope is that the hash will give better smaller bounds (if not
	 * zero bounds for the binary search on the index */
	hashIndex = RGIndexGetHashIndexFromRead(index, rg, read, readLength, 0);
	assert(hashIndex >= 0 && hashIndex < index->hashLength);

	if(index->starts[hashIndex] == UINT_MAX || 
			index->ends[hashIndex] == UINT_MAX) {
		/* Skip */
	}
	else {
		assert(index->starts[hashIndex] >=0 && index->starts[hashIndex] < index->length);
		assert(index->ends[hashIndex] >=0 && index->ends[hashIndex] < index->length);

		/* Search the index using the bounds from the hash */
		foundIndex=RGIndexGetIndex(index, 
				rg, 
				index->starts[hashIndex],  
				index->ends[hashIndex],
				read,
				&startIndex,
				&endIndex);

		if(foundIndex>0) {

			/* Check to see if we add the current number of matches, if we would go over the limit */
			if(m->numEntries + (endIndex - startIndex +1) > maxMatches) {
				RGMatchFree(m);
				m->maxReached=1;
				return;
			}

			/* Copy over */
			/* (Re)Allocate memory for the new matches */
			tmpIndex = m->numEntries;
			assert(endIndex >= startIndex);
			assert(startIndex >= 0 && startIndex < index->length);
			assert(endIndex >= 0 && endIndex < index->length);
			RGMatchReallocate(m,  m->numEntries + (endIndex-startIndex+1));

			/* Copy over */
			for(i=tmpIndex, nodeIndex = startIndex;i < m->numEntries && nodeIndex <= endIndex;i++, nodeIndex++) {
				m->positions[i] = index->positions[nodeIndex] - offset;
				m->chromosomes[i] = index->chromosomes[nodeIndex];
				m->strand[i] = direction;
			}
		}
	}
}

/* TODO */
int64_t RGIndexGetIndex(RGIndex *index,
		RGBinary *rg,
		int64_t low,
		int64_t high,
		char *read,
		int64_t *startIndex,
		int64_t *endIndex)
{
	int64_t mid=-1;
	int32_t cmp;
	int32_t cont = 1;
	int64_t tmpLow, tmpMid, tmpHigh;

	assert(low==0 || RGIndexCompareRead(index, rg, read, low-1, 0) > 0);
	assert(high==index->length-1 || RGIndexCompareRead(index, rg, read, high+1, 0) < 0); 

	while(low <= high && cont==1) {
		mid = (low+high)/2;
		cmp = RGIndexCompareRead(index, rg, read, mid, 0);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "low:%lld\tmid:%lld\thigh:%lld\tcmp:%d\n",
					(long long int)low,
					(long long int)mid,
					(long long int)high,
					cmp);
		}
		if(cmp == 0) {
			cont = 0;
		}
		else if(cmp < 0) {
			high = mid-1;
		}
		else {
			low = mid + 1;
		}
	}
	/* If we found an entry that matches, get the bounds (start and end indexes */
	if(cont == 0) {
		assert(low==0 || RGIndexCompareRead(index, rg, read, low-1, 0) > 0);
		assert(high==index->length-1 || RGIndexCompareRead(index, rg, read, high+1, 0) < 0); 
		assert(RGIndexCompareRead(index, rg, read, mid, 0) == 0);
		tmpLow = low;
		tmpMid = mid;
		tmpHigh = high;
		/*
		   fprintf(stderr, "Getting start and end:\t%lld\t%lld\t%lld\n",
		   low,
		   mid,
		   high);
		   */
		/* Get lower start Index */
		low = tmpLow;
		high = tmpMid;
		while(low < high) {
			mid = (low+high)/2;
			cmp = RGIndexCompareRead(index, rg, read, mid, 0);
			assert(cmp >= 0);
			/*
			   fprintf(stderr, "start:%lld\t%lld\t%lld\t%d\n",
			   low,
			   mid,
			   high,
			   cmp);
			   */
			if(cmp == 0) {
				high = mid;
			}
			else {
				/* mid is less than */
				low = mid+1;
			}
		}
		(*startIndex) = low;
		assert(low == high);
		assert(RGIndexCompareRead(index, rg, read, (*startIndex), 0)==0);
		assert((*startIndex) == 0 || RGIndexCompareRead(index, rg, read, (*startIndex)-1, 0)>0);
		/* Get upper start Index */
		low = tmpMid;
		high = tmpHigh;
		while(low < high) {
			mid = (low+high)/2+1;
			cmp = RGIndexCompareRead(index, rg, read, mid, 0);
			assert(cmp <= 0);
			/*
			   fprintf(stderr, "end:%lld\t%lld\t%lld\t%d\n",
			   low,
			   mid,
			   high,
			   cmp);
			   */
			if(cmp == 0) {
				low = mid;
			}
			else {
				/* mid is less than */
				high = mid-1;
			}
		}
		assert(low == high);
		/* adjust endIndex */
		(*endIndex) = low;
		assert(RGIndexCompareRead(index, rg, read, (*endIndex), 0)==0);
		assert((*endIndex) == index->length-1 || RGIndexCompareRead(index, rg, read, (*endIndex)+1, 0)<0);
		return 1;
	}
	else {
		return 0;
	}

}

/* TODO */
void RGIndexSwapAt(RGIndex *index, int64_t a, int64_t b)
{
	uint32_t tempChr, tempPos;

	tempChr = index->chromosomes[a];
	tempPos = index->positions[a];
	index->chromosomes[a] = index->chromosomes[b];
	index->positions[a] = index->positions[b];
	index->chromosomes[b] = tempChr;
	index->positions[b] = tempPos;
}

/* TODO */
int64_t RGIndexGetPivot(RGIndex *index, RGBinary *rg, int64_t low, int64_t high)
{
	int64_t pivot = (low+high)/2;
	int32_t cmp[3];
	cmp[0] = RGIndexCompareAt(index, rg, low, pivot, 0);
	cmp[1] = RGIndexCompareAt(index, rg, low, high, 0);
	cmp[2] = RGIndexCompareAt(index, rg, pivot, high, 0);

	if(cmp[0] <= 0) {
		/* low <= pivot */
		if(cmp[1] >= 0) {
			/* high <= low */
			/* so high <= low <= pivot */
			pivot = low;
		}
		else {
			/* low < high */
			if(cmp[2] <= 0) {
				/* pivot <= high */
				/* so low <= pivot <= high */
				/* choose pivot */
			}
			else {
				/* high < pivot */
				/* so low < high < pivot */
				pivot = high;
			}
		}
	}
	else {
		/* pivot < low */
		if(cmp[1] <= 0) {
			/* low <= high */
			/* so pivot < low <= high */
			pivot = low;
		}
		else {
			/* high < low */
			if(cmp[2] <= 0) {
				/* pivot <= high */
				/* so pivot <= high < low */
				pivot = high;
			}
			else {
				/* high < pivot */
				/* so high < pivot < low */
				/* choose pivot */
			}
		}
	}
	return pivot;
}

/* TODO */
int32_t RGIndexCompareAt(RGIndex *index,
		RGBinary *rg,
		int64_t a,
		int64_t b, 
		int debug)
{
	assert(a>=0 && a<index->length);
	assert(b>=0 && b<index->length);

	int32_t i, j;
	uint32_t aChr = index->chromosomes[a];
	uint32_t aPos = index->positions[a];
	uint32_t bChr = index->chromosomes[b];
	uint32_t bPos = index->positions[b];

	uint32_t aCurTilePos;
	uint32_t bCurTilePos;
	uint8_t aBase;
	uint8_t bBase;

	/* Compare base by base */
	aCurTilePos = aPos;
	bCurTilePos = bPos;

	for(i=0;i<index->numTiles;i++) { /* For each tile */
		for(j=0;j<index->tileLengths[i];j++) { /* For each position in the tile */
			aBase = ToLower(RGBinaryGetBase(rg,
						aChr,
						aCurTilePos));
			bBase = ToLower( RGBinaryGetBase(rg,
						bChr,
						bCurTilePos));

			if(debug > 0) {
				fprintf(stderr, "a[%d,%c]\tb[%d,%c]\n",
						aCurTilePos,
						aBase,
						bCurTilePos,
						bBase);
			}

			/* Compare */
			if(aBase < bBase) {
				return -1;
			}
			else if(aBase > bBase) {
				return 1;
			}
			/* Continue if the current bases are equal */

			/* Update */
			aCurTilePos++;
			bCurTilePos++;
		}
		if(i<index->numTiles-1) { /* Add gap */
			aCurTilePos += index->gaps[i];
			bCurTilePos += index->gaps[i];
		}
	}

	/* All bases were equal, return 0 */
	return 0;
}

/* TODO */
int32_t RGIndexCompareRead(RGIndex *index,
		RGBinary *rg,
		char *read,
		int64_t a,
		int debug)
{
	assert(a>=0 && a<index->length);

	int32_t i, j;
	int32_t curReadPos=0;
	int32_t readLength = strlen(read);
	uint32_t aChr = index->chromosomes[a];
	uint32_t aPos = index->positions[a];

	uint32_t aCurTilePos;
	uint8_t aBase;
	uint8_t readBase;

	/* Compare base by base */
	aCurTilePos = aPos;

	/*
	   fprintf(stderr, "%d\t%s", 
	   index->totalLength,
	   BREAK_LINE);
	   fprintf(stderr, "read[%d]:%s\n", 
	   (int)strlen(read),
	   read);
	   */
	for(i=0;i<index->numTiles && curReadPos < readLength;i++) { /* For each tile */
		for(j=0;j<index->tileLengths[i] && curReadPos < readLength;j++) { /* For each position in the tile */
			aBase=ToLower(RGBinaryGetBase(rg,
						aChr,
						aCurTilePos));
			readBase = ToLower(read[curReadPos]);

			if(debug > 0) {
				fprintf(stderr, "a[%d,%c]\tread[%d,%c]\n",
						aCurTilePos,
						aBase,
						curReadPos,
						readBase);
			}

			/* Compare */
			if(readBase < aBase) {
				return -1;
			}
			else if(readBase > aBase) {
				return 1;
			}
			/* Continue if the current bases are equal */

			/* Update */
			aCurTilePos++;
			curReadPos++;
		}
		if(i<index->numTiles-1) { /* Add gap */
			aCurTilePos += index->gaps[i];
			curReadPos += index->gaps[i];
		}
	}

	/* All bases were equal, return 0 */
	return 0;
}

/* TODO */
uint32_t RGIndexGetHashIndex(RGIndex *index,
		RGBinary *rg,
		uint32_t a,
		int debug)
{
	assert(a>=0 && a<index->length);

	int32_t i, j;
	uint32_t aChr = index->chromosomes[a];
	uint32_t aPos = index->positions[a];

	uint32_t aCurTilePos;
	uint8_t aBase;

	int32_t cur = index->hashWidth-1;
	uint32_t hashIndex = 0;

	/* Compare base by base */
	aCurTilePos = aPos;

	assert(ALPHABET_SIZE == 4);

	for(i=0;cur >= 0 && i < index->numTiles;i++) { /* For each tile */
		for(j=0;cur >= 0 && j < index->tileLengths[i];j++) { /* For each position in the tile */
			aBase = ToLower(RGBinaryGetBase(rg,
						aChr,
						aCurTilePos));
			if(debug > 0) {
				fprintf(stderr, "a[%d,%c]\tcur:%d\n",
						aCurTilePos,
						aBase,
						cur);
			}

			switch(aBase) {
				case 'a':
					/* Do nothing since a is zero base 4 */
					break;
				case 'c':
					hashIndex += pow(ALPHABET_SIZE, cur);
					break;
				case 'g':
					hashIndex += pow(ALPHABET_SIZE, cur)*2;
					break;
				case 't':
					hashIndex += pow(ALPHABET_SIZE, cur)*3;
					break;
				default:
					PrintError("RGIndexGetHashIndex",
							"aBase",
							"Could not understand base",
							Exit,
							OutOfRange);
					break;
			}
			/*
			   fprintf(stderr, "%c\t%d\t%lf\t%d\n",
			   aBase,
			   cur,
			   pow(ALPHABET_SIZE, cur),
			   hashIndex);
			   */

			cur--;

			/* Update */
			aCurTilePos++;
		}
		if(i<index->numTiles-1) { /* Add gap */
			aCurTilePos += index->gaps[i];
		}
	}

	/* All bases were equal, return 0 */
	return hashIndex;
}

/* TODO */
uint32_t RGIndexGetHashIndexFromRead(RGIndex *index,
		RGBinary *rg,
		char *read,
		int32_t readLength,
		int debug)
{
	int32_t i, j;
	int32_t curReadPos=0;
	uint8_t readBase;

	int32_t cur = index->hashWidth-1;
	uint32_t hashIndex = 0;

	for(i=0;cur >= 0 && i<index->numTiles && curReadPos < readLength;i++) { /* For each tile */
		for(j=0;cur >= 0 &&  j<index->tileLengths[i] && curReadPos < readLength;j++) { /* For each position in the tile */
			readBase = ToLower(read[curReadPos]);
			if(debug > 0) {
				fprintf(stderr, "read[%d,%c]\tcur:%d\tcurReadPos:%d\treadLength:%d\n",
						curReadPos,
						readBase,
						cur,
						curReadPos,
						readLength);
			}

			switch(readBase) {
				case 'a':
					/* Do nothing since a is zero base 4 */
					break;
				case 'c':
					hashIndex += pow(ALPHABET_SIZE, cur);
					break;
				case 'g':
					hashIndex += pow(ALPHABET_SIZE, cur)*2;
					break;
				case 't':
					hashIndex += pow(ALPHABET_SIZE, cur)*3;
					break;
				default:
					PrintError("RGIndexGetHashIndexFromRead",
							"aBase",
							"Could not understand base",
							Exit,
							OutOfRange);
					break;
			}
			cur--;

			curReadPos++;
		}
		if(i<index->numTiles-1) { /* Add gap */
			curReadPos += index->gaps[i];
		}
	}

	/* All bases were equal, return 0 */
	return hashIndex;
}

/* TODO */
/* Debug function */
void RGIndexPrintReadMasked(RGIndex *index, char *read, FILE *fp) 
{
	int curPos, curTile, curTilePos;
	for(curPos = 0, curTile = 0;curTile < index->numTiles;curTile++) {
		fprintf(fp, "[");
		for(curTilePos=0;curTilePos < index->tileLengths[curTile];curTilePos++) {
			fprintf(fp, "%c", read[curPos]);
			/* Update position */
			curPos++;
		}   
		fprintf(fp, "]\t");
		if(curTile < index->numTiles-1) {
			curPos += index->gaps[curTile];
		}   
	}   
	fprintf(fp, "\n");
}
