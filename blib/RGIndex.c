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
#include "RGIndex.h"

/* TODO */
void RGIndexCreate(RGIndex *index, RGBinary *rg, int32_t includeRepeats, int32_t includeNs) {

	/* The sort will take care of most of the work.  We just want 
	 * to make sure that we only include sequence that agrees with
	 * includeRepeats and includeNs
	 * */

	int32_t curPos=-1;
	int32_t curChr=-1;
	int32_t insert, i, j, curTilePos;
	int32_t chrIndex = 0;

	index->positions=NULL;
	index->chromosomes=NULL;
	index->length=0;
	assert(index->numTiles > 0);
	assert(index->tileLengths != NULL);

	if(VERBOSE >= 0) {
		fprintf(stderr, "Currently on [chr,pos]:\n");
		fprintf(stderr, "\r[%d,%d]",
				-1,
				-1);
	}
	/* For each chromosome */
	for(curChr=rg->startChr, chrIndex=0;curChr <= rg->endChr;curChr++, chrIndex++) { 

		/* For each position */
		for(curPos=rg->chromosomes[chrIndex].startPos;curPos<=rg->chromosomes[chrIndex].endPos;curPos++) {
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
}

/* TODO */
void RGIndexCleanUpIndex(RGIndex *index, RGBinary *rg, int32_t numThreads) 
{
	if(index->length > 0) {
		/* Sort the nodes in the index */
		fprintf(stderr, "Sorting...\n");
		RGIndexSortNodes(index, rg, numThreads);
		fprintf(stderr, "Sorted.\n");
	}
}

/* TODO */
void RGIndexSortNodes(RGIndex *index, RGBinary *rg, int32_t numThreads)
{
	int64_t i;
	ThreadRGIndexSortData *data=NULL;
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
			/*
			   fprintf(stderr, "HERE\t%Ld\t%Ld\t%Ld\n",
			   i,
			   pivots[i],
			   pivots[i+1]);
			   */
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
	for(i=1;i<index->length;i++) {
		assert(RGIndexCompareAt(index, rg, i-1, i) <= 0);
	}

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
	/*
	fprintf(stderr, "HERE: threadID:%d\tlow:%Ld\thigh:%Ld\n",
			data->threadID,
			data->low,
			data->high);
			*/
	RGIndexQuickSortNodesHelper(data->index,
			data->rg,
			data->low,
			data->high,
			data->showPercentComplete);
	if(data->showPercentComplete == 1 && VERBOSE >= 0) {
		fprintf(stderr, "\r");
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
	uint32_t tempPos;
	uint8_t tempChr;
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

		/* Proceed if we are with range */
		if(curLow < curHigh) {
			/* Choose a new pivot.  We could do this randomly (randomized quick sort)
			 * but lets just choose the middle element for now.
			 * */
			pivot = (curLow + curHigh)/2;
			assert(pivot >=0 && pivot<index->length);
			assert(curLow >=0 && curLow<index->length);
			assert(curHigh >=0 && curHigh<index->length);


			if(showPercentComplete == 1 && VERBOSE >= 0) {
				if(curPercent < 100.0*((double)(curLow - low))/total) {
					while(curPercent < 100.0*((double)(curLow - low))/total) {
						curPercent += SORT_ROTATE_INC;
					}
					fprintf(stderr, "\r%3.2lf percent complete", 100.0*((double)(curLow - low))/total);
				}
			}

			/* Swap the node at pivot with the node at curHigh */
			tempPos = index->positions[pivot];
			tempChr = index->chromosomes[pivot];
			index->positions[pivot] = index->positions[curHigh];
			index->chromosomes[pivot] = index->chromosomes[curHigh];
			index->positions[curHigh] = tempPos;
			index->chromosomes[curHigh] = tempChr;

			/* Store where the pivot should be */
			pivot = curLow;

			for(i=curLow;i<curHigh;i++) {
				assert(pivot >= 0 && pivot <= curHigh); 
				assert(i>=0 && i <= curHigh);
				if(RGIndexCompareAt(index, rg, i, curHigh) <= 0) {
					/* Swap node at i with node at the new pivot index */
					if(i!=pivot) {
						tempPos = index->positions[pivot];
						tempChr = index->chromosomes[pivot];
						index->positions[pivot] = index->positions[i];
						index->chromosomes[pivot] = index->chromosomes[i];
						index->positions[i] = tempPos;
						index->chromosomes[i] = tempChr;
					}
					/* Increment the new pivot index */
					pivot++;
				}
			}

			/* Move pivot element to correct place */
			if(pivot != curHigh) {
				tempPos = index->positions[pivot];
				tempChr = index->chromosomes[pivot];
				index->positions[pivot] = index->positions[curHigh];
				index->chromosomes[pivot] = index->chromosomes[curHigh];
				index->positions[curHigh] = tempPos;
				index->chromosomes[curHigh] = tempChr;
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
	uint32_t tempPos;
	uint8_t tempChr;

	if(low < high ) {
		/* Choose a new pivot.  We could do this randomly (randomized quick sort)
		 * but lets just choose the middle element for now.
		 * */
		pivot = (low + high)/2;
		assert(pivot >=0 && pivot<index->length);
		assert(low >=0 && low<index->length);
		assert(high >=0 && high<index->length);

		/* Partition the array.
		 * Basically, arrange everything from low to high so that everything that
		 * has value less than or equal to the pivot is on the low of the pivot, and
		 * everthing else (greater than) is on the high side. 
		 * */

		/* Swap the node at pivot with the node at high */
		tempPos = index->positions[pivot];
		tempChr = index->chromosomes[pivot];
		index->positions[pivot] = index->positions[high];
		index->chromosomes[pivot] = index->chromosomes[high];
		index->positions[high] = tempPos;
		index->chromosomes[high] = tempChr;

		/* Store where the pivot should be */
		pivot = low;

		for(i=low;i<high;i++) {
			assert(pivot >= 0 && pivot <= high); 
			assert(i>=0 && i <= high);
			if(RGIndexCompareAt(index, rg, i, high) <= 0) {
				/* Swap node at i with node at the new pivot index */
				if(i!=pivot) {
					tempPos = index->positions[pivot];
					tempChr = index->chromosomes[pivot];
					index->positions[pivot] = index->positions[i];
					index->chromosomes[pivot] = index->chromosomes[i];
					index->positions[i] = tempPos;
					index->chromosomes[i] = tempChr;
				}
				/* Increment the new pivot index */
				pivot++;
			}
		}

		/* Move pivot element to correct place */
		tempPos = index->positions[pivot];
		tempChr = index->chromosomes[pivot];
		index->positions[pivot] = index->positions[high];
		index->chromosomes[pivot] = index->chromosomes[high];
		index->positions[high] = tempPos;
		index->chromosomes[high] = tempChr;

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

	fprintf(stderr, "low:%Ld\thigh:%Ld\tlength:%Ld\n",
			low,
			high,
			length);

	while(increment > 0) {
		assert(increment < length);
		/* Perform insertion sort with jump size as increment */
		if(showPercentComplete==1 && VERBOSE >= 0) {
			fprintf(stderr, "\rincrement:%Ld\tlength:%Ld\n",
					increment,
					length);
		}
		for(i=increment+low;i<=high;i+=increment) {
			if(showPercentComplete==1 && VERBOSE >= 0 && (i-low-increment)%1000==0) {
				fprintf(stderr, "\r%Ld",
						i-low-increment);
			}
			j=i;
			while( (j>=increment+low) && RGIndexCompareAt(index, rg, j-increment, j) > 0) {
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

	index->numTiles=0;
	free(index->tileLengths);
	index->tileLengths=NULL;
	free(index->gaps);
	index->gaps=NULL;

	index->startChr=0;
	index->startPos=0;
	index->endChr=0;
	index->endPos=0;
}

/* TODO */
double RGIndexGetSize(RGIndex *index, int32_t outputSize) 
{
	double total=0;

	total += sizeof(RGIndex); /* memory used by the index base structure */
	total += sizeof(uint32_t)*index->length;/* memory used by positions */
	total += sizeof(uint8_t)*index->length;/* memory used by positions */
	total += sizeof(int32_t)*index->numTiles;/* memory used by tileLengths */
	total += sizeof(int32_t)*(index->numTiles-1);/* memory used by gaps */

	switch(outputSize) {
		case KILOBYTES:
			return (total/1024);
			break;
		case MEGABYTES:
			return (total/1048576);
			break;
		case GIGABYTES:
			return (total/1073741824);
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
		/* Print chomosomes */
		fwrite(index->chromosomes, sizeof(uint8_t), index->length, fp);

		/* Print positions */
		fwrite(index->positions, sizeof(uint32_t), index->length, fp);

		/* Print the tileLengths */
		fwrite(index->tileLengths, sizeof(int32_t), index->numTiles, fp);

		/* Print the gaps */
		fwrite(index->gaps, sizeof(int32_t), index->numTiles-1, fp);
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
		/* Read in the chromosomes.  Use the position array as temp
		 * storage */
		if(fread(index->chromosomes, sizeof(uint8_t), index->length, fp)==EOF) {
			PrintError("RGIndexRead",
					NULL,
					"Could not read in chromosomes",
					Exit,
					EndOfFile);
		}

		/* Read in positions */
		if(fread(index->positions, sizeof(uint32_t), index->length, fp)==EOF) {
			PrintError("RGIndexRead",
					NULL,
					"Could not read in positions",
					Exit,
					EndOfFile);
		}

		/* Read the tileLengths */
		if(fread(index->tileLengths, sizeof(int32_t), index->numTiles, fp)==EOF) {
			PrintError("RGIndexRead",
					NULL,
					"Could not read in tile lengths",
					Exit,
					EndOfFile);
		}

		/* Read the gaps */
		if(fread(&index->gaps, sizeof(int32_t), index->numTiles-1, fp)==EOF) {
			PrintError("RGIndexRead",
					NULL,
					"Could not read in gaps",
					Exit,
					EndOfFile);
		}
	}
}

/* TODO */
void RGIndexPrintHeader(FILE *fp, RGIndex *index, int32_t binaryOutput)
{
	if(binaryOutput == 0) {
		fprintf(fp, "%Ld\t%d\t%d\t%d\t%d\t%d\t%d\n",
				index->length,
				index->totalLength,
				index->numTiles,
				index->startChr,
				index->startPos,
				index->endChr,
				index->endPos);
	}
	else {
		/* Print Header */
		fwrite(&index->length, sizeof(int64_t), 1, fp);
		fwrite(&index->totalLength, sizeof(int32_t), 1, fp);
		fwrite(&index->numTiles, sizeof(int32_t), 1, fp);
		fwrite(&index->startChr, sizeof(int32_t), 1, fp);
		fwrite(&index->startPos, sizeof(int32_t), 1, fp);
		fwrite(&index->endChr, sizeof(int32_t), 1, fp);
		fwrite(&index->endPos, sizeof(int32_t), 1, fp);
	}
}

/* TODO */
void RGIndexReadHeader(FILE *fp, RGIndex *index, int32_t binaryInput)
{
	/* Read in header */
	if(binaryInput == 0) {
		if(fscanf(fp, "%Ld %d %d %d %d %d %d",
					&index->length,
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
		if(fread(&index->length, sizeof(int64_t), 1, fp)==EOF
				|| fread(&index->totalLength, sizeof(int32_t), 1, fp)==EOF
				|| fread(&index->numTiles, sizeof(int32_t), 1, fp)==EOF 
				|| fread(&index->startChr, sizeof(int32_t), 1, fp)==EOF
				|| fread(&index->startPos, sizeof(int32_t), 1, fp)==EOF
				|| fread(&index->endChr, sizeof(int32_t), 1, fp)==EOF
				|| fread(&index->endPos, sizeof(int32_t), 1, fp)==EOF) {
			PrintError("RGIndexReadHeader",
					NULL,
					"Could not read header",
					Exit,
					EndOfFile);
		}
	}

	/* Error checking */
	assert(index->length > 0);
	assert(index->totalLength > 0);
	assert(index->numTiles > 0);
	assert(index->startChr > 0);
	assert(index->startPos > 0);
	assert(index->endChr > 0);
	assert(index->endPos > 0);
}

/* TODO */
/* We will append the matches if matches have already been found */
void RGIndexGetMatches(RGIndex *index, RGBinary *rg, char *read, int8_t direction, int32_t offset, RGMatch *m)
{
	int64_t i;
	int64_t startIndex=-1;
	int64_t nodeIndex=-1;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "RGIndexGetMatches.  Searching for read:%s.\n",
				read);
	}

	/* Get the index of the index */
	nodeIndex = RGIndexGetFirstIndex(index, rg, read);
	if(VERBOSE >= DEBUG) {
		if(nodeIndex < 0) {
			fprintf(stderr, "Found index:%lld\n", nodeIndex);
		}
		else {
			fprintf(stderr, "Found index:%lld\tchr:%d\tpos:%d\n", 
					nodeIndex,
					(int32_t)index->chromosomes[nodeIndex],
					index->positions[nodeIndex]);
		}
	}

	/* Copy over all matches */
	while(nodeIndex >= 0 && 
			nodeIndex < index->length &&
			RGIndexCompareRead(index, rg, read, nodeIndex)==0) {

		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Will append...\n");
		}
		/* Copy over the matches */
		/* (Re)Allocate memory for the new matches */
		startIndex = m->numEntries;
		m->numEntries++;
		m->positions = realloc(m->positions, sizeof(uint32_t)*(m->numEntries)); 
		if(NULL == m->positions) {
			PrintError("RGIndexGetMatches",
					"m->positions",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		m->chromosomes = realloc(m->chromosomes, sizeof(uint8_t)*(m->numEntries)); 
		if(NULL == m->chromosomes) {
			PrintError("RGIndexGetMatches",
					"m->chromosomes",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		m->strand = realloc(m->strand, sizeof(int8_t)*(m->numEntries)); 
		if(NULL == m->strand) {
			PrintError("RGIndexGetMatches",
					"m->strand",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}

		/* Copy over */
		for(i=startIndex;i<m->numEntries;i++) {
			m->positions[i] = index->positions[nodeIndex] - offset;
			m->chromosomes[i] = index->chromosomes[nodeIndex];
			m->strand[i] = direction;
		}

		/* Update to the next node */
		nodeIndex++;
	}
}

/* TODO */
int64_t RGIndexGetFirstIndex(RGIndex *index,
		RGBinary *rg,
		char *read)
{
	int64_t low=0;
	int64_t high=index->length-1;
	int64_t mid=-1;
	int32_t cmp;
	int32_t cont = 1;

	while(low <= high && cont==1) {
		mid = (low+high)/2;
		cmp = RGIndexCompareRead(index, rg, read, mid);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "low:%lld\tmid:%lld\thigh:%lld\tcmp:%d\n",
					low,
					mid,
					high,
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
	if(1 == cont) {
		return -1;
	}
	else {
		while(mid > 0 && RGIndexCompareRead(index, rg, read, mid-1) == 0) {
			mid--;
		}
		return mid;
	}
}

/* TODO */
int32_t RGIndexCompareAt(RGIndex *index,
		RGBinary *rg,
		int64_t a,
		int64_t b)
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
		int64_t a)
{
	assert(a>=0 && a<index->length);

	int32_t i, j;
	int32_t curReadPos=0;
	int32_t readLength = strlen((char*)read);
	uint32_t aChr = index->chromosomes[a];
	uint32_t aPos = index->positions[a];

	uint32_t aCurTilePos;
	uint8_t aBase;
	uint8_t readBase;

	/* Compare base by base */
	aCurTilePos = aPos;

	for(i=0;i<index->numTiles && curReadPos < readLength;i++) { /* For each tile */
		for(j=0;j<index->tileLengths[i] && curReadPos < readLength;j++) { /* For each position in the tile */
			aBase=ToLower(RGBinaryGetBase(rg,
						aChr,
						aCurTilePos));
			readBase = ToLower(read[curReadPos]);

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
