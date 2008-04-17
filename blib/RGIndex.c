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
void RGIndexCreate(RGIndex *index, RGBinary *rg, int includeRepeats, int includeNs) {

	/* The sort will take care of most of the work.  We just want 
	 * to make sure that we only include sequence that agrees with
	 * includeRepeats and includeNs
	 * */

	unsigned int curPos, curChr, insert, i, j, curTilePos;
	unsigned int chrIndex = 0;

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

		/* Allocate a new chromosome */
		rg->chromosomes = realloc(rg->chromosomes, sizeof(unsigned char)*(chrIndex+1));
		if(NULL == rg->chromosomes) {
			PrintError("RGIndexCreate",
					"rg->chromosomes",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}

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
				index->positions = realloc(index->positions, sizeof(unsigned int)*index->length);
				if(NULL == index->positions) {
					PrintError("RGBinaryCreate",
							"index->positions",
							"Could not reallocate memory",
							Exit,
							ReallocMemory);
				}
				index->chromosomes = realloc(index->chromosomes, sizeof(unsigned char)*index->length);
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
	if(VERBOSE >= 0) {
		fprintf(stderr, "\r[%d,%d]\n",
				curChr,
				curPos);
	}
}

/* TODO */
void RGIndexCleanUpIndex(RGIndex *index, RGBinary *rg, int numThreads) 
{
	if(index->length > 0) {
		/* Sort the nodes in the index */
		fprintf(stderr, "Sorting...\n");
		RGIndexSortNodes(index, rg, numThreads);
		fprintf(stderr, "Sorted.\n");
	}
}

/* TODO */
void RGIndexSortNodes(RGIndex *index, RGBinary *rg, int numThreads)
{
	int i;
	ThreadRGIndexSortData *data=NULL;
	pthread_t *threads=NULL;
	int errCode;
	void *status=NULL;
	unsigned int *pivots;
	int max, maxIndex;

	/* Allocate memory for the thread arguments */
	data = malloc(sizeof(ThreadRGIndexSortData)*numThreads);
	if(NULL==data) {
		PrintError("RGIndexSortNodes",
				"data",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the thread pointers */
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
	pivots = malloc(sizeof(unsigned int)*(2*numThreads));
	if(NULL == pivots) {
		PrintError("RGIndexSortNodes",
				"pivots",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Get the pivots and presort */
	if(numThreads > 1) {
		RGIndexQuickSortNodesHelper(index,
				rg,
				0,
				index->length-1,
				0,
				NULL,
				0,
				index->length-1,
				1,
				pivots,
				1,
				numThreads);
		/* The last one must be less than index->length */
		pivots[2*numThreads-1]--;
	}
	else {
		assert(numThreads == 1);
		pivots[0] = 0;
		pivots[1] = index->length - 1;
	}
	for(i=0;i<2*numThreads;i+=2) {
		assert(pivots[i] >= 0 && pivots[i] < index->length);
		assert(pivots[i+1] >= 0 && pivots[i+1] < index->length);
		assert(pivots[i] <= pivots[i+1]);
		if(i==0) {
			assert(pivots[i] == 0);
		}
		if(i==2*numThreads-2) {
			assert(pivots[i+1] == index->length-1);
		}
		if(i>1) {
			assert(pivots[i] > pivots[i-1]);
		}
	}

	/* Initialize data */
	maxIndex=0;
	max = data[0].high-data[i].low;
	for(i=0;i<numThreads;i++) {
		data[i].index = index;
		data[i].rg = rg;
		data[i].threadID = i;
		data[i].low = pivots[2*i];
		data[i].high = pivots[2*i+1];
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
				&data[i]); /* data to routine */
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
		if(VERBOSE >= 0) {
			fprintf(stderr, "\rWaiting for other threads to complete...");
		}
	}

	/* Test that we sorted correctly */
	for(i=1;i<index->length;i++) {
		assert(RGIndexCompareAt(index, rg, i-1, i) <= 0);
	}

	/* Free memory */
	free(threads);
	free(data);
}

/* TODO */
void *RGIndexQuickSortNodes(void *arg)
{
	/* thread arguments */
	ThreadRGIndexSortData *data = (ThreadRGIndexSortData*)(arg);
	RGIndex *index = data->index;
	RGBinary *rg = data->rg;
	unsigned int low = data->low;
	unsigned int high = data->high;
	int showPercentComplete = data->showPercentComplete;
	double curPercent = 0.0;

	/* Call helper */
	if(showPercentComplete == 1) {
		fprintf(stderr, "0 percent complete");
	}
	RGIndexQuickSortNodesHelper(index, rg, low, high, showPercentComplete, &curPercent, low, high, 0, NULL, 0, 0);
	if(showPercentComplete == 1) {
		fprintf(stderr, "\r");
	}

	return arg;
}

void RGIndexQuickSortNodesHelper(RGIndex *index,
		RGBinary *rg,
		unsigned int low,
		unsigned int high,
		int showPercentComplete,
		double *curPercent,
		unsigned int lowTotal,
		unsigned int highTotal,
		int savePivots,
		unsigned int *pivots,
		int lowPivot,
		int highPivot)
{
	/* local variables */
	unsigned int i;
	unsigned int pivot = 0;
	unsigned int tempPos;
	unsigned char tempChr;
	unsigned int total = highTotal-lowTotal;

	if(low < high) {
		/* Choose a new pivot.  We could do this randomly (randomized quick sort)
		 * but lets just choose the middle element for now.
		 * */
		pivot = (low + high)/2;
		if(showPercentComplete == 1 && VERBOSE >= 0) {
			if((*curPercent) < 100.0*((double)low)/total) {
				while((*curPercent) < 100.0*((double)low)/total) {
					(*curPercent) += SORT_ROTATE_INC;
				}
				fprintf(stderr, "\r%3.2lf percent complete", 100.0*((double)low)/total);
			}
		}

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
			if(RGIndexCompareAt(index, rg, i, high) <= 0) {
				/* Swap node at i with node at the new pivot index */
				tempPos = index->positions[pivot];
				tempChr = index->chromosomes[pivot];
				index->positions[pivot] = index->positions[i];
				index->chromosomes[pivot] = index->chromosomes[i];
				index->positions[i] = tempPos;
				index->chromosomes[i] = tempChr;
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

		if(savePivots == 1 && lowPivot >= highPivot) {
			assert(pivots!=NULL);
			/* Save pivots if necessary */
			pivots[2*lowPivot-2] = low;
			pivots[2*lowPivot-1] = high+1;
			return;
		}
		else {

			/* Call recursively */
			if(pivot > 0) {
				RGIndexQuickSortNodesHelper(index, rg, low, pivot-1, showPercentComplete, curPercent, lowTotal, highTotal, savePivots, pivots, lowPivot, (lowPivot+highPivot)/2);
			}
			if(showPercentComplete == 1 && VERBOSE >= 0) {
				if((*curPercent) < 100.0*((double)pivot)/total) {
					while((*curPercent) < 100.0*((double)pivot)/total) {
						(*curPercent) += SORT_ROTATE_INC;
					}
					fprintf(stderr, "\r%3.2lf percent complete", 100.0*((double)pivot)/total);
				}
			}
			if(pivot < UINT_MAX) {
				RGIndexQuickSortNodesHelper(index, rg, pivot+1, high, showPercentComplete, curPercent, lowTotal, highTotal, savePivots, pivots, (lowPivot+highPivot)/2 + 1, highPivot);
			}
			if(showPercentComplete == 1 && VERBOSE >= 0) {
				if((*curPercent) < 100.0*((double)high)/total) {
					while((*curPercent) < 100.0*((double)high)/total) {
						(*curPercent) += SORT_ROTATE_INC;
					}
					fprintf(stderr, "\r%3.2lf percent complete", 100.0*((double)high)/total);
				}
			}
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
double RGIndexGetSize(RGIndex *index, int outputSize) 
{
	double total=0;

	total += sizeof(RGIndex); /* memory used by the index base structure */
	total += sizeof(unsigned int)*index->length;/* memory used by positions */
	total += sizeof(unsigned char)*index->length;/* memory used by positions */
	total += sizeof(unsigned int)*index->numTiles;/* memory used by tileLengths */
	total += sizeof(int)*(index->numTiles-1);/* memory used by gaps */

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
void RGIndexPrint(FILE *fp, RGIndex *index, int binaryOutput)
{
	unsigned int i;
	unsigned int tempInt;

	/* Print header */
	RGIndexPrintHeader(fp, index, binaryOutput);

	if(binaryOutput == 0) {

		/* Print the positions and chromosomes */
		for(i=0;i<index->length;i++) {
			fprintf(fp, "%u\t%u\n", 
					index->positions[i],
					(unsigned int)index->chromosomes[i]);
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

		/* We could just separate everything into huge arrays, but then we would double memory requirements when writing */

		/* Print the chromosomes first (better for reading since
		 * we can just use the allocated positions to read in the
		 * chromosomes.) */
		for(i=0;i<index->length;i++) {
			tempInt = index->chromosomes[i];
			tempInt = htonl(tempInt);
			fwrite(&tempInt, sizeof(unsigned int), 1, fp);
		}

		/* Print positions */
		/* Convert to network order */
		for(i=0;i<index->length;i++) {
			index->positions[i] = htonl(index->positions[i]);
		}
		/* Print the whole array */
		fwrite(index->positions, sizeof(unsigned int), index->length, fp);
		/* Convert back */
		for(i=0;i<index->length;i++) {
			index->positions[i] = htonl(index->positions[i]);
		}

		/* Print the tileLengths */
		for(i=0;i<index->numTiles;i++) {
			tempInt = index->tileLengths[i];
			tempInt = htonl(tempInt);
			fwrite(&tempInt, sizeof(unsigned int), 1, fp);
		}

		/* Print the gaps */
		for(i=0;i<index->numTiles-1;i++) {
			tempInt = index->gaps[i];
			tempInt = htonl(tempInt);
			fwrite(&tempInt, sizeof(int), 1, fp);
		}
	}
}

/* TODO */
int RGIndexRead(FILE *fp, RGIndex *index, int binaryInput)
{
	unsigned int i;
	unsigned int tempInt;

	/* Read in the header */
	RGIndexReadHeader(fp, index, binaryInput);

	assert(index->length > 0);

	/* Allocate memory for the positions */
	index->positions = malloc(sizeof(unsigned int)*index->length);
	if(NULL == index->positions) {
		PrintError("RGIndexRead",
				"index->positions",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the chromosomes */
	index->chromosomes = malloc(sizeof(unsigned int)*index->length);
	if(NULL == index->chromosomes) {
		PrintError("RGIndexRead",
				"index->chromosomes",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the tile lengths */
	index->tileLengths = malloc(sizeof(unsigned int)*index->numTiles);
	if(NULL == index->tileLengths) {
		PrintError("RGIndexRead",
				"index->tileLengths",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the gaps */
	index->gaps = malloc(sizeof(unsigned int)*(index->numTiles-1));
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
			index->chromosomes[i] = (unsigned char)tempInt;
		}

		/* Read the tileLengths */
		for(i=0;i<index->numTiles;i++) {
			if(fscanf(fp, "%u",
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
			if(fscanf(fp, "%u",
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
		if(fread(index->positions, sizeof(unsigned int), index->length, fp)==EOF) {
			PrintError("RGIndexRead",
					NULL,
					"Could not read in chromosomes",
					Exit,
					EndOfFile);
		}
		/* Convert to host order and copy over */
		for(i=0;i<index->length;i++) {
			index->positions[i] = ntohl(index->positions[i]);
			index->chromosomes[i] = (unsigned char)index->positions[i];
		}

		/* Read in positions */
		if(fread(index->positions, sizeof(unsigned int), index->length, fp)==EOF) {
			PrintError("RGIndexRead",
					NULL,
					"Could not read in positions",
					Exit,
					EndOfFile);
		}
		/* Convert to host order  */
		for(i=0;i<index->length;i++) {
			index->positions[i] = ntohl(index->positions[i]);
		}

		/* Read the tileLengths */
		if(fread(index->tileLengths, sizeof(unsigned int), index->numTiles, fp)==EOF) {
			PrintError("RGIndexRead",
					NULL,
					"Could not read in tile lengths",
					Exit,
					EndOfFile);
		}
		/* Convert to host order */
		for(i=0;i<index->numTiles;i++) {
			index->tileLengths[i] = ntohl(index->tileLengths[i]);
		}

		/* Read the gaps */
		for(i=0;i<index->numTiles-1;i++) {
			if(fread(&index->gaps[i], sizeof(unsigned int), 1, fp)==EOF) {
				PrintError("RGIndexRead",
						NULL,
						"Could not read in gap",
						Exit,
						EndOfFile);
			}
			/* Convert to host order */
			index->gaps[i] = ntohl(index->gaps[i]);
		}
	}
	return 1;
}

/* TODO */
void RGIndexPrintHeader(FILE *fp, RGIndex *index, int binaryOutput)
{
	unsigned int length=index->length;
	unsigned int totalLength=index->totalLength;
	unsigned int numTiles=index->numTiles;
	unsigned int startChr=index->startChr;
	unsigned int startPos=index->startPos;
	unsigned int endChr=index->endChr;
	unsigned int endPos=index->endPos; 

	if(binaryOutput == 0) {
		fprintf(fp, "%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
				length,
				totalLength,
				numTiles,
				startChr,
				startPos,
				endChr,
				endPos);
	}
	else {

		/* Big Endian/Little Endian conversion */
		length=htonl(length);
		totalLength=htonl(totalLength);
		numTiles=htonl(numTiles);
		startChr=htonl(startChr);
		startPos=htonl(startPos);
		endChr=htonl(endChr);
		endPos=htonl(endPos);

		/* Print Header */
		fwrite(&length, sizeof(unsigned int), 1, fp);
		fwrite(&totalLength, sizeof(unsigned int), 1, fp);
		fwrite(&numTiles, sizeof(unsigned int), 1, fp);
		fwrite(&startChr, sizeof(unsigned int), 1, fp);
		fwrite(&startPos, sizeof(unsigned int), 1, fp);
		fwrite(&endChr, sizeof(unsigned int), 1, fp);
		fwrite(&endPos, sizeof(unsigned int), 1, fp);
	}
}

/* TODO */
void RGIndexReadHeader(FILE *fp, RGIndex *index, int binaryInput)
{
	unsigned int length;
	unsigned int totalLength;
	unsigned int numTiles;
	unsigned int startChr;
	unsigned int startPos;
	unsigned int endChr;
	unsigned int endPos;

	/* Read in header */
	if(binaryInput == 0) {
		if(fscanf(fp, "%u %u %u %u %u %u %u",
					&length,
					&totalLength,
					&numTiles,
					&startChr,
					&startPos,
					&endChr,
					&endPos)==EOF) {
			PrintError("RGIndexReadHeader",
					NULL,
					"Could not read header",
					Exit,
					EndOfFile);
		}
	}
	else {
		if(fread(&length, sizeof(unsigned int), 1, fp)==EOF
				|| fread(&totalLength, sizeof(unsigned int), 1, fp)==EOF
				|| fread(&numTiles, sizeof(unsigned int), 1, fp)==EOF 
				|| fread(&startChr, sizeof(unsigned int), 1, fp)==EOF
				|| fread(&startPos, sizeof(unsigned int), 1, fp)==EOF
				|| fread(&endChr, sizeof(unsigned int), 1, fp)==EOF
				|| fread(&endPos, sizeof(unsigned int), 1, fp)==EOF) {
			PrintError("RGIndexReadHeader",
					NULL,
					"Could not read header",
					Exit,
					EndOfFile);
		}
		/* Big Endian/Little Endian conversion */
		length = ntohl(length);
		totalLength = ntohl(totalLength);
		numTiles = ntohl(numTiles);
		startChr = ntohl(startChr);
		startPos = ntohl(startPos);
		endChr = ntohl(endChr);
		endPos = ntohl(endPos);
	}

	if(VERBOSE > 0) {
		fprintf(stderr, "Read Header:%u,%u,%u,%u,%u,%u,%u\n",
				length,
				totalLength,
				numTiles,
				startChr,
				startPos,
				endChr,
				endPos);
	}

	/* Adjust header - see bpreprocess */
	index->length = length;
	index->totalLength = totalLength;
	index->numTiles = numTiles;
	index->startChr = startChr;
	index->startPos = startPos;
	index->endChr = endChr;
	index->endPos = endPos;

	/* Error checking ? */
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
int RGIndexGetMatches(RGIndex *index, RGBinary *rg, char *read, char direction, RGMatch *m)
{
	unsigned int i;
	long long int startIndex=-1;
	long long int nodeIndex=-1;

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
					(int)index->chromosomes[nodeIndex],
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
		m->positions = realloc(m->positions, sizeof(unsigned int)*(m->numEntries)); 
		if(NULL == m->positions) {
			PrintError("RGIndexGetMatches",
					"m->positions",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		m->chromosomes = realloc(m->chromosomes, sizeof(unsigned char)*(m->numEntries)); 
		if(NULL == m->chromosomes) {
			PrintError("RGIndexGetMatches",
					"m->chromosomes",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		m->strand = realloc(m->strand, sizeof(char)*(m->numEntries)); 
		if(NULL == m->strand) {
			PrintError("RGIndexGetMatches",
					"m->strand",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}

		/* Copy over */
		for(i=startIndex;i<m->numEntries;i++) {
			m->positions[i] = index->positions[nodeIndex];
			m->chromosomes[i] = index->chromosomes[nodeIndex];
			m->strand[i] = direction;
		}

		/* Update to the next node */
		nodeIndex++;
	}

	return 1;
}

/* TODO */
long long int RGIndexGetFirstIndex(RGIndex *index,
		RGBinary *rg,
		char *read)
{
	long long int low=0;;
	long long int high=index->length-1;
	long long int mid=-1;
	int cmp;
	int cont = 1;

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
int RGIndexCompareAt(RGIndex *index,
		RGBinary *rg,
		unsigned int a,
		unsigned int b)
{
	assert(a>=0 && a<index->length);
	assert(b>=0 && b<index->length);

	int i, j;
	int aChr = index->chromosomes[a];
	int aPos = index->positions[a];
	int bChr = index->chromosomes[b];
	int bPos = index->positions[b];

	int aCurTilePos;
	int bCurTilePos;
	char aBase;
	char bBase;

	/* Compare base by base */
	aCurTilePos = aPos;
	bCurTilePos = bPos;

	for(i=0;i<index->numTiles;i++) { /* For each tile */
		for(j=0;j<index->tileLengths[i];j++) { /* For each position in the tile */
			aBase = ToLower(RGBinaryGetBase(rg,
						aChr,
						aCurTilePos));
			bBase =ToLower( RGBinaryGetBase(rg,
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
int RGIndexCompareRead(RGIndex *index,
		RGBinary *rg,
		char *read,
		unsigned int a)
{
	assert(a>=0 && a<index->length);

	int i, j;
	int curReadPos=0;
	int readLength = strlen(read);
	int aChr = index->chromosomes[a];
	int aPos = index->positions[a];

	int aCurTilePos;
	char aBase;
	char readBase;

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
