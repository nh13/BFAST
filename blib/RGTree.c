#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <limits.h>
#include <string.h>
#include <pthread.h>
#include "BLibDefinitions.h"
#include "BError.h"
#include "RGMatch.h"
#include "RGTree.h"

/* TODO */
/* We assume the root has already been initialized */
/* We could insert both forward and reverse compliment of the two sequences when creating the tree,
 * but this would double the space of each tree (practically too large).  It is easier just to 
 * double the number of searches.
 * */
int RGTreeInsert(RGTree *tree, char *sequenceOne, char *sequenceTwo, unsigned int matchLength, unsigned int chromosome, unsigned int position) 
{
	assert(tree->matchLength == matchLength);

	unsigned int indexOne = RGTreeGetIndexFromSequence(sequenceOne, matchLength);
	unsigned int indexTwo = RGTreeGetIndexFromSequence(sequenceTwo, matchLength);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "RGTreeInsert: inserting indexes %d and %d.\n",
				indexOne,
				indexTwo);
	}

	/* Create a new node */
	tree->numNodes++;

	/* Allocate memory for a new node */
	tree->nodes = realloc(tree->nodes, sizeof(RGTreeNode)*tree->numNodes);
	if(NULL == tree->nodes) {
		PrintError("RGTreeInsert",
				"tree->nodes",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}

	/* Initialize node */
	tree->nodes[tree->numNodes-1].numEntries = 1;
	tree->nodes[tree->numNodes-1].indexOne = indexOne;
	tree->nodes[tree->numNodes-1].indexTwo = indexTwo;
	tree->nodes[tree->numNodes-1].positions = NULL;
	tree->nodes[tree->numNodes-1].chromosomes = NULL;

	/* Allocate memory for the node members */
	RGTreeNodeAllocate(&tree->nodes[tree->numNodes-1], 1);

	/* Copy over */
	tree->nodes[tree->numNodes-1].positions[tree->nodes[tree->numNodes-1].numEntries-1] = position;
	tree->nodes[tree->numNodes-1].chromosomes[tree->nodes[tree->numNodes-1].numEntries-1] = chromosome;

	return 1;
}

/* TODO */
void RGTreeCleanUpTree(RGTree *tree, int numThreads) 
{
	unsigned int i, j;
	unsigned int prevIndex = 0;

	if(tree->numNodes > 0) {

		/* Sort the nodes in the tree */
		if(VERBOSE >= 0) {
			fprintf(stderr, "Sorting (this will speed up):\n");
		}
		RGTreeSortNodes(tree, numThreads);
		if(VERBOSE >= 0) {
			fprintf(stderr, "\nSorting complete\n");
		}

		/* Remove duplicates */
		/* Assumes the nodes are sorted */
		if(VERBOSE >= 1) {
			fprintf(stderr, "Merging\n");
		}
		/* This will remove duplicates in-place and in O(n) time */
		for(i=1;i<tree->numNodes;i++) {
			if(RGTreeNodeCompare(&tree->nodes[i], &tree->nodes[prevIndex])==0) {
				/* Append the entry to prevIndex nodes */
				unsigned int start = tree->nodes[prevIndex].numEntries;
				/* Allocate memory for chromosomes and positions */
				tree->nodes[prevIndex].numEntries += tree->nodes[i].numEntries;
				RGTreeNodeAllocate(&tree->nodes[prevIndex], tree->nodes[prevIndex].numEntries);
				/* Copy over chromosomes and positions */
				for(j=start;j<tree->nodes[prevIndex].numEntries;j++) {
					tree->nodes[prevIndex].positions[j] = tree->nodes[i].positions[j-start];
					tree->nodes[prevIndex].chromosomes[j] = tree->nodes[i].chromosomes[j-start];
				}
				/* Free memory of the current node */
				free(tree->nodes[i].positions);
				free(tree->nodes[i].chromosomes);
				/* Nullify the current node */
				tree->nodes[i].positions = NULL;
				tree->nodes[i].chromosomes = NULL;
				tree->nodes[i].numEntries=0;
				tree->nodes[i].indexOne=-1;
				tree->nodes[i].indexTwo=-1;
			}
			else {
				prevIndex++;
				/* Move to prevIndex */
				if(prevIndex < i) {
					RGTreeNodeCopy(&tree->nodes[i], &tree->nodes[prevIndex]);
					/* Nullify the current node */
					tree->nodes[i].positions = NULL;
					tree->nodes[i].chromosomes = NULL;
					tree->nodes[i].numEntries=0;
					tree->nodes[i].indexOne=-1;
					tree->nodes[i].indexTwo=-1;
				}
			}
		}
		tree->numNodes = prevIndex+1;
		/* Reallocate memory to reflect new number of nodes */
		tree->nodes = realloc(tree->nodes, sizeof(RGTreeNode)*tree->numNodes);
		if(NULL == tree->nodes) {
			PrintError("RGTreeCleanUpTree",
					"tree->nodes",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}

		/* Sort each node */
		for(i=0;i<tree->numNodes;i++) {
			RGTreeQuickSortNode(tree, i, 0, tree->nodes[i].numEntries-1);
		}

		if(VERBOSE >= 1) {
			fprintf(stderr, "Merging complete\n");
		}
	}
}

/* TODO */
void RGTreeSortNodes(RGTree *tree, int numThreads)
{
	int i, j;
	RGTreeNode temp;
	ThreadRGTreeSortData *data=NULL;
	pthread_t *threads=NULL;
	int errCode;
	void *status=NULL;
	int splitLengths=-1;
	unsigned int *indexes=NULL;

	/* Allocate memory for the thread arguments */
	data = malloc(sizeof(ThreadRGTreeSortData)*numThreads);
	if(NULL==data) {
		PrintError("RGTreeSortNodes",                "data",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for the thread pointers */
	threads = malloc(sizeof(pthread_t)*numThreads);
	if(NULL==threads) {
		PrintError("RGTreeSortNodes",
				"threads",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Split up sort int 1/numThreads parts */
	splitLengths = tree->numNodes/numThreads;

	/* Initialize data */
	for(i=0;i<numThreads;i++) {
		data[i].tree = tree;
		data[i].threadID = i;
		data[i].low = i*splitLengths;
		if(i==numThreads-1) {
			data[i].high = tree->numNodes-1;
		}
		else {
			data[i].high = (i+1)*splitLengths - 1;
		}
		assert(data[i].low >= 0 && data[i].high < tree->numNodes);
		/*
		 *            fprintf(stderr, "i:%d\tlow:%d\thigh:%d\n",
		 *                       i,
		 *                                  data[i].low,
		 *                                             data[i].high);
		 *                                                        */
	}

	/* Check that we split correctly */
	for(i=1;i<numThreads;i++) {
		assert(data[i-1].high < data[i].low);
	}

	/* Create threads */
	for(i=0;i<numThreads;i++) {
		/* Start thread */
		errCode = pthread_create(&threads[i], /* thread struct */
				NULL, /* default thread attributes */
				RGTreeQuickSortNodes, /* start routine */
				&data[i]); /* data to routine */
		if(0!=errCode) {
			PrintError("RGTreeSortNodes",
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
			PrintError("RGTreeSortNodes",
					"pthread_join: errCode",
					"Thread returned an error",
					Exit,
					ThreadError);
		}
	}
	/* Allocate memory for indexes */
	indexes = malloc(sizeof(unsigned int)*numThreads);
	if(NULL==indexes) {
		PrintError("RGTreeSOrtNodes",
				"indexes",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Initialize indexes */
	for(i=0;i<numThreads;i++) {
		indexes[i] = data[i].low;
	}
	/* Merge the individual sorts together */
	for(i=0;i<tree->numNodes;i++) {
		int whichThread = -1;
		/* Check the result of each thread for the smallest element */
		for(j=0;j<numThreads;j++) {
			/* If the current item is valid, and the node is smaller */
			if(indexes[j] < tree->numNodes) {
				if(whichThread == -1 || RGTreeNodeCompare(&tree->nodes[indexes[j]], &tree->nodes[indexes[whichThread]]) <= 0) {
					whichThread = j;
				}
			}
		}
		assert(whichThread!=-1);
		/* Update min index */
		unsigned int minIndex = indexes[whichThread];

		/* Store the current node to swap */
		RGTreeNodeCopy(&tree->nodes[minIndex], &temp);
		/* Shift all nodes up one starting at i and ending at minIndex */
		for(j=minIndex-1;j>=i;j--) {
			RGTreeNodeCopy(&tree->nodes[j], &tree->nodes[j+1]);
		}
		/* Store the current node at i */
		RGTreeNodeCopy(&temp, &tree->nodes[i]);
		/* Update thread indexes if necessary */
		for(j=0;j<numThreads;j++) {
			if(indexes[j] <= minIndex) {
				indexes[j]++;
			}
		}
	}
	/* Test that we sorted correctly */
	for(i=1;i<tree->numNodes;i++) {
		assert( RGTreeNodeCompare(&tree->nodes[i-1], &tree->nodes[i]) <=0);
	}

	/* Free memory */
	free(indexes);
	free(threads);
	free(data);
}


/* TODO */
void *RGTreeQuickSortNodes(void *arg)
{
	/* thread arguments */
	ThreadRGTreeSortData *data = (ThreadRGTreeSortData*)(arg);
	RGTree *tree = data->tree;
	unsigned int low = data->low;
	unsigned int high = data->high;

	/* Call helper */
	RGTreeQuickSortNodesHelper(tree, low, high);

	return arg;
}


/* TODO */
void RGTreeQuickSortNodesHelper(RGTree *tree, unsigned int low, unsigned int high) 
{
	unsigned int i;
	unsigned int pivot = -1;
	RGTreeNode *temp=NULL;


	if(low < high) {
		/* Allocate temp */
		temp = malloc(sizeof(RGTreeNode));
		if(NULL == temp) {
			PrintError("RGTreeQuickSortNodes",
					"temp",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}

		/* Choose a new pivot.  We could do this randomly (randomized quick sort)
		 * but lets just choose the middle element for now.
		 * */
		pivot = (low + high)/2;

		/* Partition the array.
		 * Basically, arrange everything from low to high so that everything that
		 * has value less than or equal to the pivot is on the low of the pivot, and
		 * everthing else (greater than) is on the high side. 
		 * */

		/* Swap the node at pivot with the node at high */
		RGTreeNodeCopy(&tree->nodes[pivot], temp);
		RGTreeNodeCopy(&tree->nodes[high], &tree->nodes[pivot]);
		RGTreeNodeCopy(temp, &tree->nodes[high]);

		/* Store where the pivot should be */
		pivot = low;

		/* move elements */
		for(i=low;i<high;i++) {
			/* Compare with the pivot */
			if(RGTreeNodeCompare(&tree->nodes[i], &tree->nodes[high]) <= 0) {
				/* Swap node at i with node at the new pivot tree */
				RGTreeNodeCopy(&tree->nodes[i], temp);
				RGTreeNodeCopy(&tree->nodes[pivot], &tree->nodes[i]);
				RGTreeNodeCopy(temp, &tree->nodes[pivot]);
				/* Increment the new pivot tree */
				pivot++;
			}
		}
		/* Move pivot element to correct place */
		RGTreeNodeCopy(&tree->nodes[pivot], temp);
		RGTreeNodeCopy(&tree->nodes[high], &tree->nodes[pivot]);
		RGTreeNodeCopy(temp, &tree->nodes[high]);

		/* Free temp */
		free(temp);
		temp = NULL;

		/* Call recursively */
		if(pivot>0) {
			RGTreeQuickSortNodesHelper(tree, low, pivot-1);
		}
		if(pivot < UINT_MAX) {
			RGTreeQuickSortNodesHelper(tree, pivot+1, high);
		}
	}
}

/* TODO */
unsigned int RGTreeGetIndex(RGTree *tree,
		unsigned int indexOne,
		unsigned int indexTwo)
{
	long long int low = 0;
	long long int high = tree->numNodes-1;
	long long int mid;

	/* Binary search */
	while(low <= high) {
		mid = (low + high)/2;
		if(indexOne < tree->nodes[mid].indexOne ||
				(indexOne == tree->nodes[mid].indexOne && indexTwo < tree->nodes[mid].indexTwo)) {
			high = mid-1;
		}
		else if(indexOne == tree->nodes[mid].indexOne && indexTwo == tree->nodes[mid].indexTwo) {
			return mid;
		}
		else {
			low = mid+1;
		}
	}
	return -1;
}

/* TODO */
void RGTreeDelete(RGTree *tree)
{
	unsigned int i;

	/* Delete fields in the individual nodes */
	for(i=0;i<tree->numNodes;i++) {
		free(tree->nodes[i].positions);
		tree->nodes[i].positions=NULL;
		free(tree->nodes[i].chromosomes);
		tree->nodes[i].chromosomes=NULL;
		tree->nodes[i].numEntries=0;
		tree->nodes[i].indexOne=-1;
		tree->nodes[i].indexTwo=-1;
	}

	/* Free nodes*/
	free(tree->nodes);
	tree->nodes=NULL;
	tree->numNodes=0;
	tree->gap=0;
	tree->matchLength=0;
	tree->startChr=0;
	tree->startPos=0;
	tree->endChr=0;
	tree->endPos=0;
}

/* TODO */
double RGTreeGetSize(RGTree *tree, int outputSize) 
{
	unsigned int i;
	double total=0;

	/* Get memory used in each node */
	for(i=0;i<tree->numNodes;i++) {
		total += sizeof(RGTreeNode) + /* memory used by the node */
			sizeof(unsigned int)*tree->nodes[i].numEntries + /* memory used by positions */
			sizeof(unsigned char)*tree->nodes[i].numEntries; /* memory used by chromosomes */
	}

	/* Get memory of the RGTree */
	total += sizeof(RGTree);

	switch(outputSize) {
		case RGT_KILOBYTES:
			return (total/1000);
			break;
		case RGT_MEGABYTES:
			return (total/1000000);
			break;
		case RGT_GIGABYTES:
			return (total/1000000000);
			break;
		default:
			return total;
			break;
	}
}

/* TODO */
void RGTreePrintTree(FILE *fp, RGTree *tree, int binaryOutput)
{
	unsigned int i, j;
	unsigned int tempInt;
	unsigned int *tempIntArr=NULL;
	unsigned int tempIntArrLength=0;

	/* Print header */
	RGTreePrintHeader(fp, tree, binaryOutput);

	/* Print the nodes */
	if(binaryOutput == 0) {
		for(i=0;i<tree->numNodes;i++) {
			fprintf(fp, "\n");
			/* Print the indexes and the number of entries */
			fprintf(fp, "%d\t%d\t%d",
					tree->nodes[i].indexOne,
					tree->nodes[i].indexTwo,
					tree->nodes[i].numEntries);
			/* Print the entries */
			for(j=0;j<tree->nodes[i].numEntries;j++) {
				fprintf(fp, "\t%d\t%d",
						tree->nodes[i].chromosomes[j],
						tree->nodes[i].positions[j]);
			}
		}
		fprintf(fp, "\n");
	}
	else {
		for(i=0;i<tree->numNodes;i++) {
			/* Print the indexes and the number of entries */
			tempInt = tree->nodes[i].indexOne;
			tempInt = htonl(tempInt);
			fwrite(&tempInt, sizeof(unsigned int), 1, fp);
			tempInt = tree->nodes[i].indexTwo;
			tempInt = htonl(tempInt);
			fwrite(&tempInt, sizeof(unsigned int), 1, fp);
			tempInt = tree->nodes[i].numEntries;
			tempInt = htonl(tempInt);
			fwrite(&tempInt, sizeof(unsigned int), 1, fp);
			/* Print the entries */
			if(tree->nodes[i].numEntries > 0) {
				if(tree->nodes[i].numEntries > tempIntArrLength) {
					tempIntArrLength = tree->nodes[i].numEntries;
					tempIntArr = realloc(tempIntArr, sizeof(unsigned int)*tempIntArrLength);
					if(NULL == tempIntArr) {
						PrintError("RGTreePrintTree",
								"tempIntArr",
								"Could not reallocate memory",
								Exit,
								ReallocMemory);
					}
				}
				/* positions */
				for(j=0;j<tree->nodes[i].numEntries;j++) {
					tempIntArr[j] = (unsigned int)tree->nodes[i].positions[j];
					tempIntArr[j] = htonl(tempIntArr[j]);
				}
				fwrite(tempIntArr, sizeof(unsigned int), tree->nodes[i].numEntries, fp);
				/* chromosomes */
				for(j=0;j<tree->nodes[i].numEntries;j++) {
					tempIntArr[j] = (unsigned int)tree->nodes[i].chromosomes[j];
					tempIntArr[j] = htonl(tempIntArr[j]);
				}
				fwrite(tempIntArr,sizeof(unsigned int), tree->nodes[i].numEntries, fp);
			}
		}
		/* Free memory */
		if(tempIntArrLength > 0) {
			free(tempIntArr);
			tempIntArr=NULL;
			tempIntArrLength=0;
		}
	}
}

/* TODO */
int RGTreeReadTree(FILE *fp, RGTree *tree, int binaryInput)
{
	unsigned int i, j;
	unsigned int tempInt;
	unsigned int *tempIntArr=NULL;
	unsigned int tempIntArrLength=0;

	/* Make sure memory of the root has been allocated */
	assert(tree!=NULL);
	assert(tree->nodes==NULL);

	/* Read in the header */
	RGTreeReadHeader(fp, tree, binaryInput);

	assert(tree->numNodes > 0);

	/* Allocate memory for the nodes */
	tree->nodes = malloc(sizeof(RGTreeNode)*tree->numNodes);
	if(NULL == tree->nodes) {
		PrintError("RGTreeReadTree",
				"tree->nodes",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Read in the nodes */
	if(binaryInput == 0) {
		for(i=0;i<tree->numNodes;i++) {
			/* Read in the indexes and the number of entries */
			if(fscanf(fp, "%d %d %d",
						&tree->nodes[i].indexOne,
						&tree->nodes[i].indexTwo,
						&tree->nodes[i].numEntries)==EOF) {
				fprintf(stderr, "Error.  Could not read in indexes and numEntries.  Terminating!\n");
				exit(1);
			}

			if(tree->nodes[i].numEntries > 0) {
				/* Allocate memory for the positions and chromosomes */
				RGTreeNodeAllocate(&tree->nodes[i], tree->nodes[i].numEntries);

				/* Read in positions and chromosomes */
				for(j=0;j<tree->nodes[i].numEntries;j++) {
					if(fscanf(fp, "%d %d",
								&tempInt,
								&tree->nodes[i].positions[j])==EOF) {
						fprintf(stderr, "Error.  Could not read in position/chromosome %d.  Terminating!\n", j+1);
						exit(1);
					}
					tree->nodes[i].chromosomes[j] = tempInt;
				}
			}
			else {
				tree->nodes[i].positions = NULL;
				tree->nodes[i].chromosomes = NULL;
			}
		}
	}
	else {
		for(i=0;i<tree->numNodes;i++) {
			/* Read in the indexes and the number of entries */
			fread(&tempInt, sizeof(unsigned int), 1, fp);
			tempInt = ntohl(tempInt);
			tree->nodes[i].indexOne = tempInt;
			fread(&tempInt, sizeof(unsigned int), 1, fp);
			tempInt = ntohl(tempInt);
			tree->nodes[i].indexTwo = tempInt;
			fread(&tempInt, sizeof(unsigned int), 1, fp);
			tempInt = ntohl(tempInt);
			tree->nodes[i].numEntries = tempInt;

			if(tree->nodes[i].numEntries > 0) {
				/* Use a temp int array.  Expand when necessary. */
				if(tree->nodes[i].numEntries > tempIntArrLength) {
					/* Reallocate temp array */
					tempIntArrLength = tree->nodes[i].numEntries;
					tempIntArr = realloc(tempIntArr, sizeof(unsigned int)*tempIntArrLength);
					if(NULL == tempIntArr) {
						PrintError("RGTreeReadTree",
								"tempIntArr",
								"Could not reallocate memory",
								Exit,
								ReallocMemory);
					}
				}

				/* Allocate memory for the positions and chromosomes */
				RGTreeNodeAllocate(&tree->nodes[i], tree->nodes[i].numEntries);

				/* Read in positions */
				fread(tempIntArr, sizeof(unsigned int), tree->nodes[i].numEntries, fp);
				/* Copy over positions */
				for(j=0;j<tree->nodes[i].numEntries;j++) {
					tempIntArr[j] = ntohl(tempIntArr[j]);
					tree->nodes[i].positions[j] = tempIntArr[j];
				}
				/* Read in chromosomes */
				fread(tempIntArr, sizeof(unsigned int), tree->nodes[i].numEntries, fp);
				for(j=0;j<tree->nodes[i].numEntries;j++) {
					tempIntArr[j] = ntohl(tempIntArr[j]);
					tree->nodes[i].chromosomes[j] = (unsigned char)tempIntArr[j];
				}
			}
			else {
				tree->nodes[i].positions = NULL;
				tree->nodes[i].chromosomes = NULL;
			}
		}
		if(tempIntArrLength > 0) {
			free(tempIntArr);
			tempIntArr=NULL;
		}

	}
	return 1;
}

/* TODO */
void RGTreePrintHeader(FILE *fp, RGTree *tree, int binaryOutput)
{
	unsigned int numNodes=(unsigned int)(tree->numNodes);
	unsigned int subtractGapInformation=(unsigned int)(tree->gap<0)?(-1*tree->gap):0;
	unsigned int addGapInformation=(unsigned int)(tree->gap>=0)?(tree->gap):0;
	unsigned int matchLength=(unsigned int)tree->matchLength;
	unsigned int startChr=(unsigned int)tree->startChr;
	unsigned int startPos=(unsigned int)tree->startPos;
	unsigned int endChr=(unsigned int)tree->endChr;
	unsigned int endPos=(unsigned int)tree->endPos; 

	if(VERBOSE > 0) {
		fprintf(stderr, "Printing Header:%d,%d,%d,%d,%d,%d,%d,%d\n",
				numNodes,
				subtractGapInformation,
				addGapInformation,
				matchLength,
				startChr,
				startPos,
				endChr,
				endPos);
	}

	if(binaryOutput == 0) {
		fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
				numNodes,
				subtractGapInformation,
				addGapInformation,
				matchLength,
				startChr,
				startPos,
				endChr,
				endPos);
	}
	else {

		/* Big Endian/Little Endian conversion */
		numNodes=htonl(numNodes);
		subtractGapInformation=htonl(subtractGapInformation);
		addGapInformation=htonl(addGapInformation);
		matchLength=htonl(matchLength);
		startChr=htonl(startChr);
		startPos=htonl(startPos);
		endChr=htonl(endChr);
		endPos=htonl(endPos);

		/* Print Header */
		fwrite(&numNodes, sizeof(unsigned int), 1, fp);
		fwrite(&subtractGapInformation, sizeof(unsigned int), 1, fp);
		fwrite(&addGapInformation, sizeof(unsigned int), 1, fp);
		fwrite(&matchLength, sizeof(unsigned int), 1, fp);
		fwrite(&startChr, sizeof(unsigned int), 1, fp);
		fwrite(&startPos, sizeof(unsigned int), 1, fp);
		fwrite(&endChr, sizeof(unsigned int), 1, fp);
		fwrite(&endPos, sizeof(unsigned int), 1, fp);
	}
}

/* TODO */
void RGTreeReadHeader(FILE *fp, RGTree *tree, int binaryInput)
{
	unsigned int numNodes;
	unsigned int subtractGapInformation;
	unsigned int addGapInformation;
	unsigned int matchLength;
	unsigned int startChr;
	unsigned int startPos;
	unsigned int endChr;
	unsigned int endPos;

	/* Read in header */
	if(binaryInput == 0) {
		if(fscanf(fp, "%d %d %d %d %d %d %d %d",
					&numNodes,
					&subtractGapInformation,
					&addGapInformation,
					&matchLength,
					&startChr,
					&startPos,
					&endChr,
					&endPos)!=EOF) {
		}
		else {
			fprintf(stderr, "Error.  Could not read header file in RGTreeReadFromFile.  Terminating!\n");
			exit(1);
		}
	}
	else {
		if(fread(&numNodes, sizeof(unsigned int), 1, fp)!=EOF
				&& fread(&subtractGapInformation, sizeof(unsigned int), 1, fp)!=EOF 
				&& fread(&addGapInformation, sizeof(unsigned int), 1, fp)!=EOF 
				&& fread(&matchLength, sizeof(unsigned int), 1, fp)!=EOF 
				&& fread(&startChr, sizeof(unsigned int), 1, fp)!=EOF
				&& fread(&startPos, sizeof(unsigned int), 1, fp)!=EOF
				&& fread(&endChr, sizeof(unsigned int), 1, fp)!=EOF
				&& fread(&endPos, sizeof(unsigned int), 1, fp)!=EOF) {

			/* Big Endian/Little Endian conversion */
			numNodes = ntohl(numNodes);
			subtractGapInformation = ntohl(subtractGapInformation);
			addGapInformation = ntohl(addGapInformation);
			startChr = ntohl(startChr);
			startPos = ntohl(startPos);
			endChr = ntohl(endChr);
			endPos = ntohl(endPos);
			matchLength = ntohl(matchLength);
		}
		else {
			fprintf(stderr, "Error.  Could not read header file in RGTreeReadFromFile.  Terminating!\n");
			exit(1);
		}
	}

	if(VERBOSE > 0) {
		fprintf(stderr, "Printing Header:%d,%d,%d,%d,%d,%d,%d,%d\n",
				numNodes,
				subtractGapInformation,
				addGapInformation,
				matchLength,
				startChr,
				startPos,
				endChr,
				endPos);
	}

	/* Adjust header - see bpreprocess */
	tree->numNodes = numNodes;
	tree->gap = addGapInformation - subtractGapInformation;
	tree->matchLength = matchLength;
	tree->startChr = startChr;
	tree->startPos = startPos;
	tree->endChr = endChr;
	tree->endPos = endPos;

	/* Error checking ? */
}

/* TODO */
/* We will append the matches if matches have already been found */
int RGTreeGetMatches(RGTree *tree, unsigned int indexOne, unsigned int indexTwo, char direction, int offset, RGMatch *m) 
{
	unsigned int index, startIndex, i;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "RGTreeGetMatches.  Searching for indexes %d and %d with offset %d.\n",
				indexOne,
				indexTwo,
				offset);
	}

	/* Get the index of the tree */
	index = RGTreeGetIndex(tree, indexOne, indexTwo);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Found index:%d\n", index);
	}

	if(index >= 0 && index < tree->numNodes) {
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Will append...");
		}
		/* Copy over the matches */
		/* (Re)Allocate memory for the new matches */
		startIndex = m->numEntries;
		m->numEntries = m->numEntries + tree->nodes[index].numEntries;
		RGMatchReallocate(m, m->numEntries);

		/* Copy over */
		for(i=startIndex;i<m->numEntries;i++) {
			/* Adjust for reverse strand if necessary */
			if(direction == FORWARD) {
				/* Adjust for offset in position */
				m->positions[i] = tree->nodes[index].positions[i-startIndex] - offset;
				m->chromosomes[i] = tree->nodes[index].chromosomes[i-startIndex];
				m->strand[i] = direction;
			}
			else if(direction == REVERSE) {
				/* Adjust for offset in position */
				m->positions[i] = tree->nodes[index].positions[i-startIndex] + tree->matchLength - 1 - offset;
				m->chromosomes[i] = tree->nodes[index].chromosomes[i-startIndex];
				m->strand[i] = direction;
			}
			else {
				fprintf(stderr, "Error.  Could not understand direction [%c].  Terminating!\n", direction);
				exit(1);
			}
		}
	}

	return 1;
}

/* TODO */
unsigned int RGTreeGetIndexFromSequence(char *sequence, int matchLength)
{
	unsigned int number = 0;
	unsigned int i;
	unsigned int temp=0;

	for(i=0;i<matchLength;i++) {
		switch(sequence[i]) {
			case 'a':
				/* Same as zero */
				temp=0;
				break;
			case 'c':
				temp=1;
				break;
			case 'g':
				temp=2;
				break;
			case 't':
				temp=3;
				break;
			default:
				fprintf(stderr, "sequence:%s\n", sequence);
				fprintf(stderr, "Error: sequence not a proper character [%c] in RGTreeGetIndexFromSequence.  Terminating!\n", sequence[i]);
				exit(1);
				break;
		}
		number += temp*pow(ALPHABET_SIZE, matchLength-i-1);
	}
	return number;
}

/* TODO */
void RGTreeNodeCopy(RGTreeNode *src, RGTreeNode *dest) 
{
	dest->indexOne = src->indexOne;
	dest->indexTwo = src->indexTwo;
	dest->numEntries = src->numEntries;
	dest->positions = src->positions;
	dest->chromosomes = src->chromosomes;
}

int RGTreeNodeCompare(RGTreeNode *a, RGTreeNode *b) 
{
	if(a->indexOne < b->indexOne ||
			(a->indexOne == b->indexOne && a->indexTwo < b->indexTwo)) { 
		return -1;
	}
	else if(a->indexOne == b->indexOne && a->indexTwo == b->indexTwo) {
		return 0;
	}
	else {
		return 1;
	}
}

/* TODO */
void RGTreeQuickSortNode(RGTree *tree, unsigned int index, unsigned int low, unsigned int high)
{
	unsigned int i;
	unsigned int pivot = -1;
	unsigned int tempPos;
	unsigned char tempChr;

	if(low < high) {
		pivot = (low + high)/2;

		tempPos = tree->nodes[index].positions[pivot];
		tempChr = tree->nodes[index].chromosomes[pivot];
		tree->nodes[index].positions[pivot] = tree->nodes[index].positions[high];
		tree->nodes[index].chromosomes[pivot] = tree->nodes[index].chromosomes[high];
		tree->nodes[index].positions[high] = tempPos;
		tree->nodes[index].chromosomes[high] = tempChr;

		pivot = low;

		/* move elements */
		for(i=low;i<high;i++) {
			/* Compare with the pivot */
			if(tree->nodes[index].chromosomes[i] < tree->nodes[index].chromosomes[high]  ||
					(tree->nodes[index].chromosomes[i] == tree->nodes[index].chromosomes[high] && tree->nodes[index].positions[i] <= tree->nodes[index].positions[high])) {

				tempPos = tree->nodes[index].positions[i];
				tempChr = tree->nodes[index].chromosomes[i];
				tree->nodes[index].positions[i] = tree->nodes[index].positions[pivot];
				tree->nodes[index].chromosomes[i] = tree->nodes[index].chromosomes[pivot];
				tree->nodes[index].positions[pivot] = tempPos;
				tree->nodes[index].chromosomes[pivot] = tempChr;
				/* Increment the new pivot tree */
				pivot++;
			}
		}

		/* Move pivot element to correct place */
		tempPos = tree->nodes[index].positions[pivot];
		tempChr = tree->nodes[index].chromosomes[pivot];
		tree->nodes[index].positions[pivot] = tree->nodes[index].positions[high];
		tree->nodes[index].chromosomes[pivot] = tree->nodes[index].chromosomes[high];
		tree->nodes[index].positions[high] = tempPos;
		tree->nodes[index].chromosomes[high] = tempChr;

		/* Call recursively */
		if(pivot > 0) {
			RGTreeQuickSortNode(tree, index, low, pivot-1);
		}
		if(pivot < UINT_MAX) {
			RGTreeQuickSortNode(tree, index, pivot+1, high);
		}
	}
}

void RGTreeNodeAllocate(RGTreeNode *a, int numEntries) 
{
	a->numEntries = numEntries;
	a->positions = malloc(sizeof(unsigned int)*numEntries);
	if(NULL == a->positions) {
		PrintError("RGTreeNodeAllocate",
				"a->positions",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	a->chromosomes = malloc(sizeof(unsigned int)*numEntries);
	if(NULL == a->chromosomes) {
		PrintError("RGTreeNodeAllocate",
				"a->chromosomes",
				"Could not allocate memory",
				Exit,
				MallocMemory);

	}
}

void RGTreeNodeReallocate(RGTreeNode *a, int numEntries) 
{
	a->numEntries = numEntries;
	a->positions = realloc(a->positions, sizeof(unsigned int)*numEntries);
	if(NULL == a->positions) {
		PrintError("RGTreeNodeAllocate",
				"a->positions",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}
	a->chromosomes = realloc(a->chromosomes, sizeof(unsigned int)*numEntries);
	if(NULL == a->chromosomes) {
		PrintError("RGTreeNodeAllocate",
				"a->chromosomes",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}
}
