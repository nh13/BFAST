#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <string.h>
#include "BLibDefinitions.h"
#include "BError.h"
#include "RGTree.h"
#include "RGIndex.h"

/* TODO */
int RGIndexInsert(RGIndex *index, char *sequence, unsigned int matchLength, unsigned int chromosome, unsigned int position) 
{
	int i;
	assert(index->matchLength == matchLength);

	unsigned char *curIndex=NULL;
	int numChars = (int)ceil((2.0/8.0*matchLength)/sizeof(unsigned char));

	/* Allocate memory for the curIndex */
	curIndex = malloc(sizeof(unsigned char)*numChars);
	if(curIndex == NULL) {
		PrintError("RGIndexInsert",
				"curIndex",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	RGIndexGetIndexFromSequence(sequence, matchLength, curIndex);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "RGIndexInsert: inserting indexe:");
		for(i=0;i<numChars;i++) {
			fprintf(stderr, "\t%d", curIndex[i]);
		}
		fprintf(stderr, "\n");
	}

	/* Create a new node */
	index->numNodes++;

	/* Allocate memory for a new node */
	index->nodes = realloc(index->nodes, sizeof(RGIndexNode)*index->numNodes);
	if(NULL == index->nodes) {
		PrintError("RGIndexInsert",
				"index->nodes",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}

	/* Initialize node */
	index->nodes[index->numNodes-1].numEntries = 1;
	index->nodes[index->numNodes-1].positions = NULL;
	index->nodes[index->numNodes-1].chromosomes = NULL;
	index->nodes[index->numNodes-1].index = NULL;

	/* Allocate memory for the node members */
	index->nodes[index->numNodes-1].positions = malloc(sizeof(unsigned int));
	if(NULL == index->nodes[index->numNodes-1].positions) {
		PrintError("RGIndexInsert",
				"index->nodes[index->numNodes-1].positions",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	index->nodes[index->numNodes-1].chromosomes = malloc(sizeof(unsigned char));
	if(NULL == index->nodes[index->numNodes-1].chromosomes) {
		PrintError("RGIndexInsert",
				"index->nodes[index->numNodes-1].chromosomes",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	index->nodes[index->numNodes-1].index = malloc(sizeof(unsigned char)*numChars);
	if(NULL == index->nodes[index->numNodes-1].index) {
		PrintError("RGIndexInsert",
				"index->nodes[index->numNodes-1].index",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Copy over */
	index->nodes[index->numNodes-1].positions[0] = position;
	index->nodes[index->numNodes-1].chromosomes[0] = chromosome;
	for(i=0;i<numChars;i++) {
		index->nodes[index->numNodes-1].index[i] = curIndex[i];
	}

	/* Free memory */
	free(curIndex);

	return 1;
}

/* TODO */
void RGIndexCleanUpIndex(RGIndex *index) 
{
	unsigned int i, j;
	unsigned int prevIndex = 0;

	if(index->numNodes > 0) {

		/* Sort the nodes in the index */
		if(VERBOSE >= 0) {
			fprintf(stderr, "Sorting (this will speed up):\n");
		}
		RGIndexQuickSortNodes(index, 0, index->numNodes-1, 0);
		if(VERBOSE >= 0) {
			fprintf(stderr, "\nSorting complete\n");
		}

		/* Remove duplicates */
		/* Assumes the nodes are sorted */
		if(VERBOSE >= 1) {
			fprintf(stderr, "Merging\n");
		}
		/* This will remove duplicates in-place and in O(n) time */
		for(i=1;i<index->numNodes;i++) {
			if(RGIndexNodeCompare(&index->nodes[i], &index->nodes[prevIndex], index->matchLength)==0) {
				/* Append the entry to prevIndex nodes */
				unsigned int start = index->nodes[prevIndex].numEntries;
				/* Allocate memory for chromosomes and positions */
				index->nodes[prevIndex].numEntries += index->nodes[i].numEntries;
				index->nodes[prevIndex].positions = realloc(index->nodes[prevIndex].positions, sizeof(unsigned int)*index->nodes[prevIndex].numEntries);
				if(NULL == index->nodes[prevIndex].positions) {
					PrintError("RGIndexCleanUpIndex",
							"index->nodes[prevIndex].positions",
							"Could not reallocate memory",
							Exit,
							ReallocMemory);
				}
				index->nodes[prevIndex].chromosomes = realloc(index->nodes[prevIndex].chromosomes, sizeof(unsigned char)*index->nodes[prevIndex].numEntries);
				if(NULL == index->nodes[prevIndex].chromosomes) {
					PrintError("RGIndexCleanUpIndex",
							"index->nodes[prevIndex].chromosomes",
							"Could not reallocate memory",
							Exit,
							ReallocMemory);
				}
				/* Copy over chromosomes and positions */
				for(j=start;j<index->nodes[prevIndex].numEntries;j++) {
					index->nodes[prevIndex].positions[j] = index->nodes[i].positions[j-start];
					index->nodes[prevIndex].chromosomes[j] = index->nodes[i].chromosomes[j-start];
				}
				/* Free memory of the current node */
				free(index->nodes[i].positions);
				free(index->nodes[i].chromosomes);
				free(index->nodes[i].index);
				/* Nullify the current node */
				index->nodes[i].positions = NULL;
				index->nodes[i].chromosomes = NULL;
				index->nodes[i].index=NULL;
				index->nodes[i].numEntries=0;
			}
			else {
				prevIndex++;
				/* Move to prevIndex */
				if(prevIndex < i) {
					RGIndexNodeCopy(&index->nodes[i], &index->nodes[prevIndex], index->matchLength);
					/* Nullify the current node */
					index->nodes[i].positions = NULL;
					index->nodes[i].chromosomes = NULL;
					index->nodes[i].index = NULL;
					index->nodes[i].numEntries=0;
				}
			}
		}
		index->numNodes = prevIndex+1;
		/* Reallocate memory to reflect new number of nodes */
		index->nodes = realloc(index->nodes, sizeof(RGIndexNode)*index->numNodes);
		if(NULL == index->nodes) {
			PrintError("RGIndexCleanUpIndex",
					"index->nodes",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}

		/* Sort each node */
		for(i=0;i<index->numNodes;i++) {
			RGIndexQuickSortNode(index, i, 0, index->nodes[i].numEntries-1);
		}

		if(VERBOSE >= 1) {
			fprintf(stderr, "Merging complete\n");
		}
	}
}

/* TODO */
void RGIndexQuickSortNodes(RGIndex *index, unsigned int low, unsigned int high, unsigned int numComplete) 
{
	unsigned int i;
	unsigned int pivot = -1;
	RGIndexNode *temp=NULL;

	if(low < high) {
		/* Allocate temp */
		temp = malloc(sizeof(RGIndexNode));
		if(NULL == temp) {
			PrintError("RGIndexQuickSortNodes",
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
		RGIndexNodeCopy(&index->nodes[pivot], temp, index->matchLength);
		RGIndexNodeCopy(&index->nodes[high], &index->nodes[pivot], index->matchLength);
		RGIndexNodeCopy(temp, &index->nodes[high], index->matchLength);

		/* Store where the pivot should be */
		pivot = low;

		for(i=low;i<high;i++) {
			if(RGIndexNodeCompare(&index->nodes[i], &index->nodes[high], index->matchLength) <= 0) {
				/* Swap node at i with node at the new pivot index */
				RGIndexNodeCopy(&index->nodes[i], temp, index->matchLength);
				RGIndexNodeCopy(&index->nodes[pivot], &index->nodes[i], index->matchLength);
				RGIndexNodeCopy(temp, &index->nodes[pivot], index->matchLength);
				/* Increment the new pivot index */
				pivot++;
			}
		}
		if(VERBOSE>=0) {
			fprintf(stderr, "\r%3.1lf percent complete",
					(100.0*numComplete)/index->numNodes);
		}

		/* Move pivot element to correct place */
		RGIndexNodeCopy(&index->nodes[pivot], temp, index->matchLength);
		RGIndexNodeCopy(&index->nodes[high], &index->nodes[pivot], index->matchLength);
		RGIndexNodeCopy(temp, &index->nodes[high], index->matchLength);

		/* Free temp before recursing */
		free(temp);

		/* Call recursively */
		if(pivot > 0) {
			RGIndexQuickSortNodes(index, low, pivot-1, numComplete+1);
		}
		if(pivot < UINT_MAX) {
			RGIndexQuickSortNodes(index, pivot+1, high, pivot+1);
		}
	}

	if(VERBOSE>=0) {
		fprintf(stderr, "\r%3.1lf percent complete",
				100.0*high/index->numNodes); 
	}
}

/* TODO */
unsigned int RGIndexGetIndex(RGIndex *index,
		unsigned char *curIndex)
{
	long long int low = 0;
	long long int high = index->numNodes-1;
	long long int mid;
	int cmp;

	RGIndexNode cur;
	cur.index = curIndex;

	/* Binary search */
	while(low <= high) {
		mid = (low + high)/2;
		cmp = RGIndexNodeCompare(&cur, &index->nodes[mid], index->matchLength);
		if(cmp < 0) {
			high = mid-1;
		}
		else if(cmp == 0) {
			return mid;
		}
		else {
			low = mid+1;
		}
	}

	return -1;
}

/* TODO */
void RGIndexDelete(RGIndex *index)
{
	unsigned int i;

	/* Delete fields in the individual nodes */
	for(i=0;i<index->numNodes;i++) {
		free(index->nodes[i].positions);
		index->nodes[i].positions=NULL;
		free(index->nodes[i].chromosomes);
		index->nodes[i].chromosomes=NULL;
		free(index->nodes[i].index);
		index->nodes[i].index=NULL;
		index->nodes[i].numEntries=0;
	}

	/* Free nodes*/
	free(index->nodes);
	index->nodes=NULL;
	index->numNodes=0;
	index->matchLength=0;
	index->startChr=0;
	index->startPos=0;
	index->endChr=0;
	index->endPos=0;

}

/* TODO */
double RGIndexGetSize(RGIndex *index, int outputSize) 
{
	int i;
	double total=0;
	int numChars = (int)ceil((2.0/8.0*index->matchLength)/sizeof(unsigned char));

	/* Get memory used in each node */
	for(i=0;i<index->numNodes;i++) {
		total += sizeof(RGIndexNode) + /* memory used by the node */
			sizeof(unsigned int)*index->nodes[i].numEntries + /* memory used by positions */
			sizeof(unsigned char)*index->nodes[i].numEntries + /* memory used by chromosomes */
			sizeof(unsigned char)*numChars; /* memory used by the index */
	}

	/* Get memory of the RGIndex */
	total += sizeof(RGIndex);

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
void RGIndexPrintIndex(FILE *fp, RGIndex *index, int binaryOutput)
{
	unsigned int i, j;
	int numChars = (int)ceil((2.0/8.0*index->matchLength)/sizeof(unsigned char));
	unsigned int tempInt;
	unsigned int *tempIntArr=NULL;
	unsigned int tempIntArrLength=0;

	/* Print header */
	RGIndexPrintHeader(fp, index, binaryOutput);

	if(binaryOutput == 0) {
		/* Print the nodes */
		for(i=0;i<index->numNodes;i++) {
			fprintf(fp, "\n");
			/* Print the index */
			for(j=0;j<numChars;j++) {
				fprintf(fp, "%d\t", index->nodes[i].index[j]);
			}

			fprintf(fp, "%d",
					index->nodes[i].numEntries);
			/* Print the entries */
			for(j=0;j<index->nodes[i].numEntries;j++) {
				fprintf(fp, "\t%d\t%d",
						index->nodes[i].chromosomes[j],
						index->nodes[i].positions[j]);
			}
		}
		fprintf(fp, "\n");
	}
	else {

		/* We could just separate everything into huge arrays, but then we would double memory requirements when writing */

		/* Print the nodes */
		for(i=0;i<index->numNodes;i++) {
			/* Allocate memory if necessary */
			if(numChars > tempIntArrLength) {
				/* Reallocate memory */
				tempIntArrLength = numChars;
				tempIntArr = realloc(tempIntArr, sizeof(unsigned int)*tempIntArrLength);
				if(NULL == tempIntArr) {
					PrintError("RGIndexPrintIndex",
							"tempIntArr",
							"Could not reallocate memory",
							Exit,
							MallocMemory);
				}
			}

			/* print */
			for(j=0;j<numChars;j++) {
				/* Copy over */
				tempIntArr[j] = (unsigned int)index->nodes[i].index[j];
				tempIntArr[j] = htonl(tempIntArr[j]);
			}
			fwrite(tempIntArr, sizeof(unsigned int), numChars, fp);

			/* Print the number of entries */
			tempInt = index->nodes[i].numEntries;
			tempInt = htonl(tempInt);
			fwrite(&tempInt, sizeof(unsigned int), 1, fp);

			/* Print the entries */
			if(index->nodes[i].numEntries > 0) {
				/* Allocate memory if necessary */
				if(index->nodes[i].numEntries > tempIntArrLength) {
					/* Reallocate memory */
					tempIntArrLength = index->nodes[i].numEntries;
					tempIntArr = realloc(tempIntArr, sizeof(unsigned int)*tempIntArrLength);
					if(NULL == tempIntArr) {
						PrintError("RGIndexPrintIndex",
								"tempIntArr",
								"Could not reallocate memory",
								Exit,
								MallocMemory);
					}
				}

				/* Print the positions */ 
				for(j=0;j<index->nodes[i].numEntries;j++) {
					tempIntArr[j] = (unsigned int)index->nodes[i].positions[j];
					tempIntArr[j] = htonl(tempIntArr[j]);
				}
				fwrite(tempIntArr, sizeof(unsigned int), index->nodes[i].numEntries, fp);

				/* Print the chromosomes */
				for(j=0;j<index->nodes[i].numEntries;j++) {
					tempIntArr[j] = (unsigned int)index->nodes[i].chromosomes[j];
					tempIntArr[j] = htonl(tempIntArr[j]);
				}
				fwrite(tempIntArr, sizeof(unsigned int), index->nodes[i].numEntries, fp);
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
int RGIndexReadIndex(FILE *fp, RGIndex *index, int binaryInput)
{
	unsigned int i, j;
	unsigned int tempInt;
	unsigned int *tempIndex=NULL; /* Use this to store from file */
	unsigned int *tempIntArr=NULL;
	unsigned int tempIntArrLength=0;
	unsigned int numChars = -1;

	/* Make sure memory of the root has been allocated */
	assert(index!=NULL);
	assert(index->nodes==NULL);

	/* Read in the header */
	RGIndexReadHeader(fp, index, binaryInput);

	/* Update the structure of the indexes */
	numChars = (int)ceil((2.0/8.0*index->matchLength)/sizeof(unsigned char));

	assert(index->numNodes > 0);
	/* Allocate memory for the nodes */
	index->nodes = malloc(sizeof(RGIndexNode)*index->numNodes);
	if(NULL == index->nodes) {
		PrintError("RGIndexReadIndex",
				"index->nodes",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(index->nodes!=NULL);
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Will read %d index nodes.\n",
				index->numNodes);
	}

	if(binaryInput == 0) {
		/* Preallocate as much as possible */
		for(i=0;i<index->numNodes;i++) {
			/* Allocate memory for the index */
			index->nodes[i].index = malloc(sizeof(unsigned char)*numChars);
			if(NULL == index->nodes[i].index) {
				PrintError("RGIndexReadIndex",
						"index->nodes[i].index",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}

		/* Read in the nodes */
		for(i=0;i<index->numNodes;i++) {
			/* Read in the index */
			for(j=0;j<numChars;j++) {
				if(fscanf(fp, "%d", &tempInt)==EOF) {
					PrintError("RGIndexReadIndex",
							"tempInt",
							"Could not read in index of a node",
							Exit,
							EndOfFile);
				}
				index->nodes[i].index[j] = tempInt;
			}

			/* Read in the number of entries */
			if(fscanf(fp, "%d",
						&index->nodes[i].numEntries)==EOF) {
				PrintError("RGIndexReadIndex",
						"index->nodes[i].numEntries",
						"Could not read in the number of entries",
						Exit,
						EndOfFile);
			}

			/* Allocate memory for the positions and chromosomes */
			if(index->nodes[i].numEntries > 0) {
				index->nodes[i].positions = malloc(sizeof(unsigned int)*index->nodes[i].numEntries);
				if(NULL == index->nodes[i].positions) {
					PrintError("RGIndexReadIndex",
							"index->nodes[i].positions",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}
				index->nodes[i].chromosomes = malloc(sizeof(unsigned char)*index->nodes[i].numEntries);
				if(NULL == index->nodes[i].chromosomes) {
					PrintError("RGIndexReadIndex",
							"index->nodes[i].chromosomes",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}

				/* Read in positions and chromosomes */
				for(j=0;j<index->nodes[i].numEntries;j++) {
					if(fscanf(fp, "%d %d",
								&tempInt,
								&index->nodes[i].positions[j])==EOF) {
						PrintError("RGIndexReadIndex",
								NULL,
								"Could not read in chromosome and position",
								Exit,
								EndOfFile);
					}
					index->nodes[i].chromosomes[j] = tempInt;
				}
			}
			else {
				index->nodes[i].positions = NULL;
				index->nodes[i].chromosomes = NULL;
			}
		}
	}
	else {
		/* This will hold the index temporarily */
		tempIndex = malloc(sizeof(unsigned int)*numChars);
		if(NULL == tempIndex) {
			PrintError("RGIndexReadIndex",
					"tempIndex",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}

		/* Preallocate as much as possible */
		for(i=0;i<index->numNodes;i++) {
			/* Allocate memory for the index */
			index->nodes[i].index = malloc(sizeof(unsigned char)*numChars);
			if(NULL == index->nodes[i].index) {
				PrintError("RGIndexReadIndex",
						"index->nodes[i].index",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}

		/* Read in the nodes */
		for(i=0;i<index->numNodes;i++) {
			/* Read in the index */
			fread(tempIndex, sizeof(unsigned int), numChars, fp);
			/* Copy over index */
			for(j=0;j<numChars;j++) {
				tempIndex[j] = ntohl(tempIndex[j]);
				index->nodes[i].index[j] = (unsigned char)tempIndex[j];
				if(VERBOSE >= DEBUG) {
					fprintf(stderr, "%d\t",
							(int)index->nodes[i].index[j]);
				}
			}
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "\n");
			}

			/* Read in the number of entries */
			fread(&tempInt, sizeof(unsigned int), 1, fp);
			tempInt = ntohl(tempInt);
			index->nodes[i].numEntries = tempInt;

			if(index->nodes[i].numEntries > 0) {
				/* Use a temp int array.  Expand when necessary. */
				if(index->nodes[i].numEntries > tempIntArrLength) {
					/* Reallocate temp array */
					tempIntArrLength = index->nodes[i].numEntries;
					tempIntArr = realloc(tempIntArr, sizeof(unsigned int)*tempIntArrLength);
					if(NULL == tempIntArr) {
						PrintError("RGIndexReadIndex",
								"tempIntArr",
								"Could not reallocate memory",
								Exit,
								ReallocMemory);
					}
				}

				/* Allocate memory for the positions and chromosomes */
				index->nodes[i].positions = malloc(sizeof(unsigned int)*index->nodes[i].numEntries);
				if(NULL == index->nodes[i].positions) {
					PrintError("RGIndexReadIndex",
							"index->nodes[i].positions",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}
				index->nodes[i].chromosomes = malloc(sizeof(unsigned char)*index->nodes[i].numEntries);
				if(NULL == index->nodes[i].chromosomes) {
					PrintError("RGIndexReadIndex",
							"index->nodes[i].chromosomes",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}

				/* Read in positions */
				fread(tempIntArr, sizeof(unsigned int), index->nodes[i].numEntries, fp);
				/* Copy over positions */
				for(j=0;j<index->nodes[i].numEntries;j++) {
					tempIntArr[j] = ntohl(tempIntArr[j]);
					index->nodes[i].positions[j] = tempIntArr[j];
				}
				/* Read in chromosomes */
				fread(tempIntArr, sizeof(unsigned int), index->nodes[i].numEntries, fp);
				for(j=0;j<index->nodes[i].numEntries;j++) {
					tempIntArr[j] = ntohl(tempIntArr[j]);
					index->nodes[i].chromosomes[j] = (unsigned char)tempIntArr[j];
				}
			}
			else {
				index->nodes[i].positions = NULL;
				index->nodes[i].chromosomes = NULL;
			}

		}
		/* Free memory */
		if(tempIntArrLength > 0) {
			free(tempIntArr);
			tempIntArr=NULL;
			tempIntArrLength=0;
		}
		free(tempIndex);
	}
	return 1;
}

/* TODO */
void RGIndexPrintHeader(FILE *fp, RGIndex *index, int binaryOutput)
{
	unsigned int numNodes=(unsigned int)(index->numNodes);
	unsigned int matchLength=(unsigned int)index->matchLength;
	unsigned int startChr=(unsigned int)index->startChr;
	unsigned int startPos=(unsigned int)index->startPos;
	unsigned int endChr=(unsigned int)index->endChr;
	unsigned int endPos=(unsigned int)index->endPos; 

	if(VERBOSE > 0) {
		fprintf(stderr, "Printing Header:%d,%d,%d,%d,%d,%d\n",
				numNodes,
				matchLength,
				startChr,
				startPos,
				endChr,
				endPos);
	}
	if(binaryOutput == 0) {
		fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d",
				numNodes,
				matchLength,
				startChr,
				startPos,
				endChr,
				endPos);
	}
	else {

		/* Big Endian/Little Endian conversion */
		numNodes=htonl(numNodes);
		matchLength=htonl(matchLength);
		startChr=htonl(startChr);
		startPos=htonl(startPos);
		endChr=htonl(endChr);
		endPos=htonl(endPos);

		/* Print Header */
		fwrite(&numNodes, sizeof(unsigned int), 1, fp);
		fwrite(&matchLength, sizeof(unsigned int), 1, fp);
		fwrite(&startChr, sizeof(unsigned int), 1, fp);
		fwrite(&startPos, sizeof(unsigned int), 1, fp);
		fwrite(&endChr, sizeof(unsigned int), 1, fp);
		fwrite(&endPos, sizeof(unsigned int), 1, fp);
	}
}

/* TODO */
void RGIndexReadHeader(FILE *fp, RGIndex *index, int binaryInput)
{
	unsigned int numNodes;
	unsigned int matchLength;
	unsigned int startChr;
	unsigned int startPos;
	unsigned int endChr;
	unsigned int endPos;

	/* Read in header */
	if(binaryInput == 0) {
		if(fscanf(fp, "%d %d %d %d %d %d",
					&numNodes,
					&matchLength,
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
		if(fread(&numNodes, sizeof(unsigned int), 1, fp)==EOF
				|| fread(&matchLength, sizeof(unsigned int), 1, fp)==EOF 
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
		numNodes = ntohl(numNodes);
		startChr = ntohl(startChr);
		startPos = ntohl(startPos);
		endChr = ntohl(endChr);
		endPos = ntohl(endPos);
		matchLength = ntohl(matchLength);
	}

	if(VERBOSE > 0) {
		fprintf(stderr, "Printing Header:%d,%d,%d,%d,%d,%d\n",
				numNodes,
				matchLength,
				startChr,
				startPos,
				endChr,
				endPos);
	}

	/* Adjust header - see bpreprocess */
	index->numNodes = numNodes;
	index->matchLength = matchLength;
	index->startChr = startChr;
	index->startPos = startPos;
	index->endChr = endChr;
	index->endPos = endPos;

	/* Error checking ? */
}

/* TODO */
/* We will append the matches if matches have already been found */
int RGIndexGetMatches(RGIndex *index, unsigned char *curIndex, char direction, RGMatch *m)
{
	unsigned int i;
	unsigned int startIndex=-1;
	unsigned int nodeIndex=-1;
	unsigned int numChars = (int)ceil((2.0/8.0*index->matchLength)/sizeof(unsigned char));

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "RGIndexGetMatches.  Searching for index:");
		for(i=0;i<numChars;i++) {
			fprintf(stderr, "\t%d", curIndex[i]);
		}
		fprintf(stderr, "\n");
	}

	/* Get the index of the index */
	nodeIndex = RGIndexGetIndex(index, curIndex);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Found index:%d\n", nodeIndex);
	}

	if(nodeIndex >= 0 && nodeIndex < index->numNodes) {
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Will append...\n");
		}
		/* Copy over the matches */
		/* (Re)Allocate memory for the new matches */
		startIndex = m->numEntries;
		m->numEntries = m->numEntries + index->nodes[nodeIndex].numEntries;
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
			/* Adjust for reverse strand if necessary */
			if(direction == FORWARD) {
				m->positions[i] = index->nodes[nodeIndex].positions[i-startIndex];
				m->chromosomes[i] = index->nodes[nodeIndex].chromosomes[i-startIndex];
				m->strand[i] = direction;
			}
			else if(direction == REVERSE) {
				m->positions[i] = index->nodes[nodeIndex].positions[i-startIndex]+index->matchLength-1;
				m->chromosomes[i] = index->nodes[nodeIndex].chromosomes[i-startIndex];
				m->strand[i] = direction;
			}
			else {
				PrintError("RGIndexGetMatches",
						"direction",
						"Could not understand direction",
						Exit,
						OutOfRange);
			}

		}
	}

	return 1;
}

/* TODO */
void RGIndexNodeCopy(RGIndexNode *src, RGIndexNode *dest, int matchLength) 
{
	dest->numEntries = src->numEntries;
	dest->positions = src->positions;
	dest->chromosomes = src->chromosomes;
	dest->index = src->index;
}

int RGIndexNodeCompare(RGIndexNode *a, RGIndexNode *b, int matchLength) 
{
	int i;
	int numChars = (int)ceil((2.0/8.0*matchLength)/sizeof(unsigned char));

	/* Compare using the last unsigned char to the first unsigned char */
	for(i=numChars-1;i>=0;i--) {
		if(a->index[i] < b->index[i]) {
			return -1;
		}
		else if(a->index[i] > b->index[i]) {
			return 1;
		}
		/* otherwise they are equal */
	}
	return 0;
}

/* TODO */
void RGIndexQuickSortNode(RGIndex *index, int curNode, int low, int high)
{
	unsigned int i;
	unsigned int pivot = -1;
	unsigned int tempPos;
	unsigned char tempChr;

	if(low < high) {
		pivot = (low + high)/2;

		tempPos = index->nodes[curNode].positions[pivot];
		tempChr = index->nodes[curNode].chromosomes[pivot];
		index->nodes[curNode].positions[pivot] = index->nodes[curNode].positions[high];
		index->nodes[curNode].chromosomes[pivot] = index->nodes[curNode].chromosomes[high];
		index->nodes[curNode].positions[high] = tempPos;
		index->nodes[curNode].chromosomes[high] = tempChr;

		pivot = low;

		/* move elements */
		for(i=low;i<high;i++) {
			/* Compare with the pivot */
			if(index->nodes[curNode].chromosomes[i] < index->nodes[curNode].chromosomes[high] ||
					(index->nodes[curNode].chromosomes[i] == index->nodes[curNode].chromosomes[high] && index->nodes[curNode].positions[i] <= index->nodes[curNode].positions[high])) {

				tempPos = index->nodes[curNode].positions[i];
				tempChr = index->nodes[curNode].chromosomes[i];
				index->nodes[curNode].positions[i] = index->nodes[curNode].positions[pivot];
				index->nodes[curNode].chromosomes[i] = index->nodes[curNode].chromosomes[pivot];
				index->nodes[curNode].positions[pivot] = tempPos;
				index->nodes[curNode].chromosomes[pivot] = tempChr;
				/* Increment the new pivot index */
				pivot++;
			}
		}
		/* Move pivot element to correct place */
		tempPos = index->nodes[curNode].positions[pivot];
		tempChr = index->nodes[curNode].chromosomes[pivot];
		index->nodes[curNode].positions[pivot] = index->nodes[curNode].positions[high];
		index->nodes[curNode].chromosomes[pivot] = index->nodes[curNode].chromosomes[high];
		index->nodes[curNode].positions[high] = tempPos;
		index->nodes[curNode].chromosomes[high] = tempChr;

		/* Call recursively */
		if(pivot > 0) {
			RGIndexQuickSortNode(index, curNode, low, pivot-1);
		}
		if(pivot < UINT_MAX) {
			RGIndexQuickSortNode(index, curNode, pivot+1, high);
		}
	}
}

/* TODO */
void RGIndexGetIndexFromSequence(char *sequence, int matchLength, unsigned char *curIndex)
{
	int i, j, cur;
	int numChars = (int)ceil((2.0/8.0*matchLength)/sizeof(unsigned char));
	int temp=0;
	int number = 0;

	assert(curIndex!=NULL);
	assert(ALPHABET_SIZE == 4); /* Implemented only for a c g t */

	cur=matchLength-1; /* Start from the right most */
	for(i=numChars-1;i>=0;i--) { /* For each unsigned character/byte */
		number=0;
		/* each unsigned character is one byte (8-bits) so can hold 4 letters of sequence */
		for(j=0;cur>=0 && j<ALPHABET_SIZE;j++) {
			switch(sequence[cur]) {
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
					PrintError("RGIndexGetIndexFromSequence",
							"sequence[cur]",
							"Sequence is not a proper character",
							Exit,
							OutOfRange);
					break;
			}
			number += temp*pow(ALPHABET_SIZE, j);
			cur--;
		}
		/* Store number */
		curIndex[i]=number;
	}
}

