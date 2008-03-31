#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <string.h>
#include "BLibDefinitions.h"
#include "RGTree.h"

/* TODO */
/* We assume the root has already been initialized */
/* We could insert both forward and reverse compliment of the two sequences when creating the tree,
 * but this would double the space of each tree (practically too large).  It is easier just to 
 * double the number of searches.
 * */
int RGTreeInsert(RGTree *tree, char *sequenceOne, char *sequenceTwo, int matchLength, int chromosome, int position) 
{
	assert(tree->matchLength == matchLength);

	int indexOne = RGTreeGetIndexFromSequence(sequenceOne, matchLength);
	int indexTwo = RGTreeGetIndexFromSequence(sequenceTwo, matchLength);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "RGTreeInsert: inserting indexes %d and %d.\n",
				indexOne,
				indexTwo);
	}

	/* Create a new node */
	tree->numNodes++;

	/* Allocate memory for a new node */
	tree->nodes = (RGTreeNode*)realloc(tree->nodes, sizeof(RGTreeNode)*tree->numNodes);

	/* Initialize node */
	tree->nodes[tree->numNodes-1].numEntries = 1;
	tree->nodes[tree->numNodes-1].indexOne = indexOne;
	tree->nodes[tree->numNodes-1].indexTwo = indexTwo;
	tree->nodes[tree->numNodes-1].positions = NULL;
	tree->nodes[tree->numNodes-1].chromosomes = NULL;

	/* Allocate memory for the node members */
	tree->nodes[tree->numNodes-1].positions = (int*)malloc(sizeof(int));
	tree->nodes[tree->numNodes-1].chromosomes = (unsigned char*)malloc(sizeof(unsigned char));

	/* Copy over */
	tree->nodes[tree->numNodes-1].positions[tree->nodes[tree->numNodes-1].numEntries-1] = position;
	tree->nodes[tree->numNodes-1].chromosomes[tree->nodes[tree->numNodes-1].numEntries-1] = chromosome;

	return 1;
}

/* TODO */
void RGTreeCleanUpTree(RGTree *tree) 
{
	int i, j;
	int prevIndex = 0;

	if(tree->numNodes > 0) {

		/* Sort the nodes in the tree */
		if(VERBOSE >= 0) {
			fprintf(stderr, "Sorting nodes\n");
		}
		RGTreeQuickSortNodes(tree, 0, tree->numNodes-1, 0);
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
				int start = tree->nodes[prevIndex].numEntries;
				/* Allocate memory for chromosomes and positions */
				tree->nodes[prevIndex].numEntries += tree->nodes[i].numEntries;
				tree->nodes[prevIndex].positions = (int*)realloc(tree->nodes[prevIndex].positions, sizeof(int)*tree->nodes[prevIndex].numEntries);
				tree->nodes[prevIndex].chromosomes = (unsigned char*)realloc(tree->nodes[prevIndex].chromosomes, sizeof(unsigned char)*tree->nodes[prevIndex].numEntries);
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
		tree->nodes = (RGTreeNode*)realloc(tree->nodes, sizeof(RGTreeNode)*tree->numNodes);

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
void RGTreeQuickSortNodes(RGTree *tree, int low, int high, int numComplete) 
{
	int i;
	int pivot = -1;
	RGTreeNode temp;


	if(low < high) {
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
		RGTreeNodeCopy(&tree->nodes[pivot], &temp);
		RGTreeNodeCopy(&tree->nodes[high], &tree->nodes[pivot]);
		RGTreeNodeCopy(&temp, &tree->nodes[high]);

		/* Store where the pivot should be */
		pivot = low;

		/* move elements */
		for(i=low;i<high;i++) {
			/* Compare with the pivot */
			if(RGTreeNodeCompare(&tree->nodes[i], &tree->nodes[high]) <= 0) {
				/* Swap node at i with node at the new pivot tree */
				RGTreeNodeCopy(&tree->nodes[i], &temp);
				RGTreeNodeCopy(&tree->nodes[pivot], &tree->nodes[i]);
				RGTreeNodeCopy(&temp, &tree->nodes[pivot]);
				/* Increment the new pivot tree */
				pivot++;
			}
		}
		if(VERBOSE>=0) {
			fprintf(stderr, "\r%3.2lf percent complete",
					(100.0*numComplete)/tree->numNodes);
		}
		/* Move pivot element to correct place */
		RGTreeNodeCopy(&tree->nodes[pivot], &temp);
		RGTreeNodeCopy(&tree->nodes[high], &tree->nodes[pivot]);
		RGTreeNodeCopy(&temp, &tree->nodes[high]);

		/* Call recursively */
		RGTreeQuickSortNodes(tree, low, pivot-1, numComplete+1); 
		RGTreeQuickSortNodes(tree, pivot+1, high, pivot+1);
	}
}

/* TODO */
int RGTreeGetIndex(RGTree *tree,
		int indexOne,
		int indexTwo)
{
	int low = 0;
	int high = tree->numNodes-1;
	int mid;

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
	int i;

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
	int i;
	double total=0;

	/* Get memory used in each node */
	for(i=0;i<tree->numNodes;i++) {
		total += sizeof(RGTreeNode) + /* memory used by the node */
			sizeof(int)*tree->nodes[i].numEntries + /* memory used by positions */
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
void RGTreePrintTree(FILE *fp, RGTree *tree)
{
	int i, j;
	/*
	   unsigned int tempInt;
	   */

	/* Print header */
	RGTreePrintHeader(fp, tree);

	/* Print the nodes */
	for(i=0;i<tree->numNodes;i++) {
		fprintf(fp, "\n");
		/* Print the indexes and the number of entries */
		/*
		   tempInt = tree->nodes[i].indexOne;
		   tempInt = htonl(tempInt);
		   fwrite(&tempInt, sizeof(unsigned int), 1, fp);
		   tempInt = tree->nodes[i].indexTwo;
		   tempInt = htonl(tempInt);
		   fwrite(&tempInt, sizeof(unsigned int), 1, fp);
		   tempInt = tree->nodes[i].numEntries;
		   tempInt = htonl(tempInt);
		   fwrite(&tempInt, sizeof(unsigned int), 1, fp);
		   */
		fprintf(fp, "%d\t%d\t%d",
				tree->nodes[i].indexOne,
				tree->nodes[i].indexTwo,
				tree->nodes[i].numEntries);
		/* Print the entries */
		for(j=0;j<tree->nodes[i].numEntries;j++) {
			/*
			   tempInt = tree->nodes[i].chromosomes[j];
			   tempInt = htonl(tempInt);
			   fwrite(&tempInt, sizeof(unsigned int), 1, fp);
			   tempInt = tree->nodes[i].positions[j];
			   tempInt = htonl(tempInt);
			   fwrite(&tempInt, sizeof(unsigned int), 1, fp);
			   */
			fprintf(fp, "\t%d\t%d",
					tree->nodes[i].chromosomes[j],
					tree->nodes[i].positions[j]);
		}
	}
	fprintf(fp, "\n");
}

/* TODO */
int RGTreeReadTree(FILE *fp, RGTree *tree)
{
	int i, j;
	int tempInt;
	/*
	   unsigned int tempInt;
	   */
	/* Make sure memory of the root has been allocated */
	assert(tree!=NULL);
	assert(tree->nodes==NULL);

	/* Read in the header */
	RGTreeReadHeader(fp, tree);

	assert(tree->numNodes > 0);
	/* Allocate memory for the nodes */
	tree->nodes = (RGTreeNode*)malloc(sizeof(RGTreeNode)*tree->numNodes);

	/* Read in the nodes */
	for(i=0;i<tree->numNodes;i++) {
		/* Read in the indexes and the number of entries */
		/*
		   fread(&tempInt, sizeof(unsigned int), 1, fp);
		   tempInt = ntohl(tempInt);
		   tree->nodes[i].indexOne = tempInt;
		   fread(&tempInt, sizeof(unsigned int), 1, fp);
		   tempInt = ntohl(tempInt);
		   tree->nodes[i].indexTwo = tempInt;
		   fread(&tempInt, sizeof(unsigned int), 1, fp);
		   tempInt = ntohl(tempInt);
		   tree->nodes[i].numEntries = tempInt;
		   */
		if(fscanf(fp, "%d %d %d",
					&tree->nodes[i].indexOne,
					&tree->nodes[i].indexTwo,
					&tree->nodes[i].numEntries)==EOF) {
			fprintf(stderr, "Error.  Could not read in indexes and numEntries.  Terminating!\n");
			exit(1);
		}

		/* Allocate memory for the positions and chromosomes */
		tree->nodes[i].positions = (int*)malloc(sizeof(int)*tree->nodes[i].numEntries);
		tree->nodes[i].chromosomes = (unsigned char*)malloc(sizeof(unsigned char)*tree->nodes[i].numEntries);

		/* Read in positions and chromosomes */
		for(j=0;j<tree->nodes[i].numEntries;j++) {
			/*
			   fread(&tempInt, sizeof(unsigned int), 1, fp);
			   tempInt = ntohl(tempInt);
			   tree->nodes[i].chromosomes[j] = tempInt;
			   fread(&tempInt, sizeof(unsigned int), 1, fp);
			   tempInt = ntohl(tempInt);
			   tree->nodes[i].positions[j] = tempInt;
			   */
			if(fscanf(fp, "%d %d",
						&tempInt,
						&tree->nodes[i].positions[j])==EOF) {
				fprintf(stderr, "Error.  Could not read in position/chromosome %d.  Terminating!\n", j+1);
				exit(1);
			}
			tree->nodes[i].chromosomes[j] = tempInt;
		}
	}
	return 1;
}

/* TODO */
void RGTreePrintHeader(FILE *fp, RGTree *tree)
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

	/* Big Endian/Little Endian conversion */
	/*
	   numNodes=htonl(numNodes);
	   subtractGapInformation=htonl(subtractGapInformation);
	   addGapInformation=htonl(addGapInformation);
	   matchLength=htonl(matchLength);
	   startChr=htonl(startChr);
	   startPos=htonl(startPos);
	   endChr=htonl(endChr);
	   endPos=htonl(endPos);
	   */

	/* Print Header */
	/*
	   fwrite(&numNodes, sizeof(unsigned int), 1, fp);
	   fwrite(&subtractGapInformation, sizeof(unsigned int), 1, fp);
	   fwrite(&addGapInformation, sizeof(unsigned int), 1, fp);
	   fwrite(&matchLength, sizeof(unsigned int), 1, fp);
	   fwrite(&startChr, sizeof(unsigned int), 1, fp);
	   fwrite(&startPos, sizeof(unsigned int), 1, fp);
	   fwrite(&endChr, sizeof(unsigned int), 1, fp);
	   fwrite(&endPos, sizeof(unsigned int), 1, fp);
	   */
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

/* TODO */
void RGTreeReadHeader(FILE *fp, RGTree *tree)
{
	int numNodes;
	int subtractGapInformation;
	int addGapInformation;
	int matchLength;
	int startChr;
	int startPos;
	int endChr;
	int endPos;

	/* Read in header */
	/*
	   if(fread(&numNodes, sizeof(unsigned int), 1, fp)!=EOF
	   && fread(&subtractGapInformation, sizeof(unsigned int), 1, fp)!=EOF 
	   && fread(&addGapInformation, sizeof(unsigned int), 1, fp)!=EOF 
	   && fread(&matchLength, sizeof(unsigned int), 1, fp)!=EOF 
	   && fread(&startChr, sizeof(unsigned int), 1, fp)!=EOF
	   && fread(&startPos, sizeof(unsigned int), 1, fp)!=EOF
	   && fread(&endChr, sizeof(unsigned int), 1, fp)!=EOF
	   && fread(&endPos, sizeof(unsigned int), 1, fp)!=EOF) {
	   */
	if(fscanf(fp, "%d %d %d %d %d %d %d %d",
				&numNodes,
				&subtractGapInformation,
				&addGapInformation,
				&matchLength,
				&startChr,
				&startPos,
				&endChr,
				&endPos)!=EOF) {

		/* Big Endian/Little Endian conversion */
		/*
		   numNodes = ntohl(numNodes);
		   subtractGapInformation = ntohl(subtractGapInformation);
		   addGapInformation = ntohl(addGapInformation);
		   startChr = ntohl(startChr);
		   startPos = ntohl(startPos);
		   endChr = ntohl(endChr);
		   endPos = ntohl(endPos);
		   matchLength = ntohl(matchLength);
		   */

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
	else {
		fprintf(stderr, "Error.  Could not read header file in RGTreeReadFromFile.  Terminating!\n");
		exit(1);
	}
}

/* TODO */
/* We will append the matches if matches have already been found */
int RGTreeGetMatches(RGTree *tree, int indexOne, int indexTwo, char direction, RGMatch *m) 
{
	int index, startIndex, i;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "RGTreeGetMatches.  Searching for indexes %d and %d.\n",
				indexOne,
				indexTwo);
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
		m->positions = (int*)realloc(m->positions, sizeof(int)*(m->numEntries)); 
		m->chromosomes = (unsigned char*)realloc(m->chromosomes, sizeof(unsigned char)*(m->numEntries)); 
		m->strand = (char*)realloc(m->strand, sizeof(char)*(m->numEntries)); 

		/* Copy over */
		for(i=startIndex;i<m->numEntries;i++) {
			/* Adjust for reverse strand if necessary */
			if(direction == FORWARD) {
				m->positions[i] = tree->nodes[index].positions[i-startIndex];
				m->chromosomes[i] = tree->nodes[index].chromosomes[i-startIndex];
				m->strand[i] = direction;
			}
			else if(direction == REVERSE) {
				m->positions[i] = tree->nodes[index].positions[i-startIndex] + tree->matchLength - 1;
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
int RGTreeGetIndexFromSequence(char *sequence, int matchLength)
{
	int number = 0;
	int i;
	int temp=0;

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
void RGTreeQuickSortNode(RGTree *tree, int index, int low, int high)
{
	int i;
	int pivot = -1;
	int tempPos;
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
		RGTreeQuickSortNode(tree, index, low, pivot-1);
		RGTreeQuickSortNode(tree, index, pivot+1, high);
	}
}
