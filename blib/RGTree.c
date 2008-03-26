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
	int i, j, nextNode;
	RGNode *cur=NULL;
	RGNode *temp=NULL;
	int compOne;
	int compTwo;
	int lowerX=0;
	int midX;
	int upperX=pow(ALPHABET_SIZE, matchLength) - 1;
	int lowerY=0;
	int midY;
	int upperY=upperX;

	int sequenceOneIndex = GetIndexFromSequence(sequenceOne, matchLength);
	int sequenceTwoIndex = GetIndexFromSequence(sequenceTwo, matchLength);

	/* Initialize cur to be the root */
	if(tree->root == NULL) {
		tree->root = (void*)malloc(sizeof(RGNode));
		for(i=0;i<ALPHABET_SIZE;i++) {
			tree->root->next[i] = NULL;
		}
	}
	cur = tree->root;
	assert(cur!=NULL);

	/* Traverse the tree to insert */
	for(i=0;i<matchLength;i++) { 
		midX = (lowerX+upperX)/2;
		midY = (lowerY+upperY)/2;
		/* Four branches to new node (x,y) 
		 * B1: x'<x, y'<y
		 * B2: x'>=x, y'<y
		 * B3: x'<x, y'>=y
		 * B4: x'>=x, y'>=y
		 * */
		compOne = (sequenceOneIndex < midX)?0:1;
		compTwo = (sequenceTwoIndex < midY)?0:1;
		nextNode = compOne + 2*compTwo;

		switch(nextNode) {
			case 0:
				upperX = midX-1;
				upperY = midY-1;
				break;
			case 1:
				lowerX = midX;
				upperY = midY-1;
				break;
			case 2:
				upperX = midX-1;
				lowerY = midY;
				break;
			case 3:
				lowerX = midX;
				lowerY = midY;
				break;
			default:
				fprintf(stderr, "Error.  Could not understand nextNode [%d] in RGTree.c.  Terminating!\n", nextNode);
				exit(1);
				break;
		}

		/* Allocate a new node if it does not exist */
		if(cur->next[nextNode] == NULL) {
			/* If it is a leaf node, then allocate an array of integers */
			if(i==matchLength-1) {
				/* Create a leaf */
				cur->next[nextNode] = (void*)malloc(sizeof(RGLeaf));
				((RGLeaf*)cur->next[nextNode])->numEntries = 1;
				((RGLeaf*)cur->next[nextNode])->positions = (int*)malloc(sizeof(int));
				((RGLeaf*)cur->next[nextNode])->positions[0] = position;
				((RGLeaf*)cur->next[nextNode])->chromosomes = (char*)malloc(sizeof(char));
				((RGLeaf*)cur->next[nextNode])->chromosomes[0] = (int)chromosome;
			}
			else {
				/* Allocate the node */
				cur->next[nextNode] = malloc(sizeof(RGNode));
				assert(cur->next[nextNode] != NULL);
				/* Initialize it's pointers */
				temp = (RGNode*)cur->next[nextNode];
				for(j=0;j<ALPHABET_SIZE;j++) {
					temp->next[j] = NULL;
				}
			}
		}
		else if(i==matchLength-1) {
			/* Add to the leaf */
			((RGLeaf*)cur->next[nextNode])->numEntries++;
			((RGLeaf*)cur->next[nextNode])->positions = (int*)realloc(((RGLeaf*)cur->next[nextNode])->positions, sizeof(int)*(((RGLeaf*)cur->next[nextNode])->numEntries));
			((RGLeaf*)cur->next[nextNode])->chromosomes= (char*)realloc(((RGLeaf*)cur->next[nextNode])->chromosomes, sizeof(char)*(((RGLeaf*)cur->next[nextNode])->numEntries));
			((RGLeaf*)cur->next[nextNode])->positions[((RGLeaf*)cur->next[nextNode])->numEntries-1] = position;
			((RGLeaf*)cur->next[nextNode])->chromosomes[((RGLeaf*)cur->next[nextNode])->numEntries-1] = chromosome;
		}
		/* move to next one */
		cur = cur->next[nextNode];
	}

	return 1;
}

/* TODO */
void RGTreeDelete(RGTree *tree)
{
	RGTreeDeleteHelper((void**)&tree->root, 0, tree->depth);
}

/* TODO */
void RGTreeDeleteHelper(void **node, int curDepth, int depth)
{
	int i;
	void *temp;
	RGNode *tempRGNodePtr;
	int *tempIntPtr;
	if(NULL != (*node)) {
		/* At the leaves, we store an array of ints.  Make sure 
		 * we do not forget to delete them when we get to the 
		 * lowest depth.  But we also must be carefull not stop
		 * at the leaves.
		 * */
		if(curDepth<depth) {
			for(i=0;i<ALPHABET_SIZE;i++) {
				/* Delete the children */
				temp =(*((RGNode**)node))->next[i];
				if(curDepth == depth-1) {
					tempIntPtr= (int*)temp;
					RGTreeDeleteHelper((void**)(&tempIntPtr), curDepth+1, depth);
				}
				else {
					tempRGNodePtr = (RGNode*)temp;
					RGTreeDeleteHelper((void**)(&tempRGNodePtr), curDepth+1, depth);
				}
			}
			free((*((RGNode**)node)));
			(*((RGNode**)node)) = NULL;
		}
		else {
			/* Free the RGLeaf */
			free((*((RGLeaf**)node))->positions);
			free((*((RGLeaf**)node))->chromosomes);
			free(*((RGLeaf**)node));
		}
	}
}

/* TODO */
double RGTreeGetSize(RGTree *tree, int outputSize) 
{
	double total=RGTreeGetSizeHelper((void*)tree->root, 0, tree->depth, -1);
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
double RGTreeGetSizeHelper(void* node, int curDepth, int depth, int direction) 
{
	if(NULL!=node) {
		if(curDepth<depth) {
			return (RGTreeGetSizeHelper((void*)((RGNode*)node)->next[0], curDepth+1, depth, 0) + 
					RGTreeGetSizeHelper((void*)((RGNode*)node)->next[1], curDepth+1, depth, 1) + 
					RGTreeGetSizeHelper((void*)((RGNode*)node)->next[2], curDepth+1, depth, 2) + 
					RGTreeGetSizeHelper((void*)((RGNode*)node)->next[3], curDepth+1, depth, 3) + 
					(int)sizeof(RGNode)); /* size of the current node */
		}
		else {
			/* Get the size of the leaf */
			return ( (int)sizeof(RGLeaf*) + /* size of the current pointer */
					(int)sizeof(RGLeaf) + /* size of the current node */
					(int)sizeof(int)*(((RGLeaf*)node)->numEntries) + /* size of the int array of positions */
					(int)sizeof(char)*(((RGLeaf*)node)->numEntries) ); /* size of the char array of chromosomes */
		}
	}
	else {
		return 0;
	}
}

/* TODO */
void RGTreePrintTree(FILE *fp, RGTree *tree)
{

	RGTreePrintHeader(fp, tree);

	/* Print Tree */
	RGTreePrintTreeHelper(fp, (void*)tree->root, 0, tree->depth);
}

/* TODO */
void RGTreePrintTreeHelper(FILE *fp, void *node, int curDepth, int depth)
{
	int i;
	unsigned int tempInt;
	if(NULL!=node) {
		if(curDepth<depth) {
			for(i=0;i<ALPHABET_SIZE;i++) {
				if(((RGNode*)node)->next[i] == NULL) {
					/* The current node has no child at index i */
					tempInt = (unsigned int)(-1+RGT_ADD_TO_OUTPUT);
					if(VERBOSE >= DEBUG) {
						fprintf(stderr, "\t%d[%d]", tempInt+RGT_ADD_TO_INPUT, tempInt);
					}
					tempInt = htonl(tempInt);
					fwrite(&tempInt, sizeof(unsigned int), 1, fp);

				}
				else {
					/* The current node has a child at index i */
					fprintf(fp, "%d\t", i);
					tempInt = (unsigned int)(i+RGT_ADD_TO_OUTPUT);
					if(VERBOSE >= DEBUG) {
						fprintf(stderr, "\t%d[%d]", tempInt+RGT_ADD_TO_INPUT, tempInt);
					}
					tempInt = htonl(tempInt);
					fwrite(&tempInt, sizeof(unsigned int), 1, fp);
					RGTreePrintTreeHelper(fp, (void*)((RGNode*)node)->next[i], curDepth+1, depth);
				}
			}
		}
		else {
			/* We are at a leaf node */
			tempInt = (unsigned int)(-2+RGT_ADD_TO_OUTPUT);
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "\t%d[%d]", tempInt+RGT_ADD_TO_INPUT, tempInt);
			}
			tempInt = htonl(tempInt);
			fwrite(&tempInt, sizeof(unsigned int), 1, fp);

			/* Write contents */
			/* Write the number of entries */
			tempInt = (unsigned int)(((RGLeaf*)node)->numEntries+RGT_ADD_TO_OUTPUT);
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "\t%d[%d]", tempInt+RGT_ADD_TO_INPUT, tempInt);
			}
			tempInt = htonl(tempInt);
			fwrite(&tempInt, sizeof(unsigned int), 1, fp);
			/* Write the positions */
			for(i=0;i<((RGLeaf*)node)->numEntries;i++) {
				tempInt = (unsigned int)(((RGLeaf*)node)->positions[i]+RGT_ADD_TO_OUTPUT);
				if(VERBOSE >= DEBUG) {
					fprintf(stderr, "\t%d[%d]", tempInt+RGT_ADD_TO_INPUT, tempInt);
				}
				tempInt = htonl(tempInt);
				fwrite(&tempInt, sizeof(unsigned int), 1, fp);
			}
			/* Write the chromosomes */
			for(i=0;i<((RGLeaf*)node)->numEntries;i++) {
				tempInt = (unsigned int)(((RGLeaf*)node)->chromosomes[i]+RGT_ADD_TO_OUTPUT);
				if(VERBOSE >= DEBUG) {
					fprintf(stderr, "\t%d[%d]", tempInt+RGT_ADD_TO_INPUT, tempInt);
				}
				tempInt = htonl(tempInt);
				fwrite(&tempInt, sizeof(unsigned int), 1, fp);
			}
		}
	}
	else {
		return;
	}
}

/* TODO */
int GetIndexFromSequence(char *sequence, int matchLength)
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
				fprintf(stderr, "Error: sequence not a proper character [%c].  Terminating!\n", sequence[i]);
				exit(1);
				break;
		}
		number += temp*pow(ALPHABET_SIZE, matchLength-i-1);
	}
	return number;
}

/* TODO */
int RGTreeReadFromFile(RGTree *tree, FILE *fp)
{
	assert(tree->root==NULL);
	RGTreeReadHeader(fp, tree);

	/* Read in the tree */
	return RGTreeReadFromFileHelper((void*)tree->root, fp);
}

/* TODO */
int RGTreeReadFromFileHelper(RGNode *node, FILE *fp)
{
	int i, j;
	int tempInt;
	int numPositions;

	/* Basic idea is to read the tree from file */
	for(i=0;i<ALPHABET_SIZE;i++) { /* For each potential leaf */
		if(fread(&tempInt, sizeof(unsigned int), 1, fp)==EOF) {
			return EOF;
		}
		else {
			/* Add to input */
			tempInt += RGT_ADD_TO_INPUT;

			if(tempInt >= 0 && tempInt <=3) {
				/* Adding an internal node */
				assert(tempInt == i);
				/* Check if we have to allocate a new node */
				if(node->next[tempInt] == NULL) {
					node->next[tempInt] = malloc(sizeof(RGNode));
				}
				/* Initialize pointers */
				for(j=0;j<ALPHABET_SIZE;j++) {
					((RGNode*)node->next[tempInt])->next[j] = NULL;
				}
				/* Keep going deeper */
				if(RGTreeReadFromFileHelper((RGNode*)node->next[tempInt], fp)==EOF) {
					return EOF;
				}

			}
			else if(tempInt == -1) {
				node->next[tempInt] = NULL;
				/* Do nothing this is null */
			}
			else if(tempInt == -2) {
				/* We are at a leaf */

				/* Read the number of positions to be read */
				if(fread(&numPositions, sizeof(unsigned int), 1, fp)!=EOF) {
					numPositions += RGT_ADD_TO_INPUT;
					assert(numPositions > 0);
					/* Allocate memory for positions and chromosomes */
					assert(node->next[i]==NULL);
					node->next[i] = (void*)malloc(sizeof(RGLeaf));
					((RGLeaf*)node->next[i])->numEntries = numPositions;
					((RGLeaf*)node->next[i])->positions = (int*)malloc(sizeof(int)*numPositions);
					((RGLeaf*)node->next[i])->chromosomes = (char*)malloc(sizeof(char)*numPositions);

					/* Read in positions */
					if(fread(((RGLeaf*)node->next[i])->positions, sizeof(unsigned int), numPositions, fp)==EOF) {
						fprintf(stderr, "Error.  Trying to read %d positions.  Terminating!\n", numPositions);
						exit(1);
					}

					/* Read in chromosomes */
					if(fread(((RGLeaf*)node->next[i])->chromosomes, sizeof(unsigned int), numPositions, fp)==EOF) {
						fprintf(stderr, "Error.  Trying to read %d chromosomes.  Terminating!\n", numPositions);
						exit(1);
					}
					/* Phew */
				}
				else {
					fprintf(stderr, "Error.  Trying to read the number of positions.  Terminating!\n");
				}
			}
			else {
				fprintf(stderr, "Error.  In RGTreeReadFromFileHelper could not understand input state %d.  Terminating!\n", tempInt);
				exit(1);
			}
		}
	}
	return 1;
}

/* TODO */
void RGTreePrintHeader(FILE *fp, RGTree *tree)
{
	unsigned int subtractGapInformation=(unsigned int)(tree->gap<0)?(-1*tree->gap):0;
	unsigned int addGapInformation=(unsigned int)(tree->gap>=0)?(tree->gap):0;
	unsigned int depthOfRGTree=(unsigned int)tree->depth;
	unsigned int startChr=(unsigned int)tree->startChr;
	unsigned int startPos=(unsigned int)tree->startPos;
	unsigned int endChr=(unsigned int)tree->endChr;
	unsigned int endPos=(unsigned int)tree->endPos; 

	/* Big Endian/Little Endian conversion */
	subtractGapInformation=htonl(subtractGapInformation);
	addGapInformation=htonl(addGapInformation);
	depthOfRGTree=htonl(depthOfRGTree);
	startChr=htonl(startChr);
	startPos=htonl(startPos);
	endChr=htonl(endChr);
	endPos=htonl(endPos);

	/* Print Header */
	fwrite(&subtractGapInformation, sizeof(unsigned int), 1, fp);
	fwrite(&addGapInformation, sizeof(unsigned int), 1, fp);
	fwrite(&depthOfRGTree, sizeof(unsigned int), 1, fp);
	fwrite(&startChr, sizeof(unsigned int), 1, fp);
	fwrite(&startPos, sizeof(unsigned int), 1, fp);
	fwrite(&endChr, sizeof(unsigned int), 1, fp);
	fwrite(&endPos, sizeof(unsigned int), 1, fp);
}

/* TODO */
void RGTreeReadHeader(FILE *fp, RGTree *tree)
{
	int subtractGapInformation;
	int addGapInformation;
	int depthOfRGTree;
	int startChr;
	int startPos;
	int endChr;
	int endPos;

	/* Read in header */
	if(fread(&subtractGapInformation, sizeof(unsigned int), 1, fp)!=EOF 
			&& fread(&addGapInformation, sizeof(unsigned int), 1, fp)!=EOF 
			&& fread(&depthOfRGTree, sizeof(unsigned int), 1, fp)!=EOF 
			&& fread(&startChr, sizeof(unsigned int), 1, fp)!=EOF
			&& fread(&startPos, sizeof(unsigned int), 1, fp)!=EOF
			&& fread(&endChr, sizeof(unsigned int), 1, fp)!=EOF
			&& fread(&endPos, sizeof(unsigned int), 1, fp)!=EOF) {

		/* Big Endian/Little Endian conversion */
		subtractGapInformation = ntohl(subtractGapInformation);
		addGapInformation = ntohl(addGapInformation);
		startChr = ntohl(startChr);
		startPos = ntohl(startPos);
		endChr = ntohl(endChr);
		endPos = ntohl(endPos);
		depthOfRGTree = ntohl(depthOfRGTree);

		/* Adjust header - see bpreprocess */
		tree->gap = addGapInformation - subtractGapInformation;
		tree->depth = depthOfRGTree;
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
	int i, j, nextNode, startIndex;
	RGNode *cur=NULL;
	int compOne;
	int compTwo;
	int lowerX=0;
	int midX;
	int matchLength = tree->depth;
	int upperX=pow(ALPHABET_SIZE, matchLength) - 1;
	int lowerY=0;
	int midY;
	int upperY=upperX;

	cur = tree->root;
	assert(cur!=NULL);

	/* Traverse the tree to find the appropriate leaf */
	for(i=0;i<matchLength;i++) { 
		midX = (lowerX+upperX)/2;
		midY = (lowerY+upperY)/2;
		/* Four branches to new node (x,y) 
		 * B1: x'<x, y'<y
		 * B2: x'>=x, y'<y
		 * B3: x'<x, y'>=y
		 * B4: x'>=x, y'>=y
		 * */
		compOne = (indexOne < midX)?0:1;
		compTwo = (indexTwo < midY)?0:1;
		nextNode = compOne + 2*compTwo;

		switch(nextNode) {
			case 0:
				upperX = midX-1;
				upperY = midY-1;
				break;
			case 1:
				lowerX = midX;
				upperY = midY-1;
				break;
			case 2:
				upperX = midX-1;
				lowerY = midY;
				break;
			case 3:
				lowerX = midX;
				lowerY = midY;
				break;
			default:
				fprintf(stderr, "Error.  Could not understand nextNode [%d] in RGTree.c.  Terminating!\n", nextNode);
				exit(1);
				break;
		}

		if(cur->next[nextNode] == NULL) {
			/* If the next node does not exist, do not add any matches and return */
			return 0;
		}
		else if(i==matchLength-1) { 
			/* We are at a leaf, copy over the matches */

			/* (Re)Allocate memory for the new matches */
			startIndex = m->numEntries;
			m->numEntries = m->numEntries + ((RGLeaf*)cur->next[nextNode])->numEntries;
			m->positions = (int*)realloc(m->positions, sizeof(int)*(m->numEntries)); 
			m->chromosomes = (char*)realloc(m->chromosomes, sizeof(char)*(m->numEntries)); 
			m->strand = (char*)realloc(m->strand, sizeof(char)*(m->numEntries)); 

			/* Copy over */
			for(j=startIndex;j<m->numEntries;j++) {
				m->positions[j] = ((RGLeaf*)cur->next[nextNode])->positions[j-startIndex];
				m->chromosomes[j] = ((RGLeaf*)cur->next[nextNode])->chromosomes[j-startIndex];
				m->strand[j] = direction;
			}
		}
		/* move to next one */
		cur = cur->next[nextNode];
	}

	return 1;
}
