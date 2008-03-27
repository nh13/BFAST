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
	double midX;
	int upperX=pow(ALPHABET_SIZE, matchLength) - 1;
	int lowerY=0;
	double midY;
	int upperY=upperX;

	assert(2*matchLength == tree->depth);

	int sequenceOneIndex = GetIndexFromSequence(sequenceOne, matchLength);
	int sequenceTwoIndex = GetIndexFromSequence(sequenceTwo, matchLength);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "RGTreeInsert: inserting indexes %d and %d.\n",
				sequenceOneIndex,
				sequenceTwoIndex);
	}

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
	for(i=0;i<2*matchLength;i++) { 
		midX = (lowerX+upperX)/2.0;
		midY = (lowerY+upperY)/2.0;
		/* Four branches to new node (x,y) from search indexes (x',y') 
		 * B1: x'<x, y'<y
		 * B2: x'>x, y'<y
		 * B3: x'<x, y'>y
		 * B4: x'>x, y'>y
		 * */
		compOne = (sequenceOneIndex < midX)?0:1;
		compTwo = (sequenceTwoIndex < midY)?0:1;
		nextNode = compOne + 2*compTwo;

		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "depth:%d,lowerX:%d,midX:%.1f,upperX:%d,lowerY:%d,midY:%.1f,upperY:%d,nextNode:%d\n",
					i+1,
					lowerX,
					midX,
					upperX,
					lowerY,
					midY,
					upperY,
					nextNode);
		}

		/* Update lower and upper */
		switch(nextNode) {
			case 0:
				upperX = (int)midX;
				upperY = (int)midY;
				break;
			case 1:
				lowerX = (int)(midX+0.5);
				upperY = (int)midY;
				break;
			case 2:
				upperX = (int)midX;
				lowerY = (int)(midY+0.5);
				break;
			case 3:
				lowerX = (int)(midX+0.5);
				lowerY = (int)(midY+0.5);
				break;
			default:
				fprintf(stderr, "Error.  Could not understand nextNode [%d] in RGTree.c.  Terminating!\n", nextNode);
				exit(1);
				break;
		}

		/* Allocate a new node if it does not exist */
		if(cur->next[nextNode] == NULL) {
			/* If it is a leaf node, then allocate an array of integers */
			if(i==2*matchLength-1) {
				if(VERBOSE >= DEBUG) {
					fprintf(stderr, "Creating a new leaf\n");
				}
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
		else if(i==2*matchLength-1) {
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "Adding to the leaf\n");
			}
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

	/* Free root */
	free(tree->root);
	tree->root=NULL;
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

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Tree depth:%d\n", tree->depth);
	}
	/* Print header */
	RGTreePrintHeader(fp, tree);

	/* Print Tree */
	RGTreePrintTreeHelper(fp, (void*)tree->root, 0, tree->depth);
}

/* TODO */
void RGTreePrintTreeHelper(FILE *fp, void *node, int curDepth, int depth)
{
	int i;
	/*
	   unsigned int tempInt;
	   */
	int tempInt;
	if(NULL!=node) {
		if(curDepth == depth) {
			assert(curDepth == depth);
			/* We are at a leaf node */
			/*
			   tempInt = (unsigned int)(-2+RGT_ADD_TO_OUTPUT);
			   */
			tempInt = -2+RGT_ADD_TO_OUTPUT;
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "\t%d[%d]{%d}", tempInt+RGT_ADD_TO_INPUT, tempInt, curDepth);
			}
			/*
			   tempInt = htonl(tempInt);
			   */
			/*
			   fwrite(&tempInt, sizeof(unsigned int), 1, fp);
			   */
			fprintf(fp, "\t%d", tempInt);

			/* Write contents */
			/* Write the number of entries */
			/*
			   tempInt = (unsigned int)(((RGLeaf*)node)->numEntries+RGT_ADD_TO_OUTPUT);
			   */
			tempInt = ((RGLeaf*)node)->numEntries+RGT_ADD_TO_OUTPUT;
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "\t%d[%d]{%d}", tempInt+RGT_ADD_TO_INPUT, tempInt, curDepth);
			}
			/*
			   tempInt = htonl(tempInt);
			   */
			/*
			   fwrite(&tempInt, sizeof(unsigned int), 1, fp);
			   */
			fprintf(fp, "\t%d", tempInt);
			/* Write the positions */
			for(i=0;i<((RGLeaf*)node)->numEntries;i++) {
				/*
				   tempInt = (unsigned int)(((RGLeaf*)node)->positions[i]+RGT_ADD_TO_OUTPUT);
				   */
				tempInt = ((RGLeaf*)node)->positions[i]+RGT_ADD_TO_OUTPUT;
				if(VERBOSE >= DEBUG) {
					fprintf(stderr, "\t%d[%d]{%d}", tempInt+RGT_ADD_TO_INPUT, tempInt, curDepth);
				}
				/*
				   tempInt = htonl(tempInt);
				   */
				/*
				   fwrite(&tempInt, sizeof(unsigned int), 1, fp);
				   */
				fprintf(fp, "\t%d", tempInt);
			}
			/* Write the chromosomes */
			for(i=0;i<((RGLeaf*)node)->numEntries;i++) {
				/*
				   tempInt = (unsigned int)(((RGLeaf*)node)->chromosomes[i]+RGT_ADD_TO_OUTPUT);
				   */
				tempInt = ((int)((RGLeaf*)node)->chromosomes[i])+RGT_ADD_TO_OUTPUT;
				if(VERBOSE >= DEBUG) {
					fprintf(stderr, "\t%d[%d]{%d}", tempInt+RGT_ADD_TO_INPUT, tempInt, curDepth);
				}
				/*
				   tempInt = htonl(tempInt);
				   */
				/*
				   fwrite(&tempInt, sizeof(unsigned int), 1, fp);
				   */
				fprintf(fp, "\t%d", tempInt);
			}
		}
		else {
			assert(curDepth < depth);
			for(i=0;i<ALPHABET_SIZE;i++) {
				if(((RGNode*)node)->next[i] == NULL) {
					/* The current node has no child at index i */
					/*
					   tempInt = (unsigned int)(-1+RGT_ADD_TO_OUTPUT);
					   */
					tempInt = -1+RGT_ADD_TO_OUTPUT;
					if(VERBOSE >= DEBUG) {
						fprintf(stderr, "\t%d[%d]{%d}", tempInt+RGT_ADD_TO_INPUT, tempInt, curDepth);
					}
					/*
					   tempInt = htonl(tempInt);
					   */
					/*
					   fwrite(&tempInt, sizeof(unsigned int), 1, fp);
					   */
					fprintf(fp, "\t%d", tempInt);

				}
				else {
					/* The current node has a child at index i */
					assert(i+RGT_ADD_TO_OUTPUT>=0);
					/*
					   tempInt = (unsigned int)(i+RGT_ADD_TO_OUTPUT);
					   */
					tempInt = i+RGT_ADD_TO_OUTPUT;
					if(VERBOSE >= DEBUG) {
						fprintf(stderr, "\t%d[%d]{%d}", tempInt+RGT_ADD_TO_INPUT, tempInt, curDepth);
					}
					/*
					   tempInt = htonl(tempInt);
					   */
					/*
					   fwrite(&tempInt, sizeof(unsigned int), 1, fp);
					   */
					fprintf(fp, "\t%d", tempInt);
					RGTreePrintTreeHelper(fp, (void*)((RGNode*)node)->next[i], curDepth+1, depth);
				}
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
	int i;
	/* Make sure memory of the root has been allocated */
	assert(tree->root!=NULL);
	for(i=0;i<ALPHABET_SIZE;i++) {/*overkill*/
		assert(tree->root->next[i] == NULL);
	}

	/* Read in the header */
	RGTreeReadHeader(fp, tree);

	/* Read in the tree */
	return RGTreeReadFromFileHelper((void*)tree->root, fp, tree, 0, tree->depth);
}

/* TODO */
int RGTreeReadFromFileHelper(void *node, FILE *fp, RGTree *tree, int curDepth, int depth)
{
	int i, j;
	int tempInt;
	int numPositions;
	RGNode *curNode=NULL;
	RGLeaf *curLeaf=NULL;

	/* Basic idea is to read the tree from file */
	for(i=0;i<ALPHABET_SIZE;i++) { /* For each potential leaf */
		/*
		   if(fread(&tempInt, sizeof(unsigned int), 1, fp)==EOF) {
		   }
		   */
		if(fscanf(fp, "%d", &tempInt)==EOF) {
			return EOF;
		}
		else {
			/* Convert to proper no (Big Endian/Little Endian) */
			/*
			   tempInt = ntohl(tempInt);
			   */

			/* Add to input */
			tempInt += RGT_ADD_TO_INPUT;

			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "%d[%d]\t",
						tempInt+RGT_ADD_TO_OUTPUT,
						tempInt);
			}

			if(tempInt >= 0 && tempInt <=3) {
				assert(curDepth < depth);
				assert(tempInt == i);
				curNode = (RGNode*)node;
				/* Check if we have to allocate a new node */
				if(curNode->next[tempInt] == NULL) {
					/* Check what type of node we have to add */
					if(curDepth == depth-1) { /* Current node's child must be a leaf */
						if(VERBOSE >= DEBUG) {
							fprintf(stderr, "Child %d is a leaf (%d).\n", i, curDepth);
						}
						curNode->next[tempInt] = (void*)malloc(sizeof(RGLeaf));
						/* Initialize fields */
						curLeaf = (RGLeaf*)curNode->next[tempInt];
						curLeaf->positions = NULL;
						curLeaf->chromosomes = NULL;
						curLeaf->numEntries = 0;
					}
					else { /* Current node's child must be an internal node */
						if(VERBOSE >= DEBUG) {
							fprintf(stderr, "Child %d is an internal node (%d).\n", i, curDepth);
						}
						curNode->next[tempInt] = (void*)malloc(sizeof(RGNode));
						/* Initialize pointers */
						for(j=0;j<ALPHABET_SIZE;j++) {
							((RGNode*)curNode->next[tempInt])->next[j] = NULL;
						}
					}
				}
				else {
					if(curDepth == depth-1) { /* Current node's child must be a leaf */
						if(VERBOSE >= DEBUG) {
							fprintf(stderr, "Child %d is a leaf (%d).\n", i, curDepth);
						}
					}
					else { /* Current node's child must be an internal node */
						if(VERBOSE >= DEBUG) {
							fprintf(stderr, "Child %d is an internal node (%d).\n", i, curDepth);
						}
					}
				}
				/* Keep going deeper */
				if(RGTreeReadFromFileHelper((void*)curNode->next[tempInt], fp, tree, curDepth+1, depth)==EOF) {
					return EOF;
				}
			}
			else if(tempInt == -1) {
				assert(curDepth < depth); /* We must be at an internal node */
				/* Do nothing this is null */
				if(VERBOSE >= DEBUG) {
					fprintf(stderr, "Child %d is null (%d).\n", i, curDepth);
				}
				curNode = (RGNode*)node;
				assert(curNode->next[i] == NULL);
			}
			else if(tempInt == -2) {
				/* We are at a leaf */
				assert(curDepth == depth);
				curLeaf = (RGLeaf*)node;

				/* Read the number of positions to be read */
				if(fscanf(fp, "%d", &numPositions)==EOF) {
					fprintf(stderr, "Error.  Trying to read the number of positions.  Terminating!\n");
				}
				/*
				   if(fread(&numPositions, sizeof(unsigned int), 1, fp)!=EOF) {
				   }
				   */
				/* Convert to proper no (Big Endian/Little Endian) */
				/*
				   numPositions = ntohl(numPositions);
				   */

				/* Add offset */
				numPositions += RGT_ADD_TO_INPUT;
				assert(numPositions > 0);

				if(VERBOSE >= DEBUG) {
					fprintf(stderr, "%d[%d]\t",
							numPositions+RGT_ADD_TO_OUTPUT,
							numPositions);
				}

				/* Allocate memory for positions and chromosomes */
				assert(curLeaf != NULL);
				curLeaf->numEntries = numPositions;
				curLeaf->positions = (int*)malloc(sizeof(int)*numPositions);
				curLeaf->chromosomes = (char*)malloc(sizeof(char)*numPositions);

				/* Read in positions */
				/*
				   if(fread(curLeaf.positions, sizeof(unsigned int), numPositions, fp)==EOF) {
				   fprintf(stderr, "Error.  Trying to read %d positions.  Terminating!\n", numPositions);
				   exit(1);
				   }
				   */
				for(j=0;j<numPositions;j++) {
					if(fscanf(fp, "%d", &curLeaf->positions[j])==EOF) {
						fprintf(stderr, "Error.  Trying to read %d positions.  Terminating!\n", numPositions);
						exit(1);
					}
					/* Convert to proper no (Big Endian/Little Endian) */
					/*
					   curLeaf->positions[j] = ntohl(curLeaf->positions[j]);
					   */
					/* Add offset */
					curLeaf->positions[j] += RGT_ADD_TO_INPUT;
					if(VERBOSE >= DEBUG) {
						fprintf(stderr, "%d[%d]\t",
								curLeaf->positions[j] +RGT_ADD_TO_OUTPUT,
								curLeaf->positions[j]);
					}
				}

				/* Read in chromosomes */
				/*
				 * THIS IS WRONG SINCE chromosomes is of type char!!
				 if(fread(curLeaf->chromosomes, sizeof(unsigned int), numPositions, fp)==EOF) {
				 fprintf(stderr, "Error.  Trying to read %d chromosomes.  Terminating!\n", numPositions);
				 exit(1);
				 }
				 */
				/* Convert to proper no (Big Endian/Little Endian) */
				for(j=0;j<numPositions;j++) {
					if(fscanf(fp, "%d", &tempInt)==EOF) {
						fprintf(stderr, "Error.  Trying to read %d chromosomes.  Terminating!\n", numPositions);
						exit(1);
					}

					/*
					   curLeaf->chromosomes[j] = ntohl(curLeaf->chromosomes[j]);
					   */
					/* Add offset */
					curLeaf->chromosomes[j] = tempInt + RGT_ADD_TO_INPUT;
					if(VERBOSE >= DEBUG) {
						fprintf(stderr, "%d[%d]\t",
								curLeaf->chromosomes[j] + RGT_ADD_TO_INPUT,
								curLeaf->chromosomes[j]);
					}
				}
				if(VERBOSE >= DEBUG) {
					fprintf(stderr, "Added information to the leaf (%d).\n", curDepth);
				}
				/* No more nodes at this level */
				return 1;
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

	if(VERBOSE > 0) {
		fprintf(stderr, "Printing Header:%d,%d,%d,%d,%d,%d,%d\n",
				subtractGapInformation,
				addGapInformation,
				depthOfRGTree,
				startChr,
				startPos,
				endChr,
				endPos);
	}

	/* Big Endian/Little Endian conversion */
	/*
	   subtractGapInformation=htonl(subtractGapInformation);
	   addGapInformation=htonl(addGapInformation);
	   depthOfRGTree=htonl(depthOfRGTree);
	   startChr=htonl(startChr);
	   startPos=htonl(startPos);
	   endChr=htonl(endChr);
	   endPos=htonl(endPos);
	   */

	/* Print Header */
	/*
	   fwrite(&subtractGapInformation, sizeof(unsigned int), 1, fp);
	   fwrite(&addGapInformation, sizeof(unsigned int), 1, fp);
	   fwrite(&depthOfRGTree, sizeof(unsigned int), 1, fp);
	   fwrite(&startChr, sizeof(unsigned int), 1, fp);
	   fwrite(&startPos, sizeof(unsigned int), 1, fp);
	   fwrite(&endChr, sizeof(unsigned int), 1, fp);
	   fwrite(&endPos, sizeof(unsigned int), 1, fp);
	   */
	fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d",
			subtractGapInformation,
			addGapInformation,
			depthOfRGTree,
			startChr,
			startPos,
			endChr,
			endPos);
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
	/*
	   if(fread(&subtractGapInformation, sizeof(unsigned int), 1, fp)!=EOF 
	   && fread(&addGapInformation, sizeof(unsigned int), 1, fp)!=EOF 
	   && fread(&depthOfRGTree, sizeof(unsigned int), 1, fp)!=EOF 
	   && fread(&startChr, sizeof(unsigned int), 1, fp)!=EOF
	   && fread(&startPos, sizeof(unsigned int), 1, fp)!=EOF
	   && fread(&endChr, sizeof(unsigned int), 1, fp)!=EOF
	   && fread(&endPos, sizeof(unsigned int), 1, fp)!=EOF) {
	   */
	if(fscanf(fp, "%d %d %d %d %d %d %d",
				&subtractGapInformation,
				&addGapInformation,
				&depthOfRGTree,
				&startChr,
				&startPos,
				&endChr,
				&endPos)!=EOF) {

		/* Big Endian/Little Endian conversion */
		/*
		   subtractGapInformation = ntohl(subtractGapInformation);
		   addGapInformation = ntohl(addGapInformation);
		   startChr = ntohl(startChr);
		   startPos = ntohl(startPos);
		   endChr = ntohl(endChr);
		   endPos = ntohl(endPos);
		   depthOfRGTree = ntohl(depthOfRGTree);
		   */

		if(VERBOSE > 0) {
			fprintf(stderr, "Printing Header:%d,%d,%d,%d,%d,%d,%d\n",
					subtractGapInformation,
					addGapInformation,
					depthOfRGTree,
					startChr,
					startPos,
					endChr,
					endPos);
		}

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
	double midX;
	int upperX=pow(ALPHABET_SIZE, tree->depth/2) -1;
	int lowerY=0;
	double midY;
	int upperY=upperX;

	cur = tree->root;
	assert(cur!=NULL);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "\nRGTreeGetMatch.  Searching for indexes %d and %d.\n",
				indexOne,
				indexTwo);
	}

	/* Traverse the tree to find the appropriate leaf */
	for(i=0;i<tree->depth;i++) {
		midX = (lowerX+upperX)/2.0;
		midY = (lowerY+upperY)/2.0;
		/* Four branches to new node (x,y) from search indexes (x',y') 
		 * B1: x'<x, y'<y
		 * B2: x'>x, y'<y
		 * B3: x'<x, y'>y
		 * B4: x'>x, y'>y
		 * */
		compOne = (indexOne < midX)?0:1;
		compTwo = (indexTwo < midY)?0:1;
		nextNode = compOne + 2*compTwo;

		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "depth:%d,lowerX:%d,midX:%.1f,upperX:%d,lowerY:%d,midY:%.1f,upperY:%d,nextNode:%d\n",
					i+1,
					lowerX,
					midX,
					upperX,
					lowerY,
					midY,
					upperY,
					nextNode);
		}

		/* Update lower and upper */
		switch(nextNode) {
			case 0:
				upperX = (int)midX;
				upperY = (int)midY;
				break;
			case 1:
				lowerX = (int)(midX+0.5);
				upperY = (int)midY;
				break;
			case 2:
				upperX = (int)midX;
				lowerY = (int)(midY+0.5);
				break;
			case 3:
				lowerX = (int)(midX+0.5);
				lowerY = (int)(midY+0.5);
				break;
			default:
				fprintf(stderr, "Error.  Could not understand nextNode [%d] in RGTree.c.  Terminating!\n", nextNode);
				exit(1);
				break;
		}

		if(cur->next[nextNode] == NULL) {
			/* If the next node does not exist, do not add any matches and return */
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "NO MATCH EXISTS (at depth %d)\n", i);
			}
			return 0;
		}
		else if(i==tree->depth-1) {
			/* We are at a leaf, copy over the matches */
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "At a leaf, copying over...\n");
			}

			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "Will copy over %d entries.\n",
						((RGLeaf*)cur->next[nextNode])->numEntries);
			}

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
