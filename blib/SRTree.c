#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <string.h>
#include "BLibDefinitions.h"
#include "SRTree.h"

/* TODO */
/* We assume the root has already been initialized */
int SRTreeInsert(SRNode *root, char *sequence, int length) {
	{
		int i, j, nextNode;
		SRNode *cur=NULL;

		/* Initialize cur to be the root */
		assert(root!=NULL);
		cur = root;

		/* Traverse the tree to insert */
		for(i=0;i<length;i++) {
			switch(sequence[i]) {
				case 'a':
					nextNode = 0;
					break;
				case 'c':
					nextNode = 1;
					break;
				case 'g':
					nextNode = 2;
					break;
				case 't':
					nextNode = 3;
					break;
				default:
					fprintf(stderr, "Error.  In SRTreeInsert could not understand bp [%c].  Terminating!\n",
							sequence[i]);
					exit(1);
					break;
			}

			/* Allocate a new node if it does not exist */
			if(cur->next[nextNode] == NULL) {
				/* Allocate the node */
				cur->next[nextNode] = malloc(sizeof(SRNode));
				assert(cur->next[nextNode] != NULL);
				/* Initialize it's pointers */
				for(j=0;j<4;j++) {
					((SRNode*)cur->next[nextNode])->next[j] = NULL;
				}
			}
		}
		/* move to next one */
		cur = cur->next[nextNode];
	}

	return 1;
}

/* TODO */
void SRTreeDelete(SRNode **node, int curDepth, int depth)
{
	int i;
	if(NULL != (*node)) {
		for(i=0;i<4;i++) {
			/* Delete the children */
			SRTreeDelete((*node)->next[i], curDepth+1, depth);
		}
		free((*node));
		(*node) = NULL;
	}
}

/* TODO */
double SRTreeGetSize(SRNode *root, int depth, int outputSize) 
{
	double total=SRTreeGetSizeHelper(root, 0, depth, -1);
	switch(outputSize) {
		case SRT_KILOBYTES:
			return (total/1000);
			break;
		case SRT_MEGABYTES:
			return (total/1000000);
			break;
		case SRT_GIGABYTES:
			return (total/1000000000);
			break;
		default:
			return total;
			break;
	}
}

/* TODO */
double SRTreeGetSizeHelper(SRNode* node, int curDepth, int depth, int direction) 
{
	if(NULL!=node) {
			return (SRTreeGetSizeHelper((SRNode*)((SRNode*)node)->next[0], curDepth+1, depth, 0) + 
					SRTreeGetSizeHelper((SRNode*)((SRNode*)node)->next[1], curDepth+1, depth, 1) + 
					SRTreeGetSizeHelper((SRNode*)((SRNode*)node)->next[2], curDepth+1, depth, 2) + 
					SRTreeGetSizeHelper((SRNode*)((SRNode*)node)->next[3], curDepth+1, depth, 3) + 
					(int)sizeof(SRNode));
	}
	else {
		return 0;
	}
}

/* TODO */
/* Note: that if we have more than 4 billion sequences, the return will overflow */
int SRTreeReadFromFile(SRNode *root, FILE *fp)
{
	char sequenceName[SRT_SEQUENCE_NAME_LENGTH]="\0";
	char sequence[SRT_SEQUENCE_LENGTH]="\0";
	int length;
	int count = 0;

	while(fscanf(fp, "%s", sequenceName)!=EOF && fscanf(fp, "%s", sequence)!=EOF) {
		/* ignore the name for now */

		/* get the sequence length */
		length = strlen(sequence);
		/* move to lowercase */
		SequenceToLower(sequence, length);
		/* insert into the tree */
		SRTreeInsert(root, sequence, length);
		count = 0;
	}

	return count;
}

void SequenceToLower(char* sequence, int length) 
{
	int i;
	for(i=0;i<length;i++) {
		switch(sequence[i]) {
			case 'A':
				sequence[i] = 'a';
				break;
			case 'C':
				sequence[i] = 'c';
				break;
			case 'G':
				sequence[i] = 'g';
				break;
			case 'T':
				sequence[i] = 't';
				break;
			default:
				/* do nothing */
				break;
		}
	}
}
