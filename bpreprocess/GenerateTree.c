#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "../blib/RGTree.h"
#include "../blib/BLibDefinitions.h"
#include "Definitions.h"
#include "GenerateTree.h"

/* TODO */
void GenerateTree(RGList *rgList,
		int matchLength, 
		char *gapFileName,
		char *outputID, 
		char *outputDir)
{
	int i, j, k, l;
	char *curSequenceOne=NULL;
	char *curSequenceTwo=NULL;
	int startOne;
	int startTwo;
	char outputFileName[1024];
	FILE *fp;

	RGTree tree;

	int tempGap;
	int *gaps=NULL;
	int numGaps;
	int curGap;

	/* Open the gap file */
	if(!(fp=fopen(gapFileName, "r"))) {
		fprintf(stderr, "Error. Could not open %s for reading.  Terminating!\n",
				gapFileName);
		exit(1);
	}   

	/* Read in the gap file */
	numGaps=0;
	while(fscanf(fp, "%d", &tempGap)!=EOF) {
		numGaps++;
		gaps=realloc(gaps, sizeof(int)*numGaps);
		gaps[numGaps-1] = tempGap;
		if(numGaps>1 && gaps[numGaps-2] >= gaps[numGaps-1]) {
			fprintf(stderr, "Error.  Gaps in %s must be in ascending order.  Terminating!\n",
					gapFileName);
			exit(1);
		}
	}
	if(numGaps <= 0) {
		fprintf(stderr, "Error. Could not read any gaps from %s.  Terminating!\n",
				gapFileName);
		exit(1);
	}

	/* Allocate memory for holding the temporary sequences */
	curSequenceOne = (char*)malloc(sizeof(char)*matchLength);
	curSequenceTwo = (char*)malloc(sizeof(char)*matchLength);

	/* For every gap between the two l-mers */
	for(i=0;i<numGaps;i++) {
		curGap = gaps[i];

		/* Initialize root */
		tree.root = (RGNode*)malloc(sizeof(RGNode));
		for(j=0;j<ALPHABET_SIZE;j++) {
			tree.root->next[j] = NULL;
		}
		/* Add metadata */
		tree.depth = matchLength;
		tree.gap = curGap;
		tree.startChr = rgList->startChr;
		tree.startPos = rgList->startPos;
		tree.endChr = rgList->endChr;
		tree.endPos = rgList->endPos;

		/* Create the output file name */
		sprintf(outputFileName, "%sblatter.tree.file.%s.%d.%d.%d.%d.%d.%d.%s", 
				outputDir,
				outputID,
				tree.startChr,
				tree.startPos,
				tree.endChr,
				tree.endPos,
				curGap,
				tree.depth,
				BLATTER_TREE_FILE_EXTENSION);

		fprintf(stderr, "Generating a tree, outputing to %s.\nCurrently on [gap, chr, pos, size]:\n%d\t%d\t%d\t%d", 
				outputFileName, -1, -1, -1, -1);

		/* For every chromosome in the reference genome */
		for(j=0;j<rgList->numChrs;j++) {

			/* For every starting position of the first l-mer */
			for(k=rgList->chromosomes[j].startPos;k+2*matchLength+curGap<=rgList->chromosomes[j].endPos;k++) {
				if(k%GT_ROTATE_NUM==0) {
					fprintf(stderr, "\r%3d\t%2d\t%12d\t%10.2lfMB", 
							curGap,
							rgList->chromosomes[j].chromosome,
							k,
							RGTreeGetSize(&tree, RGT_MEGABYTES));
				}
				startOne = k;
				startTwo = k+matchLength+curGap;
				/* Copy over sequences */
				for(l=startOne;l<startOne+matchLength;l++) {
					curSequenceOne[l-startOne] = rgList->chromosomes[j].sequence[l-startOne];
				}
				for(l=startTwo;l<startTwo+matchLength;l++) {
					curSequenceOne[l-startTwo] = rgList->chromosomes[j].sequence[l-startTwo];
				}
				/* Only insert if it is a valid sequence */
				if(ValidateSequence(curSequenceOne, matchLength)==1 && ValidateSequence(curSequenceTwo, matchLength)) {
					/* Insert pair into the tree */
					RGTreeInsert(&tree, curSequenceOne, curSequenceTwo, matchLength, rgList->chromosomes[j].chromosome, startOne);
				}
			}
			fprintf(stderr, "\r%3d\t%2d\t%12d\t%10.2lfMB", 
					curGap,
					rgList->chromosomes[j].chromosome,
					k,
					RGTreeGetSize(&tree, RGT_MEGABYTES));
		}
		fprintf(stderr, "\n");

		fprintf(stderr, "Outputting tree to %s\n", outputFileName);

		/* Open the output file */
		if(!(fp=fopen(outputFileName, "wb"))) {
			fprintf(stderr, "Error. Could not open %s for writing.  Terminating!\n",
					outputFileName);
			exit(1);
		}   

		RGTreePrintTree(fp, &tree);

		/* Close the output file */
		fclose(fp);

		/* Free memory */
		RGTreeDelete(&tree);
	}
	fprintf(stderr, "\n");

	fprintf(stderr, "Generated Tree(s) successfully!\n");

	/* Free memory */
	free(curSequenceOne);
	free(curSequenceTwo);
}

int ValidateSequence(char *sequence, int length)
{
	int i;
	for(i=0;i<length;i++) {
		switch(sequence[i]) {
			case 'a':
				break;
			case 'c':
				break;
			case 'g':
				break;
			case 't':
				break;
			default:
				return 0;
		}
	}
	return 1;
}
