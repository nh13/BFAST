#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include "../blib/BError.h"
#include "../blib/RGTree.h"
#include "../blib/BLibDefinitions.h"
#include "Definitions.h"
#include "GenerateTree.h"

/* TODO */
void GenerateTree(RGList *rgList,
		int matchLength, 
		int *gaps,
		int numGaps,
		int numThreads,
		char *outputID, 
		char *outputDir,
		int binaryOutput,
		int reachedEnd)
{
	int i, j, k, l;
	char *curSequenceOne=NULL;
	char *curSequenceTwo=NULL;
	int startOne;
	int startTwo;
	char outputFileName[1024];
	FILE *fp;
	RGTree tree;
	int curGap;
	int maxGap = INT_MIN;

	assert(numGaps > 0);

	/* Allocate memory for holding the temporary sequences */
	curSequenceOne = (char*)malloc(sizeof(char)*matchLength);
	if(NULL==curSequenceOne) {
		PrintError("GenerateTree",
				"curSequenceOne",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	curSequenceTwo = (char*)malloc(sizeof(char)*matchLength);
	if(NULL==curSequenceTwo) {
		PrintError("GenerateTree",
				"curSequenceTwo",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "Will generate %d different %s files.\n%s",
				numGaps,
				BLATTER_TREE_FILE_EXTENSION,
				BREAK_LINE);
	}
	
	/* Get the max gap */
	assert(numGaps > 0);
	maxGap = gaps[0];
	for(i=1;i<numGaps;i++) {
		if(gaps[i] > maxGap) {
			maxGap = gaps[i];
		}
	}

	/* For every gap between the two l-mers */
	for(i=0;i<numGaps;i++) {
		curGap = gaps[i];

		/* Initialize the tree */
		tree.nodes = NULL;
		tree.numNodes = 0;
		tree.gap = curGap;
		tree.matchLength = matchLength;
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
				tree.endPos-((reachedEnd==1)?0:(2*matchLength-1+maxGap)),
				curGap,
				matchLength,
				BLATTER_TREE_FILE_EXTENSION);

		if(VERBOSE >=0) {
			fprintf(stderr, "Generating a tree, outputing to %s.\nCurrently on [gap, chr, pos, size]:\n%d\t%d\t%d\t%d", 
					outputFileName, -1, -1, -1, -1);
		}

		/* For every chromosome in the reference genome */
		for(j=0;j<rgList->numChrs;j++) {

			/* For every starting position of the first l-mer */
			/* Note end position is really 2*matchLength+maxGap-1 to account for max gap */
			for(k=rgList->chromosomes[j].startPos;k <= rgList->chromosomes[j].endPos - curGap + maxGap;k++) {
				if(VERBOSE >=0 && k%GT_ROTATE_NUM==0) {
					fprintf(stderr, "\r%3d\t%2d\t%12d\t%10.2lfMB", 
							curGap,
							rgList->chromosomes[j].chromosome,
							k,
							RGTreeGetSize(&tree, RGT_MEGABYTES));
				}
				startOne = k;
				startTwo = startOne+matchLength+curGap;
				if(VERBOSE>=DEBUG) {
					fprintf(stderr, "\nstartOne[%d]\tstartTwo[%d]\n",
							startOne, 
							startTwo);
				}
				/* Copy over sequences */
				for(l=startOne;l<startOne+matchLength;l++) {
					curSequenceOne[l-startOne] = rgList->chromosomes[j].sequence[l-rgList->chromosomes[j].startPos];
				}
				for(l=startTwo;l<startTwo+matchLength;l++) {
					curSequenceTwo[l-startTwo] = rgList->chromosomes[j].sequence[l-rgList->chromosomes[j].startPos];
				}
				if(VERBOSE>=DEBUG) {
					fprintf(stderr, "Inserting pair [%s] [%s] at chr%d:%d.\n",
							curSequenceOne,
							curSequenceTwo,
							rgList->chromosomes[j].chromosome,
							startOne);
				}
				/* Only insert if it is a valid sequence */
				if(ValidateSequence(curSequenceOne, matchLength)==1 && ValidateSequence(curSequenceTwo, matchLength)) {
					/* Insert pair into the tree */
					RGTreeInsert(&tree, curSequenceOne, curSequenceTwo, matchLength, rgList->chromosomes[j].chromosome, startOne);
				}
			}
			if(VERBOSE >=0) {
				fprintf(stderr, "\r%3d\t%2d\t%12d\t%10.2lfMB", 
						curGap,
						rgList->chromosomes[j].chromosome,
						k-1,
						RGTreeGetSize(&tree, RGT_MEGABYTES));
			}
		}
		/* Clean up the tree */
		if(VERBOSE >= 0) {
			fprintf(stderr, "\n");
			fprintf(stderr, "Cleaning up the tree.\n");
		}
		RGTreeCleanUpTree(&tree, numThreads);

		if(VERBOSE >= 0) {
			fprintf(stderr, "Outputting tree to %s\n", outputFileName);
		}


		/* Open the output file */
		if(!(fp=fopen(outputFileName, "wb"))) {
			PrintError("GenerateTree",
					outputFileName,
					"Could not open output file name for writing",
					Exit,
					OpenFileError);
		}   

		RGTreePrintTree(fp, &tree, binaryOutput);

		/* Close the output file */
		fclose(fp);

		/* Free memory */
		if(VERBOSE >=0) {
			fprintf(stderr, "Cleaning up.\n");
		}
		RGTreeDelete(&tree);
		if(VERBOSE >=0) {
			fprintf(stderr, "%s", BREAK_LINE);
		}
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "Generated Tree(s) successfully!\n");
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Free memory */
	free(curSequenceOne);
	free(curSequenceTwo);
}
