#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <string.h>
#include "BLibDefinitions.h"
#include "RGIndex.h"
#include "RGTree.h"
#include "RGSeqPair.h"

char ALPHABET[ALPHABET_SIZE] = "acgt";

/* TODO */
void RGSeqPairFindMatchesInIndex(RGIndex *index, 
		RGMatch *match,
		char *read)
{
	unsigned char *indexForward=NULL;
	unsigned char *indexReverse=NULL;
	char reverseSequence[SEQUENCE_LENGTH];
	int numChars = (int)ceil((2.0/8.0*index->matchLength)/sizeof(unsigned char));

	/* Allocate memory for the indexes */
	indexForward = malloc(sizeof(unsigned char)*numChars);
	assert(indexForward!=NULL);
	indexReverse = malloc(sizeof(unsigned char)*numChars);
	assert(indexReverse!=NULL);

	/* DON'T FORGET TO ALIGN TO BOTH FORWARD AND REVERSE STRANDS */
	assert(strlen(read)==index->matchLength);
	GetReverseCompliment(read, reverseSequence, strlen(read));

	/* Get indexes */
	RGIndexGetIndexFromSequence(read, index->matchLength, indexForward); 
	RGIndexGetIndexFromSequence(reverseSequence, index->matchLength, indexReverse); 

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "\nIn RGSeqPairFindMatchesInIndex.\n");
	}
	/* Get forward matches */
	RGIndexGetMatches(index, 
			indexForward,
			FORWARD,
			match);
	/* Get reverse matches */
	RGIndexGetMatches(index, 
			indexReverse,
			REVERSE,
			match);
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Found %d matches in RGSeqPairFindMatchesInIndex.\n",
				match->numEntries);
	}

	/* Remove duplicates from match */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Removing duplicates\n");
	}
	RGMatchRemoveDuplicates(match);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Removed duplicates with %d matches remaining in RGSeqPairFindMatchesInIndex.\n",
				match->numEntries);
	}

	/* Free memory */
	free(indexForward);
	free(indexReverse);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Exiting RGSeqPairFindMatchesInIndex.\n");
	}
}

/* TODO */
void RGSeqPairFindMatchesInTree(RGTree *tree, 
		RGMatch *match,
		char *read,
		int **offsets,
		int numOffsets,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions)
{
	int i;

	/* DON'T FORGET TO ALIGN TO BOTH FORWARD AND REVERSE STRANDS */
	RGSeqPair readPairs;
	readPairs.numPairs=0;
	readPairs.indexOne=NULL;
	readPairs.indexTwo=NULL;
	readPairs.strand=NULL;
	readPairs.offset=NULL;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "\nIn RGSeqPairFindMatchesInTree\n");
	}

	/* Note: we can speed this up by generating all pairs of l-mers to
	 * search for and then removing duplicated pairs. */
	RGSeqPairGeneratePairs(read,
			&readPairs,
			offsets,
			numOffsets,
			tree->matchLength,
			tree->gap,
			numMismatches,
			numInsertions,
			numDeletions,
			numGapInsertions,
			numGapDeletions);
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Generated %d pairs.\n",
				readPairs.numPairs);
	}
	if(readPairs.numPairs <= 0) { 
		/* No pairs generated, return */
		return;
	}
	if(VERBOSE >= DEBUG) {
		for(i=0;i<readPairs.numPairs;i++) { /* For each pair */
			fprintf(stderr, "%d:%d,%d,%c,%d\n",
					i+1, 
					readPairs.indexOne[i],
					readPairs.indexTwo[i],
					readPairs.strand[i],
					readPairs.offset[i]);
		}
	}

	/* All pairs have been generated and are unique.  The next step
	 * is to use each pair and search in the tree to get all matches. 
	 * Store the results in match.
	 * */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "After %d pairs, found %d matches.\n",
				0,
				match->numEntries);
	}
	for(i=0;i<readPairs.numPairs;i++) { /* For each pair */
		if(readPairs.strand[i] != FORWARD && readPairs.strand[i] != REVERSE) {
			fprintf(stderr, "Error,  Direction not recognized [%c].  Terminating!\n",
					readPairs.strand[i]);
			exit(1);
		}
		RGTreeGetMatches(tree, 
				readPairs.indexOne[i],
				readPairs.indexTwo[i],
				readPairs.strand[i],
				readPairs.offset[i],
				match);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "After %d pairs, found %d matches.\n",
					i+1,
					match->numEntries);
			int j;
			for(j=0;j<match->numEntries;j++) {
				fprintf(stderr, "%d\t%d\t%c\n",
						match->chromosomes[j],
						match->positions[j],
						match->strand[j]);
			}
		}
	}
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Found %d matches in RGSeqPairFindMatchesInTree.\n",
				match->numEntries);
	}

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "\nRemoving duplicates\n");
	}
	/* Remove duplicates from match */
	RGMatchRemoveDuplicates(match);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Removed duplicates with %d matches remaining in RGSeqPairFindMatchesInTree.\n",
				match->numEntries);
	}
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Exiting RGSeqPairFindMatchesInTree\n");
	}

}

/* TODO */
void RGSeqPairGeneratePairs(char *read,
		RGSeqPair *readPairs,
		int **offsets,
		int numOffsets,
		int matchLength,
		int gap,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions)
{
	/* DON'T FORGET TO ALIGN TO BOTH FORWARD AND REVERSE STRANDS */
	int i;
	int readLength=strlen(read);
	char reverseSequence[SEQUENCE_LENGTH];

	/* Get the reverse compliment of the read */
	GetReverseCompliment(read, reverseSequence, readLength);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Generating all possible pairs (%d offsets).\n",
				numOffsets);
		fprintf(stderr, "numMismatches:%d\tnumInsertions:%d\tnumDeletions:%d\n",
				numMismatches,
				numInsertions,
				numDeletions);
	}

	/* Go through all offsets */
	for(i=0;i<numOffsets;i++) {
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Generating pairs with offset %d.\n",
					(*offsets)[i]);
		}
		/* Go through all mismatches */
		/* Note: we allow any number (including zero) of mismatches up to
		 * numMismatches.  We also note that when mismatches is zero, this 
		 * will just insert an unaltered pair, so no need to restrict calling
		 * to when numMismatches>0.
		 * */
		/* Forward */
		RGSeqPairGenerateMismatches(read,
				readLength,
				FORWARD,
				(*offsets)[i],
				matchLength,
				gap,
				numMismatches,
				readPairs);
		/* Reverse compliment */
		RGSeqPairGenerateMismatches(reverseSequence,
				readLength,
				REVERSE,
				(*offsets)[i],
				matchLength,
				gap,
				numMismatches,
				readPairs);

		/* Go through all insertions */
		/* Note: we allow only contiguous insertions of length up to
		 * numInsertions.  We also always start from the offset.  When
		 * we model insertions, we have to remove a certain number of
		 * bases in the read, and therefore we enumerate over all
		 * possible insertions in the entire read.
		 * */
		if(numInsertions > 0) {
			/* Forward */
			RGSeqPairGenerateInsertions(read,
					readLength,
					FORWARD,
					(*offsets)[i],
					matchLength,
					gap,
					numInsertions,
					readPairs);
			/* Reverse compliment */
			RGSeqPairGenerateInsertions(reverseSequence,
					readLength,
					REVERSE,
					(*offsets)[i],
					matchLength,
					gap,
					numInsertions,
					readPairs);
		}

		/* GO through all deletions */
		/* Note: we allow only contiguous deletions of length up to 
		 * numDeletions.  We also always start from the offset.  We
		 * must add base to the reads, and therfore we enumerate
		 * over all possible deletions in the entire read.
		 * */
		if(numDeletions > 0) {
			/* Forward */
			RGSeqPairGenerateDeletions(read,
					readLength,
					FORWARD,
					(*offsets)[i],
					matchLength,
					gap,
					numDeletions,
					readPairs);
			/* Reverse compliment */
			RGSeqPairGenerateDeletions(reverseSequence,
					readLength,
					REVERSE,
					(*offsets)[i],
					matchLength,
					gap,
					numDeletions,
					readPairs);
		}

		/* Go through all possible insertions in the gap between
		 * the pair of l-mers.  If there is a gap insertion,
		 * then we will delete bases in the gap.
		 * */
		if(numGapInsertions > 0) {
			/* Forward */
			RGSeqPairGenerateGapInsertions(read,
					readLength,
					FORWARD,
					(*offsets)[i],
					matchLength,
					gap,
					numGapInsertions,
					readPairs);
			/* Reverse compliment */
			RGSeqPairGenerateGapInsertions(reverseSequence,
					readLength,
					REVERSE,
					(*offsets)[i],
					matchLength,
					gap,
					numGapInsertions,
					readPairs);
		}

		/* Go through all possible deletions in the gap between
		 * the pair of l-mers.  If there is a gap deletion, 
		 * then we will add bases to the gap.
		 * */
		if(numGapDeletions > 0) {
			/* Forward */
			RGSeqPairGenerateGapDeletions(read,
					readLength,
					FORWARD,
					(*offsets)[i],
					matchLength,
					gap,
					numGapDeletions,
					readPairs);
			/* Reverse compliment */
			RGSeqPairGenerateGapDeletions(reverseSequence,
					readLength,
					REVERSE,
					(*offsets)[i],
					matchLength,
					gap,
					numGapDeletions,
					readPairs);
		}
	}

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "All pairs generated\n");
	}

	/* Merge all pairs */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Calling RGSeqPairRemoveDuplicates\n");
	}
	RGSeqPairRemoveDuplicates(readPairs);
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Exited from RGSeqPairRemoveDuplicates\n");
	}
}

/* TODO */
void RGSeqPairGenerateMismatches(char *read,
		int readLength,
		char direction,
		int offset,
		int matchLength,
		int gap,
		int numMismatches,
		RGSeqPair *pairs)
{
	char *curOne=NULL;
	char *curTwo=NULL;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Generating indexes with %d mismatches\n",
				numMismatches);
	}

	if(offset+gap+matchLength*2 > readLength) {
		/* Out of bounds.  Don't add anything. */
		return;
	}
	else {
		/* Allocate memory */
		curOne = malloc(sizeof(char)*matchLength);
		assert(curOne!=NULL);
		curTwo = malloc(sizeof(char)*matchLength);
		assert(curTwo!=NULL);

		RGSeqPairGenerateMismatchesHelper(read,
				direction,
				offset,
				matchLength,
				gap,
				numMismatches,
				0,
				0,
				pairs,
				curOne,
				curTwo,
				0,
				0);

		/* Free memory */
		free(curOne);
		free(curTwo);
	}
}

/* TODO */
void RGSeqPairGenerateMismatchesHelper(char *read,
		char direction,
		int offset,
		int matchLength,
		int gap,
		int numMismatchesLeft,
		int numFirstPrinted,
		int numLastPrinted,
		RGSeqPair *pairs,
		char *curOne,
		char *curTwo,
		int curOneIndex,
		int curTwoIndex)
{
	int i;

	assert(direction == FORWARD || direction == REVERSE);

	if(numMismatchesLeft > 0) {
		/* No more to print */
		if(numFirstPrinted >= matchLength && numLastPrinted >= matchLength) {
			curOne[matchLength]='\0';
			curTwo[matchLength]='\0';
			/* Allocate memory */
			pairs->numPairs++;
			pairs->indexOne = realloc(pairs->indexOne, sizeof(unsigned int)*(pairs->numPairs));
			assert(pairs->indexOne!=NULL);
			pairs->indexTwo = realloc(pairs->indexTwo, sizeof(unsigned int)*(pairs->numPairs));
			assert(pairs->indexTwo!=NULL);
			pairs->strand = realloc(pairs->strand, sizeof(char)*(pairs->numPairs));
			assert(pairs->strand!=NULL);
			pairs->offset = realloc(pairs->offset, sizeof(unsigned int)*(pairs->numPairs));
			/* Copy over */
			pairs->indexOne[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curOne, matchLength);
			pairs->indexTwo[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curTwo, matchLength);
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "Pair(mistmatch):%s[%d]\t%s[%d].\n",
						curOne,
						pairs->indexOne[pairs->numPairs-1],
						curTwo,
						pairs->indexTwo[pairs->numPairs-1]);
			}
			pairs->strand[pairs->numPairs-1] = direction;
			pairs->offset[pairs->numPairs-1] = offset;
			return;
		}
		else {
			/* use mismatches */
			for(i=0;i<ALPHABET_SIZE;i++) {
				if(numFirstPrinted < matchLength) {
					curOne[curOneIndex] = ALPHABET[i];
					/* Use on first read */
					if(read[offset+numFirstPrinted] == ALPHABET[i]) {
						/* No mismatch */

						/* Keep going */
						RGSeqPairGenerateMismatchesHelper(read,
								direction,
								offset,
								matchLength,
								gap,
								numMismatchesLeft,
								numFirstPrinted+1,
								numLastPrinted,
								pairs,
								curOne,
								curTwo,
								curOneIndex+1,
								curTwoIndex);
					}
					else {
						/* Mismatch */

						/* Keep going */
						RGSeqPairGenerateMismatchesHelper(read,
								direction,
								offset,
								matchLength,
								gap,
								numMismatchesLeft-1,
								numFirstPrinted+1,
								numLastPrinted,
								pairs,
								curOne,
								curTwo,
								curOneIndex+1,
								curTwoIndex);
					}
				}
				else if(numLastPrinted < matchLength) {
					curTwo[curTwoIndex] = ALPHABET[i];

					/* Use on second read */
					if(read[offset+matchLength+gap+numLastPrinted] == ALPHABET[i]) {
						/* No mismatch */

						/* Keep going */
						RGSeqPairGenerateMismatchesHelper(read,
								direction,
								offset,
								matchLength,
								gap,
								numMismatchesLeft,
								numFirstPrinted,
								numLastPrinted+1,
								pairs,
								curOne,
								curTwo,
								curOneIndex,
								curTwoIndex+1);
					}
					else {
						/* Mismatch */

						/* Keep going */
						RGSeqPairGenerateMismatchesHelper(read,
								direction,
								offset,
								matchLength,
								gap,
								numMismatchesLeft-1,
								numFirstPrinted,
								numLastPrinted+1,
								pairs,
								curOne,
								curTwo,
								curOneIndex,
								curTwoIndex+1);
					}
				}
				else {
					fprintf(stderr, "Error.  Control should not reach here. 012345.  Terminating!\n");
					exit(1);
					return;
				}
			}
		}
	}
	else {
		/* print remaining */                                           
		if(numFirstPrinted < matchLength) {
			/* Print rest of the first read */
			for(i=offset+numFirstPrinted;i<offset+matchLength;i++) {                                    
				assert(curOneIndex<matchLength);
				curOne[curOneIndex] = read[i];
				curOneIndex++;                  
			}
		}                                       
		if(numLastPrinted < matchLength) {      
			/* Print rest of the second read */                                     
			for(i=offset+matchLength+gap+numLastPrinted;i<offset+2*matchLength+gap;i++) {                               
				assert(curTwoIndex<matchLength);
				curTwo[curTwoIndex] = read[i];                curTwoIndex++;
			}                             
		}
		curOne[matchLength]='\0';                                         
		curTwo[matchLength]='\0';
		/* Allocate memory */                                                                             
		pairs->numPairs++;
		pairs->indexOne = realloc(pairs->indexOne, sizeof(unsigned int)*(pairs->numPairs));
		assert(pairs->indexOne!=NULL);
		pairs->indexTwo = realloc(pairs->indexTwo, sizeof(unsigned int)*(pairs->numPairs));
		assert(pairs->indexTwo!=NULL);
		pairs->strand = realloc(pairs->strand, sizeof(char)*(pairs->numPairs));
		assert(pairs->strand!=NULL);
		pairs->offset = realloc(pairs->offset, sizeof(unsigned int)*(pairs->numPairs));
		/* Copy over */
		pairs->indexOne[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curOne, matchLength);
		pairs->indexTwo[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curTwo, matchLength);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Pair(mismatch):%s[%d]\t%s[%d].\n",
					curOne,
					pairs->indexOne[pairs->numPairs-1],
					curTwo,
					pairs->indexTwo[pairs->numPairs-1]);
		}
		pairs->strand[pairs->numPairs-1] = direction;
		pairs->offset[pairs->numPairs-1] = offset;
		return;
	}
}                        

/* TODO */
/* Note: Deletions have occured, so insert bases */
void RGSeqPairGenerateDeletions(char *read,
		int readLength,
		char direction,
		int offset,
		int matchLength,
		int gap,
		int numDeletions,
		RGSeqPair *pairs)
{
	char *curOne=NULL;
	char *curTwo=NULL;

	assert(direction == FORWARD || direction == REVERSE);
	if(VERBOSE>=DEBUG) {
		fprintf(stderr, "Generating pairs with %d deletions.\n",
				numDeletions);
	}

	if(offset+gap+matchLength*2 > readLength) {
		/* Out of bounds.  Don't add anything. */
		return;
	}
	else {
		/* Allocate memory */
		curOne = malloc(sizeof(char)*matchLength);
		assert(curOne!=NULL);
		curTwo = malloc(sizeof(char)*matchLength);
		assert(curTwo!=NULL);

		RGSeqPairGenerateDeletionsHelper(read,
				readLength,
				direction,
				offset,
				matchLength,
				gap, 
				numDeletions,
				numDeletions,
				0,
				pairs,
				curOne,
				curTwo,
				0,
				0);

		/* Free memory */
		free(curOne);
		free(curTwo);
	}
}

/* TODO */
/* NOTE: no error checking yet! */
void RGSeqPairGenerateDeletionsHelper(char *read,
		int readLength,
		char direction,
		int offset,
		int matchLength, 
		int gap,
		int numDeletionsLeft,
		int numDeletions,
		int deletionOffset,
		RGSeqPair *pairs,
		char *curOne,
		char *curTwo,
		int curOneIndex,
		int curTwoIndex)
{
	int i;
	int remainingOne, remainingTwo;

	remainingOne = matchLength-curOneIndex;
	remainingTwo = matchLength-curTwoIndex;

	if(numDeletionsLeft > 0) {
		/* No more to print */
		if(remainingOne <= 0 && remainingTwo <= 0 && numDeletionsLeft < numDeletions) {
			curOne[matchLength]='\0';
			curTwo[matchLength]='\0';
			/* Allocate memory */                                                                             
			pairs->numPairs++;
			pairs->indexOne = realloc(pairs->indexOne, sizeof(unsigned int)*(pairs->numPairs));
			assert(pairs->indexOne!=NULL);
			pairs->indexTwo = realloc(pairs->indexTwo, sizeof(unsigned int)*(pairs->numPairs));
			assert(pairs->indexTwo!=NULL);
			pairs->strand = realloc(pairs->strand, sizeof(char)*(pairs->numPairs));
			assert(pairs->strand!=NULL);
			pairs->offset = realloc(pairs->offset, sizeof(unsigned int)*(pairs->numPairs));
			/* Copy over */
			pairs->indexOne[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curOne, matchLength);
			pairs->indexTwo[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curTwo, matchLength);
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "Pair(deletion):%s[%d]\t%s[%d].\n",
						curOne,
						pairs->indexOne[pairs->numPairs-1],
						curTwo,
						pairs->indexTwo[pairs->numPairs-1]);
			}
			pairs->strand[pairs->numPairs-1] = direction;
			pairs->offset[pairs->numPairs-1] = offset;
			return;
		}
		else {
			/* try inserting a base */
			if( (curOneIndex > 0 && curTwoIndex == 0 && curOneIndex < matchLength) || (curOneIndex >= matchLength && curTwoIndex > 0)) {
				for(i=0;i<ALPHABET_SIZE;i++) {
					if(remainingOne > 0) {
						curOne[curOneIndex] = ALPHABET[i];
						/* Use on first read */
						RGSeqPairGenerateDeletionsHelper(read,
								readLength,
								direction,
								offset,
								matchLength,
								gap,
								numDeletionsLeft-1,
								numDeletions,
								deletionOffset+1,
								pairs,
								curOne,
								curTwo,
								curOneIndex+1,
								curTwoIndex);
					}
					else if(remainingTwo > 0) {
						curTwo[curTwoIndex] = ALPHABET[i];

						/* Use on second read */
						RGSeqPairGenerateDeletionsHelper(read,
								readLength,
								direction,
								offset,
								matchLength,
								gap,
								numDeletionsLeft-1,
								numDeletions,
								deletionOffset+1,
								pairs,
								curOne,
								curTwo,
								curOneIndex,
								curTwoIndex+1);
					}
					else {
						fprintf(stderr, "Error.  Control should not reach here. 012345.  Terminating!\n");
						exit(1);
						return;
					}
				}
			}
			if(numDeletionsLeft == numDeletions) {
				/* Try not inserting a base */
				if(remainingOne > 0) {
					curOne[curOneIndex] = read[offset+curOneIndex-deletionOffset];
					/* Use on first read */
					RGSeqPairGenerateDeletionsHelper(read,
							readLength,
							direction,
							offset,
							matchLength,
							gap,
							numDeletionsLeft,
							numDeletions,
							deletionOffset,
							pairs,
							curOne,
							curTwo,
							curOneIndex+1,
							curTwoIndex);
				}
				else if(remainingTwo > 0) {
					curTwo[curTwoIndex] = read[offset+gap+matchLength+curTwoIndex-deletionOffset];

					/* Use on second read */
					RGSeqPairGenerateDeletionsHelper(read,
							readLength,
							direction,
							offset,
							matchLength,
							gap,
							numDeletionsLeft,
							numDeletions,
							deletionOffset,
							pairs,
							curOne,
							curTwo,
							curOneIndex,
							curTwoIndex+1);
				}
			}
		}
	}
	else {
		/* print remaining */        if(remainingOne > 0) {
			/* Print rest of the first read */
			for(i=0;i<remainingOne;i++) {
				assert(curOneIndex<matchLength);
				curOne[curOneIndex] = read[offset+curOneIndex-deletionOffset];
				curOneIndex++;
			}
		}
		assert(curOneIndex == matchLength);
		if(remainingTwo > 0) {
			/* Print rest of the second read */
			for(i=0;i<remainingTwo;i++) {
				assert(curTwoIndex<matchLength);
				curTwo[curTwoIndex] = read[offset+gap+matchLength+curTwoIndex-deletionOffset];
				curTwoIndex++;
			}
		}
		assert(curTwoIndex == matchLength);
		curOne[matchLength]='\0';
		curTwo[matchLength]='\0';
		/* Allocate memory */                                                                             
		pairs->numPairs++;
		pairs->indexOne = realloc(pairs->indexOne, sizeof(unsigned int)*(pairs->numPairs));
		assert(pairs->indexOne!=NULL);
		pairs->indexTwo = realloc(pairs->indexTwo, sizeof(unsigned int)*(pairs->numPairs));
		assert(pairs->indexTwo!=NULL);
		pairs->strand = realloc(pairs->strand, sizeof(char)*(pairs->numPairs));
		assert(pairs->strand!=NULL);
		pairs->offset = realloc(pairs->offset, sizeof(unsigned int)*(pairs->numPairs));
		/* Copy over */
		pairs->indexOne[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curOne, matchLength);
		pairs->indexTwo[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curTwo, matchLength);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Pair(deletion):%s[%d]\t%s[%d].\n",
					curOne,
					pairs->indexOne[pairs->numPairs-1],
					curTwo,
					pairs->indexTwo[pairs->numPairs-1]);
		}
		pairs->strand[pairs->numPairs-1] = direction;
		pairs->offset[pairs->numPairs-1] = offset;
		return;
	}
}

/* TODO */
/* Note: Insertions have occured, so delete bases */
void RGSeqPairGenerateInsertions(char *read,
		int readLength,
		char direction,
		int offset,
		int matchLength,
		int gap,
		int numInsertions,
		RGSeqPair *pairs)
{
	char *curOne;
	char *curTwo;
	int minRemaining;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Generating pairs with %d insertions\n",
				numInsertions);
	}
	assert(direction == FORWARD || direction == REVERSE);

	/* Bounds on this will be different, since if we need 
	 * extra bases to compensate for the deletion */
	minRemaining = readLength-(offset+gap+matchLength*2);
	if(minRemaining <= 0) {
		/* Out of bounds.  Don't add anything. */
		return;
	}
	else if(minRemaining < numInsertions) {
		/* Adjust the number of insertions we can handle. We could
		 * also just use bases infront of offset, but the user
		 * specified offset for a reason. 
		 * */
		numInsertions = minRemaining;
	}

	/* Allocate memory */
	curOne = malloc(sizeof(char)*matchLength);
	assert(curOne!=NULL);
	curTwo = malloc(sizeof(char)*matchLength);
	assert(curTwo!=NULL);

	RGSeqPairGenerateInsertionsHelper(read,
			readLength,
			direction,
			offset,
			matchLength,
			gap,
			numInsertions,
			numInsertions,
			0,
			pairs,
			curOne,
			curTwo,
			0,
			0);

	/* Free memory */
	free(curOne);
	free(curTwo);
}

/* TODO */
/* NOTE: no error checking yet! */
void RGSeqPairGenerateInsertionsHelper(char *read,
		int readLength,
		char direction,
		int offset,
		int matchLength,
		int gap,
		int numInsertionsLeft,
		int numInsertions,
		int insertionOffset,
		RGSeqPair *pairs,
		char *curOne,
		char *curTwo,
		int curOneIndex,
		int curTwoIndex)
{
	int i;
	int remainingOne, remainingTwo;

	remainingOne = matchLength-curOneIndex;
	remainingTwo = matchLength-curTwoIndex;
	assert(direction == FORWARD || direction == REVERSE);

	if(numInsertionsLeft > 0) {
		/* No more to print */
		if(remainingOne <= 0 && remainingTwo <= 0 && numInsertionsLeft < numInsertions) {
			curOne[matchLength]='\0';
			curTwo[matchLength]='\0';
			/* Allocate memory */                                                                             
			pairs->numPairs++;
			pairs->indexOne = realloc(pairs->indexOne, sizeof(unsigned int)*(pairs->numPairs));
			assert(pairs->indexOne!=NULL);
			pairs->indexTwo = realloc(pairs->indexTwo, sizeof(unsigned int)*(pairs->numPairs));
			assert(pairs->indexTwo!=NULL);
			pairs->strand = realloc(pairs->strand, sizeof(char)*(pairs->numPairs));
			assert(pairs->strand!=NULL);
			pairs->offset = realloc(pairs->offset, sizeof(unsigned int)*(pairs->numPairs));
			/* Copy over */
			pairs->indexOne[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curOne, matchLength);
			pairs->indexTwo[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curTwo, matchLength);
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "Pair(insertions):%s[%d]\t%s[%d].\n",
						curOne,
						pairs->indexOne[pairs->numPairs-1],
						curTwo,
						pairs->indexTwo[pairs->numPairs-1]);
			}
			pairs->strand[pairs->numPairs-1] = direction;
			pairs->offset[pairs->numPairs-1] = offset;
			return;
		}
		else {
			/* try deleting a base */
			/* Don't delete if the previous base is the same as the one we
			 * are proposing to delete since we have already tried those permutations */
			if(remainingOne > 0) {
				if(curOneIndex == 0 || curOne[curOneIndex-1] != read[curOneIndex+insertionOffset]) { 
					RGSeqPairGenerateInsertionsHelper(read,
							readLength,
							direction,
							offset,
							matchLength,
							gap,
							numInsertionsLeft-1,
							numInsertions,
							insertionOffset+1,
							pairs,
							curOne,
							curTwo,
							curOneIndex,
							curTwoIndex);
				}
			}
			else if(remainingTwo > 0) {
				if(curTwoIndex == 0 || curTwo[curTwoIndex-1] != read[offset+gap+matchLength+curTwoIndex+insertionOffset]) {
					RGSeqPairGenerateInsertionsHelper(read,
							readLength,
							direction,
							offset,
							matchLength,
							gap,
							numInsertionsLeft-1,
							numInsertions,
							insertionOffset+1,
							pairs,
							curOne,
							curTwo,
							curOneIndex,
							curTwoIndex);
				}
			}
		}
		/* Try not deleting a base */
		/* Only do this if we haven't started deleting */ 
		if(numInsertionsLeft == numInsertions) {
			if(remainingOne > 0) {
				curOne[curOneIndex] = read[offset+curOneIndex+insertionOffset];
				/* Use on first read */
				RGSeqPairGenerateInsertionsHelper(read,
						readLength,
						direction,
						offset,
						matchLength,
						gap,
						numInsertionsLeft,
						numInsertions,
						insertionOffset,
						pairs,
						curOne,
						curTwo,
						curOneIndex+1,
						curTwoIndex);
			}
			else if(remainingTwo > 0) {
				curTwo[curTwoIndex] = read[offset+gap+matchLength+curTwoIndex+insertionOffset];

				/* Use on second read */            
				RGSeqPairGenerateInsertionsHelper(read,
						readLength,
						direction,
						offset,
						matchLength,
						gap,
						numInsertionsLeft,
						numInsertions,
						insertionOffset,
						pairs,
						curOne,
						curTwo,
						curOneIndex,
						curTwoIndex+1);
			}
		}
	}
	else {
		/* print remaining */

		if(remainingOne > 0) {
			/* Print rest of the first read */
			for(i=0;i<remainingOne;i++) {
				assert(curOneIndex<matchLength);
				curOne[curOneIndex] = read[offset+curOneIndex+insertionOffset];
				curOneIndex++;
			}
		}
		assert(curOneIndex == matchLength);
		if(remainingTwo > 0) {
			/* Print rest of the second read */
			for(i=0;i<remainingTwo;i++) {
				assert(curTwoIndex<matchLength);
				curTwo[curTwoIndex] = read[offset+gap+matchLength+curTwoIndex+insertionOffset];
				curTwoIndex++;
			}
		}
		assert(curTwoIndex == matchLength);
		curOne[matchLength]='\0';
		curTwo[matchLength]='\0';
		/* Allocate memory */                                                                             
		pairs->numPairs++;
		pairs->indexOne = realloc(pairs->indexOne, sizeof(unsigned int)*(pairs->numPairs));
		assert(pairs->indexOne!=NULL);
		pairs->indexTwo = realloc(pairs->indexTwo, sizeof(unsigned int)*(pairs->numPairs));
		assert(pairs->indexTwo!=NULL);
		pairs->strand = realloc(pairs->strand, sizeof(char)*(pairs->numPairs));
		assert(pairs->strand!=NULL);
		pairs->offset = realloc(pairs->offset, sizeof(unsigned int)*(pairs->numPairs));
		/* Copy over */
		pairs->indexOne[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curOne, matchLength);
		pairs->indexTwo[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curTwo, matchLength);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Pair(insertion):%s[%d]\t%s[%d].\n",
					curOne,
					pairs->indexOne[pairs->numPairs-1],
					curTwo,
					pairs->indexTwo[pairs->numPairs-1]);
		}
		pairs->strand[pairs->numPairs-1] = direction;
		pairs->offset[pairs->numPairs-1] = offset;
		return;
	}
}

/* TODO */
/* Note: Deletions have occured, so insert bases in the gap */
void RGSeqPairGenerateGapDeletions(char *read,
		int readLength,
		char direction,
		int offset,
		int matchLength,
		int gap,
		int numGapDeletions,
		RGSeqPair *pairs)
{
	char curOne[SEQUENCE_LENGTH];
	char curTwo[SEQUENCE_LENGTH];

	assert(direction == FORWARD || direction == REVERSE);
	if(VERBOSE>=DEBUG) {
		fprintf(stderr, "Generating pairs with %d gap deletions.\n",
				numGapDeletions);
	}

	if(gap <= 0) {
		/* Out of bounds.  Don't add anything. */
		return;
	}
	else {
		RGSeqPairGenerateGapDeletionsHelper(read,
				readLength,
				direction,
				offset,
				matchLength,
				gap, 
				numGapDeletions,
				pairs,
				curOne,
				curTwo);
	}
}

/* TODO */
/* NOTE: no error checking yet! */
/* We assume that all insertions in the gap are grouped together */
void RGSeqPairGenerateGapDeletionsHelper(char *read,
		int readLength,
		char direction,
		int offset,
		int matchLength, 
		int gap,
		int numGapDeletions,
		RGSeqPair *pairs,
		char *curOne,
		char *curTwo)
{
	int i, j, k;
	assert(gap > 0);

	/* Choose the starting position of the gap insertion */ 
	for(i=0;i<gap;i++) {
		for(j=1;j<=numGapDeletions;j++) { /* For the total # of gaps we will insert */
			if(numGapDeletions > gap) {
				/* We will insert more than there is found in the gap.  This will cause a combinatorial
				 * explosion since we can do a simplifying shift. */
				/* Ignore for now */
			}
			else {
				/* Copy over first read */
				for(k=offset+i;k<offset+i+matchLength;k++) {
					curOne[k-offset-i] = read[k];
				}
				curOne[matchLength]='\0';

				/* Copy over second read */
				for(k=offset+matchLength+gap-j;k<offset+2*matchLength+gap-j;k++) {
					curTwo[k-offset-matchLength-gap+j] = read[k];
				}
				curTwo[matchLength]='\0';

				/* Allocate memory */                                                                             
				pairs->numPairs++;
				pairs->indexOne = realloc(pairs->indexOne, sizeof(unsigned int)*(pairs->numPairs));
				assert(pairs->indexOne!=NULL);
				pairs->indexTwo = realloc(pairs->indexTwo, sizeof(unsigned int)*(pairs->numPairs));
				assert(pairs->indexTwo!=NULL);
				pairs->strand = realloc(pairs->strand, sizeof(char)*(pairs->numPairs));
				assert(pairs->strand!=NULL);
				pairs->offset = realloc(pairs->offset, sizeof(unsigned int)*(pairs->numPairs));
				/* Copy over */
				pairs->indexOne[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curOne, matchLength);
				pairs->indexTwo[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curTwo, matchLength);
				if(VERBOSE >= DEBUG) {
					fprintf(stderr, "Pair(gap deletion):%s[%d]\t%s[%d].\n",
							curOne,
							pairs->indexOne[pairs->numPairs-1],
							curTwo,
							pairs->indexTwo[pairs->numPairs-1]);
				}
				pairs->strand[pairs->numPairs-1] = direction;
				pairs->offset[pairs->numPairs-1] = offset+i;
			}
		}
	}
}

/* TODO */
/* Note: Insertions have occured, so delete bases */
void RGSeqPairGenerateGapInsertions(char *read,
		int readLength,
		char direction,
		int offset,
		int matchLength,
		int gap,
		int numGapInsertions,
		RGSeqPair *pairs)
{
	char *curOne;
	char *curTwo;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Generating pairs with %d insertions\n",
				numGapInsertions);
	}
	assert(direction == FORWARD || direction == REVERSE);

	/* Allocate memory */
	curOne = malloc(sizeof(char)*matchLength);
	assert(curOne!=NULL);
	curTwo = malloc(sizeof(char)*matchLength);
	assert(curTwo!=NULL);

	RGSeqPairGenerateGapInsertionsHelper(read,
			readLength,
			direction,
			offset,
			matchLength,
			gap,
			numGapInsertions,
			pairs,
			curOne,
			curTwo);

	/* Free memory */
	free(curOne);
	free(curTwo);
}

/* TODO */
/* NOTE: no error checking yet! */
void RGSeqPairGenerateGapInsertionsHelper(char *read,
		int readLength,
		char direction,
		int offset,
		int matchLength,
		int gap,
		int numGapInsertions,
		RGSeqPair *pairs,
		char *curOne,
		char *curTwo)
{
	int i, j;
	int startOffset = offset; /* The number of bases we can shift at the start of the read */
	int endOffset = readLength - offset - 2*matchLength - gap; /* The number of bases we can shift at the end of the read */
	assert(direction == FORWARD || direction == REVERSE);

	/* Choose the number of insertions */
	for(i=1;i<=numGapInsertions;i++) {
		/* Always choose to towards the end */
		int endShift = -1;
		int startShift = -1;
		if(endOffset >= i) {
			/* Shift the second sequence all end */
			startShift = 0;
			endShift = i;
		}
		else {
			startShift = endOffset - i;
			endShift = endOffset;
		}
		if(startShift <= startOffset && endShift <= endOffset) {
			/* Copy over first read */
			for(j=offset-startShift;j<offset-startShift+matchLength;j++) {
				assert(j>=0 && j<readLength);
				curOne[j-offset+startShift] = read[j];
			}
			curOne[matchLength]='\0';

			/* Copy over second read */
			for(j=offset+matchLength+gap+endShift;j<offset+2*matchLength+gap+endShift;j++) {
				assert(j>=0 && j<readLength);
				curTwo[j-offset-matchLength-gap-endShift] = read[j];
			}
			curTwo[matchLength]='\0';

			/* Allocate memory */                                                                             
			pairs->numPairs++;
			pairs->indexOne = realloc(pairs->indexOne, sizeof(unsigned int)*(pairs->numPairs));
			assert(pairs->indexOne!=NULL);
			pairs->indexTwo = realloc(pairs->indexTwo, sizeof(unsigned int)*(pairs->numPairs));
			assert(pairs->indexTwo!=NULL);
			pairs->strand = realloc(pairs->strand, sizeof(char)*(pairs->numPairs));
			assert(pairs->strand!=NULL);
			pairs->offset = realloc(pairs->offset, sizeof(unsigned int)*(pairs->numPairs));
			/* Copy over */
			pairs->indexOne[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curOne, matchLength);
			pairs->indexTwo[pairs->numPairs-1] = RGTreeGetIndexFromSequence(curTwo, matchLength);
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "Pair(gap deletion):%s[%d]\t%s[%d].\n",
						curOne,
						pairs->indexOne[pairs->numPairs-1],
						curTwo,
						pairs->indexTwo[pairs->numPairs-1]);
			}
			pairs->strand[pairs->numPairs-1] = direction;
			pairs->offset[pairs->numPairs-1] = offset-startShift;
		}
	}
}

/* TODO */
void RGSeqPairRemoveDuplicates(RGSeqPair *s)
{
	unsigned int i;
	unsigned int prevIndex=0;

	if(s->numPairs <= 0) {
		return;
	}

	/* Sort the data structure */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Sorting\n");
	}
	RGSeqPairQuickSort(s, 0, s->numPairs-1);
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Sorted!\n");
	}

	/* Remove duplicates */
	prevIndex=0;
	for(i=1;i<s->numPairs;i++) {
		if(RGSeqPairCompareAtIndex(s, prevIndex, s, i)==0) { 
			/* Ignore */
		}
		else {
			prevIndex++;
			/* Copy over to temporary pair */
			RGSeqPairCopyAtIndex(s, i, s, prevIndex);
		}
	}

	/* Reallocate pair */
	s->indexOne = realloc(s->indexOne, sizeof(unsigned int)*(prevIndex+1));
	assert(s->indexOne!=NULL);
	s->indexTwo = realloc(s->indexTwo, sizeof(unsigned int)*(prevIndex+1));
	assert(s->indexTwo!=NULL);
	s->strand = realloc(s->strand, sizeof(char)*(prevIndex+1));
	assert(s->strand!=NULL);
	s->offset = realloc(s->offset, sizeof(unsigned int)*(prevIndex+1));
	assert(s->offset!=NULL);
	s->numPairs = prevIndex+1;
}

/* TO DO */
void RGSeqPairQuickSort(RGSeqPair *s, int low, int high)
{
	unsigned int i;
	unsigned int pivot=-1;
	RGSeqPair temp;

	if(low < high) {
		/* Allocate memory for the temp RGSeqPair indexes and strand */
		temp.indexOne = malloc(sizeof(unsigned int));
		temp.indexTwo = malloc(sizeof(unsigned int));
		temp.strand = malloc(sizeof(char));
		temp.offset = malloc(sizeof(unsigned int));

		pivot = (low + high)/2;

		RGSeqPairCopyAtIndex(s, pivot, &temp, 0);
		RGSeqPairCopyAtIndex(s, high, s, pivot);
		RGSeqPairCopyAtIndex(&temp, 0, s, high);

		pivot = low;

		for(i=low;i<high;i++) {
			if(RGSeqPairCompareAtIndex(s, i, s, high) <= 0) {
				RGSeqPairCopyAtIndex(s, i, &temp, 0);
				RGSeqPairCopyAtIndex(s, pivot, s, i);
				RGSeqPairCopyAtIndex(&temp, 0, s, pivot);
				pivot++;
			}
		}
		RGSeqPairCopyAtIndex(s, pivot, &temp, 0);
		RGSeqPairCopyAtIndex(s, high, s, pivot);
		RGSeqPairCopyAtIndex(&temp, 0, s, high);

		/* Free memory before recursive call */
		free(temp.indexOne);
		free(temp.indexTwo);
		free(temp.strand);
		free(temp.offset);

		RGSeqPairQuickSort(s, low, pivot-1);
		RGSeqPairQuickSort(s, pivot+1, high);

	}
}

/* TODO */
void GetReverseCompliment(char *s,
		char *r,
		int length) 
{
	int i;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "read:");
		for(i=0;i<length;i++) {
			fprintf(stderr, "%c", s[i]);
		}
		fprintf(stderr, "\n");
	}

	/* Get reverse compliment read */
	for(i=length-1;i>=0;i--) {
		switch(s[length-1-i]) {
			case 'a':
				r[i] = 't';
				break;
			case 'c':
				r[i] = 'g';
				break;
			case 'g':
				r[i] = 'c';
				break;
			case 't':
				r[i] = 'a';
				break;
			default:
				fprintf(stderr, "Error.  In GetReverseCompliment, could not understand %c.  Terminating!\n", s[length-1-i]);
				exit(1);
				break;
		}
	}
	r[length]='\0';
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "reverse compliment:");
		for(i=0;i<length;i++) {
			fprintf(stderr, "%c", r[i]);
		}
		fprintf(stderr, "\n");
	}
}

int RGSeqPairCompareAtIndex(RGSeqPair *pOne, int iOne, RGSeqPair *pTwo, int iTwo) 
{
	if(pOne->indexOne[iOne] < pTwo->indexOne[iTwo]  ||
			(pOne->indexOne[iOne] == pTwo->indexOne[iTwo] && pOne->indexTwo[iOne] < pTwo->indexTwo[iTwo]) ||  
			(pOne->indexOne[iOne] == pTwo->indexOne[iTwo] && pOne->indexTwo[iOne] ==  pTwo->indexTwo[iTwo] && pOne->strand[iOne] <  pTwo->strand[iTwo]) ||
			(pOne->indexOne[iOne] == pTwo->indexOne[iTwo] && pOne->indexTwo[iOne] ==  pTwo->indexTwo[iTwo] && pOne->strand[iOne] <  pTwo->strand[iTwo] && pOne->offset[iOne] < pTwo->offset[iTwo])) {
		return -1;
	}
	else if(pOne->indexOne[iOne] == pTwo->indexOne[iTwo] && pOne->indexTwo[iOne] ==  pTwo->indexTwo[iTwo] && pOne->strand[iOne] == pTwo->strand[iTwo] && pOne->offset[iOne] == pTwo->offset[iTwo]) {
		return 0;
	}
	else {
		return 1;
	}
}

void RGSeqPairCopyAtIndex(RGSeqPair *src, int srcIndex, RGSeqPair *dest, int destIndex)
{
	dest->indexOne[destIndex] = src->indexOne[srcIndex];
	dest->indexTwo[destIndex] = src->indexTwo[srcIndex];
	dest->strand[destIndex] = src->strand[srcIndex];
	dest->offset[destIndex] = src->offset[srcIndex];
}
