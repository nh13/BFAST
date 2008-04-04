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
		char *sequence)
{
	unsigned char *indexForward=NULL;
	unsigned char *indexReverse=NULL;
	char reverseSequence[SEQUENCE_LENGTH];
	int numChars = (int)ceil((2.0/8.0*index->matchLength)/sizeof(unsigned char));

	/* Allocate memory for the indexes */
	indexForward = (unsigned char*)malloc(sizeof(unsigned char)*numChars);
	indexReverse = (unsigned char*)malloc(sizeof(unsigned char)*numChars);

	/* DON'T FORGET TO ALIGN TO BOTH FORWARD AND REVERSE STRANDS */
	assert(strlen(sequence)==index->matchLength);
	GetReverseCompliment(sequence, reverseSequence, strlen(sequence));

	/* Get indexes */
	RGIndexGetIndexFromSequence(sequence, index->matchLength, indexForward); 
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
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Exiting RGSeqPairFindMatchesInIndex.\n");
	}
}

/* TODO */
void RGSeqPairFindMatchesInTree(RGTree *tree, 
		RGMatch *match,
		char *sequence,
		int **offsets,
		int numOffsets,
		int numMismatches,
		int numInsertions,
		int numDeletions)
{
	int i;

	/* DON'T FORGET TO ALIGN TO BOTH FORWARD AND REVERSE STRANDS */
	RGSeqPair seqPairs;
	seqPairs.numPairs=0;
	seqPairs.indexOne=NULL;
	seqPairs.indexTwo=NULL;
	seqPairs.strand=NULL;
	seqPairs.offset=NULL;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "\nIn RGSeqPairFindMatchesInTree\n");
	}

	/* Note: we can speed this up by generating all pairs of l-mers to
	 * search for and then removing duplicated pairs. */
	RGSeqPairGeneratePairs(sequence,
			&seqPairs,
			offsets,
			numOffsets,
			tree->matchLength,
			tree->gap,
			numMismatches,
			numInsertions,
			numDeletions);
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Generated %d pairs.\n",
				seqPairs.numPairs);
	}
	if(seqPairs.numPairs <= 0) { 
		/* No pairs generated, return */
		return;
	}
	if(VERBOSE >= DEBUG) {
		for(i=0;i<seqPairs.numPairs;i++) { /* For each pair */
			fprintf(stderr, "%d:%d,%d,%c,%d\n",
					i+1, 
					seqPairs.indexOne[i],
					seqPairs.indexTwo[i],
					seqPairs.strand[i],
					seqPairs.offset[i]);
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
	for(i=0;i<seqPairs.numPairs;i++) { /* For each pair */
		if(seqPairs.strand[i] != FORWARD && seqPairs.strand[i] != REVERSE) {
			fprintf(stderr, "Error,  Direction not recognized [%c].  Terminating!\n",
					seqPairs.strand[i]);
			exit(1);
		}
		RGTreeGetMatches(tree, 
				seqPairs.indexOne[i],
				seqPairs.indexTwo[i],
				seqPairs.strand[i],
				seqPairs.offset[i],
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
void RGSeqPairGeneratePairs(char *sequence,
		RGSeqPair *seqPairs,
		int **offsets,
		int numOffsets,
		int matchLength,
		int gap,
		int numMismatches,
		int numInsertions,
		int numDeletions)
{
	/* DON'T FORGET TO ALIGN TO BOTH FORWARD AND REVERSE STRANDS */
	int i;
	int sequenceLength=strlen(sequence);
	char reverseSequence[SEQUENCE_LENGTH];

	/* Get the reverse compliment of the sequence */
	GetReverseCompliment(sequence, reverseSequence, sequenceLength);

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
		RGSeqPairGenerateMismatches(sequence,
				sequenceLength,
				FORWARD,
				(*offsets)[i],
				matchLength,
				gap,
				numMismatches,
				seqPairs);
		/* Reverse compliment */
		RGSeqPairGenerateMismatches(reverseSequence,
				sequenceLength,
				REVERSE,
				(*offsets)[i],
				matchLength,
				gap,
				numMismatches,
				seqPairs);

		/* Go through all insertions */
		/* Note: we allow only contiguous insertions of length up to
		 * numInsertions.  We also always start from the offset.  When
		 * we model insertions, we have to remove a certain number of
		 * bases in the sequence, and therefore we enumerate over all
		 * possible insertions in the entire sequence.
		 * */
		if(numInsertions > 0) {
			/* Forward */
			RGSeqPairGenerateInsertions(sequence,
					sequenceLength,
					FORWARD,
					(*offsets)[i],
					matchLength,
					gap,
					numInsertions,
					seqPairs);
			/* Reverse compliment */
			RGSeqPairGenerateInsertions(reverseSequence,
					sequenceLength,
					REVERSE,
					(*offsets)[i],
					matchLength,
					gap,
					numInsertions,
					seqPairs);
		}

		/* GO through all deletions */
		/* Note: we allow only contiguous deletions of length up to 
		 * numDeletions.  We also always start from the offset.  We
		 * must add base to the sequences, and therfore we enumerate
		 * over all possible deletions in the entire sequence.
		 * */
		if(numDeletions > 0) {
			/* Forward */
			RGSeqPairGenerateDeletions(sequence,
					sequenceLength,
					FORWARD,
					(*offsets)[i],
					matchLength,
					gap,
					numDeletions,
					seqPairs);
			/* Reverse compliment */
			RGSeqPairGenerateDeletions(reverseSequence,
					sequenceLength,
					REVERSE,
					(*offsets)[i],
					matchLength,
					gap,
					numDeletions,
					seqPairs);
		}
	}

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "All pairs generated\n");
	}

	/* Merge all pairs */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Calling RGSeqPairRemoveDuplicates\n");
	}
	RGSeqPairRemoveDuplicates(seqPairs);
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Exited from RGSeqPairRemoveDuplicates\n");
	}
}

/* TODO */
void RGSeqPairGenerateMismatches(char *seq,
		int seqLength,
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

	if(offset+gap+matchLength*2 > seqLength) {
		/* Out of bounds.  Don't add anything. */
		return;
	}
	else {
		/* Allocate memory */
		curOne = (char*)malloc(sizeof(char)*matchLength);
		curTwo = (char*)malloc(sizeof(char)*matchLength);

		RGSeqPairGenerateMismatchesHelper(seq,
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
void RGSeqPairGenerateMismatchesHelper(char *seq,
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
			pairs->indexOne = realloc(pairs->indexOne, sizeof(int)*(pairs->numPairs));
			pairs->indexTwo = realloc(pairs->indexTwo, sizeof(int)*(pairs->numPairs));
			pairs->strand = realloc(pairs->strand, sizeof(char)*(pairs->numPairs));
			pairs->offset = realloc(pairs->offset, sizeof(int)*(pairs->numPairs));
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
					/* Use on first sequence */
					if(seq[offset+numFirstPrinted] == ALPHABET[i]) {
						/* No mismatch */

						/* Keep going */
						RGSeqPairGenerateMismatchesHelper(seq,
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
						RGSeqPairGenerateMismatchesHelper(seq,
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

					/* Use on second sequence */
					if(seq[offset+matchLength+gap+numLastPrinted] == ALPHABET[i]) {
						/* No mismatch */

						/* Keep going */
						RGSeqPairGenerateMismatchesHelper(seq,
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
						RGSeqPairGenerateMismatchesHelper(seq,
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
			/* Print rest of the first sequence */
			for(i=offset+numFirstPrinted;i<offset+matchLength;i++) {                                    
				assert(curOneIndex<matchLength);
				curOne[curOneIndex] = seq[i];
				curOneIndex++;                  
			}
		}                                       
		if(numLastPrinted < matchLength) {      
			/* Print rest of the second sequence */                                     
			for(i=offset+matchLength+gap+numLastPrinted;i<offset+2*matchLength+gap;i++) {                               
				assert(curTwoIndex<matchLength);
				curTwo[curTwoIndex] = seq[i];                curTwoIndex++;
			}                             
		}
		curOne[matchLength]='\0';                                         
		curTwo[matchLength]='\0';
		/* Allocate memory */                                                                             
		pairs->numPairs++;
		pairs->indexOne = realloc(pairs->indexOne, sizeof(int)*(pairs->numPairs));
		pairs->indexTwo = realloc(pairs->indexTwo, sizeof(int)*(pairs->numPairs));
		pairs->strand = realloc(pairs->strand, sizeof(char)*(pairs->numPairs));
		pairs->offset = realloc(pairs->offset, sizeof(int)*(pairs->numPairs));
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
void RGSeqPairGenerateDeletions(char *seq,
		int seqLength,
		char direction,
		int offset,
		int matchLength,
		int gap,
		int numDeletions,
		RGSeqPair *pairs)
{
	char curOne[SEQUENCE_LENGTH];
	char curTwo[SEQUENCE_LENGTH];

	assert(direction == FORWARD || direction == REVERSE);
	if(VERBOSE>=DEBUG) {
		fprintf(stderr, "Generating pairs with %d deletions.\n",
				numDeletions);
	}

	if(offset+gap+matchLength*2 > seqLength) {
		/* Out of bounds.  Don't add anything. */
		return;
	}
	else {

		RGSeqPairGenerateDeletionsHelper(seq,
				seqLength,
				direction,
				offset,
				matchLength,
				gap, 
				numDeletions,
				0,
				pairs,
				curOne,
				curTwo,
				0,
				0);
	}
}

/* TODO */
/* NOTE: no error checking yet! */
void RGSeqPairGenerateDeletionsHelper(char *seq,
		int seqLength,
		char direction,
		int offset,
		int matchLength, 
		int gap,
		int numDeletionsLeft,
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
		if(remainingOne <= 0 && remainingTwo <= 0) {
			curOne[matchLength]='\0';
			curTwo[matchLength]='\0';
			/* Allocate memory */                                                                             
			pairs->numPairs++;
			pairs->indexOne = realloc(pairs->indexOne, sizeof(int)*(pairs->numPairs));
			pairs->indexTwo = realloc(pairs->indexTwo, sizeof(int)*(pairs->numPairs));
			pairs->strand = realloc(pairs->strand, sizeof(char)*(pairs->numPairs));
			pairs->offset = realloc(pairs->offset, sizeof(int)*(pairs->numPairs));
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
			for(i=0;i<ALPHABET_SIZE;i++) {
				if(remainingOne > 0) {
					curOne[curOneIndex] = ALPHABET[i];
					/* Use on first sequence */
					RGSeqPairGenerateDeletionsHelper(seq,
							seqLength,
							direction,
							offset,
							matchLength,
							gap,
							numDeletionsLeft-1,
							deletionOffset+1,
							pairs,
							curOne,
							curTwo,
							curOneIndex+1,
							curTwoIndex);
				}
				else if(remainingTwo > 0) {
					curTwo[curTwoIndex] = ALPHABET[i];

					/* Use on second sequence */
					RGSeqPairGenerateDeletionsHelper(seq,
							seqLength,
							direction,
							offset,
							matchLength,
							gap,
							numDeletionsLeft-1,
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
			/* Try not inserting a base */
			if(remainingOne > 0) {
				curOne[curOneIndex] = seq[offset+curOneIndex-deletionOffset];
				/* Use on first sequence */
				RGSeqPairGenerateDeletionsHelper(seq,
						seqLength,
						direction,
						offset,
						matchLength,
						gap,
						numDeletionsLeft,
						deletionOffset,
						pairs,
						curOne,
						curTwo,
						curOneIndex+1,
						curTwoIndex);
			}
			else if(remainingTwo > 0) {
				curTwo[curTwoIndex] = seq[offset+gap+matchLength+curTwoIndex-deletionOffset];

				/* Use on second sequence */
				RGSeqPairGenerateDeletionsHelper(seq,
						seqLength,
						direction,
						offset,
						matchLength,
						gap,
						numDeletionsLeft,
						deletionOffset,
						pairs,
						curOne,
						curTwo,
						curOneIndex,
						curTwoIndex+1);
			}
		}
	}
	else {
		/* print remaining */        if(remainingOne > 0) {
			/* Print rest of the first sequence */
			for(i=0;i<remainingOne;i++) {
				assert(curOneIndex<matchLength);
				curOne[curOneIndex] = seq[offset+curOneIndex-deletionOffset];
				curOneIndex++;
			}
		}
		assert(curOneIndex == matchLength);
		if(remainingTwo > 0) {
			/* Print rest of the second sequence */
			for(i=0;i<remainingTwo;i++) {
				assert(curTwoIndex<matchLength);
				curTwo[curTwoIndex] = seq[offset+gap+matchLength+curTwoIndex-deletionOffset];
				curTwoIndex++;
			}
		}
		assert(curTwoIndex == matchLength);
		curOne[matchLength]='\0';
		curTwo[matchLength]='\0';
		/* Allocate memory */                                                                             
		pairs->numPairs++;
		pairs->indexOne = realloc(pairs->indexOne, sizeof(int)*(pairs->numPairs));
		pairs->indexTwo = realloc(pairs->indexTwo, sizeof(int)*(pairs->numPairs));
		pairs->strand = realloc(pairs->strand, sizeof(char)*(pairs->numPairs));
		pairs->offset = realloc(pairs->offset, sizeof(int)*(pairs->numPairs));
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
void RGSeqPairGenerateInsertions(char *seq,
		int seqLength,
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
	minRemaining = seqLength-(offset+gap+matchLength*2);
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
	curOne = (char*)malloc(sizeof(char)*matchLength);
	curTwo = (char*)malloc(sizeof(char)*matchLength);

	RGSeqPairGenerateInsertionsHelper(seq,
			seqLength,
			direction,
			offset,
			matchLength,
			gap,
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
void RGSeqPairGenerateInsertionsHelper(char *seq,
		int seqLength,
		char direction,
		int offset,
		int matchLength,
		int gap,
		int numInsertionsLeft,
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
		if(remainingOne <= 0 && remainingTwo <= 0) {
			curOne[matchLength]='\0';
			curTwo[matchLength]='\0';
			/* Allocate memory */                                                                             
			pairs->numPairs++;
			pairs->indexOne = realloc(pairs->indexOne, sizeof(int)*(pairs->numPairs));
			pairs->indexTwo = realloc(pairs->indexTwo, sizeof(int)*(pairs->numPairs));
			pairs->strand = realloc(pairs->strand, sizeof(char)*(pairs->numPairs));
			pairs->offset = realloc(pairs->offset, sizeof(int)*(pairs->numPairs));
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
			RGSeqPairGenerateInsertionsHelper(seq,
					seqLength,
					direction,
					offset,
					matchLength,
					gap,
					numInsertionsLeft-1,
					insertionOffset+1,
					pairs,
					curOne,
					curTwo,
					curOneIndex,
					curTwoIndex);
		}
		/* Try not deleting a base */
		if(remainingOne > 0) {
			curOne[curOneIndex] = seq[offset+curOneIndex+insertionOffset];
			/* Use on first sequence */
			RGSeqPairGenerateInsertionsHelper(seq,
					seqLength,
					direction,
					offset,
					matchLength,
					gap,
					numInsertionsLeft,
					insertionOffset,
					pairs,
					curOne,
					curTwo,
					curOneIndex+1,
					curTwoIndex);
		}
		else if(remainingTwo > 0) {
			curTwo[curTwoIndex] = seq[offset+gap+matchLength+curTwoIndex+insertionOffset];

			/* Use on second sequence */            
			RGSeqPairGenerateInsertionsHelper(seq,
					seqLength,
					direction,
					offset,
					matchLength,
					gap,
					numInsertionsLeft,
					insertionOffset,
					pairs,
					curOne,
					curTwo,
					curOneIndex,
					curTwoIndex+1);
		}
	}
	else {
		/* print remaining */

		if(remainingOne > 0) {
			/* Print rest of the first sequence */
			for(i=0;i<remainingOne;i++) {
				assert(curOneIndex<matchLength);
				curOne[curOneIndex] = seq[offset+curOneIndex+insertionOffset];
				curOneIndex++;
			}
		}
		assert(curOneIndex == matchLength);
		if(remainingTwo > 0) {
			/* Print rest of the second sequence */
			for(i=0;i<remainingTwo;i++) {
				assert(curTwoIndex<matchLength);
				curTwo[curTwoIndex] = seq[offset+gap+matchLength+curTwoIndex+insertionOffset];
				curTwoIndex++;
			}
		}
		assert(curTwoIndex == matchLength);
		curOne[matchLength]='\0';
		curTwo[matchLength]='\0';
		/* Allocate memory */                                                                             
		pairs->numPairs++;
		pairs->indexOne = realloc(pairs->indexOne, sizeof(int)*(pairs->numPairs));
		pairs->indexTwo = realloc(pairs->indexTwo, sizeof(int)*(pairs->numPairs));
		pairs->strand = realloc(pairs->strand, sizeof(char)*(pairs->numPairs));
		pairs->offset = realloc(pairs->offset, sizeof(int)*(pairs->numPairs));
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
void RGSeqPairRemoveDuplicates(RGSeqPair *s)
{
	int i;
	RGSeqPair t;
	RGSeqPair prev;
	int curIndex=0;

	/* Sort the data structure */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Sorting\n");
	}
	RGSeqPairQuickSort(s, 0, s->numPairs-1);
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Sorted!\n");
	}

	/* Allocate prev */
	prev.indexOne = (int*)malloc(sizeof(int));
	prev.indexTwo = (int*)malloc(sizeof(int));
	prev.strand = (char*)malloc(sizeof(char));
	prev.offset = (int*)malloc(sizeof(int));
	prev.numPairs = 1;
	/* Allocate temporary */
	t.indexOne = (int*)malloc(sizeof(int)*(s->numPairs));
	t.indexTwo = (int*)malloc(sizeof(int)*(s->numPairs));
	t.strand = (char*)malloc(sizeof(char)*(s->numPairs));
	t.offset = (int*)malloc(sizeof(int)*(s->numPairs));
	t.numPairs = s->numPairs;

	/* Initialize prev */
	prev.indexOne[0] = -1;
	prev.indexTwo[0] = -1;
	prev.strand[0] = 'z';
	prev.offset[0]= 0;

	/* Remove duplicates */
	for(i=0;i<s->numPairs;i++) {
		if(RGSeqPairCompareAtIndex(&prev, 0, s, i)==0) { 
			/* Ignore */
		}
		else {
			/* Copy over to temporary pair */
			RGSeqPairCopyAtIndex(s, i, &t, curIndex);
			curIndex++;

			/* Save previous */
			RGSeqPairCopyAtIndex(s, i, &prev, 0);
		}
	}

	/* Reallocate pair */
	s->indexOne = (int*)realloc(s->indexOne, sizeof(int)*curIndex);
	s->indexTwo = (int*)realloc(s->indexTwo, sizeof(int)*curIndex);
	s->strand = (char*)realloc(s->strand, sizeof(char)*curIndex);
	s->offset = (int*)realloc(s->offset, sizeof(int)*curIndex);
	s->numPairs = curIndex;

	/* Copy over */
	for(i=0;i<curIndex;i++) {
		RGSeqPairCopyAtIndex(&t, i, s, i);
	}

	/* Free prev memory */
	free(prev.indexOne);
	free(prev.indexTwo);
	free(prev.strand);
	free(prev.offset);
	/* Free temporary memory */
	free(t.indexOne);
	free(t.indexTwo);
	free(t.strand);
	free(t.offset);
}

/* TO DO */
void RGSeqPairQuickSort(RGSeqPair *s, int low, int high)
{
	int i;
	int pivot=-1;
	RGSeqPair temp;

	if(low < high) {
		/* Allocate memory for the temp RGSeqPair indexes and strand */
		temp.indexOne = (int*)malloc(sizeof(int));
		temp.indexTwo = (int*)malloc(sizeof(int));
		temp.strand = (char*)malloc(sizeof(char));
		temp.offset = (int*)malloc(sizeof(int));

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
		fprintf(stderr, "sequence:");
		for(i=0;i<length;i++) {
			fprintf(stderr, "%c", s[i]);
		}
		fprintf(stderr, "\n");
	}

	/* Get reverse compliment sequence */
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
