#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <string.h>
#include "BLibDefinitions.h"
#include "RGSeqPair.h"
#include "RGTree.h"

char ALPHABET[ALPHABET_SIZE] = "ACGT";

/* TODO */
void RGSeqPairFindMatches(RGTree *tree, 
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

	/* Note: we can speed this up by generating all pairs of l-mers to
	 * search for and then removing duplicated pairs. */
	RGSeqPairGeneratePairs(sequence,
			&seqPairs,
			offsets,
			tree->depth,
			tree->gap,
			numOffsets,
			numMismatches,
			numInsertions,
			numDeletions);

	/* All pairs have been generated and are unique.  The next step
	 * is to use each pair and search in the tree to get all matches. 
	 * Store the results in match.
	 * */
	for(i=0;i<seqPairs.numPairs;i++) { /* For each pair */
		RGTreeGetMatches(tree, 
				seqPairs.indexOne[i],
				seqPairs.indexTwo[i],
				seqPairs.strand[i],
				match);
	}
	/* Remove duplicates from match */
	RGMatchRemoveDuplicates(match);

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

	GetReverseCompliment(sequence, reverseSequence, sequenceLength);

	/* Go through all offsets */
	for(i=0;i<numOffsets;i++) {
		/* Go through all mismatches */
		/* Note: we allow any number (including zero) of mismatches up to
		 * numMismatches.  We also note that when mismatches is zero, this 
		 * will just insert an unaltered pair, so no need to restrict calling
		 * to when numMismatches>0.
		 * */
		/* Forward */
		RGSeqPairGenerateMismatches(sequence,
				sequenceLength,
				'+',
				(*offsets)[i],
				matchLength,
				gap,
				numMismatches,
				seqPairs);
		/* Reverse compliment */
		RGSeqPairGenerateMismatches(reverseSequence,
				sequenceLength,
				'-',
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
					'-',
				(*offsets)[i],
					matchLength,
					gap,
					numInsertions,
					seqPairs);
			/* Reverse compliment */
			RGSeqPairGenerateInsertions(reverseSequence,
					sequenceLength,
					'+',
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
					'+',
				(*offsets)[i],
					matchLength,
					gap,
					numDeletions,
					seqPairs);
			/* Reverse compliment */
			RGSeqPairGenerateDeletions(reverseSequence,
					sequenceLength,
					'-',
				(*offsets)[i],
					matchLength,
					gap,
					numDeletions,
					seqPairs);
		}
	}

	/* Merge all pairs */
	RGSeqPairRemoveDuplicates(seqPairs);
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
			/* Copy over */
			pairs->indexOne[pairs->numPairs-1] = GetIndexFromSequence(curOne, matchLength);
			pairs->indexTwo[pairs->numPairs-1] = GetIndexFromSequence(curTwo, matchLength);
			pairs->strand[pairs->numPairs-1] = direction;
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
		/* Copy over */
		pairs->indexOne[pairs->numPairs-1] = GetIndexFromSequence(curOne, matchLength);
		pairs->indexTwo[pairs->numPairs-1] = GetIndexFromSequence(curTwo, matchLength);
		pairs->strand[pairs->numPairs-1] = direction;
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
			/* Copy over */
			pairs->indexOne[pairs->numPairs-1] = GetIndexFromSequence(curOne, matchLength);
			pairs->indexTwo[pairs->numPairs-1] = GetIndexFromSequence(curTwo, matchLength);
			pairs->strand[pairs->numPairs-1] = direction;
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
		/* Copy over */
		pairs->indexOne[pairs->numPairs-1] = GetIndexFromSequence(curOne, matchLength);
		pairs->indexTwo[pairs->numPairs-1] = GetIndexFromSequence(curTwo, matchLength);
		pairs->strand[pairs->numPairs-1] = direction;
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

	/* Bounds on this will be different, since if we need 
	 * extra bases to compensate for the deletion */
	minRemaining = seqLength-(offset+gap+matchLength*2+numInsertions);
	if(minRemaining < 0) {
		if(numInsertions - minRemaining <= 0) {
			/* Out of bounds.  Don't add anything. */
			return;
		}
		else {
			/* Adjust the number of insertions we can handle. We could
			 * also just use bases infront of offset, but the user
			 * specified offset for a reason. 
			 * */
			numInsertions -= minRemaining;
		}
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
			/* Copy over */
			pairs->indexOne[pairs->numPairs-1] = GetIndexFromSequence(curOne, matchLength);
			pairs->indexTwo[pairs->numPairs-1] = GetIndexFromSequence(curTwo, matchLength);
			pairs->strand[pairs->numPairs-1] = direction;
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
		/* Copy over */
		pairs->indexOne[pairs->numPairs-1] = GetIndexFromSequence(curOne, matchLength);
		pairs->indexTwo[pairs->numPairs-1] = GetIndexFromSequence(curTwo, matchLength);
		pairs->strand[pairs->numPairs-1] = direction;
		return;
	}
}

/* TODO */
void RGSeqPairRemoveDuplicates(RGSeqPair *s)
{
	int i;
	int prevIndexOne = -1;
	int prevIndexTwo = -1;
	char prevStrand = 'z';
	RGSeqPair t;
	int curIndex=0;

	/* Merge sort the data structure */
	RGSeqPairMergeSort(s, 0, s->numPairs-1);

	/* Allocate temporary */
	t.indexOne = (int*)malloc(sizeof(int)*(s->numPairs));
	t.indexTwo = (int*)malloc(sizeof(int)*(s->numPairs));
	t.strand = (char*)malloc(sizeof(char)*(s->numPairs));

	/* Remove duplicates */
	for(i=0;i<s->numPairs;i++) {
		if(s->indexOne[i] == prevIndexOne &&
				s->indexTwo[i] == prevIndexTwo &&
				s->strand[i] == prevStrand) {
			/* Ignore */
		}
		else {
			/* Copy over to temporary pair */
			t.indexOne[curIndex] = s->indexOne[i];
			t.indexTwo[curIndex] = s->indexTwo[i];
			t.strand[curIndex] = s->strand[i];
			curIndex++;

			/* Save previous */
			prevIndexOne = s->indexOne[i];
			prevIndexTwo = s->indexTwo[i];
			prevStrand = s->strand[i];
		}
	}

	/* Reallocate pair */
	s->indexOne = (int*)realloc(s->indexOne, sizeof(int)*curIndex);
	s->indexTwo = (int*)realloc(s->indexTwo, sizeof(int)*curIndex);
	s->strand = (char*)realloc(s->strand, sizeof(char)*curIndex);

	/* Copy over */
	for(i=0;i<curIndex;i++) {
		s->indexOne[i] = t.indexOne[i];
		s->indexTwo[i] = t.indexTwo[i];
		s->strand[i] = t.strand[i];
	}

	/* Free temporary memory */
	free(t.indexOne);
	free(t.indexTwo);
	free(t.strand);
}

/* TO DO */
void RGSeqPairMergeSort(RGSeqPair *s, int low, int high)
{
	/* NOTE: when high-low < 20 we could use selection sort since it is faster 
	 * on smaller lengths. */

	/* MergeSort! */
	int mid = (low + high)/2;
	int start_upper = mid + 1;
	int end_upper = high;
	int start_lower = low;
	int end_lower = mid;
	int ctr, i;

	int *tempIndexOne;
	int *tempIndexTwo;
	char *tempStrand;

	if(low >= high) {
		return;
	}

	/* Partition the list into two lists and then sort them recursively */
	RGSeqPairMergeSort(s, low, mid);
	RGSeqPairMergeSort(s, mid+1, high);

	/* Allocate temporary memory */
	tempIndexOne = (int*)malloc((high-low+1)*sizeof(int));
	tempIndexTwo = (int*)malloc((high-low+1)*sizeof(int));
	tempStrand = (char*)malloc((high-low+1)*sizeof(char));

	/* Merge the two lists */
	ctr = 0;
	while( (start_lower<=end_lower) && (start_upper<=end_upper) )
	{
		if(s->indexOne[start_lower] < s->indexOne[start_upper] ||
				(s->indexOne[start_lower] == s->indexOne[start_upper] && s->indexTwo[start_lower] < s->indexTwo[start_upper]) ||
				(s->indexOne[start_lower] == s->indexOne[start_upper] && s->indexTwo[start_lower] == s->indexTwo[start_upper] && s->strand[start_lower] <= s->strand[start_upper])) { 
			tempIndexOne[ctr] = s->indexOne[start_lower];
			tempIndexTwo[ctr] = s->indexTwo[start_lower];
			tempStrand[ctr] = s->strand[start_lower];
			start_lower++;
		}
		else {
			tempIndexOne[ctr] = s->indexOne[start_upper];
			tempIndexTwo[ctr] = s->indexTwo[start_upper];
			tempStrand[ctr] = s->strand[start_upper];
			start_upper++;
		}
		ctr++;
		assert(ctr<high-low+1);
	}
	if(start_lower<=end_lower) {
		while(start_lower<=end_lower) {
			tempIndexOne[ctr] = s->indexOne[start_lower];
			tempIndexTwo[ctr] = s->indexTwo[start_lower];
			tempStrand[ctr] = s->strand[start_lower];
			ctr++;
			start_lower++;
		}
	}
	else {
		while(start_upper<=end_upper) {
			tempIndexOne[ctr] = s->indexOne[start_upper];
			tempIndexTwo[ctr] = s->indexTwo[start_upper];
			tempStrand[ctr] = s->strand[start_upper];
			ctr++;
			start_upper++;
		}
	}
	for(i=low, ctr=0;i<=high;i++, ctr++) {
		s->indexOne[i] = tempIndexOne[ctr];
		s->indexTwo[i] = tempIndexTwo[ctr];
		s->strand[i] = tempStrand[ctr];
	}
}


/* TODO */
void GetReverseCompliment(char *s,
		char *r,
		int length) 
{
	int i;

	/* Get reverse compliment sequence */
	for(i=length-1;i>=0;i--) {
		switch(s[length-1-i]) {
			case 'A':
				r[i] = 'T';
				break;
			case 'C':
				r[i] = 'G';
				break;
			case 'G':
				r[i] = 'C';
				break;
			case 'T':
				r[i] = 'A';
				break;
			default:
				fprintf(stderr, "Error.  In GetReverseCompliment, could not understand %c.  Terminating!\n", s[length-1-i]);
				exit(1);
				break;
		}
	}
	r[length]='\0';
}
