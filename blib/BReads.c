#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <string.h>
#include <assert.h>
#include "BLibDefinitions.h"
#include "BError.h"
#include "BLib.h"
#include "BIndex.h"
#include "BMatch.h"
#include "BRanges.h"
#include "BReads.h"

/* TODO 
 * remove directionality in BReads
 * remove directionality in generating reads.
 * */

char ALPHABET[ALPHABET_SIZE] = "acgt";

/* TODO */
void BReadsFindMatches(BIndex *index, 
		BBinary *rg,
		BMatch *match,
		int *offsets,
		int numOffsets,
		int space,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions,
		int maxKeyMatches,
		int maxNumMatches)
{
	int64_t i;
	int64_t numEntries = 0;
	int read->length=0;
	BString read, reverseRead;
	BReads reads;
	BRanges ranges;

	/* Initialize */
	BReadsInitialize(&reads);
	BRangesInitialize(&ranges);

	BStringInitialize(&read);
	BStringInitialize(&reverseRead);

	/* Remove adaptor and first color if we are in color space */
	if(space==ColorSpace) {
		/* Allocate */
		BStringAllocate(&read, match->read->length-2);
		BStringAllocate(&reverseRead, match->read->length-2);

		/* First letter is adapter, second letter is the color (unusable) */
		for(i=2;i<match->read.length;i++) {
			read.string[i-2] = match->read.string[i];
		}
		read->length = match->read->length-2;
		read.string[read->length] = '\0';

		/* In color space, the reverse compliment is just the reverse of the colors */
		BStringReverse(&reverseRead, &read);
		
		/* Update the colors in the read */
		ConvertColorsToStorage(&read, &read->length);
		ConvertColorsToStorage(&reverseRead, &read->length);
	}
	else {
		assert(space==NTSpace);
		/* Copy over */
		BStringCopy(&read, &match->read);
		/* Get the reverse compliment */
		BStringReverseCompliment(&reverseRead, &read);
	}

	/* Generate reads */
	BReadsGenerateReads(read,
			index,
			&reads,
			FORWARD,
			offsets,
			numOffsets,
			space,
			numMismatches,
			numInsertions,
			numDeletions,
			numGapInsertions,
			numGapDeletions);

	/* Generate reads */
	BReadsGenerateReads(reverseRead,
			index,
			&reads,
			REVERSE,
			offsets,
			numOffsets,
			space,
			numMismatches,
			numInsertions,
			numDeletions,
			numGapInsertions,
			numGapDeletions);

	/* Merge all reads */
	/* This may be necessary for a large number of generated reads, but omit for now */
	/*
	   if(numMismatches > 0 || 
	   numInsertions > 0 ||
	   numDeletions > 0 ||
	   numGapInsertions > 0 ||
	   numGapDeletions > 0) {
	   BReadsRemoveDuplicates(reads);
	   }
	   */

	/* Get the matches */
	for(i=0;i<reads.numReads && match->maxReached == 0;i++) {
		if(1==WillGenerateValidKey(index, &reads.reads[i])) {
			BIndexGetRanges(index, 
					rg,
					&reads.reads[i],
					reads.strand[i],
					reads.offset[i],
					maxKeyMatches,
					&ranges);
		}
	}

	/* Remove duplicate ranges */
	/* This exploits the fact that the ranges are non-overlapping */
	BRangesRemoveDuplicates(&ranges);

	/* Get the total number of matches we wish to copy over */
	for(i=0;i<ranges.numEntries;i++) {
		numEntries += ranges.endIndex[i] - ranges.startIndex[i] + 1;
	}

	/* Only copy to rgmatches if there are matches and fewer matches than
	 * max Matches */
	if(numEntries > maxNumMatches) {
		match->maxReached = 1;
	}
	else if(0 < numEntries) {
		/* Allocate memory for the matches */
		BMatchAllocate(match, numEntries);

		/* Transfer ranges to matches */
		BRangesCopyToBMatch(&ranges,
				index,
				match);
	}

	/* Remove duplicates */
	BMatchRemoveDuplicates(match,
			maxNumMatches);

	/* In color space we removed the first base/color so we need to 
	 * decrement the positions by one.
	 * */
	if(space == ColorSpace) {
		for(i=0;i<match->numEntries;i++) {
			match->positions[i]--;
		}
	}

	/* Free memory */
	BRangesFree(&ranges);
	BReadsFree(&reads);
	BStringFree(&read);
	BStringFree(&reverseRead);
}

/* TODO */
/* We may want to include enumeration of SNPs in color space */
void BReadsGenerateReads(BString *read,
		BIndex *index,
		BReads *reads,
		char direction,
		int *offsets,
		int numOffsets,
		int space,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions)
{
	int i;

	/* Go through all offsets */
	for(i=0;i<numOffsets && (read->length - offsets[i]) >= index->width;i++) {

		/* Insert the perfect match */
		BReadsGeneratePerfectMatch(read,
				direction,
				offsets[i],
				index,
				reads);

		/* Go through all mismatches */
		/* Note: we allow any number (including zero) of mismatches up to
		 * numMismatches.  
		 * */
		if(numMismatches > 0) {
			BReadsGenerateMismatches(read,
					direction,
					offsets[i],
					numMismatches,
					index,
					reads);
		}

		/* Go through all deletions */
		/* Note: we allow only contiguous deletions of length up to 
		 * numDeletions.  We also always start from the offset.  We
		 * must add base to the reads, and therfore we enumerate
		 * over all possible deletions in the entire read.
		 * */
		if(numDeletions > 0) {
			BReadsGenerateDeletions(read,
					direction,
					offsets[i],
					numDeletions,
					index,
					reads);
		}

		/* Go through all insertions */
		/* Note: we allow only contiguous insertions of length up to
		 * numInsertions.  We also always start from the offset.  When
		 * we model insertions, we have to remove a certain number of
		 * bases in the read, and therefore we enumerate over all
		 * possible insertions in the entire read.
		 * */
		if(numInsertions > 0) {
			/* Forward */
			BReadsGenerateInsertions(read,
					direction,
					offsets[i],
					numInsertions,
					index,
					reads);
		}

		/* Go through all possible insertions in the gap between
		 * the pair of l-mers.  If there is a gap insertion,
		 * then we will delete bases in the gap.
		 * */
		if(numGapInsertions > 0) {
			BReadsGenerateGapInsertions(read,
					direction,
					offsets[i],
					numGapInsertions,
					index,
					reads);
		}

		/* Go through all possible deletions in the gap between
		 * the pair of l-mers.  If there is a gap deletion, 
		 * then we will add bases to the gap.
		 * */
		if(numGapDeletions > 0) {
			/* Forward */
			BReadsGenerateGapDeletions(read,
					direction,
					offsets[i],
					numGapDeletions,
					index,
					reads);
		}
	}

}

/* TODO */
void BReadsGeneratePerfectMatch(char *read,
		char direction,
		int offset,
		BIndex *index,
		BReads *reads)
{
	int32_t i;
	char *curRead;

	/* Check bounds */
	if(read->length < index->width + offset) {
		return;
	}

	curRead = malloc(sizeof(char)*(index->width+1));
	if(NULL == curRead) {
		PrintError("BReadsPerfectMatchesHelper",
				"curRead",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Copy over */
	for(i=offset;i<index->width + offset;i++) {
		curRead[i-offset] = read->string[i];
	}
	curRead[index->width]='\0';

	/* Append */
	BReadsAppend(reads, curRead, index->width, direction, offset);

	free(curRead);
	curRead = NULL;
}

/* TODO */
void BReadsGenerateMismatches(char *read,
		int read->length,
		char direction,
		int offset,
		int numMismatches,
		BIndex *index,
		BReads *reads)
{
	BString curRead;

	BStringInitialize(&curRead);

	/* Check bounds */
	if(read->length < index->width+offset) {
		return;
	}

	/* Allocate memory */
	curRead = malloc(sizeof(char)*(index->width+1));
	if(NULL == curRead) {
		PrintError("BReadsGenerateMismatches",
				"curRead",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	BReadsGenerateMismatchesHelper(read,
			direction,
			offset,
			numMismatches,
			&curRead,
			0,
			index,
			reads);

	/* Free memory */
	free(curRead);
}

/* TODO */
void BReadsGenerateMismatchesHelper(char *read,
		char direction,
		int offset,
		int numMismatchesLeft,
		BString *curRead,
		int curReadIndex,
		BIndex *index,
		BReads *reads)
{
	int i;

	if(curReadIndex > index->width) {
		return;
	}

	if(numMismatchesLeft > 0) {
		assert(curReadIndex <= index->width);
		/* No more to print */
		if(curReadIndex == index->width) {
			curRead[index->width]='\0';
			BReadsAppend(reads, curRead, index->width, direction, offset);
			return;
		}
		else {
			assert(curReadIndex < index->width);
			/* use mismatches */
			for(i=0;i<ALPHABET_SIZE;i++) {
				int tempReadIndex = curReadIndex;
				while(index->mask[tempReadIndex] == 0 &&
						tempReadIndex < read->length) {
					curRead[tempReadIndex] = read->string[offset+tempReadIndex];
					tempReadIndex++;
				}
				assert(tempReadIndex < read->length);
				if(index->mask[tempReadIndex] == 0) {
					return;
				}
				curRead[tempReadIndex] = ALPHABET[i];
				if(read->string[offset+tempReadIndex] == ALPHABET[i]) {
					/* No mismatch */
					/* Keep going */
					BReadsGenerateMismatchesHelper(read,
							direction,
							offset,
							numMismatchesLeft,
							curRead,
							tempReadIndex+1,
							index,
							reads);
				}
				else {
					/* Mismatch */

					/* Keep going */
					BReadsGenerateMismatchesHelper(read,
							direction,
							offset,
							numMismatchesLeft-1,
							curRead,
							tempReadIndex+1,
							index,
							reads);
				}
			}
		}
	}
	else {
		/* print remaining */                                           
		while(curReadIndex < index->width) {
			curRead[curReadIndex] = read->string[curReadIndex+offset];
			curReadIndex++;

		}
		assert(curReadIndex == index->width);
		curRead[index->width]='\0';
		/* Append */
		BReadsAppend(reads, curRead, index->width, direction, offset);
	}
}

/* TODO */
/* Note: Deletions have occured, so insert bases */
void BReadsGenerateDeletions(char *read,
		char direction,
		int offset,
		int numDeletions,
		BIndex *index,
		BReads *reads)
{
	BString curRead;
	BStringInitialize(&curRead);
	int i;

	/* Allocate memory */
	curRead = malloc(sizeof(char)*(index->width+1));
	if(NULL == curRead) {
		PrintError("BReadsGenerateDeletions",
				"curRead",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	for(i=1;i<=numDeletions;i++) {
		BReadsGenerateDeletionsHelper(read,
				direction,
				offset,
				i,
				i,
				0,
				curRead,
				0,
				index,
				reads);
	}

	/* Free memory */
	free(curRead);
}

/* TODO */
/* NOTE: no error checking yet! */
/* Deletion occured, so insert bases */
void BReadsGenerateDeletionsHelper(char *read,
		char direction,
		int offset,
		int numDeletionsLeft,
		int numDeletions,
		int deletionOffset,
		char *curRead,
		int curReadIndex,
		BIndex *index,
		BReads *reads)
{
	int i;

	if(curReadIndex > index->width) {
		return;
	}

	if(numDeletionsLeft > 0) {
		/* No more to print */
		if(curReadIndex == index->width) {
			if(numDeletionsLeft != numDeletions) {
				curRead[index->width]='\0';
				/* Append */
				BReadsAppend(reads, curRead, index->width, direction, offset);
			}
			return;
		}
		else {
			/* Update curReadIndex etc. based on current tile */
			int tempReadIndex = curReadIndex;
			while(index->mask[tempReadIndex] == 0 &&
					tempReadIndex < read->length) {
				curRead[tempReadIndex] = read->string[offset+tempReadIndex-deletionOffset];
				tempReadIndex++;
			}
			assert(tempReadIndex < read->length);
			if(index->mask[tempReadIndex] == 0) {
				return;
			}
			/* try inserting a base - do not insert at the beginning or the end of a read */
			if(curReadIndex > 0 && curReadIndex < read->length-1) {
				for(i=0;i<ALPHABET_SIZE;i++) {
					curRead[curReadIndex] = ALPHABET[i];
					/* Use on first read */
					BReadsGenerateDeletionsHelper(read,
							direction,
							offset,
							numDeletionsLeft-1,
							numDeletions,
							deletionOffset+1,
							curRead,
							tempReadIndex+1,
							index,
							reads);
				}
			}
			/* This will enforce that insertions occur together */
			if(numDeletionsLeft == numDeletions) {
				/* Try not inserting a base */
				curRead[curReadIndex] = read->string[offset+curReadIndex-deletionOffset];
				BReadsGenerateDeletionsHelper(read,
						direction,
						offset,
						numDeletionsLeft,
						numDeletions,
						deletionOffset,
						curRead,
						tempReadIndex+1,
						index,
						reads);
			}
		}
	}
	else {
		/* print remaining */                                           
		while(curReadIndex < index->width) {
			curRead[curReadIndex] = read->string[curReadIndex+offset-deletionOffset];
			curReadIndex++;
		}
		curRead[index->width]='\0';
		assert(curReadIndex == index->width);
		/* Append */
		BReadsAppend(reads, curRead, index->width, direction, offset);
		return;

	}
}

/* TODO */
/* Note: Insertions have occured, so delete bases */
void BReadsGenerateInsertions(char *read,
		char direction,
		int offset,
		int numInsertions,
		BIndex *index,
		BReads *reads)
{
	BString curRead;
	int maxNumInsertions = 0;
	int i;

	BStringInitialize(&curRead);

	/* Get the total number of insertions (delete bases) possible */
	maxNumInsertions = read->length - index->width;

	if(maxNumInsertions <= 0) {
		/* Cannot delete any bases */
		return;
	}
	else if(maxNumInsertions < numInsertions) {
		/* Can only delete a certain # of bases */
		numInsertions = maxNumInsertions;
	}

	/* Allocate memory */
	curRead = malloc(sizeof(char)*(index->width+1));
	if(NULL == curRead) {
		PrintError("BReadsGenerateInsertions",
				"curRead",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Try up to the number of insertions */
	for(i=1;i<=numInsertions;i++) {
		BReadsGenerateInsertionsHelper(read,
				direction,
				offset,
				i,
				i,
				0,
				curRead,
				0,
				index,
				reads);
	}

	/* Free memory */
	free(curRead);
}

/* TODO */
/* NOTE: no error checking yet! */
/* Try deleting bases from the read */
void BReadsGenerateInsertionsHelper(char *read,
		char direction,
		int offset,
		int numInsertionsLeft,
		int numInsertions,
		int insertionOffset,
		char *curRead,
		int curReadIndex,
		BIndex *index,
		BReads *reads)
{
	if(curReadIndex > index->width) {
		return;
	}

	if(numInsertionsLeft > 0) {
		/* No more to print */
		if(curReadIndex >= index->width) {
			if(numInsertionsLeft != numInsertions) {
				curRead[index->width]='\0';
				/* Append */
				BReadsAppend(reads, curRead, index->width, direction, offset);
			}
			return;
		}
		/* try deleting a base */
		if(curReadIndex == 0 || 
				read->string[curReadIndex-1] != read->string[curReadIndex]) {
			BReadsGenerateInsertionsHelper(read,
					direction,
					offset,
					numInsertionsLeft-1,
					numInsertions,
					insertionOffset+1,
					curRead,
					curReadIndex,
					index,
					reads);
		}
		/* Try not deleting a base */
		/* Only do this if we haven't started deleting */ 
		if(numInsertionsLeft == numInsertions) {
			int tempReadIndex = curReadIndex;
			while(index->mask[tempReadIndex] == 0 &&
					tempReadIndex < read->length) {
				curRead[tempReadIndex] = read->string[offset+tempReadIndex+insertionOffset];
				tempReadIndex++;
			}
			assert(tempReadIndex < read->length);
			if(index->mask[tempReadIndex] == 0) {
				return;
			}
			curRead[tempReadIndex] = read->string[offset+tempReadIndex+insertionOffset];
			BReadsGenerateInsertionsHelper(read,
					direction,
					offset,
					numInsertionsLeft,
					numInsertions,
					insertionOffset,
					curRead,
					tempReadIndex+1,
					index,
					reads);
		}
	}
	else {
		/* print remaining */                                           
		while(curReadIndex<index->width) {
			curRead[curReadIndex] = read->string[curReadIndex+offset+insertionOffset];
			curReadIndex++;
		}
		curRead[index->width]='\0';
		assert(curReadIndex == index->width);
		/* Append */
		BReadsAppend(reads, curRead, index->width, direction, offset);
		return;
	}
}

/* TODO */
/* Note: Deletions have occured, so insert bases in the gaps */
void BReadsGenerateGapDeletions(char *read,
		char direction,
		int offset,
		int numGapDeletions,
		BIndex *index,
		BReads *reads)
{
	BString curRead;
	BStringInitialize(&curRead);

	if(numGapDeletions <= 0) {
		return;
	}
	/* Allocate memory */
	curRead = malloc(sizeof(char)*(read->length+1));
	if(NULL == curRead) {
		PrintError("BReadsGenerateGapDeletions",
				"curRead",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	BReadsGenerateGapDeletionsHelper(read,
			direction,
			offset,
			numGapDeletions,
			curRead,
			index,
			reads);

	/* Free memory */
	free(curRead);
}

/* TODO */
/* NOTE: no error checking yet! */
/* We assume that all insertions in the gap are grouped together */
void BReadsGenerateGapDeletionsHelper(char *read,
		char direction,
		int offset,
		int numGapDeletions,
		char *curRead,
		BIndex *index,
		BReads *reads)
{
	int i, j, k;
	int readPos;
	int curReadPos;

	int prevType=0;
	int gapLength=0;
	int numToInsert=0;

	/* Choose a gap to insert bases */
	readPos = offset;
	curReadPos = 0;
	for(i=0;i<index->width;i++) {
		/* Previous base was a 1 and we are not at a 0 */
		if(prevType == 1 &&
				index->mask[i] == 0) {

			/* Get the gap length */
			gapLength = 0;
			for(j=i;index->mask[j]==0;j++) {
				gapLength++;
			}
			/* Insert the min(gapLength, numGapDeletions) bases into the gap, since
			 * we require the bases in the gap */
			numToInsert = (gapLength < numGapDeletions)?gapLength:numGapDeletions;

			/* j is the current number of bases we are inserting */
			for(j=1;j<=numToInsert;j++) {
				int tempCurReadPos = curReadPos;
				int tempReadPos = readPos;

				/* Insert bases into the gap */
				for(k=0;k<j;k++) {
					curRead[tempCurReadPos] = NULL_LETTER;
					tempCurReadPos++;
				}
				/* Copy over the bases after the inserted bases */
				while(tempCurReadPos < index->width) {
					curRead[tempCurReadPos] = read->string[tempReadPos];
					tempCurReadPos++;
					tempReadPos++;
				}
				assert(tempCurReadPos == index->width);
				curRead[index->width]='\0';
				/* Append */
				BReadsAppend(reads, curRead, index->width, direction, offset);
			}
		}
		/* Copy base */
		curRead[curReadPos] = read->string[readPos];
		curReadPos++;
		readPos++;
		prevType = index->mask[i];
	}
}

/* TODO */
/* Note: Insertions have occured, so delete bases */
void BReadsGenerateGapInsertions(char *read,
		char direction,
		int offset,
		int numGapInsertions,
		BIndex *index,
		BReads *reads)
{
	BString curRead;
	BStringInitialize(&curRead);

	if(numGapInsertions <= 0) {
		/* Out of bounds.  Don't add anything. */
		return;
	}

	/* Allocate memory */
	curRead = malloc(sizeof(char)*(read->length+1));
	if(NULL == curRead) {
		PrintError("BReadsGenerateGapInsertions",
				"curRead",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	BReadsGenerateGapInsertionsHelper(read,
			direction,
			offset,
			numGapInsertions,
			curRead,
			index,
			reads);

	/* Free memory */
	free(curRead);
}

/* TODO */
/* NOTE: no error checking yet! */
/* Delete bases in the gaps */
void BReadsGenerateGapInsertionsHelper(char *read,
		char direction,
		int offset,
		int numGapInsertions,
		char *curRead,
		BIndex *index,
		BReads *reads)
{
	int i, j;
	int readPos;
	int curReadPos;

	int gapLength = 0;
	int prevType = 1;
	int numToDelete = 0;

	/* Choose a gap to try to remove bases */
	readPos = offset;
	curReadPos = 0;
	for(i=0;i<index->width;i++) {
		/* Previous base was a 1 and we are not at a 0 */
		if(prevType == 1 &&
				index->mask[i] == 0) {

			/* Get the gap length */
			gapLength = 0;
			for(j=i;index->mask[j]==0;j++) {
				gapLength++;
			}

			/* Delete min(gapLength, (read->length - readPos) - (index->width - curReadPos)).  We can only 
			 * remove as many bases as we can shift into the gap. 
			 */
			numToDelete = (gapLength < (read->length - readPos) - (index->width - curReadPos))?gapLength:(read->length - readPos) - (index->width - curReadPos);
			/* j is the current number of bases we are deleting */
			for(j=1;j<=numToDelete;j++) {
				int tempCurReadPos = curReadPos;
				int tempReadPos = readPos+j;

				/* Copy over the bases after the deleted bases */
				while(tempCurReadPos < index->width) {
					assert(tempReadPos < read->length);
					curRead[tempCurReadPos] = read->string[tempReadPos];
					tempCurReadPos++;
					tempReadPos++;
				}
				assert(tempCurReadPos == index->width);
				curRead[index->width]='\0';
				/* Append */
				BReadsAppend(reads, curRead, index->width, direction, offset);
			}
		}
		/* Copy base */
		curRead[curReadPos] = read->string[readPos];
		curReadPos++;
		readPos++;
		prevType = index->mask[i];
	}
}

/* TODO */
void BReadsRemoveDuplicates(BReads *s)
{
	int32_t i;
	int32_t prevIndex=0;

	if(s->numReads <= 0) {
		return;
	}

	/* Sort the data structure */
	BReadsQuickSort(s, 0, s->numReads-1);

	/* Remove duplicates */
	prevIndex=0;
	for(i=1;i<s->numReads;i++) {
		if(BReadsCompareAtIndex(s, prevIndex, s, i)==0) { 
			/* Ignore */
		}
		else {
			prevIndex++;
			/* Copy over to temporary pair */
			BReadsCopyAtIndex(s, i, s, prevIndex);
		}
	}

	/* Reallocate pair */
	BReadsReallocate(s, prevIndex+1);

}

/* TO DO */
void BReadsQuickSort(BReads *s, int low, int high)
{
	int32_t i;
	int32_t pivot=-1;
	BReads *temp;

	if(low < high) {
		/* Allocate memory for the temp BReads indexes */
		temp = malloc(sizeof(BReads));
		if(NULL == temp) {
			PrintError("BReadsQuickSort",
					"temp",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		BReadsInitialize(temp);
		BReadsAllocate(temp, 1);
		assert(temp->numReads == 1);

		pivot = (low + high)/2;

		BReadsCopyAtIndex(s, pivot, temp, 0);
		BReadsCopyAtIndex(s, high, s, pivot);
		BReadsCopyAtIndex(temp, 0, s, high);

		pivot = low;

		for(i=low;i<high;i++) {
			if(BReadsCompareAtIndex(s, i, s, high) <= 0) {
				if(i!=pivot) {
					BReadsCopyAtIndex(s, i, temp, 0);
					BReadsCopyAtIndex(s, pivot, s, i);
					BReadsCopyAtIndex(temp, 0, s, pivot);
				}
				pivot++;
			}
		}
		BReadsCopyAtIndex(s, pivot, temp, 0);
		BReadsCopyAtIndex(s, high, s, pivot);
		BReadsCopyAtIndex(temp, 0, s, high);

		/* Free memory before recursive call */
		assert(temp->numReads == 1);
		BReadsFree(temp);
		free(temp);
		temp=NULL;

		BReadsQuickSort(s, low, pivot-1);
		BReadsQuickSort(s, pivot+1, high);
	}
}

int BReadsCompareAtIndex(BReads *pOne, int iOne, BReads *pTwo, int iTwo) 
{
	int cmp;

	cmp = BStringCompare(&pOne->reads[iOne], &pTwo->reads[iTwo])
	if(cmp < 0 ||
			(cmp == 0 && pOne->offset[iOne] < pTwo->offset[iTwo]) || 
			(cmp == 0 && pOne->offset[iOne] < pTwo->offset[iTwo] && pOne->strand[iOne] < pTwo->strand[iTwo])) {
		return -1;
	}
	else if(cmp == 0 && pOne->offset[iOne] == pTwo->offset[iTwo] && pOne->strand[iOne] == pTwo->strand[iTwo]) {
		return 0;
	}
	else {
		return 1;
	}
}

void BReadsCopyAtIndex(BReads *src, int srcIndex, BReads *dest, int destIndex)
{
	if(dest != src || srcIndex != destIndex) {
		BStrincCopy(&dest->reads[destIndex], &src->reads[srcIndex]);
		dest->offset[destIndex] = src->offset[srcIndex];
		dest->strand[destIndex] = src->strand[srcIndex];
	}
}

void BReadsAllocate(BReads *reads, int numReads)
{
	assert(reads->numReads == 0);
	reads->numReads = numReads;
	reads->reads = malloc(sizeof(BString*)*reads->numReads);
	if(NULL == reads->reads) {
		PrintError("BReadsAllocate",
				"reads->reads",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	reads->offset = malloc(sizeof(int32_t)*(reads->numReads));
	if(NULL == reads->offset) {
		PrintError("BReadsAllocate",
				"reads->offset",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	reads->strand = malloc(sizeof(int8_t)*(reads->numReads));
	if(NULL == reads->strand) {
		PrintError("BReadsAllocate",
				"reads->strand",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
}

void BReadsReallocate(BReads *reads, int numReads) 
{
	int i;
	if(numReads > 0) {
		/* Remember to free the reads that will be reallocated if we go to less */
		if(numReads < reads->numReads) {
			for(i=numReads;i<reads->numReads;i++) {
				BStringFree(&reads->reads[i]);
			}
		}
		reads->numReads = numReads;
		reads->reads = realloc(reads->reads, sizeof(BString*)*(reads->numReads));
		if(NULL == reads->reads) {
			PrintError("BReadsReallocate",
					"reads->reads",
					"Could not reallocate memory",
					Exit,
					MallocMemory);
		}
		reads->offset = realloc(reads->offset, sizeof(int32_t)*(reads->numReads));
		if(NULL == reads->offset) {
			PrintError("BReadsReallocate",
					"reads->offset",
					"Could not reallocate memory",
					Exit,
					MallocMemory);
		}
		reads->strand = realloc(reads->strand, sizeof(int8_t)*(reads->numReads));
		if(NULL == reads->strand) {
			PrintError("BReadsReallocate",
					"reads->strand",
					"Could not reallocate memory",
					Exit,
					MallocMemory);
		}
	}
	else {
		BReadsFree(reads);
	}
}

void BReadsFree(BReads *reads) 
{
	int i;

	/* Free memory from reads */
	for(i=0;i<reads->numReads;i++) {
		BStringFree(&reads->reads[i]);
	}
	free(reads->reads);
	free(reads->read->length);
	free(reads->strand);
	free(reads->offset);
	BReadsInitialize(reads);
}

void BReadsInitialize(BReads *reads) 
{
	reads->reads=NULL;
	reads->read->length=NULL;
	reads->strand=NULL;
	reads->offset=NULL;
	reads->numReads=0;
}

void BReadsAppend(BReads *reads, 
		BString *read,
		int8_t direction,
		int32_t offset) 
{
	/* Copy over */
	BStringCopy(&reads->reads[reads->numReads-1], read);
	reads->read->length[reads->numReads-1] = read->length;
	reads->offset[reads->numReads-1] = offset;
	reads->strand[reads->numReads-1] = direction;
}

/* TODO */
/* Debugging procedure */
void BReadsPrint(BReads *reads, BIndex *index) 
{
	int i;
	for(i=0;i<reads->numReads;i++) {
		BIndexPrintReadMasked(index, reads->reads[i], 0, stderr);
		fprintf(stderr, "%s\t%d\t%c\t%d\n",
				reads->reads.string[i],
				reads->read.length[i],
				reads->strand[i],
				reads->offset[i]);
	}
}
