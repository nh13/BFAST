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
#include "RGIndex.h"
#include "RGReads.h"

/* TODO 
 * remove directionality in RGReads
 * remove directionality in generating reads.
 * */

char ALPHABET[ALPHABET_SIZE] = "acgt";

/* TODO */
void RGReadsFindMatches(RGIndex *index, 
		RGBinary *rg,
		RGMatch *match,
		char *read,
		int readLength,
		int *offsets,
		int numOffsets,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions,
		int maxMatches)
{
	int i;
	char reverseRead[SEQUENCE_LENGTH]="\0";
	RGReads reads;

	/* Initialize reads */ 
	RGReadsInitialize(&reads);

	/* Get the reverse compliment */
	GetReverseComplimentAnyCase(read, reverseRead, readLength);

	/* Generate reads */
	RGReadsGenerateReads(read,
			readLength,
			index,
			&reads,
			FORWARD,
			offsets,
			numOffsets,
			numMismatches,
			numInsertions,
			numDeletions,
			numGapInsertions,
			numGapDeletions);

	/* Generate reads */
	RGReadsGenerateReads(reverseRead,
			readLength,
			index,
			&reads,
			REVERSE,
			offsets,
			numOffsets,
			numMismatches,
			numInsertions,
			numDeletions,
			numGapInsertions,
			numGapDeletions);

	/* Merge all reads */
	/*
	   if(numMismatches > 0 || 
	   numInsertions > 0 ||
	   numDeletions > 0 ||
	   numGapInsertions > 0 ||
	   numGapDeletions > 0) {
	   RGReadsRemoveDuplicates(reads);
	   }
	   */

	/* Get the matches */
	for(i=0;i<reads.numReads && match->maxReached == 0;i++) {
		RGIndexGetMatches(index, 
				rg,
				reads.reads[i],
				reads.readLength[i],
				reads.strand[i],
				reads.offset[i],
				match,
				maxMatches);
	}

	/* Remove duplicates from match */
	RGMatchRemoveDuplicates(match, maxMatches);

	/* Free memory */
	RGReadsFree(&reads);
}

/* TODO */
void RGReadsGenerateReads(char *read,
		int readLength,
		RGIndex *index,
		RGReads *reads,
		char direction,
		int *offsets,
		int numOffsets,
		int numMismatches,
		int numInsertions,
		int numDeletions,
		int numGapInsertions,
		int numGapDeletions)
{
	int i;

	/* Go through all offsets */
	for(i=0;i<numOffsets && (readLength - offsets[i]) >= index->totalLength;i++) {

		/* Insert the perfect match */
		RGReadsGeneratePerfectMatch(read,
				readLength,
				direction,
				offsets[i],
				index->numTiles,
				index->tileLengths,
				index->gaps,
				index->totalLength,
				numMismatches,
				reads);

		/* Go through all mismatches */
		/* Note: we allow any number (including zero) of mismatches up to
		 * numMismatches.  
		 * */
		if(numMismatches > 0) {
			RGReadsGenerateMismatches(read,
					readLength,
					direction,
					offsets[i],
					index->numTiles,
					index->tileLengths,
					index->gaps,
					index->totalLength,
					numMismatches,
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
			RGReadsGenerateInsertions(read,
					readLength,
					direction,
					offsets[i],
					index->numTiles,
					index->tileLengths,
					index->gaps,
					index->totalLength,
					numInsertions,
					reads);
		}

		/* GO through all deletions */
		/* Note: we allow only contiguous deletions of length up to 
		 * numDeletions.  We also always start from the offset.  We
		 * must add base to the reads, and therfore we enumerate
		 * over all possible deletions in the entire read.
		 * */
		if(numDeletions > 0) {
			RGReadsGenerateDeletions(read,
					readLength,
					direction,
					offsets[i],
					index->numTiles,
					index->tileLengths,
					index->gaps,
					index->totalLength,
					numDeletions,
					reads);
		}

		/* Go through all possible insertions in the gap between
		 * the pair of l-mers.  If there is a gap insertion,
		 * then we will delete bases in the gap.
		 * */
		if(numGapInsertions > 0) {
			RGReadsGenerateGapInsertions(read,
					readLength,
					direction,
					offsets[i],
					index->numTiles,
					index->tileLengths,
					index->gaps,
					index->totalLength,
					numGapInsertions,
					reads);
		}

		/* Go through all possible deletions in the gap between
		 * the pair of l-mers.  If there is a gap deletion, 
		 * then we will add bases to the gap.
		 * */
		if(numGapDeletions > 0) {
			/* Forward */
			RGReadsGenerateGapDeletions(read,
					readLength,
					direction,
					offsets[i],
					index->numTiles,
					index->tileLengths,
					index->gaps,
					index->totalLength,
					numGapDeletions,
					reads);
		}
	}

}

void RGReadsGeneratePerfectMatch(char *read,
		int readLength,
		char direction,
		int offset,
		int32_t numTiles,
		int32_t *tileLengths,
		int32_t *gaps,
		int64_t totalLength,
		int numMismatches,
		RGReads *reads)
{
	int32_t i;
	char *curRead;

	/* Check bounds */
	if(readLength < totalLength+offset) {
		return;
	}

	curRead = malloc(sizeof(char)*(totalLength+1));
	if(NULL == curRead) {
		PrintError("RGReadsGenerateMismatchesHelper",
				"curRead",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Copy over */
	for(i=offset;i<totalLength+offset;i++) {
		curRead[i-offset] = read[i];
	}
	curRead[totalLength]='\0';

	/* Append */
	RGReadsAppend(reads, curRead, totalLength, direction, offset);

	free(curRead);
}

/* TODO */
void RGReadsGenerateMismatches(char *read,
		int readLength,
		char direction,
		int offset,
		int32_t numTiles,
		int32_t *tileLengths,
		int32_t *gaps,
		int64_t totalLength,
		int numMismatches,
		RGReads *reads)
{
	char *curRead=NULL;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Generating indexes with %d mismatches\n",
				numMismatches);
	}

	/* Check bounds */
	if(readLength < totalLength+offset) {
		return;
	}

	assert(readLength >= totalLength+offset);
	/* Allocate memory */
	curRead = malloc(sizeof(char)*(totalLength+1));
	if(NULL == curRead) {
		PrintError("RGReadsGenerateMismatches",
				"curRead",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	RGReadsGenerateMismatchesHelper(read,
			readLength,
			direction,
			offset,
			numTiles,
			tileLengths,
			gaps,
			totalLength,
			numMismatches,
			0,
			reads,
			curRead,
			0,
			0);

	/* Free memory */
	free(curRead);
}

/* TODO */
void RGReadsGenerateMismatchesHelper(char *read,
		int readLength,
		char direction,
		int offset,
		int32_t numTiles,
		int32_t *tileLengths,
		int32_t *gaps,
		int64_t totalLength,
		int numMismatchesLeft,
		int readIndex,
		RGReads *reads,
		char *curRead,
		int curTile,
		int curTileIndex)
{
	int i, j;

	if(readIndex > totalLength) {
		return;
	}

	if(numMismatchesLeft > 0) {
		/* No more to print */
		assert(readIndex <= totalLength);
		if(readIndex == totalLength) {
			curRead[totalLength]='\0';
			RGReadsAppend(reads, curRead, totalLength, direction, offset);
			return;
		}
		else {
			assert(readIndex < totalLength);
			/* use mismatches */
			for(i=0;i<ALPHABET_SIZE;i++) {
				int tempReadIndex = readIndex;
				int tempCurTile = curTile;
				int tempCurTileIndex = curTileIndex;
				/* Check if we are at the end of the tile */
				if(curTileIndex == tileLengths[curTile]) {
					/* Move to the next tile */
					for(j=0;j<gaps[tempCurTile];j++) {
						curRead[tempReadIndex] = read[offset+tempReadIndex];
						tempReadIndex++;
						if(tempReadIndex > readLength) {
							return;
						}
					}
					tempCurTile++;
					tempCurTileIndex=0;
				}
				else {
					/* Increment position in the tile */
					tempCurTileIndex++;
				}
				if(tempReadIndex > readLength) {
					return;
				}
				assert(tempReadIndex <= readLength);
				curRead[tempReadIndex] = ALPHABET[i];
				if(read[offset+tempReadIndex] == ALPHABET[i]) {
					/* No mismatch */
					/* Keep going */
					RGReadsGenerateMismatchesHelper(read,
							readLength,
							direction,
							offset,
							numTiles,
							tileLengths,
							gaps,
							totalLength,
							numMismatchesLeft,
							tempReadIndex+1,
							reads,
							curRead,
							tempCurTile,
							tempCurTileIndex);
				}
				else {
					/* Mismatch */

					/* Keep going */
					RGReadsGenerateMismatchesHelper(read,
							readLength,
							direction,
							offset,
							numTiles,
							tileLengths,
							gaps,
							totalLength,
							numMismatchesLeft-1,
							tempReadIndex+1,
							reads,
							curRead,
							tempCurTile,
							tempCurTileIndex);
				}
			}
		}
	}
	else {
		/* print remaining */                                           
		while(readIndex < totalLength) {
			curRead[readIndex] = read[readIndex+offset];
			readIndex++;

		}
		assert(readIndex == totalLength);
		curRead[totalLength]='\0';
		/* Append */
		RGReadsAppend(reads, curRead, totalLength, direction, offset);
	}
}

/* TODO */
/* Note: Deletions have occured, so insert bases */
void RGReadsGenerateDeletions(char *read,
		int readLength,
		char direction,
		int offset,
		int32_t numTiles,
		int32_t *tileLengths,
		int32_t *gaps,
		int64_t totalLength,
		int numDeletions,
		RGReads *reads)
{
	char *curRead=NULL;

	if(VERBOSE>=DEBUG) {
		fprintf(stderr, "Generating reads with %d deletions.\n",
				numDeletions);
	}

	/* Allocate memory */
	curRead = malloc(sizeof(char)*(totalLength+1));
	if(NULL == curRead) {
		PrintError("RGReadsGenerateDeletions",
				"curRead",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	RGReadsGenerateDeletionsHelper(read,
			readLength,
			direction,
			offset,
			numTiles,
			tileLengths,
			gaps,
			totalLength,
			numDeletions,
			numDeletions,
			0,
			0,
			reads,
			curRead,
			0,
			0);

	/* Free memory */
	free(curRead);
}

/* TODO */
/* NOTE: no error checking yet! */
/* Deletion occured, so insert bases */
void RGReadsGenerateDeletionsHelper(char *read,
		int readLength,
		char direction,
		int offset,
		int32_t numTiles,
		int32_t *tileLengths,
		int32_t *gaps,
		int64_t totalLength,
		int numDeletionsLeft,
		int numDeletions,
		int deletionOffset,
		int readIndex,
		RGReads *reads,
		char *curRead,
		int curTile,
		int curTileIndex)
{
	int i, j;

	if(readIndex > totalLength) {
		return;
	}

	if(numDeletionsLeft > 0) {
		/* No more to print */
		if(readIndex == totalLength) {
			if(numDeletionsLeft != numDeletions) {
				curRead[totalLength]='\0';
				/* Append */
				RGReadsAppend(reads, curRead, totalLength, direction, offset);
			}
			return;
		}
		else {
			/* Update readIndex etc. based on current tile */
			int tempReadIndex = readIndex;
			int tempCurTile = curTile;
			int tempCurTileIndex = curTileIndex;
			/* Check if we are at the end of the tile */
			if(curTileIndex == tileLengths[curTile]) {
				/* Move to the next tile */
				for(j=0;j<gaps[tempCurTile];j++) {
					curRead[tempReadIndex] = read[offset+tempReadIndex-deletionOffset];
					tempReadIndex++;
				}
				tempCurTile++;
				tempCurTileIndex=0;
			}
			else {
				/* Increment position in the tile */
				tempCurTileIndex++;
			}
			/* try inserting a base - do not insert at the beginning or the end of a read */
			if(readIndex > 0 && readIndex < readLength-1) {
				for(i=0;i<ALPHABET_SIZE;i++) {
					curRead[readIndex] = ALPHABET[i];
					/* Use on first read */
					RGReadsGenerateDeletionsHelper(read,
							readLength,
							direction,
							offset,
							numTiles,
							tileLengths,
							gaps,
							totalLength,
							numDeletionsLeft-1,
							numDeletions,
							deletionOffset+1,
							tempReadIndex+1,
							reads,
							curRead,
							tempCurTile,
							tempCurTileIndex);
				}
			}
			/* This will enforce that insertions occur together */
			if(numDeletionsLeft == numDeletions) {
				/* Try not inserting a base */
				curRead[readIndex] = read[offset+readIndex-deletionOffset];
				RGReadsGenerateDeletionsHelper(read,
						readLength,
						direction,
						offset,
						numTiles,
						tileLengths,
						gaps,
						totalLength,
						numDeletionsLeft,
						numDeletions,
						deletionOffset,
						tempReadIndex+1,
						reads,
						curRead,
						tempCurTile,
						tempCurTileIndex);
			}
		}
	}
	else {
		/* print remaining */                                           
		while(readIndex < totalLength) {
			curRead[readIndex] = read[readIndex+offset-deletionOffset];
			readIndex++;
		}
		curRead[totalLength]='\0';
		assert(readIndex == totalLength);
		/* Append */
		RGReadsAppend(reads, curRead, totalLength, direction, offset);
		return;

	}
}

/* TODO */
/* Note: Insertions have occured, so delete bases */
void RGReadsGenerateInsertions(char *read,
		int readLength,
		char direction,
		int offset,
		int32_t numTiles,
		int32_t *tileLengths,
		int32_t *gaps,
		int64_t totalLength,
		int numInsertions,
		RGReads *reads)
{
	char *curRead=NULL;
	int maxNumInsertions = 0;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Generating reads with %d insertions\n",
				numInsertions);
	}

	/* Get the total number of insertions (delete bases) possible */
	maxNumInsertions = readLength - totalLength;

	if(maxNumInsertions <= 0) {
		/* Cannot delete any bases */
		return;
	}
	else if(maxNumInsertions < numInsertions) {
		/* Can only delete a certain # of bases */
		numInsertions = maxNumInsertions;
	}

	/* Allocate memory */
	curRead = malloc(sizeof(char)*(totalLength+1));
	if(NULL == curRead) {
		PrintError("RGReadsGenerateInsertions",
				"curRead",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	RGReadsGenerateInsertionsHelper(read,
			readLength,
			direction,
			offset,
			numTiles,
			tileLengths,
			gaps,
			totalLength,
			numInsertions,
			numInsertions,
			0,
			0,
			reads,
			curRead,
			0,
			0);

	/* Free memory */
	free(curRead);
}

/* TODO */
/* NOTE: no error checking yet! */
void RGReadsGenerateInsertionsHelper(char *read,
		int readLength,
		char direction,
		int offset,
		int32_t numTiles,
		int32_t *tileLengths,
		int32_t *gaps,
		int64_t totalLength,
		int numInsertionsLeft,
		int numInsertions,
		int insertionOffset,
		int readIndex,
		RGReads *reads,
		char *curRead,
		int curTile,
		int curTileIndex)
{
	int i, j;

	if(readIndex > totalLength) {
		return;
	}

	if(numInsertionsLeft > 0) {
		/* No more to print */
		if(readIndex >= totalLength) {
			if(numInsertionsLeft != numInsertions) {
				curRead[totalLength]='\0';
				/* Append */
				RGReadsAppend(reads, curRead, totalLength, direction, offset);
			}
			return;
		}
		/* try deleting a base */
		/* Don't delete if the previous base is the same as the one we
		 * are proposing to delete since we have already tried those permutations */
		if(readIndex == 0 || read[readIndex-1] != read[readIndex]) {
			RGReadsGenerateInsertionsHelper(read,
					readLength,
					direction,
					offset,
					numTiles,
					tileLengths,
					gaps,
					totalLength,
					numInsertionsLeft-1,
					numInsertions,
					insertionOffset+1,
					readIndex,
					reads,
					curRead,
					curTile,
					curTileIndex);
		}
		/* Try not deleting a base */
		/* Only do this if we haven't started deleting */ 
		if(numInsertionsLeft == numInsertions) {
			/* Update readIndex based on current tile etc */
			int tempReadIndex = readIndex;
			int tempCurTile = curTile;
			int tempCurTileIndex = curTileIndex;
			/* Check if we are at the end of the tile */
			if(curTileIndex == tileLengths[curTile]) {
				/* Move to the next tile */
				/* Save sequence in the gap to test for redundancy when deleting */
				for(i=readIndex;i<readIndex+gaps[tempCurTile];i++) {
					curRead[tempReadIndex] = read[offset+tempReadIndex+insertionOffset];
					tempReadIndex++;
				}
				for(j=0;j<gaps[tempCurTile];j++) {
					curRead[tempReadIndex] = read[offset+tempReadIndex+insertionOffset];
					tempReadIndex++;
				}
				tempCurTile++;
				tempCurTileIndex=0;
			}
			else {
				/* Increment position in the tile */
				tempCurTileIndex++;
			}
			curRead[tempReadIndex] = read[offset+tempReadIndex+insertionOffset];
			RGReadsGenerateInsertionsHelper(read,
					readLength,
					direction,
					offset,
					numTiles,
					tileLengths,
					gaps,
					totalLength,
					numInsertionsLeft,
					numInsertions,
					insertionOffset,
					tempReadIndex+1,
					reads,
					curRead,
					tempCurTile,
					tempCurTileIndex);
		}
	}
	else {
		/* print remaining */                                           
		while(readIndex<totalLength) {
			curRead[readIndex] = read[readIndex+offset+insertionOffset];
			readIndex++;
		}
		curRead[totalLength]='\0';
		assert(readIndex == totalLength);
		/* Append */
		RGReadsAppend(reads, curRead, totalLength, direction, offset);
		return;
	}
}

/* TODO */
/* Note: Deletions have occured, so insert bases in the gaps */
void RGReadsGenerateGapDeletions(char *read,
		int readLength,
		char direction,
		int offset,
		int32_t numTiles,
		int32_t *tileLengths,
		int32_t *gaps,
		int64_t totalLength,
		int numGapDeletions,
		RGReads *reads)
{
	char *curRead = NULL;

	if(VERBOSE>=DEBUG) {
		fprintf(stderr, "Generating reads with %d gap deletions.\n",
				numGapDeletions);
	}

	/* Get the total number of insertions (delete bases) possible */
	if(numGapDeletions <= 0) {
		/* Out of bounds.  Don't add anything. */
		return;
	}
	else {
		/* Allocate memory */
		curRead = malloc(sizeof(char)*(readLength+1));
		if(NULL == curRead) {
			PrintError("RGReadsGenerateGapDeletions",
					"curRead",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		RGReadsGenerateGapDeletionsHelper(read,
				readLength,
				direction,
				offset,
				numTiles,
				tileLengths,
				gaps,
				totalLength,
				numGapDeletions,
				reads,
				curRead);

		/* Free memory */
		free(curRead);
	}
}

/* TODO */
/* NOTE: no error checking yet! */
/* We assume that all insertions in the gap are grouped together */
void RGReadsGenerateGapDeletionsHelper(char *read,
		int readLength,
		char direction,
		int offset,
		int32_t numTiles,
		int32_t *tileLengths,
		int32_t *gaps,
		int64_t totalLength,
		int numGapDeletions,
		RGReads *reads,
		char *curRead)
{
	int i, j, k;
	int readPos;
	int curReadPos;

	/* Choose a gap to insert bases */
	readPos = offset;
	curReadPos = 0;
	for(i=0;i<numTiles-1 && readPos < readLength;i++) {
		/* Copy over tile before gap */
		for(j=0;j<tileLengths[i] && readPos < readLength;j++) {
			curRead[curReadPos] = read[readPos];
			curReadPos++;
			readPos++;
		}

		/* Insert the bases at the beginning of the gap */
		/* Only insert min(numGapDeletions, gaps[i]) bases into the gap, since we 
		 * wish to use the bases in the gap */
		for(j=1;j<=numGapDeletions && 
				j <= gaps[i] && 
				j <= (readLength - totalLength - readPos);j++) { /* The number of bases to insert */
			int tempCurReadPos = curReadPos;
			int tempReadPos = readPos;

			/* Insert bases into the gap */
			for(k=0;k<j;k++) {
				curRead[tempCurReadPos] = NULL_LETTER;
				tempCurReadPos++;
			}
			/* Copy over the bases after the current tile */
			while(tempCurReadPos < totalLength) {
				curRead[tempCurReadPos] = read[tempReadPos];
				tempCurReadPos++;
				tempReadPos++;
			}
			assert(tempCurReadPos == totalLength);
			curRead[totalLength]='\0';
			/* Append */
			RGReadsAppend(reads, curRead, totalLength, direction, offset);
		}

		/* Add the gap to our current position */
		for(j=0;j<gaps[i] && readPos < readLength;j++) {
			curRead[curReadPos] = read[readPos];
			curReadPos++;
			readPos++;
		}
	}
}

/* TODO */
/* Note: Insertions have occured, so delete bases */
void RGReadsGenerateGapInsertions(char *read,
		int readLength,
		char direction,
		int offset,
		int32_t numTiles,
		int32_t *tileLengths,
		int32_t *gaps,
		int64_t totalLength,
		int numGapInsertions,
		RGReads *reads)
{
	char *curRead=NULL;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Generating reads with %d gap insertions\n",
				numGapInsertions);
	}

	/* Get the total number of insertions (delete bases) possible */
	if(numGapInsertions <= 0) { 
		/* Out of bounds.  Don't add anything. */
		return;
	}
	else {

		/* Allocate memory */
		curRead = malloc(sizeof(char)*(readLength+1));
		if(NULL == curRead) {
			PrintError("RGReadsGenerateGapInsertions",
					"curRead",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}

		RGReadsGenerateGapInsertionsHelper(read,
				readLength,
				direction,
				offset,
				numTiles,
				tileLengths,
				gaps,
				totalLength,
				numGapInsertions,
				reads,
				curRead);

		/* Free memory */
		free(curRead);
	}
}

/* TODO */
/* NOTE: no error checking yet! */
/* Delete bases in the gaps */
void RGReadsGenerateGapInsertionsHelper(char *read,
		int readLength,
		char direction,
		int offset,
		int32_t numTiles,
		int32_t *tileLengths,
		int32_t *gaps,
		int64_t totalLength,
		int numGapInsertions,
		RGReads *reads,
		char *curRead)
{
	int i, j, k;
	int readPos;
	int curReadPos;

	/* Choose the number of bases to remove */
	for(i=1;i<=numGapInsertions &&
			i <= readLength - totalLength;i++) {
		/* Find a start position that will handle the number of bases we will remove */

		/* Initialize the start position to be at the offset */
		int startPos = offset;
		/* We first try to use the bases at the end of the read.  So get the number of
		 * bases we would need to shift the start position (the end can't handle all of 
		 * them in this case).
		 * */
		int frontRemove = (totalLength + i) - (readLength - offset);
		if(frontRemove > 0) {
			startPos -= frontRemove;
		}
		/* Only generate a read if we have enough bases */
		if(startPos >= 0) {
			assert( (readLength - startPos) >= totalLength);
			/* Choose a gap to remove bases */
			readPos = startPos;
			curReadPos = 0;
			for(j=0;j<numTiles-1;j++) {
				/* Get the bases before the gap in which we will insert bases */ 
				for(k=0;k<tileLengths[j];k++) {
					curRead[curReadPos] = read[readPos];
					curReadPos++;
					readPos++;
				}

				/* Delete the bases in the gap */
				int tempCurReadPos = curReadPos;
				int tempReadPos = readPos;
				/* Skip over j bases */
				tempReadPos += i;
				while(tempCurReadPos < totalLength) {
					curRead[tempCurReadPos] = read[tempReadPos];
					tempCurReadPos++;
					tempReadPos++;
				}
				curRead[totalLength]='\0';
				/* Append - not the different offset, which corresponds to the
				 * (possibly) new start position. */
				RGReadsAppend(reads, curRead, totalLength, direction, startPos);

				/* Add the bases in the gap */
				for(k=0;k<gaps[j];k++) {
					curRead[curReadPos] = read[readPos];
					curReadPos++;
					readPos++;
				}
			}
		}
	}
}

/* TODO */
void RGReadsRemoveDuplicates(RGReads *s)
{
	int32_t i;
	int32_t prevIndex=0;

	if(s->numReads <= 0) {
		return;
	}

	/* Sort the data structure */
	RGReadsQuickSort(s, 0, s->numReads-1);

	/* Remove duplicates */
	prevIndex=0;
	for(i=1;i<s->numReads;i++) {
		if(RGReadsCompareAtIndex(s, prevIndex, s, i)==0) { 
			/* Ignore */
		}
		else {
			prevIndex++;
			/* Copy over to temporary pair */
			RGReadsCopyAtIndex(s, i, s, prevIndex);
		}
	}

	/* Reallocate pair */
	RGReadsReallocate(s, prevIndex+1);

}

/* TO DO */
void RGReadsQuickSort(RGReads *s, int low, int high)
{
	int32_t i;
	int32_t pivot=-1;
	RGReads *temp;

	if(low < high) {
		/* Allocate memory for the temp RGReads indexes */
		temp = malloc(sizeof(RGReads));
		if(NULL == temp) {
			PrintError("RGReadsQuickSort",
					"temp",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		RGReadsInitialize(temp);
		RGReadsAllocate(temp, 1);
		temp->reads[0] = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == temp->reads[0]) {
			PrintError("RGReadsQuickSort",
					"temp->reads[0]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		temp->reads[0][0]='\0';
		assert(temp->numReads == 1);

		pivot = (low + high)/2;

		RGReadsCopyAtIndex(s, pivot, temp, 0);
		RGReadsCopyAtIndex(s, high, s, pivot);
		RGReadsCopyAtIndex(temp, 0, s, high);

		pivot = low;

		for(i=low;i<high;i++) {
			if(RGReadsCompareAtIndex(s, i, s, high) <= 0) {
				if(i!=pivot) {
					RGReadsCopyAtIndex(s, i, temp, 0);
					RGReadsCopyAtIndex(s, pivot, s, i);
					RGReadsCopyAtIndex(temp, 0, s, pivot);
				}
				pivot++;
			}
		}
		RGReadsCopyAtIndex(s, pivot, temp, 0);
		RGReadsCopyAtIndex(s, high, s, pivot);
		RGReadsCopyAtIndex(temp, 0, s, high);

		/* Free memory before recursive call */
		assert(temp->numReads == 1);
		RGReadsFree(temp);
		free(temp);
		temp=NULL;

		RGReadsQuickSort(s, low, pivot-1);
		RGReadsQuickSort(s, pivot+1, high);
	}
}

int RGReadsCompareAtIndex(RGReads *pOne, int iOne, RGReads *pTwo, int iTwo) 
{
	int cmp;

	cmp = strcmp(pOne->reads[iOne], pTwo->reads[iTwo]);
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

void RGReadsCopyAtIndex(RGReads *src, int srcIndex, RGReads *dest, int destIndex)
{
	if(dest != src || srcIndex != destIndex) {
		strcpy(dest->reads[destIndex], src->reads[srcIndex]);
		dest->readLength[destIndex] = src->readLength[srcIndex];
		dest->offset[destIndex] = src->offset[srcIndex];
		dest->strand[destIndex] = src->strand[srcIndex];
	}
}

void RGReadsAllocate(RGReads *reads, int numReads)
{
	assert(reads->numReads == 0);
	reads->numReads = numReads;
	reads->reads = malloc(sizeof(char*)*reads->numReads);
	if(NULL == reads->reads) {
		PrintError("RGReadsAllocate",
				"reads->reads",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	reads->readLength= malloc(sizeof(int32_t)*(reads->numReads));
	if(NULL == reads->readLength) {
		PrintError("RGReadsAllocate",
				"reads->readLength",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	reads->offset = malloc(sizeof(int32_t)*(reads->numReads));
	if(NULL == reads->offset) {
		PrintError("RGReadsAllocate",
				"reads->offset",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	reads->strand = malloc(sizeof(int8_t)*(reads->numReads));
	if(NULL == reads->strand) {
		PrintError("RGReadsAllocate",
				"reads->strand",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
}

void RGReadsReallocate(RGReads *reads, int numReads) 
{
	int i;
	if(numReads > 0) {
		/* Remember to free the reads that will be reallocated if we go to less */
		if(numReads < reads->numReads) {
			for(i=numReads;i<reads->numReads;i++) {
				free(reads->reads[i]);
			}
		}
		reads->numReads = numReads;
		reads->reads = realloc(reads->reads, sizeof(char*)*(reads->numReads));
		if(NULL == reads->reads) {
			PrintError("RGReadsReallocate",
					"reads->reads",
					"Could not reallocate memory",
					Exit,
					MallocMemory);
		}
		reads->readLength = realloc(reads->readLength, sizeof(int32_t)*(reads->numReads));
		if(NULL == reads->readLength) {
			PrintError("RGReadsReallocate",
					"reads->readLength",
					"Could not reallocate memory",
					Exit,
					MallocMemory);
		}
		reads->offset = realloc(reads->offset, sizeof(int32_t)*(reads->numReads));
		if(NULL == reads->offset) {
			PrintError("RGReadsReallocate",
					"reads->offset",
					"Could not reallocate memory",
					Exit,
					MallocMemory);
		}
		reads->strand = realloc(reads->strand, sizeof(int8_t)*(reads->numReads));
		if(NULL == reads->strand) {
			PrintError("RGReadsReallocate",
					"reads->strand",
					"Could not reallocate memory",
					Exit,
					MallocMemory);
		}
	}
	else {
		RGReadsFree(reads);
	}
}

void RGReadsFree(RGReads *reads) 
{
	int i;

	/* Free memory from reads */
	for(i=0;i<reads->numReads;i++) {
		free(reads->reads[i]);
		reads->reads[i] = NULL;
	}
	free(reads->reads);
	free(reads->readLength);
	free(reads->strand);
	free(reads->offset);
	RGReadsInitialize(reads);
}

void RGReadsInitialize(RGReads *reads) 
{
	reads->reads=NULL;
	reads->readLength=NULL;
	reads->strand=NULL;
	reads->offset=NULL;
	reads->numReads=0;
}

void RGReadsAppend(RGReads *reads, 
		char *read,
		int32_t readLength,
		int8_t direction,
		int32_t offset) 
{
	char *FnName = "RGReadsAppend";

	/* Allocate memory */
	RGReadsReallocate(reads, reads->numReads+1);
	/* Allocate memory for read */
	reads->reads[reads->numReads-1] = malloc(sizeof(char)*(readLength+1));
	if(NULL == reads->reads[reads->numReads-1]) {
		PrintError(FnName,
				"reads->reads[reads->numReads-1]",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Copy over */
	strcpy(reads->reads[reads->numReads-1], read);
	reads->readLength[reads->numReads-1] = readLength;
	reads->offset[reads->numReads-1] = offset;
	reads->strand[reads->numReads-1] = direction;
}

/* TODO */
/* Debugging procedure */
void RGReadsPrint(RGReads *reads, RGIndex *index) 
{
	int i;
	for(i=0;i<reads->numReads;i++) {
		RGIndexPrintReadMasked(index, reads->reads[i], 0, stderr);
		fprintf(stderr, "%s\t%d\t%c\t%d\n",
				reads->reads[i],
				reads->readLength[i],
				reads->strand[i],
				reads->offset[i]);
	}
}
