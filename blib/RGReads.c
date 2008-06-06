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
	reads.numReads = 0;
	reads.reads = NULL;
	reads.strand = NULL;
	reads.offset = NULL;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "\nIn RGReadsFindMatchesInIndex.\n");
	}

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

	/* Get the matches */
	for(i=0;i<reads.numReads && match->maxReached == 0;i++) {
		RGIndexGetMatches(index, 
				rg,
				reads.reads[i],
				reads.strand[i],
				reads.offset[i],
				match,
				maxMatches);
	}

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Found %d matches in RGReadsFindMatchesInIndex.\n",
				match->numEntries);
	}

	/* Remove duplicates from match */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Removing duplicates\n");
	}
	RGMatchRemoveDuplicates(match, maxMatches);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Removed duplicates with %d matches remaining in RGReadsFindMatchesInIndex.\n",
				match->numEntries);
	}

	/* Free memory */
	RGReadsFree(&reads);

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Exiting RGReadsFindMatchesInIndex.\n");
	}
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

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Generating all possible reads (%d offsets).\n",
				numOffsets);
	}

	/* Go through all offsets */
	for(i=0;i<numOffsets;i++) {
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Generating reads with offset %d.\n",
					offsets[i]);
		}

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

	/* Go through all the reads and concatenate to total length */
	/*
	   for(i=0;i<reads->numReads;i++) {
	   reads->reads[i][index->totalLength]='\0';
	   }
	   */

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "All reads generated\n");
	}

	/* Merge all reads */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Calling RGReadsRemoveDuplicates\n");
	}
	if(numMismatches > 0 || 
			numInsertions > 0 ||
			numDeletions > 0 ||
			numGapInsertions > 0 ||
			numGapDeletions > 0) {
		RGReadsRemoveDuplicates(reads);
	}
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Exited from RGReadsRemoveDuplicates\n");
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

	/* Check bounds */
	if(readLength < totalLength+offset) {
		return;
	}
	/* Update the number of reads */
	reads->numReads++;
	/* Allocate memory */
	RGReadsReallocate(reads, reads->numReads);
	reads->reads[reads->numReads-1] = malloc(sizeof(char)*(totalLength+1));
	if(NULL == reads->reads[reads->numReads-1]) {
		PrintError("RGReadsGenerateMismatchesHelper",
				"reads->reads[reads->numReads-1]",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Copy over */
	for(i=offset;i<totalLength+offset;i++) {
		reads->reads[reads->numReads-1][i-offset] = read[i];
	}
	reads->reads[reads->numReads-1][totalLength] = '\0';
	reads->offset[reads->numReads-1] = offset;
	reads->strand[reads->numReads-1] = direction;
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
			/* Update the number of reads */
			reads->numReads++;
			/* Allocate memory */
			RGReadsReallocate(reads, reads->numReads);
			reads->reads[reads->numReads-1] = malloc(sizeof(char)*(totalLength+1));
			if(NULL == reads->reads[reads->numReads-1]) {
				PrintError("RGReadsGenerateMismatchesHelper",
						"reads->reads[reads->numReads-1]",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
			/* Copy over */
			strcpy(reads->reads[reads->numReads-1], curRead);
			reads->offset[reads->numReads-1] = offset;
			reads->strand[reads->numReads-1] = direction;
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
		/* Update the number of reads */
		reads->numReads++;
		/* Allocate memory */
		RGReadsReallocate(reads, reads->numReads);
		reads->reads[reads->numReads-1] = malloc(sizeof(char)*(totalLength+1));
		if(NULL == reads->reads[reads->numReads-1]) {
			PrintError("RGReadsGenerateMismatchesHelper",
					"reads->reads[reads->numReads-1]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Copy over */
		strcpy(reads->reads[reads->numReads-1], curRead);
		reads->offset[reads->numReads-1] = offset;
		reads->strand[reads->numReads-1] = direction;
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
				reads->numReads++;
				/* Allocate memory */
				RGReadsReallocate(reads, reads->numReads);
				reads->reads[reads->numReads-1] = malloc(sizeof(char)*(totalLength+1));
				if(NULL == reads->reads[reads->numReads-1]) {
					PrintError("RGReadsGenerateDeletionsHelper",
							"reads->reads[reads->numReads-1]",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}
				/* Copy over */
				strcpy(reads->reads[reads->numReads-1], curRead);
				reads->offset[reads->numReads-1] = offset;
				reads->strand[reads->numReads-1] = direction;
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
		/* Update the number of reads */
		reads->numReads++;
		/* Allocate memory */
		RGReadsReallocate(reads, reads->numReads);
		reads->reads[reads->numReads-1] = malloc(sizeof(char)*(totalLength+1));
		if(NULL == reads->reads[reads->numReads-1]) {
			PrintError("RGReadsGenerateDeletionsHelper",
					"reads->reads[reads->numReads-1]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Copy over */
		strcpy(reads->reads[reads->numReads-1], curRead);
		reads->offset[reads->numReads-1] = offset;
		reads->strand[reads->numReads-1] = direction;
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
				reads->numReads++;
				/* Allocate memory */
				RGReadsReallocate(reads, reads->numReads);
				reads->reads[reads->numReads-1] = malloc(sizeof(char)*(totalLength+1));
				if(NULL == reads->reads[reads->numReads-1]) {
					PrintError("RGReadsGenerateInsertionsHelper",
							"reads->reads[reads->numReads-1]",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}
				/* Copy over */
				strcpy(reads->reads[reads->numReads-1], curRead);
				reads->offset[reads->numReads-1] = offset;
				reads->strand[reads->numReads-1] = direction;
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
		/* Update the number of reads */
		reads->numReads++;
		/* Allocate memory */
		RGReadsReallocate(reads, reads->numReads);
		reads->reads[reads->numReads-1] = malloc(sizeof(char)*(totalLength+1));
		if(NULL == reads->reads[reads->numReads-1]) {
			PrintError("RGReadsGenerateInsertionsHelper",
					"reads->reads[reads->numReads-1]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Copy over */
		strcpy(reads->reads[reads->numReads-1], curRead);
		reads->offset[reads->numReads-1] = offset;
		reads->strand[reads->numReads-1] = direction;
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
	int i;
	int maxNumGapDeletions=0;

	if(VERBOSE>=DEBUG) {
		fprintf(stderr, "Generating reads with %d gap deletions.\n",
				numGapDeletions);
	}

	/* Check bounds */
	/* If we are to insert n bases, then we can only insert into gaps
	 * that are at least as big.
	 * */
	for(i=0;i<numTiles-1;i++) {
		if(gaps[i] > maxNumGapDeletions) {
			maxNumGapDeletions = gaps[i];
		}
	}

	if(numGapDeletions > maxNumGapDeletions) {
		numGapDeletions = maxNumGapDeletions;
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
	int i, j, k, l, m;
	int curPos = 0;
	int readPos;
	int curReadPos;

	/* Choose a gap to insert bases */
	for(i=0;i<numTiles-1;i++) {
		/* Move to the current gap */
		curPos += tileLengths[i];
		/* Choose a starting position in the gap for inserting bases */
		for(j=0;j<gaps[i];j++) {
			/* For the total number of bases we will insert */
			/* If k > gaps[i], we will insert more than there 
			 * is found in the gap.  This will cause a combinatorial
			 * explosion since we can do a simplifying shift. */
			for(k=1;k<=numGapDeletions && k <= gaps[i];k++) {
				/* Copy over */

				curReadPos = 0;
				readPos = offset + k - 1;

				/* First copy over the bases up to the current insertion point */
				for(l=0;l<=i;l++) { /* For the given tile */
					for(m=0;m<tileLengths[l];m++) { /* For each base in the tile */
						curRead[curReadPos] = read[readPos];
						curReadPos++;
						readPos++;
					}
					if(l<i) { /* For each gap up to (not including) the current one */
						for(m=0;m<gaps[l];m++) {
							curRead[curReadPos] = read[readPos];
							curReadPos++;
							readPos++;
						}
					}
				}
				assert(curReadPos == curPos);
				for(m=0;m<j;m++) {
					curRead[curReadPos] = read[readPos];
					curReadPos++;
					readPos++;
				}
				/* Move up to insertion point */
				/* Insert bases into the gap - let's choose 'a' since it will not 
				 * be used anyway */
				for(m=0;m<k;m++) {
					curRead[curReadPos] = 'a';
					curReadPos++;
				}
				/* Copy over the bases after the current tile */
				for(m=curReadPos;m<totalLength;m++) {
					curRead[curReadPos] = read[readPos];
					curReadPos++;
					readPos++;
				}
				curRead[curReadPos]='\0';
				if(curReadPos >= totalLength) {
					/* Update the number of reads */
					reads->numReads++;
					/* Allocate memory */
					RGReadsReallocate(reads, reads->numReads);
					reads->reads[reads->numReads-1] = malloc(sizeof(char)*(readLength+1));
					if(NULL == reads->reads[reads->numReads-1]) {
						PrintError("RGReadsGenerateGapDeletionsHelper",
								"reads->reads[reads->numReads-1]",
								"Could not allocate memory",
								Exit,
								MallocMemory);
					}
					/* Copy over */
					strcpy(reads->reads[reads->numReads-1], curRead);
					reads->offset[reads->numReads-1] = offset;
					reads->strand[reads->numReads-1] = direction;
				}
			}
		}
		/* Add the gap to our current position */
		curPos += gaps[i-1];
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
	int i;
	char *curRead=NULL;
	int maxNumGapInsertions = 0;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Generating reads with %d gap insertions\n",
				numGapInsertions);
	}

	/* Check bounds */
	/* If we are to delete n bases, then we can only delete the 
	 * max gap and also we need at least n extra bases to fill in.
	 * */
	for(i=0;i<numTiles-1;i++) {
		if(gaps[i] > maxNumGapInsertions) {
			maxNumGapInsertions = gaps[i];
		}
	}

	if(maxNumGapInsertions < numGapInsertions) {
		numGapInsertions = maxNumGapInsertions;
	}

	if(totalLength + numGapInsertions > readLength) {
		numGapInsertions = totalLength - readLength;
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
	int i, j, l, m;
	int curReadPos=0;

	/* Choose which gap to remove bases */
	for(i=0;i<numTiles-1;i++) { /* Remove from the ith gap */
		/* Only delete the maximum amount in the gap */
		int tempNumGapInsertions = numGapInsertions;
		if(tempNumGapInsertions > gaps[i]) {
			tempNumGapInsertions = gaps[i];
		}
		for(j=1;j<=tempNumGapInsertions;j++) { /* Choose the number of bases to delete from the gap */
			curReadPos=0;
			/* Copy over bases before the insertion */
			for(l=0;l<=i;l++) { /* For each previous tile */
				for(m=0;m<tileLengths[l];m++) { /* For each base in the tile */
					curRead[curReadPos] = read[offset+curReadPos];
					curReadPos++;
				}
				for(m=0;m<gaps[l] && l!=i;m++) { /* For each base in the gaps */
					curRead[curReadPos] = read[offset+curReadPos];
					curReadPos++;
				}
			}
			/* For the first j bases in the gap */
			for(m=0;m<j-1;m++) {
				curRead[curReadPos] = read[offset+curReadPos];
				curReadPos++;
			}
			/* Copy over bases after the insertion */
			for(l=i+1;l<numTiles;l++) { /* For each tile */
				for(m=0;m<tileLengths[l];m++) { /* For each base in the tile */
					curRead[curReadPos] = read[offset+curReadPos+j];
					curReadPos++;
				}
				for(m=0;l!=numTiles-1 && m<gaps[l];m++) { /* For each base in the gaps */
					curRead[curReadPos] = read[offset+curReadPos+j];
					curReadPos++;
				}
			}
			curRead[curReadPos]='\0';
			if(curReadPos >= totalLength) {
				/* Update the number of reads */
				reads->numReads++;
				/* Allocate memory */
				RGReadsReallocate(reads, reads->numReads);
				reads->reads[reads->numReads-1] = malloc(sizeof(char)*(readLength+1));
				if(NULL == reads->reads[reads->numReads-1]) {
					PrintError("RGReadsGenerateGapDeletionsHelper",
							"reads->reads[reads->numReads-1]",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}
				/* Copy over */
				strcpy(reads->reads[reads->numReads-1], curRead);
				reads->offset[reads->numReads-1] = offset;
				reads->strand[reads->numReads-1] = direction;
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
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Sorting\n");
	}
	RGReadsQuickSort(s, 0, s->numReads-1);
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Sorted!\n");
	}

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
		RGReadsAllocate(temp, 1);
		temp->reads[0] = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == temp->reads[0]) {
			PrintError("RGReadsQuickSort",
					"temp->reads[0]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}

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
		for(i=0;i<temp->numReads;i++) {
			free(temp->reads[i]);
		}
		free(temp->reads);
		free(temp->offset);
		free(temp->strand);
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
	strcpy(dest->reads[destIndex], src->reads[srcIndex]);
	dest->offset[destIndex] = src->offset[srcIndex];
	dest->strand[destIndex] = src->strand[srcIndex];
}

void RGReadsAllocate(RGReads *reads, int numReads)
{
	reads->numReads = numReads;
	reads->reads = malloc(sizeof(char*)*reads->numReads);
	if(NULL == reads->reads) {
		PrintError("RGReadsAllocate",
				"reads->reads",
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
	reads->numReads = numReads;
	reads->reads = realloc(reads->reads, sizeof(char*)*(reads->numReads));
	if(NULL == reads->reads) {
		PrintError("RGReadsReallocate",
				"reads->reads",
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

void RGReadsFree(RGReads *reads) 
{
	int i;

	/* Free memory from reads */
	for(i=0;i<reads->numReads;i++) {
		free(reads->reads[i]);
		reads->reads[i] = NULL;
	}
	free(reads->reads);
	reads->reads=NULL;
	free(reads->strand);
	reads->strand=NULL;
	free(reads->offset);
	reads->offset=NULL;
}
