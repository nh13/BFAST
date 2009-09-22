#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "BLib.h"
#include "BLibDefinitions.h"
#include "BError.h"
#include "RGMatches.h"
#include "AlignedRead.h"
#include "AlignedEnd.h"
#include "AlignedEntry.h"
#include "ScoringMatrix.h"
#include "AlignNTSpace.h"
#include "AlignColorSpace.h"
#include "RGMatch.h"
#include "Align.h"

int AlignRGMatches(RGMatches *m,
		RGBinary *rg,
		AlignedRead *a,
		int32_t space,
		int32_t offset,
		ScoringMatrix *sm,
		int32_t ungapped,
		int32_t unconstrained,
		int32_t bestOnly,
		int32_t usePairedEndLength,
		int32_t pairedEndLength,
		int32_t mirroringType,
		int32_t forceMirroring,
		AlignMatrixNT ***matrixNT,
		AlignMatrixCS ***matrixCS,
		int32_t *maxReadLength,
		int32_t *maxReferenceLength
		)
{
	double bestScore;
	int32_t i;
	int32_t numLocalAlignments = 0;
	int32_t numAligned=0;

	/* Check to see if we should try to align one read with no candidate
	 * locations if the other one has candidate locations.
	 *
	 * This assumes that the the first read is 5'->3' before the second read.
	 * */
	if(2 == m->numEnds && usePairedEndLength == 1) {
		RGMatchesMirrorPairedEnd(m,
				rg,
				pairedEndLength,
				mirroringType,
				forceMirroring);
	}

	AlignedReadAllocate(a,
			m->readName,
			m->numEnds,
			space);

	/* Align each end individually */
	for(i=0;i<m->numEnds;i++) {
		/* Align an end */
		AlignRGMatchesOneEnd(&m->ends[i],
				rg,
				&a->ends[i],
				space,
				offset,
				sm,
				ungapped,
				unconstrained,
				bestOnly,
				&bestScore,
				&numAligned,
				matrixNT,
				matrixCS,
				maxReadLength,
				maxReferenceLength
				);
		if(BestOnly == bestOnly) {
			numLocalAlignments += AlignRGMatchesKeepBestScore(&a->ends[i],
					bestScore);
		}
		else {
			numLocalAlignments += numAligned;
		}
	}
	return numLocalAlignments;
}

/* TODO */
void AlignRGMatchesOneEnd(RGMatch *m,
		RGBinary *rg,
		AlignedEnd *end,
		int32_t space,
		int32_t offset,
		ScoringMatrix *sm,
		int32_t ungapped,
		int32_t unconstrained,
		int32_t bestOnly,
		double *bestScore,
		int32_t *numAligned,
		AlignMatrixNT ***matrixNT,
		AlignMatrixCS ***matrixCS,
		int32_t *maxReadLength,
		int32_t *maxReferenceLength)
{
	char *FnName="AlignRGMatchOneEnd";
	int32_t i;
	char **references=NULL;
	char **masks=NULL;
	int32_t *referenceLengths=NULL;
	int32_t *referencePositions=NULL;
	char read[SEQUENCE_LENGTH]="\0";
	int32_t readLength, updateMatrix=0;
	int32_t ctr=0;

	(*bestScore)=DBL_MIN;

	strcpy(read, m->read);

	if(NTSpace == space) {
		readLength = m->readLength;
	}
	else {
		readLength = ConvertReadFromColorSpace(read, m->readLength);
	}
	if((*maxReadLength) < readLength) {
		(*maxReadLength) = readLength;
		updateMatrix = 1;
	}

	/* Allocate */
	AlignedEndAllocate(end,
			m->read,
			m->qual,
			m->numEntries);

	/* Get all the references */
	references = malloc(sizeof(char*)*m->numEntries);
	if(NULL==references) {
		PrintError(FnName, "references", "Could not allocate memory", Exit, MallocMemory);
	}
	masks = malloc(sizeof(char*)*m->numEntries);
	if(NULL==masks) {
		PrintError(FnName, "masks", "Could not allocate memory", Exit, MallocMemory);
	}
	referenceLengths = malloc(sizeof(int32_t)*m->numEntries);
	if(NULL==referenceLengths) {
		PrintError(FnName, "referenceLengths", "Could not allocate memory", Exit, MallocMemory);
	}
	referencePositions = malloc(sizeof(int32_t)*m->numEntries);
	if(NULL==referencePositions) {
		PrintError(FnName, "referencePositions", "Could not allocate memory", Exit, MallocMemory);
	}
	for((*numAligned)=0,i=0,ctr=0;i<m->numEntries;i++) {
		references[ctr]=NULL; /* This is needed for RGBinaryGetReference */
		/* Get references */
		RGBinaryGetReference(rg,
				m->contigs[i],
				m->positions[i],
				m->strands[i], 
				offset,
				&references[ctr],
				readLength,
				&referenceLengths[ctr],
				&referencePositions[ctr]);
		assert(referenceLengths[ctr] > 0);
		/* Initialize entries */
		if(readLength <= referenceLengths[ctr]) {
			/* Copy over mask */
			masks[i] = RGMatchMaskToString(m->masks[i], m->readLength);
			/* Update contig name and strand */
			end->entries[ctr].contigNameLength = rg->contigs[m->contigs[i]-1].contigNameLength;
			end->entries[ctr].contigName = malloc(sizeof(char)*(end->entries[ctr].contigNameLength+1));
			if(NULL==end->entries[ctr].contigName) {
				PrintError(FnName, "end->entries[ctr].contigName", "Could not allocate memory", Exit, MallocMemory);
			}
			strcpy(end->entries[ctr].contigName, rg->contigs[m->contigs[i]-1].contigName);
			end->entries[ctr].contig = m->contigs[i];
			end->entries[ctr].strand = m->strands[i];
			/* The rest should be filled in later */
			end->entries[ctr].position = -1; 
			end->entries[ctr].score=DBL_MIN;
			end->entries[ctr].length = 0;
			end->entries[ctr].referenceLength = 0;
			end->entries[ctr].read = end->entries[ctr].reference = end->entries[ctr].colorError = NULL;

			if((*maxReferenceLength) < referenceLengths[ctr]) {
				(*maxReferenceLength) = referenceLengths[ctr];
				updateMatrix = 1;
			}

			(*numAligned)++;
			ctr++;
		}
		else {
			/* Free retrieved reference sequence */
			free(references[ctr]);
			references[ctr]=NULL;
			masks[ctr]=NULL;
		}
	}

	/* Reallocate entries if necessary */
	if(ctr < end->numEntries) {
		AlignedEndReallocate(end,
				ctr);
	}

	/* Reallocate matrix */
	if(1 == updateMatrix) {
		if(NTSpace == space) {
			(*matrixNT) = realloc((*matrixNT), sizeof(AlignMatrixNT*)*((*maxReadLength)+1));
			if(NULL==(*matrixNT)) {
				PrintError(FnName, "(*matrixNT)", "Could not allocate memory", Exit, MallocMemory);
			}
			for(i=0;i<(*maxReadLength)+1;i++) {
				(*matrixNT)[i] = realloc((*matrixNT)[i], sizeof(AlignMatrixNT)*((*maxReferenceLength)+1));
				if(NULL==(*matrixNT)[i]) {
					PrintError(FnName, "(*matrixNT)[i]", "Could not allocate memory", Exit, MallocMemory);
				}
			}
		}
		else {
			(*matrixCS) = realloc((*matrixCS), sizeof(AlignMatrixCS*)*((*maxReadLength)+1));
			if(NULL==(*matrixCS)) {
				PrintError(FnName, "(*matrixCS)", "Could not allocate memory", Exit, MallocMemory);
			}
			for(i=0;i<(*maxReadLength)+1;i++) {
				(*matrixCS)[i] = realloc((*matrixCS)[i], sizeof(AlignMatrixCS)*((*maxReferenceLength)+1));
				if(NULL==(*matrixCS)[i]) {
					PrintError(FnName, "(*matrix)[i]", "Could not allocate memory", Exit, MallocMemory);
				}
			}
		}
	}

#ifndef UNOPTIMIZED_SMITH_WATERMAN
	int32_t foundExact;
	foundExact = 0;
	/* Try exact alignment */
	for(i=0;i<end->numEntries;i++) {
		if(1==AlignExact(read, 
					readLength, 
					references[i], 
					referenceLengths[i],
					sm,
					&end->entries[i],
					referencePositions[i],
					end->entries[i].strand,
					space,
					offset)) {
			foundExact=1;
			if((*bestScore) < end->entries[i].score) {
				(*bestScore) = end->entries[i].score;
			}
		}
	}

	/* If we are to only output the best alignments and we have found an exact alignment, return */
	if(1==foundExact && bestOnly == BestOnly) {
		for(i=0;i<end->numEntries;i++) {
			free(references[i]);
			free(masks[i]);
		}
		free(references);
		free(masks);
		free(referenceLengths);
		free(referencePositions);
		return;
	}
#endif

#ifdef UNOPTIMIZED_SMITH_WATERMAN
	if(Ungapped == ungapped) {
#endif
		for(i=0;i<end->numEntries;i++) {
			if(!(DBL_MIN < end->entries[i].score)) {
				AlignUngapped(read,
						masks[i],
						readLength,
						references[i],
						referenceLengths[i],
						unconstrained,
						sm,
						&end->entries[i],
						space,
						offset,
						referencePositions[i],
						end->entries[i].strand);
				if((*bestScore) < end->entries[i].score) {
					(*bestScore) = end->entries[i].score;
				}
			}
		}

		/* Return if we are only to be searching for mismatches */
#ifndef UNOPTIMIZED_SMITH_WATERMAN
		if(Ungapped == ungapped) {
#endif
			for(i=0;i<end->numEntries;i++) {
				free(references[i]);
				free(masks[i]);
			}
			free(references);
			free(masks);
			free(referenceLengths);
			free(referencePositions);

			return;
			/* These compiler commands aren't necessary, but are here for vim tab indenting */
#ifdef UNOPTIMIZED_SMITH_WATERMAN
		}
#else 
	}
#endif

	/* Run Full */
	for(i=0;i<end->numEntries;i++) {
		AlignFullWithBound(read,
				masks[i],
				readLength,
				references[i],
				referenceLengths[i],
				sm,
				&end->entries[i],
				space,
				referencePositions[i],
				end->entries[i].strand,
				(BestOnly == bestOnly)?(*bestScore):end->entries[i].score,
				matrixNT,
				matrixCS);
		if((*bestScore) < end->entries[i].score) {
			(*bestScore) = end->entries[i].score;
		}
	}

	for(i=0;i<end->numEntries;i++) {
		free(references[i]);
		free(masks[i]);
	}
	free(references);
	free(masks);
	free(referenceLengths);
	free(referencePositions);
}

/* TODO */
int32_t AlignExact(char *read,
		int32_t readLength,
		char *reference,
		int32_t referenceLength,
		ScoringMatrix *sm,
		AlignedEntry *a,
		int32_t position,
		char strand,
		int32_t space,
		int32_t offset)
{
	char *FnName="AlignExact";
	int32_t i;
	int32_t foundOffset = -1;

	// Use this when we want to perform exact substring matching 
	//foundOffset = KnuthMorrisPratt(read, readLength, reference, referenceLength);

	// Use this when we want to perform exact string matching
	for(i=0, foundOffset=offset;i<readLength;i++) {
		if(ToUpper(read[i]) != ToUpper(reference[offset+i])) {
			foundOffset = -1;
			break;
		}
	}

	if(foundOffset < 0) {
		return -1;
	}
	else {
		char prevReadBase = COLOR_SPACE_START_NT;
		char referenceAligned[SEQUENCE_LENGTH]="\0";
		char colorErrorAligned[SEQUENCE_LENGTH]="\0";
		int32_t score=0;

		for(i=0;i<readLength;i++) {
			if(ColorSpace == space) {
				char curColor='X';
				if(0 == ConvertBaseToColorSpace(prevReadBase, read[i], &curColor)) {
					PrintError(FnName, "curColor", "Could not convert base to color space", Exit, OutOfRange);
				}
				/* Add score for color error, if any */
				score += ScoringMatrixGetColorScore(curColor,
						curColor,
						sm);

				colorErrorAligned[i] = GAP;
			}
			//assert(ToLower(read[i]) == ToLower(reference[i+foundOffset]));
			score += ScoringMatrixGetNTScore(read[i], read[i], sm);
			referenceAligned[i] = reference[i+foundOffset];
		}
		referenceAligned[readLength]='\0';
		colorErrorAligned[readLength]='\0';

		AlignedEntryUpdateAlignment(a,
				(FORWARD==strand) ? (position + foundOffset) : (position + referenceLength - readLength - foundOffset),
				score,
				readLength,
				readLength,
				read,
				referenceAligned,
				(NTSpace == space) ? NULL : colorErrorAligned);
	}

	return 1;
}

void AlignUngapped(char *read,
		char *mask,
		int32_t readLength,
		char *reference,
		int32_t referenceLength,
		int32_t unconstrained,
		ScoringMatrix *sm,
		AlignedEntry *a,
		int32_t space,
		int32_t offset,
		uint32_t position,
		char strand)
{
	// Note: we do not allow any wiggle room in ungapped alignment
	char *FnName="AlignUngapped";
	switch(space) {
		case NTSpace:
			AlignNTSpaceUngapped(read,
					mask,
					readLength,
					reference,
					referenceLength,
					unconstrained,
					sm,
					a,
					offset, //(Unconstrained == unconstrained) ? 0 : offset,
					position,
					strand);
			break;
		case ColorSpace:
			AlignColorSpaceUngapped(read,
					mask,
					readLength,
					reference,
					referenceLength,
					unconstrained,
					sm,
					a,
					offset, //(Unconstrained == unconstrained) ? 0 : offset,
					position,
					strand);
			break;
		default:
			PrintError(FnName, "space", "Could not understand space", Exit, OutOfRange);
			break;
	}
}

void AlignFullWithBound(char *read,
		char *mask,
		int32_t readLength,
		char *reference,
		int32_t referenceLength,
		ScoringMatrix *sm,
		AlignedEntry *a,
		int32_t space,
		int32_t position,
		char strand,
		double lowerBound,
		AlignMatrixNT ***matrixNT,
		AlignMatrixCS ***matrixCS) 
{
	char *FnName="AlignFullWithBound";
	int64_t maxH, maxV;

	maxV = maxH = 0;
	/* Get the maximum number of vertical and horizontal moves allowed */
	if(sm->gapOpenPenalty < sm->gapExtensionPenalty) {
		if(NTSpace == space) {
			/* b = nt match score */
			/* p = gap open */
			/* e = gap extend */
			/* N = read length */
			/* Find x such that b(N) + p + e(x - 1) < Bound */
			/* x < (e + Bound - p - b(N)) / e */
			maxH = GETMAX(0, (int32_t)ceil((lowerBound - sm->ntMatch*readLength  - sm->gapOpenPenalty + sm->gapExtensionPenalty) / sm->gapExtensionPenalty));
			/* Find x such that b(N - x) + p + e(x - 1) < lowerBound */
			/* b(N) - x(b) + p + e(x) - e < lowerBound */
			/* - x(b) + e(x) < lowerBound + e -p - b(N) */
			/* x < (lowerBound - e - p - b(N) ) / (e - b)*/
			maxV = GETMAX(0, ceil((lowerBound - sm->ntMatch*readLength  - sm->gapOpenPenalty + sm->gapExtensionPenalty) / (sm->gapExtensionPenalty - sm->ntMatch)));
		}
		else {
			/* c = color match score */
			/* b = nt match score */
			/* p = gap open */
			/* e = gap extend */
			/* N = read length */
			/* Find x such that (c + b)N + p + e(x - 1) < Bound */
			maxH = GETMAX(0, (int32_t)ceil((lowerBound - (sm->colorMatch + sm->ntMatch)*readLength  - sm->gapOpenPenalty + sm->gapExtensionPenalty) / sm->gapExtensionPenalty));
			/* Find x such that (c + b)(N - x) + p + e(x - 1) < lowerBound */
			maxV = GETMAX(0, ceil((lowerBound - (sm->colorMatch + sm->ntMatch)*readLength  - sm->gapOpenPenalty + sm->gapExtensionPenalty) / (sm->gapExtensionPenalty - sm->colorMatch - sm->ntMatch)));
		}
		assert(maxH >= 0 && maxV >= 0);
	}
	else {
		// Default to maximums
		maxH = GETMIN(maxH, readLength);
		maxV = GETMIN(maxV, readLength);
		/*
		   PrintError(FnName, PACKAGE_BUGREPORT, "This is currently not implemented, please report", Exit, OutOfRange);
		   */
	}
	if(maxH == 0 && maxV == 0) {
		/* Use result from searching only mismatches */
		return;
	}

	/* Do full alignment */

	/* Free relevant entries */
	free(a->read);
	free(a->reference);
	free(a->colorError);
	a->read = a->reference = a->colorError = NULL;

	/* Get the maximum number of vertical and horizontal moves */
	maxH = GETMIN(maxH, readLength);
	maxV = GETMIN(maxV, readLength);

	switch(space) {
		case NTSpace:
			AlignNTSpaceFullWithBound(read,
					readLength,
					reference,
					referenceLength,
					sm,
					a,
					maxH,
					maxV,
					(*matrixNT),
					position,
					strand);
			break;
		case ColorSpace:
			AlignColorSpaceFullWithBound(read,
					readLength,
					reference,
					referenceLength,
					sm,
					a,
					maxH,
					maxV,
					(*matrixCS),
					position,
					strand);
			break;
		default:
			PrintError(FnName, "space", "Could not understand space", Exit, OutOfRange);
			break;
	}
}

int32_t AlignRGMatchesKeepBestScore(AlignedEnd *end,
		double bestScore)
{
	char *FnName="AlignRGMatchesKeepBestScore";
	int32_t curIndex, i;
	int32_t numLocalAlignments = 0;

	for(curIndex=0, i=0;
			i<end->numEntries;
			i++) {
		if(DBL_MIN < end->entries[i].score) {
			numLocalAlignments++;
		}
		if(bestScore < end->entries[i].score) {
			PrintError(FnName, "bestScore", "Best score is incorrect", Exit, OutOfRange);
		}
		else if(!(end->entries[i].score < bestScore)) {
			/* Free */
			AlignedEntryFree(&end->entries[i]);
		}
		else {
			/* Copy over to cur index */
			AlignedEntryCopyAtIndex(end->entries, curIndex, end->entries, i);
			curIndex++;
		}
	}
	assert(curIndex > 0);

	end->numEntries = curIndex;
	end->entries = realloc(end->entries, sizeof(AlignedEntry*)*end->numEntries);
	if(NULL == end->entries) {
		PrintError(FnName, "end->entries", "Could not reallocate memory", Exit, MallocMemory);
	}

	return numLocalAlignments;
}
