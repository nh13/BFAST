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
#include "Align.h"

int AlignRGMatches(RGMatches *m,
		RGBinary *rg,
		AlignedRead *a,
		int32_t space,
		int32_t offsetLength,
		ScoringMatrix *sm,
		int32_t alignmentType,
		int32_t bestOnly,
		int32_t usePairedEndLength,
		int32_t pairedEndLength,
		int32_t mirroringType,
		int32_t forceMirroring)
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
				offsetLength,
				sm,
				alignmentType,
				bestOnly,
				&bestScore,
				&numAligned
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
		int32_t offsetLength,
		ScoringMatrix *sm,
		int32_t alignmentType,
		int32_t bestOnly,
		double *bestScore,
		int32_t *numAligned)
{
	char *FnName="AlignRGMatchOneEnd";
	int32_t i;
	char **references=NULL;
	int32_t *referenceLengths=NULL;
	int32_t *referencePositions=NULL;
	char read[SEQUENCE_LENGTH]="\0";
	int32_t readLength;
	int32_t ctr=0;

	(*bestScore)=DBL_MIN;

	strcpy(read, m->read);

	if(NTSpace == space) {
		readLength = m->readLength;
	}
	else {
		readLength = ConvertReadFromColorSpace(read, m->readLength);
	}

	/* Allocate */
	AlignedEndAllocate(end,
			m->read,
			m->qual,
			m->numEntries);

	references = malloc(sizeof(char*)*m->numEntries);
	if(NULL==references) {
		PrintError(FnName,
				"references",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	referenceLengths = malloc(sizeof(int32_t)*m->numEntries);
	if(NULL==referenceLengths) {
		PrintError(FnName,
				"referenceLengths",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	referencePositions = malloc(sizeof(int32_t)*m->numEntries);
	if(NULL==referencePositions) {
		PrintError(FnName,
				"referencePositions",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	for((*numAligned)=0,i=0,ctr=0;i<m->numEntries;i++) {
		references[ctr]=NULL; /* This is needed for RGBinaryGetReference */
		/* Get references */
		RGBinaryGetReference(rg,
				m->contigs[i],
				m->positions[i],
				m->strands[i], 
				offsetLength,
				&references[ctr],
				readLength,
				&referenceLengths[ctr],
				&referencePositions[ctr]);
		assert(referenceLengths[ctr] > 0);
		/* Initialize entries */
		if(readLength <= referenceLengths[ctr]) {
			end->entries[ctr].contigNameLength = rg->contigs[m->contigs[i]-1].contigNameLength;
			end->entries[ctr].contigName = malloc(sizeof(char)*(end->entries[ctr].contigNameLength+1));
			if(NULL==end->entries[ctr].contigName) {
				PrintError(FnName,
						"end->entries[ctr].contigName",
						"Could not allocate memory",
						Exit,
						MallocMemory);
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
			(*numAligned)++;
			ctr++;
		}
		else {
			/* Free retrieved reference sequence */
			free(references[ctr]);
			references[ctr]=NULL;
		}
	}

	/* Reallocate entries if necessary */
	if(ctr < end->numEntries) {
		AlignedEndReallocate(end,
				ctr);
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
					end->entries[i].strand,
					referencePositions[i],
					space)) {
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
		}
		free(references);
		free(referenceLengths);
		free(referencePositions);
		return;
	}
#endif

#ifdef UNOPTIMIZED_SMITH_WATERMAN
	if(MismatchesOnly == alignmentType) {
#endif
		for(i=0;i<end->numEntries;i++) {
			if(!(DBL_MIN < end->entries[i].score)) {
				AlignMismatchesOnly(read,
						readLength,
						references[i],
						referenceLengths[i],
						sm,
						&end->entries[i],
						space,
						end->entries[i].strand,
						referencePositions[i]);
				if((*bestScore) < end->entries[i].score) {
					(*bestScore) = end->entries[i].score;
				}
			}
		}

		/* Return if we are only to be searching for mismatches */
#ifndef UNOPTIMIZED_SMITH_WATERMAN
		if(MismatchesOnly == alignmentType) {
#endif
			for(i=0;i<end->numEntries;i++) {
				free(references[i]);
			}
			free(references);
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
				readLength,
				references[i],
				referenceLengths[i],
				sm,
				&end->entries[i],
				space,
				end->entries[i].strand,
				referencePositions[i],
				(BestOnly == bestOnly)?(*bestScore):end->entries[i].score);
		if((*bestScore) < end->entries[i].score) {
			(*bestScore) = end->entries[i].score;
		}
	}

	for(i=0;i<end->numEntries;i++) {
		free(references[i]);
	}
	free(references);
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
		char strand,
		int32_t position,
		int32_t space) 
{
	char *FnName="AlignExact";
	int32_t i;
	int32_t offset;

	offset = KnuthMorrisPratt(read, readLength, reference, referenceLength);

	if(offset < 0) {
		return -1;
	}
	else {
		char prevReadBase = COLOR_SPACE_START_NT;
		a->strand = strand;
		a->position = (FORWARD==strand)?(position + offset):(position + referenceLength - readLength - offset);
		a->score = 0;
		a->length = readLength;
		a->referenceLength = readLength;

		/* Allocate memory */
		assert(NULL==a->read);
		a->read = malloc(sizeof(char)*(a->length+1));
		if(NULL==a->read) {
			PrintError(FnName,
					"a->read",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		assert(NULL==a->reference);
		a->reference = malloc(sizeof(char)*(a->length+1));
		if(NULL==a->reference) {
			PrintError(FnName,
					"a->reference",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		assert(NULL==a->colorError);
		if(ColorSpace == space) {
			a->colorError = malloc(sizeof(char)*SEQUENCE_LENGTH);
			if(NULL==a->colorError) {
				PrintError(FnName,
						"a->colorError",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}

		for(i=0;i<readLength;i++) {
			if(ColorSpace == space) {
				char curColor='X';
				if(0 == ConvertBaseToColorSpace(prevReadBase, read[i], &curColor)) {
					PrintError(FnName,
							"curColor",
							"Could not convert base to color space",
							Exit,
							OutOfRange);
				}
				/* Add score for color error, if any */
				a->score += ScoringMatrixGetColorScore(curColor,
						curColor,
						sm);

				a->colorError[i] = GAP;
			}
			assert(ToLower(read[i]) == ToLower(reference[i+offset]));
			a->score += ScoringMatrixGetNTScore(read[i], read[i], sm);
			a->read[i] = read[i];
			a->reference[i] = reference[i+offset];
		}
	}
	a->read[a->length] = '\0';
	a->reference[a->length] = '\0';
	if(ColorSpace == space) {
		a->colorError[a->length] = '\0';
	}

	return 1;
}

void AlignMismatchesOnly(char *read,
		int32_t readLength,
		char *reference,
		int32_t referenceLength,
		ScoringMatrix *sm,
		AlignedEntry *a,
		int32_t space,
		char strand,
		int32_t position)
{
	char *FnName="AlignMismatchesOnly";
	switch(space) {
		case NTSpace:
			AlignNTSpaceMismatchesOnly(read,
					readLength,
					reference,
					referenceLength,
					sm,
					a,
					strand,
					position);
			break;
		case ColorSpace:
			AlignColorSpaceMismatchesOnly(read,
					readLength,
					reference,
					referenceLength,
					sm,
					a,
					strand,
					position);
			break;
		default:
			PrintError(FnName,
					"space",
					"Could not understand space",
					Exit,
					OutOfRange);
			break;
	}
}

void AlignFullWithBound(char *read,
		int32_t readLength,
		char *reference,
		int32_t referenceLength,
		ScoringMatrix *sm,
		AlignedEntry *a,
		int32_t space,
		char strand,
		int32_t position,
		double lowerBound)
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
		   PrintError(FnName,
		   PACKAGE_BUGREPORT,
		   "This is currently not implemented, please report",
		   Exit,
		   OutOfRange);
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
					strand,
					position,
					maxH,
					maxV);
			break;
		case ColorSpace:
			AlignColorSpaceFullWithBound(read,
					readLength,
					reference,
					referenceLength,
					sm,
					a,
					strand,
					position,
					maxH,
					maxV);
			break;
		default:
			PrintError(FnName,
					"space",
					"Could not understand space",
					Exit,
					OutOfRange);
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
			PrintError(FnName,
					"bestScore",
					"Best score is incorrect",
					Exit,
					OutOfRange);
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
		PrintError(FnName,
				"end->entries",
				"Could not reallocate memory",
				Exit,
				MallocMemory);
	}

	return numLocalAlignments;
}
