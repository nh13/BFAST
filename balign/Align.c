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
#include "../blib/BLib.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/RGMatches.h"
#include "../blib/AlignEntries.h"
#include "../blib/AlignEntry.h"
#include "ScoringMatrix.h"
#include "AlignNTSpace.h"
#include "AlignColorSpace.h"
#include "Align.h"

int AlignRGMatches(RGMatches *m,
		RGBinary *rg,
		AlignEntries *a,
		int32_t space,
		int32_t pairedEnd,
		int32_t scoringType,
		int32_t offsetLength,
		ScoringMatrix *sm,
		int32_t alignmentType,
		int32_t bestOnly,
		int32_t usePairedEndLength,
		int32_t pairedEndLength,
		int32_t forceMirroring)
{
	double bestScore;
	int32_t numLocalAlignments = 0;

	/* Check to see if we should try to align one read with no candidate
	 * locations if the other one has candidate locations.
	 *
	 * This assumes that the the first read is 5'->3' before the second read.
	 * */
	if(PairedEnd == pairedEnd && usePairedEndLength == 1) {
		RGMatchesMirrorPairedEnd(m,
				pairedEndLength,
				forceMirroring);
	}

	AlignEntriesAllocate(a,
			(char*)m->readName,
			m->matchOne.numEntries,
			m->matchTwo.numEntries,
			pairedEnd,
			space);

	/* Align each individually */
	bestScore = AlignRGMatchesOneEnd(&m->matchOne,
			rg,
			a->entriesOne,
			space,
			scoringType,
			offsetLength,
			sm,
			alignmentType,
			bestOnly);
	if(BestOnly == bestOnly) {
		numLocalAlignments += AlignRGMatchesKeepBestScore(&a->entriesOne,
				&a->numEntriesOne,
				bestScore);
	}
	else {
		numLocalAlignments += a->numEntriesOne;
	}
	if(PairedEnd == pairedEnd) {
		bestScore = AlignRGMatchesOneEnd(&m->matchTwo,
				rg,
				a->entriesTwo,
				space,
				scoringType,
				offsetLength,
				sm,
				alignmentType,
				bestOnly);
		if(BestOnly == bestOnly) {
			numLocalAlignments += AlignRGMatchesKeepBestScore(&a->entriesTwo,
					&a->numEntriesTwo,
					bestScore);
		}
		else {
			numLocalAlignments += a->numEntriesTwo;
		}
	}
	return numLocalAlignments;
}

/* TODO */
double AlignRGMatchesOneEnd(RGMatch *m,
		RGBinary *rg,
		AlignEntry *entries,
		int32_t space,
		int32_t scoringType,
		int32_t offsetLength,
		ScoringMatrix *sm,
		int32_t alignmentType,
		int32_t bestOnly)
{
	char *FnName="AlignRGMatchOneEnd";
	int32_t i;
	char **references=NULL;
	int32_t *referenceLengths=NULL;
	int32_t *referencePositions=NULL;
	double bestScore=DBL_MIN;
	int32_t foundExact;
	char read[SEQUENCE_LENGTH]="\0";
	int32_t readLength;

	strcpy(read, (char*)m->read);
	if(NTSpace == space) {
		readLength = m->readLength;
	}
	else {
		readLength = ConvertReadFromColorSpace(read, m->readLength);
	}

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
	for(i=0;i<m->numEntries;i++) {
		references[i]=NULL; /* This is needed for RGBinaryGetReference */
		/* Get references */
		RGBinaryGetReference(rg,
				m->contigs[i],
				m->positions[i],
				m->strands[i], 
				offsetLength,
				&references[i],
				readLength,
				&referenceLengths[i],
				&referencePositions[i]);
		assert(referenceLengths[i] > 0);
		/* Initialize entries */
		entries[i].contigNameLength = rg->contigs[m->contigs[i]-1].contigNameLength;
		entries[i].contigName = malloc(sizeof(char)*(entries[i].contigNameLength+1));
		if(NULL==entries[i].contigName) {
			PrintError(FnName,
					"entries[i].contigName",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		strcpy(entries[i].contigName, rg->contigs[m->contigs[i]-1].contigName);
		entries[i].contig = m->contigs[i];
		entries[i].strand = m->strands[i];
		/* The rest should be filled in later */
		entries[i].position = -1; 
		entries[i].score=DBL_MIN;
		entries[i].length = 0;
		entries[i].referenceLength = 0;
		entries[i].read = entries[i].reference = entries[i].colorError = NULL;
	}

#ifndef UNOPTIMIZED_SMITH_WATERMAN
	foundExact = 0;
	/* Try exact alignment */
	for(i=0;i<m->numEntries;i++) {
		if(1==AlignExact(read, 
					readLength, 
					references[i], 
					referenceLengths[i],
					scoringType,
					sm,
					&entries[i],
					m->strands[i],
					referencePositions[i],
					space)) {
			foundExact=1;
			if(bestScore < entries[i].score) {
				bestScore = entries[i].score;
			}
		}
	}

	/* If we are to only output the best alignments and we have found an exact alignment, return */
	if(1==foundExact && bestOnly == BestOnly) {
		for(i=0;i<m->numEntries;i++) {
			free(references[i]);
		}
		free(references);
		free(referenceLengths);
		free(referencePositions);

		return bestScore;
	}
#endif

	for(i=0;i<m->numEntries;i++) {
		if(!(DBL_MIN < entries[i].score)) {
			AlignMismatchesOnly(read,
					readLength,
					references[i],
					referenceLengths[i],
					scoringType,
					sm,
					&entries[i],
					space,
					m->strands[i],
					referencePositions[i]);
			if(bestScore < entries[i].score) {
				bestScore = entries[i].score;
			}
		}
	}

	/* Return if we are only to be searching for mismatches */
	if(MismatchesOnly == alignmentType) {
		for(i=0;i<m->numEntries;i++) {
			free(references[i]);
		}
		free(references);
		free(referenceLengths);
		free(referencePositions);

		return bestScore;
	}

	/* Run Full */
	for(i=0;i<m->numEntries;i++) {
		AlignFullWithBound(read,
				readLength,
				references[i],
				referenceLengths[i],
				scoringType,
				sm,
				&entries[i],
				space,
				m->strands[i],
				referencePositions[i],
				(BestOnly == bestOnly)?bestScore:entries[i].score);
		if(bestScore < entries[i].score) {
			bestScore = entries[i].score;
		}
	}

	for(i=0;i<m->numEntries;i++) {
		free(references[i]);
	}
	free(references);
	free(referenceLengths);
	free(referencePositions);

	return bestScore;
}

/* TODO */
int32_t AlignExact(char *read,
		int32_t readLength,
		char *reference,
		int32_t referenceLength,
		int32_t scoringType,
		ScoringMatrix *sm,
		AlignEntry *a,
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
		a->colorError = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL==a->colorError) {
			PrintError(FnName,
					"a->colorError",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}

		for(i=0;i<readLength;i++) {
			if(ColorSpace == space && 
					ColorSpace == scoringType) {
				uint8_t curColor='X';
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

				a->colorError[i] = '0';
			}
			assert(ToLower(read[i]) == ToLower(reference[i+offset]));
			a->score += ScoringMatrixGetNTScore(read[i], read[i], sm);
			a->read[i] = a->reference[i] = read[i];
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
		int32_t scoringType,
		ScoringMatrix *sm,
		AlignEntry *a,
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
					scoringType,
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
		int32_t scoringType,
		ScoringMatrix *sm,
		AlignEntry *a,
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
		/* c = max color sub score */
		/* b = max nt sub score */
		/* p = gap open */
		/* e = gap extend */
		/* Find x such that (c + b)N + p + e(x - 1) < Bound */
		maxH = MAX(0, (int32_t)ceil((lowerBound - (sm->maxColorScore + sm->maxNTScore)*readLength  - sm->gapOpenPenalty + sm->gapExtensionPenalty) / sm->gapExtensionPenalty));
		/* Find x such that (c + b)(N - x) + p + e(x - 1) < lowerBound */
		maxV = MAX(0, ceil((lowerBound - (sm->maxColorScore + sm->maxNTScore)*readLength  - sm->gapOpenPenalty + sm->gapExtensionPenalty) / (sm->gapExtensionPenalty - sm->maxColorScore - sm->maxNTScore)));
		assert(maxH >= 0 && maxV >= 0);
	}
	else {
		PrintError(FnName,
				PACKAGE_BUGREPORT,
				"This is currently not implemented, please report",
				Exit,
				OutOfRange);
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
	maxH = MIN(maxH, readLength);
	maxV = MIN(maxV, readLength);

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
					scoringType,
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

int32_t AlignRGMatchesKeepBestScore(AlignEntry **entries,
		int32_t *numEntries,
		double bestScore)
{
	char *FnName="AlignRGMatchesKeepBestScore";
	int32_t curIndex, i;
	int32_t numLocalAlignments = 0;

	for(curIndex=0, i=0;
			i<(*numEntries);
			i++) {
		if(DBL_MIN < (*entries)[i].score) {
			numLocalAlignments++;
		}
		if(bestScore < (*entries)[i].score) {
			PrintError(FnName,
					"bestScore",
					"Best score is incorrect",
					Exit,
					OutOfRange);
		}
		else if(!((*entries)[i].score < bestScore)) {
			/* Free */
			AlignEntryFree(&(*entries)[i]);
		}
		else {
			/* Copy over */
			AlignEntryCopyAtIndex((*entries), i, (*entries), curIndex);
			curIndex++;
		}
	}
	assert(curIndex > 0);

	(*numEntries) = curIndex;
	(*entries) = realloc((*entries), sizeof(AlignEntry*)*(*numEntries));
	if(NULL == (*entries)) {
		PrintError(FnName,
				"(*entries)",
				"Could not reallocate memory",
				Exit,
				MallocMemory);
	}

	return numLocalAlignments;
}
