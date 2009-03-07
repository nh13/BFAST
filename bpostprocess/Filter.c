#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include "../blib/BError.h"
#include "../blib/AlignedEntry.h"
#include "../blib/AlignedEnd.h"
#include "../blib/AlignedRead.h"
#include "../blib/BLib.h"
#include "Definitions.h"
#include "Filter.h"

/* TODO */
int FilterAlignedRead(AlignedRead *a,
		int algorithmReads,
		int mappingQuality,
		int minMappingQuality,
		int minScoreReads,
		int startContig,
		int startPos,
		int endContig,
		int endPos,
		int maxMismatches,
		int maxColorErrors,
		int minDistancePaired,
		int maxDistancePaired,
		int useDistancePaired,
		int contigAbPaired,
		int inversionsPaired,
		int unpaired)
{
	char *FnName="FilterAlignedRead";
	int foundType;
	int32_t *foundTypes=NULL;
	AlignedRead tmpA;
	AlignedEnd tmpEnd;
	int32_t i, j;
	int32_t best, bestIndex, numBest;
	int32_t curContigDistance, curPositionDistance;
	double curUniquenessScore;

	AlignedReadInitialize(&tmpA);

	/* We should only modify "a" if it is going to be reported */ 
	/* Copy in case we do not find anything to report */
	AlignedReadCopy(&tmpA, a);

	if(NoFiltering != algorithmReads) {
		/* Filter each alignment individually */
		for(i=0;i<tmpA.numEnds;i++) {
			for(j=0;j<tmpA.ends[i].numEntries;j++) {
				/* Check if we should filter */
				if(0<FilterAlignedEntry(&tmpA.ends[i].entries[j],
							tmpA.space,
							minScoreReads,
							startContig,
							startPos,
							endContig,
							endPos,
							maxMismatches,
							maxColorErrors)) {
					/* Copy end here */
					if(j < tmpA.ends[i].numEntries-1) {
						AlignedEntryCopy(&tmpA.ends[i].entries[j], 
								&tmpA.ends[i].entries[tmpA.ends[i].numEntries-1]);
					}
					AlignedEndReallocate(&tmpA.ends[i], tmpA.ends[i].numEntries-1);
					/* Since we removed, do not incrememt */
					j--;
				}
			}
		}
	}

	foundType=NoneFound;
	foundTypes=malloc(sizeof(int32_t)*tmpA.numEnds);
	if(NULL == foundTypes) {
		PrintError(FnName,
				"foundTypes",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	for(i=0;i<tmpA.numEnds;i++) {
		foundTypes[i]=NoneFound;
	}

	/* Pick alignment for each end individually (is this a good idea?) */
	for(i=0;i<tmpA.numEnds;i++) {
		/* Choose each end */
		switch(algorithmReads) {
			case NoFiltering:
			case AllNotFiltered:
				foundTypes[i] = (0<tmpA.ends[i].numEntries)?Found:NoneFound;
				if(1==mappingQuality) {
					AlignedEndInitialize(&tmpEnd);
					AlignedEndCopy(&tmpEnd, &tmpA.ends[i]);
					for(j=0;j<tmpEnd.numEntries;j++) {
						tmpA.ends[i].entries[j].score=GetUniquenessScore(&tmpEnd,
								j);
					}
					AlignedEndFree(&tmpEnd);
				}
				break;
			case Unique:
				foundTypes[i]=(1==tmpA.ends[i].numEntries)?Found:NoneFound;
				if(1==mappingQuality) {
					if(Found == foundTypes[i]) {
						tmpA.ends[i].entries[0].score = 1;
					}
				}
				break;
			case BestMappingQuality:
				/* Replace score with uniqueness score if necessary */
				assert(0 == mappingQuality); 
				for(i=0;i<tmpA.numEnds;i++) {
					AlignedEndInitialize(&tmpEnd);
					AlignedEndCopy(&tmpEnd, &tmpA.ends[i]);
					for(j=0;j<tmpEnd.numEntries;j++) {
						tmpA.ends[i].entries[j].score=GetUniquenessScore(&tmpEnd,
								j);
					}
					AlignedEndFree(&tmpEnd);
				}
				/* Remove those that have to small of a uniqueness score */
				for(j=0;j<tmpA.ends[i].numEntries;j++) {
					if(tmA.ends[i].entries[j].score < minMappingQuality) {
						/* Copy end here */
						if(j < tmpA.ends[i].numEntries-1) {
							AlignedEntryCopy(&tmpA.ends[i].entries[j], 
									&tmpA.ends[i].entries[tmpA.ends[i].numEntries-1]);
						}
						AlignedEndReallocate(&tmpA.ends[i], tmpA.ends[i].numEntries-1);
						/* Since we removed, do not incrememt */
						j--;
					}
				}
				/* Continue to find the best score, since the score
				 * is now our mapping quality.  Note that we should have
				 * mappingQuality == 0 */
				assert(0 == mappingQuality);
				/* NO BREAK */
			case BestScore:
				best = INT_MIN;
				bestIndex = -1;
				numBest = 0;
				for(j=0;j<tmpA.ends[i].numEntries;j++) {
					if(best < tmpA.ends[i].entries[j].score) {
						best = tmpA.ends[i].entries[j].score;
						bestIndex = j;
						numBest = 1;
					}
					else if(best == tmpA.ends[i].entries[j].score) {
						numBest++;
					}
				}
				if(1==numBest) {
					foundTypes[i] = Found;
					if(1==mappingQuality) {
						/* Calculate uniqueness score */
						curUniquenessScore=GetUniquenessScore(&tmpA.ends[i],
								bestIndex);
						if(curUniquenessScore < minMappingQuality) {
							foundTypes[i] = NoneFound;
						}
						else {
							tmpA.ends[i].entries[bestIndex].score = curUniquenessScore;
						}
					}
					if(Found == foundTypes[i]) {
						/* Copy to front */
						AlignedEntryCopy(&tmpA.ends[i].entries[0], 
								&tmpA.ends[i].entries[bestIndex]);
						AlignedEndReallocate(&tmpA.ends[i], 1);
					}
				}
				break;
			default:
				PrintError(FnName,
						"algorithmReads",
						"Could not understand algorithmReads",
						Exit,
						OutOfRange);
				break;
		}
	}

	if(1 == tmpA.numEnds) {
		foundType=foundTypes[0];
	}
	else if(2 == tmpA.numEnds) {
		foundType=(Found==foundTypes[0] && Found==foundTypes[1])?Found:NoneFound;

		if(Found == foundType) {
			if(1 == useDistancePaired) {
				curContigDistance = a->ends[1].entries[0].contig - a->ends[0].entries[0].contig;
				curPositionDistance = a->ends[1].entries[0].position - a->ends[0].entries[0].position;

				if(0 != curContigDistance ||
						curPositionDistance < minDistancePaired ||
						maxDistancePaired < curPositionDistance) {
					if(a->ends[0].entries[0].strand == a->ends[1].entries[0].strand) {
						if(1 == contigAbPaired) {
							/* Same strand - contig abnormality */
							foundType = ContigAb;
						}
						else {
							foundType = NoneFound;
						}
					}
					else {
						if(1 == inversionsPaired) { 
							/* Different strand - inversion */
							foundType = Inversion;
						}
						else {
							foundType = NoneFound;
						}
					}
				}
			}
		}
	}
	else {
		/* Call found if all have been found */
		foundType=Found;
		for(i=0;Found==foundType && i<tmpA.numEnds;i++) {
			if(NoneFound == foundTypes[i]) {
				foundType=NoneFound;
			}
		}
	}
	/* See if we should check for unpaired alignments */
	if(1 == unpaired && NoneFound == foundType) {
		/* Call unpaired if at least one is found */
		foundType=NoneFound;
		for(i=0;NoneFound==foundType && i<tmpA.numEnds;i++) {
			if(Found == foundTypes[i]) {
				foundType=Found;
			}
		}
		if(Found == foundType) {
			foundType = Unpaired;
		}
	}

	/* If we found, then copy back */
	if(NoneFound != foundType) {
		AlignedReadFree(a);
		AlignedReadCopy(a, &tmpA);
	}
	AlignedReadFree(&tmpA);
	free(foundTypes);

	return foundType;
}

/* TODO */
int FilterAlignedEntry(AlignedEntry *a,
		int space,
		int minScoreReads,
		int startContig,
		int startPos,
		int endContig,
		int endPos,
		int maxMismatches,
		int maxColorErrors)
{
	int32_t numMismatches, numColorErrors;
	numMismatches=numColorErrors=0;
	/* Check if the alignment has at least the minimum score */
	if(a->score < minScoreReads) {
		return 1;
	}
	/* Check the genomic location */
	if((startContig != 0 && a->contig < startContig) || 
			(startContig != 0 && a->contig == startContig && startPos != 0 && a->position < startPos) ||
			(endContig != 0 && a->contig == endContig && endPos != 0 && a->position > endPos) ||
			(endContig != 0 && a->contig > endContig)) {
		return 2;
	}
	/* Check mismatches */
	numMismatches = GetNumMismatchesInAlignedEntry(a);
	/* Check color errors */
	numColorErrors = GetNumColorErrorsInAlignedEntry(a, space);
	if(maxMismatches < numMismatches) {
		return 3;
	}
	if(ColorSpace == space) {
		if(maxColorErrors < numColorErrors) {
			return 4;
		}
	}
	return 0;
}

/* TODO */
double GetMappingQuality(AlignedEnd *end,
		int32_t bestIndex) 
{
	int32_t nextBest, numNextBest;
	int32_t i;

	nextBest = INT_MIN;
	numNextBest = 0;

	for(i=0;i<end->numEntries;i++) {
		if(i != bestIndex) {
			if(nextBest < end->entries[i].score) {
				nextBest = end->entries[i].score;
				numNextBest = 1;
			}
			else if(nextBest == end->entries[i].score) {
				numNextBest++;
			}
		}
	}

	if(0 < numNextBest && 0 < nextBest) {
		return (end->entries[bestIndex].score
				/ (end->entries[bestIndex].score + numNextBest*nextBest));
	}
	else {
		return 1.0;
	}
}
