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
#include "../blib/RGMatch.h"
#include "../blib/RGMatches.h"
#include "../blib/AlignedRead.h"
#include "../blib/AlignedEnd.h"
#include "../blib/AlignedEntry.h"
#include "../blib/ScoringMatrix.h"
#include "../balign/AlignNTSpace.h"
#include "../balign/AlignColorSpace.h"
#include "../balign/Align.h"
#include "Repair.h"

int Repair(AlignedRead *src,
		RGBinary *rg,
		AlignedRead *dest,
		int32_t space,
		ScoringMatrix *sm,
		int32_t alignmentType,
		int32_t minPairedEndDistance,
		int32_t maxPairedEndDistance,
		int32_t mirroringType,
		int32_t strandedness)
{
	char *FnName="Repair";
	double bestScore;
	int32_t i;
	int32_t numLocalAlignments = 0;
	int32_t numAligned=0;
	int32_t offsetLength;
	RGMatch m;

	if(2 != src->numEnds) {
		PrintError(FnName,
				"2 != src->numEnds",
				"Only unpaired reads may be processed",
				Exit,
				OutOfRange);
	}

	/* Copy over */
	AlignedReadCopy(dest, src);

	for(i=0;i<dest->numEnds;i++) {
		if(1 == dest->ends[i].numEntries) {
			RGMatchInitialize(&m);
			/* Create an artificial RGMatches */
			offsetLength = (int32_t)( (maxPairedEndDistance - minPairedEndDistance)/2.0 + 0.5);
			RepairMirrorFromAlignedEndToRGMatch(&m,
					i+1,
					dest->ends,
					space,
					offsetLength,
					mirroringType,
					strandedness);
			assert(0 <= offsetLength);
			/* Run */
			AlignRGMatchesOneEnd(&m,
					rg,
					&dest->ends[(i+1)%2],
					space,
					offsetLength,
					sm,
					alignmentType,
					AllAlignments,
					&bestScore,
					&numAligned
					);
			numLocalAlignments += numAligned;
			RGMatchFree(&m);
		}
	}
	return numLocalAlignments;
}

/* TODO */

void RepairMirrorFromAlignedEndToRGMatch(RGMatch *dest,
		int32_t whichEnd,
		AlignedEnd *ends,
		int32_t space,
		int32_t offsetLength,
		int32_t mirroringType,
		int32_t strandedness)
{
	char *FnName="RepairMirrorFromAlignedEndToRGMatch";
	int32_t i, ctr, numEntries;

	numEntries=0;
	switch(mirroringType) {
		case MirrorForward:
		case MirrorReverse:
			numEntries = 1;
			break;
		case MirrorBoth:
			numEntries = 2;
		default:
			break;
			PrintError(FnName,
					"mirroringType",
					"Could not understand mirroringType",
					Exit,
					OutOfRange);
	}
	switch(strandedness) {
		case StrandSame:
		case StrandOpposite:
			break;
		case StrandBoth:
			numEntries *= 2;
			break;
		default:
			PrintError(FnName,
					"strandedness",
					"Could not understand strandedness",
					Exit,
					OutOfRange);
	}
	/* Reallocate */
	RGMatchAllocate(dest, numEntries);
	/* Copy over read and qual */
	dest->readLength = ends[(2==whichEnd)?0:1].readLength;
	dest->read = malloc(sizeof(char)*(1+dest->readLength));
	if(NULL == dest->read) {
		PrintError(FnName,
				"dest->read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	strcpy(dest->read, ends[(2==whichEnd)?0:1].read);
	dest->qualLength = ends[(2==whichEnd)?0:1].qualLength;
	dest->qual = malloc(sizeof(char)*(1+dest->qualLength));
	if(NULL == dest->qual) {
		PrintError(FnName,
				"dest->qual",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	strcpy(dest->qual, ends[(2==whichEnd)?0:1].qual);


	/* Create "fake" matches */
	ctr=0;

	/* Based on mirroring type */
	if(FORWARD == ends[whichEnd-1].entries[0].strand) {
		if(MirrorForward == mirroringType ||
				MirrorBoth == mirroringType) {
			dest->contigs[ctr] = ends[whichEnd-1].entries[0].contig;
			dest->positions[ctr] = ends[whichEnd-1].entries[0].position + 
				ends[whichEnd-1].entries[0].referenceLength +
				offsetLength;
			dest->strands[ctr] = FORWARD;
			ctr++;
		}
		if(MirrorReverse == mirroringType ||
				MirrorBoth == mirroringType) {
			dest->contigs[ctr] = ends[whichEnd-1].entries[0].contig;
			dest->positions[ctr] = ends[whichEnd-1].entries[0].position -
				offsetLength - 
				ends[(2==whichEnd)?0:1].readLength +
				(ColorSpace==space)?1:0;
			dest->strands[ctr] = FORWARD;
			ctr++;
		}
	}
	else {
		if(MirrorForward == mirroringType ||
				MirrorBoth == mirroringType) {
			dest->contigs[ctr] = ends[whichEnd-1].entries[0].contig;
			dest->positions[ctr] = ends[whichEnd-1].entries[0].position -
				ends[whichEnd-1].entries[0].referenceLength -
				offsetLength;
			dest->strands[ctr] = REVERSE;
			ctr++;
		}
		if(MirrorReverse == mirroringType ||
			MirrorBoth == mirroringType) {
				dest->contigs[ctr] = ends[whichEnd-1].entries[0].contig;
				dest->positions[ctr] = ends[whichEnd-1].entries[0].position +
					offsetLength +
					ends[(2==whichEnd)?0:1].readLength -
					(ColorSpace==space)?1:0;
				dest->strands[ctr] = REVERSE;
				ctr++;
			}
	}
	/* Adjust/add strandedness */
	for(i=0;i<ctr;i++) {
		switch(strandedness) {
			case StrandSame:
				/* Ignore - this was the default. */
				break;
			case StrandOpposite:
				dest->strands[i] = (FORWARD==dest->strands[i])?REVERSE:FORWARD;
				break;
			case StrandBoth:
				dest->contigs[i+ctr] = dest->contigs[i];
				dest->positions[i+ctr] = dest->positions[i];
				dest->strands[i+ctr] = (FORWARD==dest->strands[i])?REVERSE:FORWARD;
				break;
			default:
				PrintError(FnName,
						"strandedness",
						"Could not understand strandedness",
						Exit,
						OutOfRange);
		}
	}
}
