#ifndef REPAIR_H_
#define REPAIR_H_

#include "../blib/BLibDefinitions.h"
#include "../blib/ScoringMatrix.h"

int Repair(AlignedRead *src,
		RGBinary *rg,
		AlignedRead *dest,
		int32_t space,
		ScoringMatrix *sm,
		int32_t alignmentType,
		int32_t minPairedEndDistance,
		int32_t maxPairedEndDistance,
		int32_t mirroringType,
		int32_t strandedness);
void RepairMirrorFromAlignedEndToRGMatch(RGMatch *dest,
		int32_t whichEnd,
		AlignedEnd *ends,
		int32_t space,
		int32_t offsetLength,
		int32_t mirroringType,
		int32_t strandedness);
#endif
