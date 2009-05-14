#ifndef INPUTOUTPUTTOFILES_H_
#define INPUTOUTPUTTOFILES_H_

#include "../blib/AlignedRead.h"

void ReadInputFilterAndOutput(RGBinary *rg,
		char *inputFileName,
		int binaryInput,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int algorithm,
		int minScore,
		int minQual,
		int maxMismatches,
		int maxColorErrors,
		int minDistancePaired,
		int maxDistancePaired,
		int useDistancePaired,
		int contigAbPaired,
		int inversionsPaired,
		int unpaired,
		char *outputID,
		char *outputDir,
		int outputFormat);

#endif
