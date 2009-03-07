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
		int algorithmReads,
		int mappingQuality,
		int minMappingQuality,
		int minScoreReads,
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

void PrintHeader(FILE *fp,
		int outputFormat);

#endif
