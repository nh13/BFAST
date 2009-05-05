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
		int minScores,
		int maxMismatches,
		int maxColorErrors,
		int minDistancePaired,
		int maxDistancePaired,
		int useDistancePaired,
		int contigAbPaired,
		int inversionsPaired,
		int unpaired,
		char *scoringMatrixFileName,
		int avgMismatchQuality,
		int space,
		char *outputID,
		char *outputDir,
		int outputFormat);

#endif
