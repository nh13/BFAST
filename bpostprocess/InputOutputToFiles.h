#ifndef INPUTOUTPUTTOFILES_H_
#define INPUTOUTPUTTOFILES_H_

#include "../blib/AlignEntries.h"

void ReadInputFilterAndOutput(RGBinary *rg,
		char *inputFileName,
		int binaryInput,
		int pairedEnd,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int algorithmReads,
		int uniquenessScore,
		int minUniquenessScore,
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

void PrintAlignEntriesToOutputFormat(AlignEntries *a,
		RGBinary *rg,
		FILE *fp,
		int outputFormat,
		int outputBinary);

void PrintAlignEntriesToMAF(AlignEntries *a,
		RGBinary *rg,
		FILE *fp);

void PrintAlignEntryToMAF(AlignEntry *a,
		RGBinary *rg,
		char *readName,
		int pairedEnd,
		int colorSpace,
		int readNum,
		FILE *fp);

#endif
