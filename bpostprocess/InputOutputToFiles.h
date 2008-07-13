#ifndef INPUTOUTPUTTOFILES_H_
#define INPUTOUTPUTTOFILES_H_

#include "../blib/AlignEntries.h"

void ReadInputFilterAndOutput(char *inputFileName,
		int binaryInput,
		int pairedEnd,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int algorithmReads,
		int minScoreReads,
		int algorithmReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int meanDistancePaired,
		int chrAbPaired,
		int inversionsPaired,
		char *outputID,
		char *outputDir,
		int outputFormat);

void PrintHeader(FILE *fp,
		int outputFormat);

void PrintAlignEntriesToOutputFormat(AlignEntries *a,
		FILE *fp,
		int outputFormat);

void PrintAlignEntriesToMAF(AlignEntries *a,
		FILE *fp);

void PrintAlignEntryToMAF(AlignEntry *a,
		char *readName,
		int pairedEnd,
		int readNum,
		FILE *fp);

#endif
