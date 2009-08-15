#ifndef INPUTOUTPUTTOFILES_H_
#define INPUTOUTPUTTOFILES_H_

#include "../blib/AlignedRead.h"

void ReadInputFilterAndOutput(RGBinary *rg,
		char *inputFileName,
		int algorithm,
		int queueLength,
		char *outputID,
		char *outputDir,
		int outputFormat);

int32_t GetAlignedReads(gzFile, AlignedRead*, int32_t);
#endif
