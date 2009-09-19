#ifndef INPUTOUTPUTTOFILES_H_
#define INPUTOUTPUTTOFILES_H_

#include "AlignedRead.h"

void ReadInputFilterAndOutput(RGBinary *rg,
		char *inputFileName,
		int algorithm,
		int queueLength,
		int outputFormat,
		char *outputID,
		char *unmappedFileName);

int32_t GetAlignedReads(gzFile, AlignedRead*, int32_t);
#endif
