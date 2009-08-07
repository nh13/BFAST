#ifndef INPUTOUTPUTTOFILES_H_
#define INPUTOUTPUTTOFILES_H_

#include "../blib/AlignedRead.h"

void ReadInputFilterAndOutput(RGBinary *rg,
		char *inputFileName,
		int algorithm,
		char *outputID,
		char *outputDir,
		int outputFormat);

#endif
