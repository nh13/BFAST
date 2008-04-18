#ifndef PRINTOUTPUTFILES_H_
#define PRINTOUTPUTFILES_H_

#include <stdio.h>
#include "../blib/AlignEntry.h"

void ConvertAndPrint(void *,
		int,
		int,
		int,
		char*,
		char*,
		int);

void ConvertAndPrintAlignEntryToWig(AlignEntry*,
		int,
		FILE*);

#endif
