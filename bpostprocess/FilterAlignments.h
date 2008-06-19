#ifndef FILTERALIGNMENTS_H_
#define FILTERALIGNMENTS_H_

#include <stdio.h>
#include "../blib/AlignEntry.h"

void FilterAlignments(char*,
		int,
		int,
		int, 
		int,
		int,
		int,
		int,
		int,
		char*,
		char*,
		char*,
		int);

int FilterEntries(AlignEntry**,
		int,
		int,
		int,
		int,
		int,
		int,
		int,
		int);


#endif
