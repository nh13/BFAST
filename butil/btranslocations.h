#ifndef BTRANSLOCATIONS_H_
#define BTRANSLOCATIONS_H_

#include "../blib/AlignEntries.h"

typedef struct {
	int32_t contigStart;
	int32_t contigEnd;
	int32_t positionStart;
	int32_t positionEnd;
} Range;

#define BTRANSLOCATIONS_ROTATE_NUM 100000

void ParseRange(Range*, char*);
int32_t CheckRange(Range*, int32_t, int32_t);

#endif
