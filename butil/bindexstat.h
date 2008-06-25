#ifndef BINDEXSTAT_H_
#define BINDEXSTAT_H_

#include "../blib/RGIndex.h"
#include "../blib/RGBinary.h"

void PrintSummary(RGIndex*, RGBinary*);
void PrintHistogram(RGIndex*, RGBinary*, int, int, FILE*);

#endif
