#ifndef ALIGNENTRIES_H_
#define ALIGNENTRIES_H_

#include <stdio.h>
#include "BLibDefinitions.h"
#include "RGBinary.h"

void AlignEntriesPrint(AlignEntries*, FILE*, int);
int AlignEntriesRead(AlignEntries*, FILE*, int, int, int);
void AlignEntriesRemoveDuplicates(AlignEntries*, int);
void AlignEntriesQuickSort(AlignEntries*, int, int);
void AlignEntriesMergeSort(AlignEntries*, int, int);
void AlignEntriesAllocate(AlignEntries*, char*, int, int, int, int);
void AlignEntriesReallocate(AlignEntries*, int, int, int, int);
void AlignEntriesFree(AlignEntries*);
void AlignEntriesInitialize(AlignEntries*);
void AlignEntriesKeepOnly(AlignEntries*, int, int, int, int);
void AlignEntriesCopy(AlignEntries*, AlignEntries*);
#endif

