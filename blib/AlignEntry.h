#ifndef ALIGNENTRY_H_
#define ALIGNENTRY_H_

#include <stdio.h>
#include "BLibDefinitions.h"
#include "RGBinary.h"

enum {AlignEntrySortByAll, AlignEntrySortByChrPos};

typedef struct {
	char *read; /* The read */
	char *reference;
	char *colorError;
	unsigned int length; /* The length of the alignment */
	int32_t chromosome;
	uint32_t position;
	char strand;
	double score;
} AlignEntry;

int AlignEntryPrint(AlignEntry*, FILE*, int);
int AlignEntryRead(AlignEntry*, FILE*, int);
int AlignEntryRemoveDuplicates(AlignEntry**, int, int);
void AlignEntryQuickSort(AlignEntry**, int, int, int, int, double*, int);
void AlignEntryMergeSort(AlignEntry**, int, int, int, int, double*, int);
void AlignEntryCopyAtIndex(AlignEntry*, int, AlignEntry*, int);
int AlignEntryCompareAtIndex(AlignEntry*, int, AlignEntry*, int, int);
int AlignEntryGetOneRead(AlignEntry**, FILE*);
int AlignEntryGetAll(AlignEntry**, FILE*);
void AlignEntryCopy(AlignEntry*, AlignEntry*);
void AlignEntryFree(AlignEntry*);
void AlignEntryInitialize(AlignEntry*);
void AlignEntryCheckReference(AlignEntry*, RGBinary*, int);
int AlignEntryGetPivot(AlignEntry*, int, int, int);
#endif
