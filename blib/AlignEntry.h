#ifndef ALIGNENTRY_H_
#define ALIGNENTRY_H_

#include <stdio.h>
#include "BLibDefinitions.h"


enum {AlignEntrySortByAll, AlignEntrySortByChrPos};

/* TODO */
typedef struct {
	char readName[SEQUENCE_NAME_LENGTH];
	char *read; /* The read */
	char *reference;
	unsigned int length; /* The length of the alignment */
	unsigned int referenceLength; /* The number of bases excluding gaps in the reference string */
	unsigned int chromosome;
	unsigned int position;
	char strand;
	double score;
} AlignEntry;

void AlignEntryPrint(AlignEntry*, FILE*);
int AlignEntryRead(AlignEntry*, FILE*);
int AlignEntryRemoveDuplicates(AlignEntry**, int, int);
void AlignEntryQuickSort(AlignEntry**, int, int, int, int, double*, int);
void AlignEntryCopyAtIndex(AlignEntry*, int, AlignEntry*, int);
int AlignEntryCompareAtIndex(AlignEntry*, int, AlignEntry*, int, int);
int AlignEntryGetOneRead(AlignEntry**, FILE*);
int AlignEntryGetAll(AlignEntry**, FILE*);
void AlignEntryCopy(AlignEntry*, AlignEntry*);
#endif
