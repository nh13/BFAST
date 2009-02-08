#ifndef ALIGNENTRIESCONVERT_H_
#define ALIGNENTRIESCONVERT_H_

#include <stdlib.h>
#include <stdio.h>

#include "AlignEntries.h"
#include "AlignEntry.h"
#include "BError.h"

void AlignEntriesConvertPrintHeader(FILE*, int);
void AlignEntriesConvertPrintOutputFormat(AlignEntries*, RGBinary*, FILE*, int, int);
void AlignEntriesConvertPrintMAF(AlignEntries*, RGBinary*, FILE*);
void AlignEntriesConvertPrintAlignEntryToMAF(AlignEntry*, RGBinary*, char*, int, int, int, FILE*);
void AlignEntriesConvertPrintGFF(AlignEntries*, FILE*);
void AlignEntriesConvertPrintAlignEntryToGFF(AlignEntry*, char*, int, int, int, FILE*);
#endif
