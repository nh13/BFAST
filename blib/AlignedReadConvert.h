#ifndef ALIGNEDREADCONVERT_H_
#define ALIGNEDREADCONVERT_H_

#include <stdlib.h>
#include <stdio.h>

#include "AlignedRead.h"
#include "AlignedEntry.h"
#include "BError.h"

void AlignedReadConvertPrintHeader(FILE*, int);
void AlignedReadConvertPrintOutputFormat(AlignedRead*, RGBinary*, FILE*, int, int);
void AlignedReadConvertPrintMAF(AlignedRead*, RGBinary*, FILE*);
void AlignedReadConvertPrintAlignedEntryToMAF(AlignedEntry*, RGBinary*, char*, char*, int, int, int, FILE*);
void AlignedReadConvertPrintGFF(AlignedRead*, FILE*);
void AlignedReadConvertPrintAlignedEntryToGFF(AlignedEntry*, char*, char*, int, int, int, FILE*);
#endif
