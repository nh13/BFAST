#ifndef ALIGNEDREADCONVERT_H_
#define ALIGNEDREADCONVERT_H_

#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>

#include "AlignedRead.h"
#include "AlignedEntry.h"
#include "BError.h"

void AlignedReadConvertPrintHeader(FILE*, RGBinary*, int);
void AlignedReadConvertPrintOutputFormat(AlignedRead*, RGBinary*, FILE*, gzFile, char*, int, int);
void AlignedReadConvertPrintMAF(AlignedRead*, RGBinary*, FILE*);
void AlignedReadConvertPrintAlignedEntryToMAF(AlignedEntry*, RGBinary*, char*, char*, int, int, int, FILE*);
void AlignedReadConvertPrintGFF(AlignedRead*, FILE*);
void AlignedReadConvertPrintAlignedEntryToGFF(AlignedEntry*, char*, char*, int, int, int, FILE*);
void AlignedReadConvertPrintSAM(AlignedRead*, char*, FILE*);
void AlignedReadConvertPrintAlignedEntryToSAM(AlignedRead*, int32_t, int32_t, char*, FILE*);
void AlignedReadConvertPrintAlignedEntryToCIGAR(AlignedEntry*, int32_t, FILE*);

#endif
