#ifndef BFIXBAF_H_
#define BFIXBAF_H_

void ConvertBAF(char*, char*, int32_t, int32_t, int32_t);
int32_t AlignedReadReadOld(AlignedRead*, FILE*, int32_t);
int32_t AlignedEndReadOld(AlignedEnd*, FILE*, int32_t, int32_t);
int32_t AlignedEntryReadOld(AlignedEntry*, FILE*, int32_t, int32_t);

#endif
