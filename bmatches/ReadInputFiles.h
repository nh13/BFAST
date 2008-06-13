#ifndef READINPUTFILES_H_
#define READINPUTFILES_H_

#include <stdio.h>
#include "../blib/RGIndex.h"

int ReadSequencesToTempFile(FILE*, FILE***, char***, int, int, int, int, char*, int);
int ReadNextSequence(FILE*, char**, int*, char**, int*, char**, int);
void ReadRGIndex(char*, RGIndex*, int);
int ReadFileNames(char*, char***);
int ReadOffsets(char*, int**);
void ReadTempSequencesAndOutput(FILE*, FILE*, FILE*, int, int);
void SequenceToLower(char*, int);

#endif
