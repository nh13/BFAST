#ifndef READINPUTFILES_H_
#define READINPUTFILES_H_

#include <stdio.h>
#include "../blib/RGTree.h"
#include "../blib/RGIndex.h"

int ReadSequencesToTempFile(FILE*, FILE***, int, int, int, int, int);
int ReadNextSequence(FILE*, char**, char**, char**, int);
void ReadRGIndex(char*, RGIndex*, int);
void ReadRGTree(char*, RGTree*, int);
int ReadFileNames(char*, char***);
int ReadOffsets(char*, int**);
int ReadRGTreeFileNames(char*, char***, int**, int*);
void ReadTempSequencesAndOutput(FILE*, FILE*, FILE*, int);
void SequenceToLower(char*, int);

#endif
