#ifndef RGSEQPAIR_H_
#define RGSEQPAIR_H_

#include <stdio.h>
#include "BLibDefinitions.h"
#include "RGMatch.h"
#include "RGTree.h"

void RGSeqPairFindMatches(RGTree*, RGMatch*, char*, int**, int, int, int, int);
void RGSeqPairGeneratePairs(char*, RGSeqPair*, int**, int, int, int, int, int, int);
void RGSeqPairGenerateMismatches(char*, int, char, int, int, int, int, RGSeqPair*);
void RGSeqPairGenerateMismatchesHelper(char*, char, int, int, int, int, int, int, RGSeqPair*, char*, char*, int, int);
void RGSeqPairGenerateDeletions(char*, int, char, int, int, int, int, RGSeqPair*);
void RGSeqPairGenerateDeletionsHelper(char*, int, char, int, int, int, int, int, RGSeqPair*, char*, char*, int, int);
void RGSeqPairGenerateInsertions(char*, int, char, int, int, int, int, RGSeqPair*);
void RGSeqPairGenerateInsertionsHelper(char*, int, char, int, int, int, int, int, RGSeqPair*, char*, char*, int, int);

void RGSeqPairRemoveDuplicates(RGSeqPair*);
void RGSeqPairMergeSort(RGSeqPair*, int, int);

void GetReverseCompliment(char*, char*, int);

int RGSeqPairCompareAtIndex(RGSeqPair*, int, RGSeqPair*, int);
void RGSeqPairCopyAtIndex(RGSeqPair*, int, RGSeqPair*, int);

#endif

