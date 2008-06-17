#ifndef ALIGN_H_
#define ALIGN_H_
#include "../blib/AlignEntry.h"
#include "Definitions.h"

int AlignmentGetScore(char*, int, char*, int, ScoringMatrix*, AlignEntry*);
double GetScoreFromMatrix(char, char, ScoringMatrix*);
double GetMaximumOfTwoDoubles(double, double);
void ReverseSequence(char*, int);
void GetGapScore(double, double, double*, int*);
#endif
