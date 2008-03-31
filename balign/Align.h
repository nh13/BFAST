#ifndef ALIGN_H_
#define ALIGN_H_

#include "Definitions.h"

int AlignmentGetScore(char*, int, char*, int, ScoringMatrix*, AlignOutput*);
double GetScoreFromMatrix(char, char, ScoringMatrix*);
double GetMaximumOfTwoDoubles(double, double);
#endif
