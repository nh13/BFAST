#ifndef ALIGN_H_
#define ALIGN_H_

#include "Definitions.h"

double AlignmentGetScore(char*, int, char*, int, ScoringMatrix*);
double GetScoreFromMatrix(char, char, ScoringMatrix*);
double GetMaximumOfTwoDoubles(double, double);
#endif
