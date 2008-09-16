#ifndef SCORINGMATRIX_H_
#define SCORINGMATRIX_H_

#include "Definitions.h"

int ScoringMatrixRead(char*, ScoringMatrix*, int);
void ScoringMatrixInitialize(ScoringMatrix*);
void ScoringMatrixFree(ScoringMatrix*);
double ScoringMatrixGetNTScore(char, char, ScoringMatrix*);
double ScoringMatrixGetColorScore(uint8_t, uint8_t, ScoringMatrix*);

#endif
