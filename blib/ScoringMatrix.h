#ifndef SCORINGMATRIX_H_
#define SCORINGMATRIX_H_

#include "BLibDefinitions.h"

int ScoringMatrixRead(char*, ScoringMatrix*, int);
void ScoringMatrixInitialize(ScoringMatrix*);
void ScoringMatrixFree(ScoringMatrix*);
int32_t ScoringMatrixGetNTScore(char, char, ScoringMatrix*);
int32_t ScoringMatrixGetColorScore(char, char, ScoringMatrix*);
int32_t ScoringMatrixCheck(ScoringMatrix*, int32_t);

#endif
