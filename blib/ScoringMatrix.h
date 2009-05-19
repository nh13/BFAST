#ifndef SCORINGMATRIX_H_
#define SCORINGMATRIX_H_

#include "BLibDefinitions.h"

#define ScoringMatrixGetNTScore(a, b, sm) (ToUpper(a) == ToUpper(b)) ? sm->ntMatch : sm->ntMismatch
#define ScoringMatrixGetColorScore(a, b, sm) (a == b)  ? sm->colorMatch : sm->colorMismatch

int ScoringMatrixRead(char*, ScoringMatrix*, int);
void ScoringMatrixInitialize(ScoringMatrix*);
int32_t ScoringMatrixCheck(ScoringMatrix*, int32_t);

#endif
