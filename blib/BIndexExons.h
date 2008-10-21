#ifndef BINDEXEXONS_H_
#define BINDEXEXONS_H_

#include "BLibDefinitions.h"

void BIndexExonsRead(char*, BIndexExons*);
int BIndexExonsWithin(BIndexExons*, uint32_t, uint32_t, uint32_t, uint32_t);
void BIndexExonsInitialize(BIndexExons*);
void BIndexExonsDelete(BIndexExons*);

#endif
