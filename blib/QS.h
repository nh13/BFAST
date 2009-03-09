#ifndef QS_H_
#define QS_H_

#include "BLibDefinitions.h"

void QSInitialize(QS*, int32_t, int32_t);
void QSAdd(QS*, int32_t);
void QSMerge(QS*, QS*);
int32_t QSGet(QS*, int32_t);
void QSFree(QS*);

#endif
