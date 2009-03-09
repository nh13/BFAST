#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#include "BLibDefinitions.h"
#include "BError.h"
#include "QS.h"

void QSInitialize(QS *qs, int32_t minScore, int32_t maxScore) 
{
	char *FnName="QSInitialize";
	int32_t i;
	assert(minScore <= maxScore);
	qs->scores=NULL;
	qs->maxDiff = maxScore - minScore;
	assert(0 <= qs->maxDiff);
	qs->total=0;

	qs->scores = calloc(sizeof(int32_t), qs->maxDiff+1);
	if(NULL == qs->scores) {
		PrintError(FnName,
				"qs->scores",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	for(i=0;i<=qs->maxDiff;i++) {
		qs->scores[i] = 1; /* Psuedo counting, sp prob > 0 */
	}
}

void QSAdd(QS *qs, int32_t scoreDiff)
{
	assert(0 <= scoreDiff && 
			scoreDiff <= qs->maxDiff);
	qs->scores[scoreDiff]++;
	qs->total++;
}

void QSMerge(QS *dest, QS *src)
{
	assert(dest->maxDiff == src->maxDiff);
	int32_t i;

	for(i=0;i<=dest->maxDiff;i++) {
		dest->scores[i] += src->scores[i];
	}
	dest->total += src->total;
}

int32_t QSGet(QS *qs, int32_t scoreDiff)
{
	if(scoreDiff < 0 || qs->maxDiff < scoreDiff) {
		return 0;
	}
	/* c = 10/log(10) ~ 3.434 */
	/* return -c log (pr) */
	assert(0 <= scoreDiff && 
			scoreDiff <= qs->maxDiff);
	return (int32_t)(-3.434*(log(qs->scores[scoreDiff]) - log(qs->total))); 
}

void QSFree(QS *qs)
{
	free(qs->scores);
	qs->scores=NULL;
	qs->maxDiff=-1;
}
