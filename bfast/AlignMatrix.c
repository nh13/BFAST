#include <stdlib.h>
#include "BLibDefinitions.h"
#include "BError.h"
#include "AlignMatrix.h"

void AlignMatrixReallocate(AlignMatrix *m, int32_t nrow, int32_t ncol)
{
	char *FnName="AlignMatrixReallocate";
	int32_t i, prevNRow, prevNCol;

	prevNRow = m->nrow;
	prevNCol = m->ncol;

	/* Shrink rows */
	if(nrow < prevNRow) {
		for(i=nrow;i<prevNRow;i++) {
			free(m->cells[i]);
		}
		m->cells = realloc(m->cells, sizeof(AlignMatrixCell*)*nrow);
		if(NULL == m->cells) {
			PrintError(FnName, "m->cells", "Could not reallocate memory", Exit, ReallocMemory);
		}
		m->nrow = nrow;
	}
	/* Shrink cols */
	if(ncol < prevNCol) {
		for(i=0;i<m->nrow;i++) {
			m->cells[i] = realloc(m->cells[i], sizeof(AlignMatrixCell)*ncol);
			if(NULL == m->cells[i]) {
				PrintError(FnName, "m->cells[i]", "Could not reallocate memory", Exit, ReallocMemory);
			}
		}
		m->ncol = ncol;
	}
	/* Expand rows */
	if(prevNRow < nrow) {
		m->cells = realloc(m->cells, sizeof(AlignMatrixCell*)*nrow);
		if(NULL == m->cells) {
			PrintError(FnName, "m->cells", "Could not reallocate memory", Exit, ReallocMemory);
		}
		for(i=prevNRow;i<nrow;i++) {
			m->cells[i] = malloc(sizeof(AlignMatrixCell)*ncol); // Allocate to the new number of columns ^^
			if(NULL == m->cells[i]) {
				PrintError(FnName, "m->cells[i]", "Could not allocate memory", Exit, MallocMemory);
			}
		}
		m->nrow = nrow;
	}
	/* Expand cols */
	if(prevNCol < ncol) {
		for(i=0;i<prevNRow;i++) { // if prevNRow != nrow, then we have already allocated ^^
			m->cells[i] = realloc(m->cells[i], sizeof(AlignMatrixCell)*ncol); 
			if(NULL == m->cells[i]) {
				PrintError(FnName, "m->cells[i]", "Could not reallocate memory", Exit, ReallocMemory);
			}
		}
		m->ncol = ncol;
	}
}

void AlignMatrixInitialize(AlignMatrix *m)
{
	m->cells=NULL;
	m->nrow=m->ncol=0;
}

void AlignMatrixFree(AlignMatrix *m)
{
	int32_t i;
	for(i=0;i<m->nrow;i++) {
		free(m->cells[i]);
	}
	free(m->cells);
	AlignMatrixInitialize(m);
}
