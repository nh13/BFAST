#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <string.h>
#include "BLibDefinitions.h"
#include "BError.h"
#include "RGMatch.h"

/* TODO */
int32_t RGMatchRead(FILE *fp,
		RGMatch *m,
		int32_t binaryInput)
{
	char *FnName = "RGMatchRead";
	int32_t i;
	int32_t tempInt;

	/* Read the matches from the input file */
	if(binaryInput == 0) {
		/* Read the read length */
		if(fscanf(fp, "%d", &m->readLength)==EOF) {
			return EOF;
		}
		assert(m->readLength > 0);
		assert(m->readLength < SEQUENCE_LENGTH);

		/* Allocate memory for the read */
		m->read = malloc(sizeof(int8_t)*(m->readLength+1));
		if(NULL==m->read) {
			PrintError(FnName,
					"read",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Read in the read */
		if(fscanf(fp, "%s", m->read) ==EOF) {
			PrintError(FnName,
					"m->read",
					"Could not read in the read",
					Exit,
					EndOfFile);
		}

		/* Read in if we have reached the maximum number of matches */
		if(fscanf(fp, "%d", &m->maxReached)==EOF) {
			PrintError(FnName,
					"m->maxReached",
					"Could not read in m->maxReached",
					Exit,
					EndOfFile);
		}

		/* Read in the number of matches */
		if(fscanf(fp, "%d", &m->numEntries)==EOF) {
			PrintError(FnName,
					"m->numEntries",
					"Could not read in m->numEntries",
					Exit,
					EndOfFile);
		}
		assert(m->numEntries >= 0);

		/* Allocate memory for the matches */
		RGMatchReallocate(m, m->numEntries);

		/* Read first sequence matches */
		for(i=0;i<m->numEntries;i++) {
			if(fscanf(fp, "%d %d %c", 
						&tempInt,
						&m->positions[i],
						&m->strand[i])==EOF) {
				PrintError(FnName,
						NULL,
						"Could not read in match",
						Exit,
						EndOfFile);
			}
			m->chromosomes[i] = tempInt;
		}
	}
	else {
		/* Read in the read length */
		if(fread(&m->readLength, sizeof(int32_t), 1, fp)!=1) {
			if(feof(fp) != 0) {
				return EOF;
			}
			else {
				PrintError(FnName,
						"m->readLength",
						"Could not read in read length",
						Exit,
						ReadFileError);
			}
		}
		assert(m->readLength < SEQUENCE_LENGTH);
		assert(m->readLength > 0);

		/* Allocate memory for the read */
		m->read = malloc(sizeof(int8_t)*(m->readLength+1));
		if(NULL==m->read) {
			PrintError(FnName,
					"read",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}

		/* Read in the read */
		if(fread(m->read, sizeof(int8_t), m->readLength, fp)!=m->readLength) {
			PrintError(FnName,
					"m->read",
					"Could not read in the read",
					Exit,
					ReadFileError);
		}
		m->read[m->readLength]='\0';

		/* Read in if we have reached the maximum number of matches */
		if(fread(&m->maxReached, sizeof(int32_t), 1, fp)!=1) {
			PrintError(FnName,
					"m->maxReached",
					"Could not read in m->maxReached",
					Exit,
					ReadFileError);
		}
		assert(m->maxReached == 0 || m->maxReached == 1);

		/* Read in the number of matches */
		if(fread(&m->numEntries, sizeof(int32_t), 1, fp)!=1) {
			PrintError(FnName,
					"m->numEntries",
					"Could not read in m->numEntries",
					Exit,
					ReadFileError);
		}
		assert(m->numEntries >= 0);

		/* Allocate memory for the matches */
		RGMatchReallocate(m, m->numEntries);

		/* Read first sequence matches */
		if(fread(m->chromosomes, sizeof(uint8_t), m->numEntries, fp)!=m->numEntries) {
			PrintError(FnName,
					"m->chromosomes",
					"Could not read in chromosomes",
					Exit,
					ReadFileError);
		}
		if(fread(m->positions, sizeof(uint32_t), m->numEntries, fp)!=m->numEntries) {
			PrintError(FnName,
					"m->positions",
					"Could not read in positions",
					Exit,
					ReadFileError);
		}
		if(fread(m->strand, sizeof(int8_t), m->numEntries, fp)!=m->numEntries) {
			PrintError(FnName,
					"m->strand",
					"Could not read in strand",
					Exit,
					ReadFileError);
		}
	}

	return 1;
}

/* TODO */
void RGMatchPrint(FILE *fp,
		RGMatch *m,
		int32_t binaryOutput)
{
	char *FnName = "RGMatchPrint";
	int32_t i;
	assert(fp!=NULL);
	assert(m->readLength > 0);

	/* Print the matches to the output file */
	if(binaryOutput == 0) {
		if(0 > fprintf(fp, "%d\t%s\t%d\t%d",
					m->readLength,
					m->read,
					m->maxReached,
					m->numEntries)) {
			PrintError(FnName,
					NULL,
					"Could not write m->readLength, m->read, m->maxReached, and m->numEntries",
					Exit,
					WriteFileError);
		}

		for(i=0;i<m->numEntries;i++) {
			if(0 > fprintf(fp, "\t%d\t%d\t%c", 
						m->chromosomes[i],
						m->positions[i],
						m->strand[i])) {
				PrintError(FnName,
						NULL,
						"Could not write m->chromosomes[i], m->positions[i], and m->strand[i]",
						Exit,
						WriteFileError);
			}
		}
		if(0 > fprintf(fp, "\n")) {
			PrintError(FnName,
					NULL,
					"Could not write newline",
					Exit,
					WriteFileError);
		}
	}
	else {
		/* Print read length, read, maximum reached, and number of entries. */
		if(fwrite(&m->readLength, sizeof(int32_t), 1, fp) != 1 ||
				fwrite(m->read, sizeof(int8_t), m->readLength, fp) != m->readLength ||
				fwrite(&m->maxReached, sizeof(int32_t), 1, fp) != 1 ||
				fwrite(&m->numEntries, sizeof(int32_t), 1, fp) != 1) {
			PrintError(FnName,
					NULL,
					"Could not write m->readLength, m->read, m->maxReached, and m->numEntries",
					Exit,
					WriteFileError);
		}

		/* Print the chromosomes, positions, and strands */
		if(fwrite(m->chromosomes, sizeof(uint8_t), m->numEntries, fp) != m->numEntries ||
				fwrite(m->positions, sizeof(uint32_t), m->numEntries, fp) != m->numEntries ||
				fwrite(m->strand, sizeof(int8_t), m->numEntries, fp) != m->numEntries) {
			PrintError(FnName,
					NULL,
					"Could not write chromosomes, positions and strands",
					Exit,
					WriteFileError);
		}
	}
}

/* TODO */
void RGMatchRemoveDuplicates(RGMatch *m,
		int32_t maxNumMatches)
{
	int32_t i;
	int32_t prevIndex=0;

	/* Check to see if the max has been reached.  If so free all matches and return.
	 * We should remove duplicates before checking against maxNumMatches. */
	if(m->maxReached == 1) {
		/* Clear the matches but don't free the read name */
		RGMatchClearMatches(m);
		m->maxReached=1;
		return;
	}

	if(m->numEntries > 0) {
		/* Quick sort the data structure */
		RGMatchQuickSort(m, 0, m->numEntries-1);

		/* Remove duplicates */
		prevIndex=0;
		for(i=1;i<m->numEntries;i++) {
			if(RGMatchCompareAtIndex(m, prevIndex, m, i)==0) {
				/* ignore */
			}
			else {
				prevIndex++;
				/* Copy to prevIndex (incremented) */
				RGMatchCopyAtIndex(m, i, m, prevIndex);
			}
		}

		/* Reallocate pair */
		/* does not make sense if there are no entries */
		RGMatchReallocate(m, prevIndex+1);

		/* Check to see if we have too many matches */
		if(m->numEntries > maxNumMatches) {
			/* Clear the entries but don't free the read */
			RGMatchClearMatches(m);
			m->maxReached=1;
		}
		else { 
			m->maxReached = 0;
		}
	}
}

/* TODO */
void RGMatchQuickSort(RGMatch *m, int32_t low, int32_t high)
{
	int32_t i;
	int32_t pivot=-1;
	RGMatch *temp=NULL;

	if(low < high) {
		/* Allocate memory for the temp used for swapping */
		temp=malloc(sizeof(RGMatch));
		RGMatchInitialize(temp);
		if(NULL == temp) {
			PrintError("RGMatchQuickSort",
					"temp",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		RGMatchAllocate(temp, 1);

		pivot = (low+high)/2;

		RGMatchCopyAtIndex(m, pivot, temp, 0);
		RGMatchCopyAtIndex(m, high, m, pivot);
		RGMatchCopyAtIndex(temp, 0, m, high);

		pivot = low;

		for(i=low;i<high;i++) {
			if(RGMatchCompareAtIndex(m, i, m, high) <= 0) {
				if(i!=pivot) {
					RGMatchCopyAtIndex(m, i, temp, 0);
					RGMatchCopyAtIndex(m, pivot, m, i);
					RGMatchCopyAtIndex(temp, 0, m, pivot);
				}
				pivot++;
			}
		}
		RGMatchCopyAtIndex(m, pivot, temp, 0);
		RGMatchCopyAtIndex(m, high, m, pivot);
		RGMatchCopyAtIndex(temp, 0, m, high);

		/* Free temp before the recursive call, otherwise we have a worst
		 * case of O(n) space (NOT IN PLACE) 
		 * */
		RGMatchFree(temp);
		free(temp);
		temp=NULL;

		RGMatchQuickSort(m, low, pivot-1);
		RGMatchQuickSort(m, pivot+1, high);
	}
}

/* TODO */
int32_t RGMatchCompareAtIndex(RGMatch *mOne, int32_t indexOne, RGMatch *mTwo, int32_t indexTwo) 
{
	assert(indexOne >= 0 && indexOne < mOne->numEntries);
	assert(indexTwo >= 0 && indexTwo < mTwo->numEntries);
	if(mOne->chromosomes[indexOne] < mTwo->chromosomes[indexTwo] ||
			(mOne->chromosomes[indexOne] == mTwo->chromosomes[indexTwo] && mOne->positions[indexOne] < mTwo->positions[indexTwo]) ||
			(mOne->chromosomes[indexOne] == mTwo->chromosomes[indexTwo] && mOne->positions[indexOne] == mTwo->positions[indexTwo] && mOne->strand[indexOne] < mTwo->strand[indexTwo])) {
		return -1;
	}
	else if(mOne->chromosomes[indexOne] ==  mTwo->chromosomes[indexTwo] && mOne->positions[indexOne] == mTwo->positions[indexTwo] && mOne->strand[indexOne] == mTwo->strand[indexTwo]) {
		return 0;
	}
	else {
		return 1;
	}
}

/* TODO */
void RGMatchAppend(RGMatch *src, RGMatch *dest)
{
	char *FnName = "RGMatchAppend";
	int32_t i, start;

	/* Make sure we are not appending to ourselves */
	assert(src != dest);

	/* Check to see if we need to copy over the read as well */
	if(dest->readLength <= 0) {
		assert(dest->read == NULL);
		dest->readLength = src->readLength;
		/* Allocate memory */
		dest->read = malloc(sizeof(int8_t)*(dest->readLength+1));
		if(NULL==dest->read) {
			PrintError(FnName,
					"dest->read",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}   
		strcpy((char*)dest->read, (char*)src->read);
	}

	start = dest->numEntries;
	/* Allocate memory for the entires */
	RGMatchReallocate(dest, dest->numEntries + src->numEntries);

	assert(dest->numEntries == start + src->numEntries);
	assert(start <= dest->numEntries);

	/* Copy over the entries */
	for(i=start;i<dest->numEntries;i++) {
		RGMatchCopyAtIndex(src, i-start, dest, i);
	}
}

/* TODO */
void RGMatchCopyAtIndex(RGMatch *src, int32_t srcIndex, RGMatch *dest, int32_t destIndex)
{
	assert(srcIndex >= 0 && srcIndex < src->numEntries);
	assert(destIndex >= 0 && destIndex < dest->numEntries);

	if(src != dest || srcIndex != destIndex) {
		dest->positions[destIndex] = src->positions[srcIndex];
		dest->chromosomes[destIndex] = src->chromosomes[srcIndex];
		dest->strand[destIndex] = src->strand[srcIndex];
	}
}

/* TODO */
void RGMatchAllocate(RGMatch *m, int32_t numEntries)
{
	char *FnName = "RGMatchAllocate";
	assert(m->numEntries==0);
	m->numEntries = numEntries;
	assert(m->positions==NULL);
	m->positions = malloc(sizeof(uint32_t)*numEntries); 
	if(NULL == m->positions) {
		PrintError(FnName,
				"m->positions",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(m->chromosomes==NULL);
	m->chromosomes = malloc(sizeof(uint8_t)*numEntries); 
	if(NULL == m->chromosomes) {
		PrintError(FnName,
				"m->chromosomes",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(m->strand==NULL);
	m->strand = malloc(sizeof(int8_t)*numEntries); 
	if(NULL == m->strand) {
		PrintError(FnName,
				"m->strand",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
}

/* TODO */
void RGMatchReallocate(RGMatch *m, int32_t numEntries)
{
	char *FnName = "RGMatchReallocate";
	if(numEntries > 0) {
		m->numEntries = numEntries;
		m->positions = realloc(m->positions, sizeof(uint32_t)*numEntries); 
		if(numEntries > 0 && NULL == m->positions) {
			fprintf(stderr, "numEntries:%d\n", numEntries);
			PrintError(FnName,
					"m->positions",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		m->chromosomes = realloc(m->chromosomes, sizeof(uint8_t)*numEntries); 
		if(numEntries > 0 && NULL == m->chromosomes) {
			PrintError(FnName,
					"m->chromosomes",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		m->strand = realloc(m->strand, sizeof(int8_t)*numEntries); 
		if(numEntries > 0 && NULL == m->strand) {
			PrintError(FnName,
					"m->strand",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
	}
	else {
		/* Free just the matches part, not the meta-data */
		RGMatchClearMatches(m);
	}
}

/* TODO */
/* Does not free read */
void RGMatchClearMatches(RGMatch *m) 
{
	m->maxReached=0;
	m->numEntries=0;
	/* Free */
	free(m->chromosomes);
	free(m->positions);
	free(m->strand);
	m->chromosomes=NULL;
	m->positions=NULL;
	m->strand=NULL;
}

/* TODO */
void RGMatchFree(RGMatch *m) 
{
	free(m->read);
	free(m->chromosomes);
	free(m->positions);
	free(m->strand);
	RGMatchInitialize(m);
}

/* TODO */
void RGMatchInitialize(RGMatch *m)
{
	m->readLength=0;
	m->read=NULL;
	m->maxReached=0;
	m->numEntries=0;
	m->chromosomes=NULL;
	m->positions=NULL;
	m->strand=NULL;
}

/* TODO */
void RGMatchCheck(RGMatch *m)
{
	char *FnName="RGMatchCheck";
	/* Basic asserts */
	assert(m->readLength >= 0);
	assert(m->maxReached == 0 || m->maxReached == 1);
	assert(m->maxReached == 0 || m->numEntries == 0);
	assert(m->numEntries >= 0);
	/* Check that if the read length is greater than zero the read is not null */
	if(m->readLength > 0 && m->read == NULL) {
		PrintError(FnName,
				NULL,
				"m->readLength > 0 && m->read == NULL",
				Exit,
				OutOfRange);
	}
	/* Check that the read length matches the read */
	if(((int)strlen((char*)m->read)) != m->readLength) {
		PrintError(FnName,
				NULL,
				"m->readLength and strlen(m->read) do not match",
				Exit,
				OutOfRange);
	}
	/* Check that if the max has been reached then there are no entries */
	if(1==m->maxReached && m->numEntries > 0) {
		PrintError(FnName,
				NULL,
				"1==m->maxReached and m->numEntries>0",
				Exit,
				OutOfRange);
	}
	/* Check that if the number of entries is greater than zero that the entries are not null */
	if(m->numEntries > 0 && (m->chromosomes == NULL || m->positions == NULL || m->strand == NULL)) {
		PrintError(FnName,
				NULL,
				"m->numEntries > 0 && (m->chromosomes == NULL || m->positions == NULL || m->strand == NULL)",
				Exit,
				OutOfRange);
	}
}

/* TODO */
void RGMatchFilterOutOfRange(RGMatch *m,
		int32_t startChr,
		int32_t startPos,
		int32_t endChr,
		int32_t endPos,
		int32_t maxNumMatches)
{
	int32_t i, prevIndex;

	/* Filter chr/pos */
	/* Remove duplicates */
	prevIndex = -1;
	int filter;
	for(i=0;i<m->numEntries;i++) {
		filter = 0;
		if(m->chromosomes[i] < startChr || 
				(m->chromosomes[i] == startChr && (m->positions[i] + m->readLength - 1) < startPos) ||
				(m->chromosomes[i] == endChr && m->positions[i] > endPos) ||
				(m->chromosomes[i] > endChr)) {
			/* ignore */
		}
		else {
			/* Do not filter */
			prevIndex++;
			/* Copy chr/pos at i to chr/pos at prevIndex */
			RGMatchCopyAtIndex(m, i, m, prevIndex);
		}
	}

	/* Reallocate pair */
	RGMatchReallocate(m, prevIndex+1);

	/* Filter based on the maximum number of matches */
	if(maxNumMatches != 0 && m->numEntries > maxNumMatches) {
		/* Do not align this one */
		RGMatchClearMatches(m);
		assert(m->readLength > 0);
	}
}
