#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <string.h>

#include "BLib.h"
#include "BLibDefinitions.h"
#include "BError.h"
#include "RGMatch.h"

/* TODO */
int32_t RGMatchRead(gzFile fp,
		RGMatch *m)
{

	char *FnName = "RGMatchRead";

	/* Read in the read length */
	if(gzread64(fp, &m->readLength, sizeof(int32_t))!=sizeof(int32_t)||
			gzread64(fp, &m->qualLength, sizeof(int32_t))!=sizeof(int32_t)) {
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
	m->read = malloc(sizeof(char)*(m->readLength+1));
	if(NULL==m->read) {
		PrintError(FnName,
				"read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	m->qual = malloc(sizeof(char)*(m->qualLength+1));
	if(NULL==m->qual) {
		PrintError(FnName,
				"qual",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Read in the read */
	if(gzread64(fp, m->read, sizeof(char)*m->readLength)!=sizeof(char)*m->readLength||
			gzread64(fp, m->qual, sizeof(char)*m->qualLength)!=sizeof(char)*m->qualLength) {
		PrintError(FnName,
				"m->read",
				"Could not read in the read and qual",
				Exit,
				ReadFileError);
	}
	m->read[m->readLength]='\0';
	m->qual[m->qualLength]='\0';

	/* Read in if we have reached the maximum number of matches */
	if(gzread64(fp, &m->maxReached, sizeof(int32_t))!=sizeof(int32_t)) {
		PrintError(FnName,
				"m->maxReached",
				"Could not read in m->maxReached",
				Exit,
				ReadFileError);
	}
	assert(0 == m->maxReached || 1 == m->maxReached);

	/* Read in the number of matches */
	if(gzread64(fp, &m->numEntries, sizeof(int32_t))!=sizeof(int32_t)) {
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
	if(gzread64(fp, m->contigs, sizeof(uint32_t)*m->numEntries)!=sizeof(uint32_t)*m->numEntries) {
		PrintError(FnName,
				"m->contigs",
				"Could not read in contigs",
				Exit,
				ReadFileError);
	}
	if(gzread64(fp, m->positions, sizeof(uint32_t)*m->numEntries)!=sizeof(uint32_t)*m->numEntries) {
		PrintError(FnName,
				"m->positions",
				"Could not read in positions",
				Exit,
				ReadFileError);
	}
	if(gzread64(fp, m->strands, sizeof(char)*m->numEntries)!=sizeof(char)*m->numEntries) {
		PrintError(FnName,
				"m->strands",
				"Could not read in strand",
				Exit,
				ReadFileError);
	}

	return 1;
}

/* TODO */
int32_t RGMatchReadText(FILE *fp,
		RGMatch *m)
{

	char *FnName = "RGMatchRead";
	int32_t i;
	char read[SEQUENCE_LENGTH]="\0";
	char qual[SEQUENCE_LENGTH]="\0";

	/* Read the read and qual */
	if(fscanf(fp, "%s %s",
				read,
				qual)==EOF) {
		return EOF;
	}
	m->readLength = strlen(read);
	m->qualLength = strlen(qual);
	assert(m->readLength > 0);
	assert(m->readLength < SEQUENCE_LENGTH);

	/* Allocate memory for the read */
	m->read = malloc(sizeof(char)*(m->readLength+1));
	if(NULL==m->read) {
		PrintError(FnName,
				"read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	m->qual = malloc(sizeof(char)*(m->qualLength+1));
	if(NULL==m->qual) {
		PrintError(FnName,
				"qual",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	strcpy(m->read, read);
	strcpy(m->qual, qual);

	/* Read in if we have reached the maximum number of matches */
	if(fscanf(fp, "%d", &m->maxReached)==EOF) {
		PrintError(FnName,
				"m->maxReached",
				"Could not read in m->maxReached",
				Exit,
				EndOfFile);
	}
	assert(1==m->maxReached || 0 == m->maxReached);

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
		if(fscanf(fp, "%u %d %c", 
					&m->contigs[i],
					&m->positions[i],
					&m->strands[i])==EOF) {
			PrintError(FnName,
					NULL,
					"Could not read in match",
					Exit,
					EndOfFile);
		}
	}

	return 1;
}

/* TODO */
void RGMatchPrint(gzFile fp,
		RGMatch *m)
{
	char *FnName = "RGMatchPrint";
	assert(fp!=NULL);
	assert(m->readLength > 0);
	assert(m->qualLength > 0);

	/* Print the matches to the output file */
	/* Print read length, read, maximum reached, and number of entries. */
	if(gzwrite64(fp, &m->readLength, sizeof(int32_t))!=sizeof(int32_t) ||
			gzwrite64(fp, &m->qualLength, sizeof(int32_t))!=sizeof(int32_t) ||
			gzwrite64(fp, m->read, sizeof(char)*m->readLength)!=sizeof(char)*m->readLength ||
			gzwrite64(fp, m->qual, sizeof(char)*m->qualLength)!=sizeof(char)*m->qualLength ||
			gzwrite64(fp, &m->maxReached, sizeof(int32_t))!=sizeof(int32_t) ||
			gzwrite64(fp, &m->numEntries, sizeof(int32_t))!=sizeof(int32_t)) {
		PrintError(FnName,
				NULL,
				"Could not write m->readLength, m->qualLength, m->read, m->qual, m->maxReached, and m->numEntries",
				Exit,
				WriteFileError);
	}

	/* Print the contigs, positions, and strands */
	if(gzwrite64(fp, m->contigs, sizeof(uint32_t)*m->numEntries)!=sizeof(uint32_t)*m->numEntries ||
			gzwrite64(fp, m->positions, sizeof(uint32_t)*m->numEntries)!=sizeof(uint32_t)*m->numEntries ||
			gzwrite64(fp, m->strands, sizeof(char)*m->numEntries)!=sizeof(char)*m->numEntries) {
		PrintError(FnName,
				NULL,
				"Could not write contigs, positions and strands",
				Exit,
				WriteFileError);
	}
}

/* TODO */
void RGMatchPrintText(FILE *fp,
		RGMatch *m)
{
	char *FnName = "RGMatchPrint";
	int32_t i;
	assert(fp!=NULL);
	assert(m->readLength > 0);
	assert(m->qualLength > 0);

	/* Print the matches to the output file */
	if(0 > fprintf(fp, "%s\t%s\t%d\t%d",
				m->read,
				m->qual,
				m->maxReached,
				m->numEntries)) {
		PrintError(FnName,
				NULL,
				"Could not write m->read, m->qual, m->maxReached, and m->numEntries",
				Exit,
				WriteFileError);
	}

	for(i=0;i<m->numEntries;i++) {
		assert(m->contigs[i] > 0);
		if(0 > fprintf(fp, "\t%u\t%d\t%c",
					m->contigs[i],
					m->positions[i],
					m->strands[i])) {
			PrintError(FnName,
					NULL,
					"Could not write m->contigs[i], m->positions[i], and m->strands[i]",
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

/* TODO */
void RGMatchPrintFastq(FILE *fp,
		char *readName,
		RGMatch *m)
{
	char *FnName = "RGMatchPrintFastq";
	assert(fp!=NULL);
	assert(m->readLength > 0);
	assert(m->qualLength > 0);

	if(0 > fprintf(fp, "@%s\n%s\n+\n%s\n",
				readName,
				m->read,
				m->qual)) {
		PrintError(FnName,
				NULL,
				"Could not to file",
				Exit,
				WriteFileError);
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
	if(1 == m->maxReached) {
		/* Clear the matches but don't free the read name */
		RGMatchClearMatches(m);
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
				RGMatchCopyAtIndex(m, prevIndex, m, i);
			}
		}

		/* Reallocate pair */
		/* does not make sense if there are no entries */
		RGMatchReallocate(m, prevIndex+1);

		/* Check to see if we have too many matches */
		if(maxNumMatches < m->numEntries) {
			/* Clear the entries but don't free the read */
			RGMatchClearMatches(m);
			m->maxReached = 1;
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

		RGMatchCopyAtIndex(temp, 0, m, pivot);
		RGMatchCopyAtIndex(m, pivot, m, high);
		RGMatchCopyAtIndex(m, high, temp, 0);

		pivot = low;

		for(i=low;i<high;i++) {
			if(RGMatchCompareAtIndex(m, i, m, high) <= 0) {
				if(i!=pivot) {
					RGMatchCopyAtIndex(temp, 0, m, i);
					RGMatchCopyAtIndex(m, i, m, pivot);
					RGMatchCopyAtIndex(m, pivot, temp, 0);
				}
				pivot++;
			}
		}
		RGMatchCopyAtIndex(temp, 0, m, pivot);
		RGMatchCopyAtIndex(m, pivot, m, high);
		RGMatchCopyAtIndex(m, high, temp, 0);

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
	if(mOne->contigs[indexOne] < mTwo->contigs[indexTwo] ||
			(mOne->contigs[indexOne] == mTwo->contigs[indexTwo] && mOne->positions[indexOne] < mTwo->positions[indexTwo]) ||
			(mOne->contigs[indexOne] == mTwo->contigs[indexTwo] && mOne->positions[indexOne] == mTwo->positions[indexTwo] && mOne->strands[indexOne] < mTwo->strands[indexTwo])) {
		return -1;
	}
	else if(mOne->contigs[indexOne] ==  mTwo->contigs[indexTwo] && mOne->positions[indexOne] == mTwo->positions[indexTwo] && mOne->strands[indexOne] == mTwo->strands[indexTwo]) {
		return 0;
	}
	else {
		return 1;
	}
}

/* TODO */
void RGMatchAppend(RGMatch *dest, RGMatch *src)
{
	char *FnName = "RGMatchAppend";
	int32_t i, start;

	/* Make sure we are not appending to ourselves */
	assert(src != dest);
	assert(NULL != dest);
	assert(NULL != src);

	/* Check to see if we need to copy over the read as well */
	if(dest->readLength <= 0) {
		assert(dest->read == NULL);
		dest->readLength = src->readLength;
		dest->qualLength = src->qualLength;
		/* Allocate memory */
		dest->read = malloc(sizeof(char)*(dest->readLength+1));
		if(NULL==dest->read) {
			PrintError(FnName,
					"dest->read",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}   
		assert(dest->qual == NULL);
		dest->qual = malloc(sizeof(char)*(dest->qualLength+1));
		if(NULL==dest->qual) {
			PrintError(FnName,
					"dest->qual",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Copy over */
		strcpy(dest->read, src->read);
		strcpy(dest->qual, src->qual);
	}

	/* if the max has been reached by the start or dest, then ignore */
	if(1 != dest->maxReached && 1 != src->maxReached) { 
		/* Allocate memory for the entries */
		start = dest->numEntries;
		RGMatchReallocate(dest, dest->numEntries + src->numEntries);

		assert(dest->numEntries == start + src->numEntries);
		assert(start <= dest->numEntries);

		/* Copy over the entries */
		for(i=start;i<dest->numEntries;i++) {
			RGMatchCopyAtIndex(dest, i, src, i-start);
		}
	}
	else {
		/* Clear matches and set max reached flag */
		RGMatchClearMatches(dest);
	}
}

/* TODO */
void RGMatchCopyAtIndex(RGMatch *dest, int32_t destIndex, RGMatch *src, int32_t srcIndex)
{
	assert(srcIndex >= 0 && srcIndex < src->numEntries);
	assert(destIndex >= 0 && destIndex < dest->numEntries);

	if(src != dest || srcIndex != destIndex) {
		dest->positions[destIndex] = src->positions[srcIndex];
		dest->contigs[destIndex] = src->contigs[srcIndex];
		dest->strands[destIndex] = src->strands[srcIndex];
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
	assert(m->contigs==NULL);
	m->contigs = malloc(sizeof(uint32_t)*numEntries); 
	if(NULL == m->contigs) {
		PrintError(FnName,
				"m->contigs",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(m->strands==NULL);
	m->strands = malloc(sizeof(char)*numEntries); 
	if(NULL == m->strands) {
		PrintError(FnName,
				"m->strands",
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
			/*
			   fprintf(stderr, "numEntries:%d\n", numEntries);
			   */
			PrintError(FnName,
					"m->positions",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		m->contigs = realloc(m->contigs, sizeof(uint32_t)*numEntries); 
		if(numEntries > 0 && NULL == m->contigs) {
			PrintError(FnName,
					"m->contigs",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		m->strands = realloc(m->strands, sizeof(char)*numEntries); 
		if(numEntries > 0 && NULL == m->strands) {
			PrintError(FnName,
					"m->strands",
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
	m->numEntries=0;
	/* Free */
	free(m->contigs);
	free(m->positions);
	free(m->strands);
	m->contigs=NULL;
	m->positions=NULL;
	m->strands=NULL;
}

/* TODO */
void RGMatchFree(RGMatch *m) 
{
	free(m->read);
	free(m->qual);
	free(m->contigs);
	free(m->positions);
	free(m->strands);
	RGMatchInitialize(m);
}

/* TODO */
void RGMatchInitialize(RGMatch *m)
{
	m->readLength=0;
	m->qualLength=0;
	m->read=NULL;
	m->qual=NULL;
	m->maxReached=0;
	m->numEntries=0;
	m->contigs=NULL;
	m->positions=NULL;
	m->strands=NULL;
}

/* TODO */
void RGMatchCheck(RGMatch *m)
{
	char *FnName="RGMatchCheck";
	/* Basic asserts */
	assert(m->readLength >= 0);
	assert(m->qualLength >= 0);
	assert(m->maxReached == 0 || m->maxReached == 1);
	assert(m->maxReached == 0 || m->numEntries == 0);
	assert(m->numEntries >= 0);
	/* Check that if the read length is greater than zero the read is not null */
	if(m->readLength > 0 && m->read == NULL && m->qual == NULL) {
		PrintError(FnName,
				NULL,
				"m->readLength > 0 && m->read == NULL && m->qual == NULL",
				Exit,
				OutOfRange);
	}
	/* Check that the read length matches the read */
	if(((int)strlen(m->read)) != m->readLength) {
		PrintError(FnName,
				NULL,
				"m->readLength and strlen(m->read) do not match",
				Exit,
				OutOfRange);
	}
	/* Check that the qual length matches the qual */
	if(((int)strlen(m->qual)) != m->qualLength) {
		PrintError(FnName,
				NULL,
				"m->qualLength and strlen(m->qual) do not match",
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
	if(m->numEntries > 0 && (m->contigs == NULL || m->positions == NULL || m->strands == NULL)) {
		PrintError(FnName,
				NULL,
				"m->numEntries > 0 && (m->contigs == NULL || m->positions == NULL || m->strands == NULL)",
				Exit,
				OutOfRange);
	}
}

/* TODO */
void RGMatchFilterOutOfRange(RGMatch *m,
		int32_t maxNumMatches)
{
	/* Filter based on the maximum number of matches */
	if(maxNumMatches != 0 && m->numEntries > maxNumMatches) {
		/* Do not align this one */
		RGMatchClearMatches(m);
		m->maxReached=1;
		assert(m->readLength > 0);
	}
}
