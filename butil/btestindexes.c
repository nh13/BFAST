#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "../blib/BError.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BLib.h"
#include "btestindexes.h"

#define GENERATE_HELPER_ROTATE_NUM 1000000
#define RANDOM_ROTATE_NUM 100
#define ROTATE_NUM 10000

#define NUM_MISMATCHES_START 1
#define NUM_MISMATCHES_END 6 
#define NUM_INSERTIONS_START 1
#define NUM_INSERTIONS_END 1
#define NUM_DELETIONS_START 1
#define NUM_DELETIONS_END 1

/* Is a utility that tests, searches for, and compares layouts for indexes against certain events,
 * such as errors, mismatches and insertions.
 * */

char Colors[5] = "01234";

/* TODO */
void RunSearchForIndexes(int readLength,
		int numEventsToSample,
		int numIndexesToSample,
		int keySize,
		int maxKeyWidth,
		int maxIndexSetSize,
		int accuracyThreshold,
		int space,
		int maxNumMismatches,
		int maxNumColorErrors)
{
	int i, j;
	IndexSet set;
	Index curBest, curIndex;
	AccuracyProfile curBestP, curIndexP;

	/* Seed random number */
	srand(time(NULL));

	/* Initialize index set */
	IndexSetInitialize(&set);

	/* Will always seed with contiguous 1s mask */
	IndexSetSeed(&set,
			keySize);

	for(i=2;i<=maxIndexSetSize;i++) { /* Add one index to the set */
		/* Initialize */
		IndexInitialize(&curBest); 
		AccuracyProfileInitialize(&curBestP);
		for(j=0;j<numIndexesToSample;j++) { /* Sample the space of possible indexes */
			/* Initialze cur */
			IndexInitialize(&curIndex);
			AccuracyProfileInitialize(&curIndexP);

			/* Get random index */
			do {
				IndexFree(&curIndex);
				IndexGetRandom(&curIndex,
						keySize,
						maxKeyWidth);
			}
			while(1==IndexSetContains(&set, &curIndex));

			/* Push random index */
			IndexSetPush(&set,
					&curIndex);

			/* Get accuracy profile */
			AccuracyProfileUpdate(&set,
					&curIndexP,
					readLength,
					numEventsToSample,
					space,
					maxNumMismatches,
					maxNumColorErrors);
			assert(curIndexP.numReads == numEventsToSample);

			/* Compare accuracy profile */
			if(AccuracyProfileCompare(&curBestP, &curIndexP, accuracyThreshold) < 0) {
				/* Copy index over to the current best */
				IndexCopy(&curBest, &curIndex);
				AccuracyProfileCopy(&curBestP, &curIndexP);
			}
			/* Pop random index */
			IndexSetPop(&set);

			/* Free cur index */
			IndexFree(&curIndex);
			AccuracyProfileFree(&curIndexP);
		}
		/* Add the best index to the set */
		IndexSetPush(&set,
				&curBest);
		/* Save accuracy profile */
		/* Free best index */
		IndexFree(&curBest);
	}

	/* Print */
	IndexSetPrint(&set, stderr);

	/* Free */
	IndexSetFree(&set);
}

/* TODO */
void RunEvaluateIndexes(char *inputFileName,
		int readLength,
		int numEventsToSample,
		int space,
		int maxNumMismatches,
		int maxInsertionLength,
		int maxNumColorErrors)
{
	char *FnName="RunEvaluateIndexes";
	int setSize;
	IndexSet set, curSet;

	assert(space == 1 || maxNumColorErrors == 0);

	/* Seed random number */
	srand(time(NULL));

	IndexSetInitialize(&curSet);

	/* Read in */
	IndexSetRead(&set, inputFileName);

	for(setSize=1;setSize<=set.numIndexes;setSize++) { /* For increasing set size */

		/* Add an index to the set */
		IndexSetPush(&curSet, &set.indexes[setSize-1]); 

		switch(space) {
			case NTSpace:
				RunEvaluateIndexesNTSpace(&curSet,
						readLength,
						numEventsToSample,
						maxNumMismatches,
						maxInsertionLength);
				break;
			case ColorSpace:
				RunEvaluateIndexesColorSpace(&curSet,
						readLength,
						numEventsToSample,
						maxNumMismatches,
						maxInsertionLength,
						maxNumColorErrors);
				break;
			default:
				PrintError(FnName,
						"space",
						"Could not understand space",
						Exit,
						OutOfRange);
		}
	}

	/* Free memory */
	IndexSetFree(&set);
	IndexSetFree(&curSet);
}

void RunEvaluateIndexesNTSpace(IndexSet *set,
		int readLength,
		int numEventsToSample,
		int maxNumMismatches,
		int maxInsertionLength)
{
	int i, j;

	/* Print the header */
	fprintf(stderr, "%d\n", set->numIndexes);
	fprintf(stderr, "%-5s\t", "MM"); /* # of Mismatches */
	fprintf(stderr, "%-5s\t", "BP:0"); /* Mismatches accuracy */
	fprintf(stderr, "%-5s\t", "1=del"); /* Deletions */
	for(j=1;j<=maxInsertionLength;j++) {
		fprintf(stderr, "%2s%-3d\t", "2i", j); /* Insertions */
	}
	fprintf(stderr, "\n");

	/* SNPs - include no SNPs */
	/* Mismatches including zero */
	for(i=0;i<=maxNumMismatches;i++) {
		fprintf(stderr, "%-5d\t", i);
		fprintf(stderr, "%1.3lf\t",
				GetNumCorrect(set,
					readLength,
					numEventsToSample,
					i,
					0,
					NoIndelType,
					0,
					NTSpace)/((double)numEventsToSample));
		/* Deletion with Mismatches */
		fprintf(stderr, "%1.3lf\t",
				GetNumCorrect(set,
					readLength,
					numEventsToSample,
					i,
					0,
					DeletionType,
					0,
					NTSpace)/((double)numEventsToSample));
		/* Insertions with Mismatches */
		for(j=1;j<=maxInsertionLength;j++) {
			fprintf(stderr, "%1.3lf\t",
					GetNumCorrect(set,
						readLength,
						numEventsToSample,
						i,
						0,
						InsertionType,
						j,
						NTSpace)/((double)numEventsToSample));
		}
		fprintf(stderr, "\n");
	}
}

void RunEvaluateIndexesColorSpace(IndexSet *set,
		int readLength,
		int numEventsToSample,
		int maxNumMismatches,
		int maxInsertionLength,
		int maxNumColorErrors)
{
	int i, j;
	/* Print the header */
	fprintf(stderr, "%d\n", set->numIndexes);
	/* TODO */

	/* Get accuracy and print out */
	for(i=0;i<=maxNumColorErrors;i++) {
		/* SNPs with color errors - include no SNPs */
		for(j=0;j<=maxNumMismatches;j++) {
			fprintf(stderr, "%1.3lf\t",
					GetNumCorrect(set,
						readLength,
						numEventsToSample,
						j,
						i,
						NoIndelType,
						0,
						ColorSpace)/((double)numEventsToSample));
		}
		/* Deletion with color errors */
		fprintf(stderr, "%1.3lf\t",
				GetNumCorrect(set,
					readLength,
					numEventsToSample,
					0,
					i,
					DeletionType,
					0,
					ColorSpace)/((double)numEventsToSample));
		/* Insertions with color errors */
		for(j=1;j<=maxInsertionLength;j++) {
			fprintf(stderr, "%1.3lf\t",
					GetNumCorrect(set,
						readLength,
						numEventsToSample,
						0,
						i,
						InsertionType,
						j,
						ColorSpace)/((double)numEventsToSample));
		}
		fprintf(stderr, "\n");
	}
}

int32_t GetNumCorrect(IndexSet *set,
		int readLength,
		int numEventsToSample,
		int numSNPs,
		int numColorErrors,
		int indelType,
		int insertionLength,
		int space)
{
	char *FnName="GetNumCorrect";
	assert(space == 1 || numColorErrors == 0);
	assert(insertionLength <= 0 || indelType == InsertionType);

	int32_t i;
	int32_t numCorrect = 0;
	int32_t breakpoint;
	Read curRead, r1, r2;

	ReadInitialize(&curRead);

	for(i=0;i<numEventsToSample;i++) {
		/* Get random read with SNPs and ColorErrors */
		ReadGetRandom(&curRead,
				readLength,
				numSNPs,
				numColorErrors,
				space);
		/* Get the breakpoint:
		 * SNPs - no breakpoint (0)
		 * Deletion - breakpoint within the read 
		 * Insertion - breakpoint within the read, including start
		 * */
		switch(indelType) {
			case NoIndelType:
				/* Only SNPs and color errors */
				assert(insertionLength == 0);
				/* Check read */
				numCorrect += IndexSetCheckRead(set, &curRead);
				break;
			case DeletionType:
				assert(insertionLength == 0);
				/* Get where the break point occured for the deletion */
				breakpoint = ( rand()%(readLength - 1) ) + 1;
				assert(breakpoint > 0);
				assert(readLength - breakpoint > 0);
				/* Split read into two reads based on the breakpoint */
				ReadSplit(&curRead, &r1, &r2, breakpoint, 0);
				/* Check either end of the read after split */
				if(1==IndexSetCheckRead(set, &r1) ||
						1==IndexSetCheckRead(set, &r2)) {
					numCorrect++;
				}
				/* Free read */
				ReadFree(&r1);
				ReadFree(&r2);
				break;
			case InsertionType:
				assert(insertionLength > 0);
				/* Get where the insertion occured relative to the start of hte read */
				breakpoint = (rand()%(readLength-insertionLength));
				/* Split read into two reads */
				ReadSplit(&curRead, &r1, &r2, breakpoint, insertionLength);
				/* Check either end of the read after split, substracting the insertion */
				if(1==IndexSetCheckRead(set, &r1) ||
						1==IndexSetCheckRead(set, &r2)) {
					numCorrect++;
				}
				/* Free read */
				ReadFree(&r1);
				ReadFree(&r2);
				break;
			default:
				PrintError(FnName,
						"indelType",
						"Could not understand indel type",
						Exit,
						OutOfRange);
		}
		/* Free read */
		ReadFree(&curRead);
	}
	return numCorrect;
}

void IndexSetRead(IndexSet *set,
		char *inputFileName)
{
	char *FnName="IndexSetRead";
	FILE *fp;
	Index index;

	if(!(fp=fopen(inputFileName, "rb"))) {
		PrintError(FnName,
				inputFileName,
				"Could not open file for reading",
				Exit,
				OpenFileError);
	}

	IndexSetInitialize(set);
	IndexInitialize(&index);

	while(EOF!=IndexRead(&index, fp)) {
		IndexSetPush(set, &index);
		IndexFree(&index);
	}

	fclose(fp);
}

int IndexSetContains(IndexSet *set,
		Index *index)
{
	int i;
	for(i=0;i<set->numIndexes;i++) {
		/* If they are the same, return 1 */
		if(IndexCompare(index, &set->indexes[i])==0) {
			return 1;
		}
	}
	/* Return zero */
	return 0;
}

int32_t IndexSetCheckRead(IndexSet *set,
		Read *r)
{
	int i;
	for(i=0;i<set->numIndexes;i++) {
		if(IndexCheckRead(&set->indexes[i],
					r) == 1) {
			return 1;
		}
	}
	return 0;
}

void IndexSetPush(IndexSet *set, 
		Index *index)
{
	char *FnName="IndexSetPush";
	set->numIndexes++;
	set->indexes = realloc(set->indexes, sizeof(Index)*set->numIndexes);
	if(NULL == set->indexes) {
		PrintError(FnName,
				"set->indexes",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}
	IndexInitialize(&set->indexes[set->numIndexes-1]);
	IndexCopy(&set->indexes[set->numIndexes-1], index);
}

void IndexSetPop(IndexSet *set)
{
	char *FnName="IndexSetPop";
	IndexFree(&set->indexes[set->numIndexes-1]);
	set->numIndexes--;
	set->indexes = realloc(set->indexes, sizeof(Index)*set->numIndexes);
	if(NULL == set->indexes) {
		PrintError(FnName,
				"set->indexes",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}
}

/* Seed index set with one index with a contiguous mask */
void IndexSetSeed(IndexSet *set,
		int keySize)
{
	char *FnName="IndexSetSeed";
	int i;

	/* Allocate for the index set */
	set->numIndexes=1;
	set->indexes = malloc(sizeof(Index)*set->numIndexes);
	if(NULL == set->indexes) {
		PrintError(FnName,
				"set->indexes",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate index */
	IndexAllocate(&set->indexes[set->numIndexes-1],
			keySize,
			keySize);
	/* Initialize the mask to ones */
	for(i=0;i<set->indexes[set->numIndexes-1].keySize;i++) {
		set->indexes[set->numIndexes-1].mask[i] = 1;
	}
}

void IndexSetInitialize(IndexSet *set) 
{
	set->indexes=NULL;
	set->numIndexes=0;
}

void IndexSetFree(IndexSet *set) 
{
	int i;
	for(i=0;i<set->numIndexes;i++) {
		IndexFree(&set->indexes[i]);
	}
	free(set->indexes);
	IndexSetInitialize(set);
}

void IndexSetPrint(IndexSet *set,
		FILE *fp)
{
	int i;
	for(i=0;i<set->numIndexes;i++) {
		IndexPrint(&set->indexes[i],
				fp);
	}
}

int IndexRead(Index *index, 
		FILE *fp)
{
	char *FnName = "IndexRead";
	char tempMask[2056]="\0";
	int32_t i;

	/* Read */
	if(EOF == fscanf(fp, "%s", tempMask)) {
		return EOF;
	}
	IndexAllocate(index,
			0,
			(int)strlen(tempMask));
	index->keySize = 0;
	for(i=0;i<index->keyWidth;i++) {
		switch(tempMask[i]) {
			case '0':
				index->mask[i] = 0;
				break;
			case '1':
				index->mask[i] = 1;
				index->keySize++;
				break;
			default:
				PrintError(FnName,
						"mask",
						"Could not read mask",
						Exit,
						OutOfRange);
		}
	}

	return 1;
}

int IndexCompare(Index *a, Index *b) 
{
	int i;
	if(a->keySize == b->keySize &&
			a->keyWidth == b->keyWidth) {
		/* key size and key width are the same */

		/* Compare masks */
		for(i=0;i<a->keyWidth;i++) {
			if(a->mask[i] != b->mask[i]) {
				/* Different */
				return 1;
			}
		}
		/* They must be the same */
		return 0;
	}
	else {
		return 1;
	}
}

int32_t IndexCheckRead(Index *index,
		Read *r)
{
	int i, j;
	int success;

	if(index->keyWidth > r->length) {
		return 0;
	}

	for(i=0;i<r->length - index->keyWidth + 1;i++) { /* For all possible offsets */
		success = 1;
		for(j=0;1==success && j<index->keyWidth;j++) { /* Go over the index mask */
			if(index->mask[j] == 1 && r->profile[j+i] == 1) {
				success = 0;
			}
		}
		if(1==success) {
			return 1;
		}
	}
	return 0;
}

void IndexCopy(Index *dest, Index *src)
{
	int i;

	if(NULL != dest->mask) {
		IndexFree(dest);
	}

	IndexAllocate(dest,
			src->keySize,
			src->keyWidth);
	for(i=0;i<src->keyWidth;i++) {
		dest->mask[i] = src->mask[i];
	}
}

void IndexGetRandom(Index *index,
		int keySize,
		int maxKeyWidth)
{
	char *FnName="IndexGetRandom";
	int i, j, k;
	int numLeft;
	int32_t *bins=NULL;
	int numBins = keySize-1;

	/* Generate random masks by having bins inbetween each "1".  We 
	 * then randomly insert zeros into the bins */

	/* Allocate memory for the bins */
	bins = malloc(sizeof(int32_t)*numBins);
	if(NULL == bins) {
		PrintError(FnName,
				"bins",
				"Could not allocate memory",
				Exit,
				OutOfRange);
	}
	/* Initialize bins */
	for(i=0;i<numBins;i++) {
		bins[i] = 0;
	}

	/* Choose a number of zeros to insert into the bins */
	numLeft = rand()%(maxKeyWidth - keySize + 1);
	assert(numLeft >=0 && numLeft <= maxKeyWidth - keySize);

	/* Allocate memory for the index */
	IndexAllocate(index,
			keySize,
			keySize + numLeft);

	/* Insert into bins */
	while(numLeft > 0) {
		/* choose a bin between 1 and keySize-1 */
		i = (rand()%numBins); /* Note: this is not truly inform, but a good approximation */
		assert(i>=0 && i<numBins);
		bins[i]++;
		numLeft--;
	}

	/* Generate index based off the bins */
	/* First base is always a 1 */ 
	for(i=0, j=1, index->mask[0] = 1;
			i<index->keySize-1;
			i++, j++) {
		/* Insert zero based on the bin size */
		for(k=0;
				k<bins[i];
				k++, j++) {
			assert(j<index->keyWidth);
			index->mask[j] = 0;
		}
		/* Add a one */
		assert(j<index->keyWidth);
		index->mask[j] = 1;
	}
	assert(index->keyWidth == j);

	/* Free memory */
	free(bins);
	bins=NULL;
}

void IndexAllocate(Index *index,
		int keySize,
		int keyWidth)
{
	char *FnName = "IndexAllocate";
	index->keySize = keySize;
	index->keyWidth = keyWidth;
	index->mask = malloc(sizeof(int32_t)*index->keyWidth);
	if(NULL == index->mask) {
		PrintError(FnName,
				"index->mask",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
}

void IndexInitialize(Index *index)
{
	index->mask = NULL;
	index->keySize = 0;
	index->keyWidth = 0;
}

void IndexFree(Index *index)
{
	free(index->mask);
	IndexInitialize(index);
}

void IndexPrint(Index *index,
		FILE *fp)
{
	int i;
	fprintf(fp, "%d\t%d\t",
			index->keyWidth,
			index->keySize);
	for(i=0;i<index->keyWidth;i++) {
		fprintf(fp, "%1d", index->mask[i]);
	}
	fprintf(fp, "\n");
}

int AccuracyProfileCompare(AccuracyProfile *a, 
		AccuracyProfile *b,
		int accuracyThreshold)
{
	int i;
	assert(a->length > 0 || b->length > 0);

	if(a->numReads <= 0) {
		return -1;
	}
	if(b->numReads <= 0) {
		return 1;
	}
	assert(a->numReads > 0);
	assert(b->numReads > 0);
	assert(a->length == b->length);
	for(i=0;i<a->length;i++) {
		if((a->correct[i]/a->numReads) < accuracyThreshold || (b->correct[i]/b->numReads) < accuracyThreshold) {
			if(a->correct[i]/a->numReads < b->correct[i]/b->numReads) {
				return -1;
			}
			else if(a->correct[i]/a->numReads > b->correct[i]/b->numReads) {
				return 1;
			}
		}
	}

	return 0;
}

void AccuracyProfileCopy(AccuracyProfile *dest, AccuracyProfile *src) 
{
	int i;
	if(dest->numReads > 0) {
		AccuracyProfileFree(dest);
	}
	AccuracyProfileAllocate(dest, src->numSNPs, src->numColorErrors);
	for(i=0;i<dest->length;i++) {
		dest->correct[i] = src->correct[i];
	}
}

void AccuracyProfileUpdate(IndexSet *set, 
		AccuracyProfile *p,
		int readLength,
		int numEventsToSample,
		int space,
		int maxNumMismatches,
		int maxNumColorErrors)
{
	int i, j, ctr;
	/* Get the estimated accuracy profile for the index set */

	/* Allocate memory for the profile */
	AccuracyProfileAllocate(p, maxNumMismatches, maxNumColorErrors);
	p->numReads = numEventsToSample;

	/* Go through */
	for(i=0,ctr=0;i<p->numColorErrors;i++) { /* color errors are prioritized */
		for(j=0;j<p->numSNPs;j++) { /* SNPs are secondary */
			p->correct[ctr] = GetNumCorrect(set,
					readLength,
					numEventsToSample,
					j,
					i,
					NoIndelType,
					0,
					space);
			ctr++;
		}
	}
}

void AccuracyProfileAllocate(AccuracyProfile *a,
		int numSNPs,
		int numColorErrors)
{
	char *FnName="AccuracyProfileAllocate";
	AccuracyProfileInitialize(a);
	a->numSNPs = numSNPs;
	a->numColorErrors = numColorErrors;
	a->length = (a->numSNPs + 1)*(a->numColorErrors + 1);
	a->correct = malloc(sizeof(int32_t)*a->length);
	if(NULL == a->correct) {
		PrintError(FnName,
				"a->correct",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
}

void AccuracyProfileInitialize(AccuracyProfile *a) 
{
	a->numReads = 0;
	a->correct = NULL;
	a->length = 0;
	a->numSNPs = -1;
	a->numColorErrors = -1;
}

void AccuracyProfileFree(AccuracyProfile *a)
{
	free(a->correct);
	AccuracyProfileInitialize(a);
}

void ReadSplit(Read *curRead,
		Read *r1,
		Read *r2,
		int breakpoint,
		int insertionLength)
{
	int i, ctr;

	/* Read 1 */
	if(breakpoint > 0) {
		ReadAllocate(r1, breakpoint);
	}
	for(i=0;i<breakpoint;i++) {
		assert(i < r1->length);
		r1->profile[i] = curRead->profile[i];
	}
	/* Read 2 */
	if(curRead->length - breakpoint - insertionLength > 0) {
		ReadAllocate(r2, curRead->length - breakpoint - insertionLength);
	}
	for(i=breakpoint+insertionLength, ctr=0;
			i<curRead->length;
			i++,ctr++) {
		assert(ctr < r2->length);
		r2->profile[ctr] = curRead->profile[i];
	}
}

void ReadGetRandom(Read *r, 
		int readLength,
		int numSNPs,
		int numColorErrors,
		int space)
{
	char *FnName="ReadGetRandom";
	int i;
	int numSNPsLeft = numSNPs;
	int numColorErrorsLeft = numColorErrors;
	int index;
	char original;
	char *read=NULL;
	char *originalColorSpace = NULL;
	int tmpReadLength = 0;

	assert(numSNPs <= readLength);
	assert(numColorErrors <= readLength);

	/* Allocate memory for a read */
	ReadAllocate(r, readLength);

	/* Initialize to no SNPs or color errors */
	for(i=0;i<r->length;i++) {
		r->profile[i] = 0;
	}

	assert(space == 1 || numColorErrors == 0);
	if(space == 0) {
		/* Insert random SNPS */
		while(numSNPsLeft > 0) {
			/* Pick a position to convert */
			index = rand()%(r->length);

			if(r->profile[index] == 0) {
				r->profile[index] = 1;
				numSNPsLeft--;
			}
		}
	}
	else {

		read = malloc(sizeof(char)*(readLength+1));
		if(NULL == read) {
			PrintError(FnName,
					"read",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		originalColorSpace = malloc(sizeof(char)*(readLength+1));
		if(NULL == originalColorSpace) {
			PrintError(FnName,
					"originalColorSpace",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}

		/* Get a random NT read */
		for(i=0;i<readLength;i++) {
			read[i] = DNA[rand()%4];
		}
		read[readLength]='\0';

		/* Get the color space of the original read */
		strcpy(originalColorSpace, read);
		tmpReadLength = readLength;
		ConvertReadToColorSpace(&originalColorSpace,
				&tmpReadLength);

		/* Insert random SNPs */
		while(numSNPsLeft > 0) {
			/* Pick a position to convert */
			index = rand()%(r->length);

			if(r->profile[index] == 0) {
				r->profile[index] = 1;
				numSNPsLeft--;
				/* Modify base to a new base */
				for(original = read[index];
						original == read[index];
						read[index] = DNA[rand()%4]) {
				}
			}
		}
		/* Convert to color space */
		tmpReadLength = readLength;
		ConvertReadToColorSpace(&read,
				&tmpReadLength);
		/* Insert color errors */
		while(numColorErrorsLeft > 0) {
			/* Pick a position to convert */
			index = rand()%(r->length);

			if(r->profile[index] != 2) {
				r->profile[index] = 2;
				numColorErrorsLeft--;
				/* Modify base to a new color */
				for(original = read[index];
						original == read[index];
						read[index] = Colors[rand()%4]) {
				}
			}
		}
		/* Compare the two profiles to get an end profile */
		for(i=0;i<r->length;i++) {
			if(originalColorSpace[i+1] == read[i+1]) {
				r->profile[i] = 0;
			}
			else {
				r->profile[i] = 1;
			}
		}

		free(read);
		read = NULL;
		free(originalColorSpace);
		originalColorSpace = NULL;
	}
}

void ReadInitialize(Read *r)
{
	r->length = 0;
	r->profile = NULL;
}

void ReadAllocate(Read *r, 
		int readLength)
{
	char *FnName = "ReadAllocate";
	r->length = readLength;
	r->profile = malloc(sizeof(int32_t)*r->length);
	if(NULL == r->profile) {
		PrintError(FnName,
				"r->profile",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
}

void ReadFree(Read *r)
{
	free(r->profile);
	ReadInitialize(r);
}

void ReadPrint(Read *r, FILE *fp) 
{
	int i;
	for(i=0;i<r->length;i++) {
		fprintf(fp, "%1d", r->profile[i]);
	}
	fprintf(fp, "\n");
}

void PrintUsage()
{
	fprintf(stderr, "Usage: test.indexes [OPTIONS]...\n");
	fprintf(stderr, "******************************* Algorithm Options (no defaults) *******************************\n");
	fprintf(stderr, "\t-a\tINT\talgorithm\n\t\t\t\t0: search for indexes\n\t\t\t\t1: evaluate indexes\n");
	fprintf(stderr, "\t-r\tINT\tread length (for all) \n");
	fprintf(stderr, "\t-S\tINT\tnumber of events to sample\n");
	fprintf(stderr, "\t-A\tINT\tspace 0: nucleotide space 1: color space\n");
	fprintf(stderr, "******************************* Search Options (for -a 0) *************************************\n");
	fprintf(stderr, "\t-s\tINT\tnumber of indexes to sample\n");
	fprintf(stderr, "\t-l\tINT\tkey size\n");
	fprintf(stderr, "\t-w\tINT\tmaximum key width\n");
	fprintf(stderr, "\t-n\tINT\tmaximum index set size\n");
	fprintf(stderr, "\t-t\tINT\taccuracy percent threshold (0-100)\n");
	fprintf(stderr, "******************************* Evaluate Options (for -a 1) ***********************************\n");
	fprintf(stderr, "\t-f\tSTRING\tinput file name\n");
	fprintf(stderr, "\t-I\tINT\tmaximum insertion length (-a 1)\n");
	fprintf(stderr, "******************************* Event Options (default =0 ) ***********************************\n");
	fprintf(stderr, "\t-M\tINT\tmaximum number of mismatches\n");
	fprintf(stderr, "\t-e\tINT\tnumber of color errors (-A 1)\n");
	fprintf(stderr, "******************************* Miscellaneous Options  ****************************************\n");
	fprintf(stderr, "\t-p\tNULL\tprints the program parameters\n");
	fprintf(stderr, "\t-h\tNULL\tprints this message\n");
}

void PrintProgramParameters(arguments *args)
{
	/* Print program parameters */
	fprintf(stdout, "%s", BREAK_LINE);
	fprintf(stdout, "Printing program parameters:\n");
	fprintf(stdout, "algorithm:\t\t\t%d\t[%s]\n", args->algorithm, Algorithm[args->algorithm]); 
	fprintf(stdout, "read length:\t\t\t%d\n", args->readLength);
	fprintf(stdout, "number of events to sample:\t%d\n", args->numEventsToSample);
	fprintf(stdout, "space:\t\t\t\t%d\n", args->space);
	fprintf(stdout, "number of indexes to sample:\t%d\n", args->numIndexesToSample);
	fprintf(stdout, "key size:\t\t\t%d\n", args->keySize);
	fprintf(stdout, "key width:\t\t\t%d\n", args->maxKeyWidth);
	fprintf(stdout, "max index set size:\t\t%d\n", args->maxIndexSetSize);
	fprintf(stdout, "accuracy percent threshold:\t%d\n", args->accuracyThreshold);
	fprintf(stdout, "input file name:\t\t%s\n", args->inputFileName);
	fprintf(stdout, "maximum insertion length:\t%d\n", args->maxInsertionLength);
	fprintf(stdout, "maximum number of mismatches:\t%d\n", args->maxNumMismatches);
	fprintf(stdout, "maximum number of color errors:\t%d\n", args->maxNumColorErrors);
	fprintf(stdout, "%s", BREAK_LINE);
}

void AssignDefaultValues(arguments *args) 
{
	args->algorithm=0;
	strcpy(args->inputFileName, "\0");
	args->readLength=0;
	args->numEventsToSample=0;
	args->numIndexesToSample=0;
	args->keySize=0;
	args->maxKeyWidth=0;
	args->maxIndexSetSize=0;
	args->accuracyThreshold=0;
	args->space=0;
	args->maxNumMismatches=0;
	args->maxInsertionLength=0;
	args->maxNumColorErrors=0;
}

void ValidateArguments(arguments *args)
{
	char *FnName="ValidateArguments";

	if(args->algorithm < 0 || args->algorithm > 1) {
		PrintError(FnName, "Command line argument", "algorithm", Exit, OutOfRange);
	}
	if(args->readLength <= 0) {
		PrintError(FnName, "Command line argument", "readLength", Exit, OutOfRange);
	}
	if(args->numEventsToSample <= 0) {
		PrintError(FnName, "Command line argument", "numEventsToSample", Exit, OutOfRange);
	}
	if(args->numIndexesToSample < 0 ||
			(args->algorithm == 0 && args->numIndexesToSample <= 0)) {
		PrintError(FnName, "Command line argument", "numIndexesToSample", Exit, OutOfRange);
	}
	if(args->algorithm == 0) {
		if(args->keySize <= 0) {
			PrintError(FnName, "Command line argument", "keySize", Exit, OutOfRange);
		}
		if(args->maxKeyWidth <= 0) {
			PrintError(FnName, "Command line argument", "maxKeyWidth", Exit, OutOfRange);
		}
		if(args->maxIndexSetSize <= 0) {
			PrintError(FnName, "Command line argument", "maxIndexSetSize", Exit, OutOfRange);
		}
		if(args->accuracyThreshold < 0) {
			PrintError(FnName, "Command line argument", "accuracyThreshold", Exit, OutOfRange);
		}
	}
	if(args->space < 0 || 1 < args->space) {
		PrintError(FnName, "Command line argument", "space", Exit, OutOfRange);
	}
	if(args->maxNumMismatches < 0) {
		PrintError(FnName, "Command line argument", "maxNumMismatches", Exit, OutOfRange);
	}
	if(args->maxInsertionLength < 0) {
		PrintError(FnName, "Command line argument", "maxInsertionLength", Exit, OutOfRange);
	}
	if(args->space == 1) {
		if(args->maxNumColorErrors < 0) {
			PrintError(FnName, "Command line argument", "maxNumColorErrors", Exit, OutOfRange);
		}
	}
	else {
		if(args->maxNumColorErrors > 0) {
			PrintError(FnName, "Command line argument", "maxNumColorErrors", Exit, OutOfRange);
		}
	}
}

void ParseCommandLineArguments(int argc, char *argv[], arguments *args) 
{
	int i;
	if(argc==1) {
		PrintUsage();
		exit(1);
	}
	for(i=1;i<argc;i+=2) {
		if(argv[i][0] != '-' ||
				strlen(argv[i]) != 2) {
			fprintf(stderr, "*** Error.  Could not understand command line option %s.  Terminating! ***\n",
					argv[i]);
			exit(1);
		}
		switch(argv[i][1]) {
			case 'a':
				args->algorithm = atoi(argv[i+1]);
				break;
			case 'A':
				args->space = atoi(argv[i+1]);
				break;
			case 'e':
				args->maxNumColorErrors = atoi(argv[i+1]);
				break;
			case 'f':
				strcpy(args->inputFileName, argv[i+1]);
				break;
			case 'I':
				args->maxInsertionLength = atoi(argv[i+1]);
				break;
			case 'h':
				PrintUsage();
				exit(1);
				break;
			case 'l':
				args->keySize = atoi(argv[i+1]);
				break;
			case 'M':
				args->maxNumMismatches = atoi(argv[i+1]);
				break;
			case 'n':
				args->maxIndexSetSize = atoi(argv[i+1]);
				break;
			case 'p':
				args->algorithm = ProgramParameters;
				break;
			case 'r':
				args->readLength = atoi(argv[i+1]);
				break;
			case 's':
				args->numIndexesToSample = atoi(argv[i+1]);
				break;
			case 'S':
				args->numEventsToSample = atoi(argv[i+1]);
				break;
			case 't':
				args->accuracyThreshold = atoi(argv[i+1]);
				break;
			case 'w':
				args->maxKeyWidth = atoi(argv[i+1]);
				break;
			default:
				fprintf(stderr, "*** Error.  Could not understand command line option %s.  Terminating! ***\n",
						argv[i]);
				exit(1);
				break;

		}
	}
}

int main(int argc, char *argv[])
{
	/* Command line arguments */
	arguments args;

	/* Assign default values */
	AssignDefaultValues(&args);

	/* Parse command line arguments */
	ParseCommandLineArguments(argc, argv, &args);

	/* Validate command line arguments */
	ValidateArguments(&args);

	/* Print program parameters */
	PrintProgramParameters(&args);

	switch(args.algorithm) {
		case SearchForIndexes:
			RunSearchForIndexes(args.readLength,
					args.numEventsToSample,
					args.numIndexesToSample,
					args.keySize,
					args.maxKeyWidth,
					args.maxIndexSetSize,
					args.accuracyThreshold,
					args.space,
					args.maxNumMismatches,
					args.maxNumColorErrors);
			break;
		case EvaluateIndexes:
			RunEvaluateIndexes(args.inputFileName,
					args.readLength,
					args.numEventsToSample,
					args.space,
					args.maxNumMismatches,
					args.maxInsertionLength,
					args.maxNumColorErrors);
			break;
		case ProgramParameters:
			/* Do nothing */
			break;
		default:
			fprintf(stderr, "Error.  Could not understand program mode [%d].  Terminating!\n", args.algorithm);
			exit(1);

	}
	return 0;
}
