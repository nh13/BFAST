#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <float.h>

#include "BLibDefinitions.h"
#include "ScoringMatrix.h"
#include "QS.h"
#include "BError.h"
#include "BLib.h"
#include "AlignedEntry.h"
#include "AlignedEnd.h"

/* 
   4.343*log(1:256)
   */
static int32_t logN[256] = {
	0.000000, 3.010338, 4.771273, 6.020676, 6.989789, 7.781611,
	8.451088, 9.031015, 9.542546, 10.000127, 10.414059, 10.791950, 11.139575,
	11.461426, 11.761062, 12.041353, 12.304646, 12.552885, 12.787698, 13.010465,
	13.222361, 13.424397, 13.617451, 13.802288, 13.979578, 14.149913, 14.313820,
	14.471764, 14.624166, 14.771400, 14.913806, 15.051691, 15.185332, 15.314984,
	15.440877, 15.563223, 15.682216, 15.798037, 15.910848, 16.020803, 16.128043,
	16.232699, 16.334892, 16.434736, 16.532335, 16.627790, 16.721191, 16.812626,
	16.902176, 16.989916, 17.075919, 17.160251, 17.242978, 17.324158, 17.403848,
	17.482102, 17.558972, 17.634504, 17.708745, 17.781738, 17.853525, 17.924145,
	17.993634, 18.062029, 18.129364, 18.195671, 18.260980, 18.325322, 18.388725,
	18.451215, 18.512819, 18.573561, 18.633465, 18.692555, 18.750851, 18.808375,
	18.865147, 18.921186, 18.976512, 19.031142, 19.085093, 19.138382, 19.191025,
	19.243037, 19.294434, 19.345230, 19.395439, 19.445074, 19.494148, 19.542673,
	19.590663, 19.638128, 19.685080, 19.731529, 19.777487, 19.822964, 19.867970,
	19.912514, 19.956606, 20.000254, 20.043468, 20.086257, 20.128628, 20.170590,
	20.212150, 20.253316, 20.294096, 20.334496, 20.374524, 20.414186, 20.453490,
	20.492441, 20.531045, 20.569310, 20.607240, 20.644842, 20.682121, 20.719083,
	20.755733, 20.792077, 20.828118, 20.863863, 20.899317, 20.934483, 20.969367,
	21.003972, 21.038305, 21.072367, 21.106165, 21.139702, 21.172982, 21.206009,
	21.238786, 21.271318, 21.303608, 21.335660, 21.367477, 21.399063, 21.430420,
	21.461553, 21.492464, 21.523157, 21.553634, 21.583899, 21.613955, 21.643804,
	21.673449, 21.702893, 21.732139, 21.761189, 21.790046, 21.818713, 21.847192,
	21.875485, 21.903595, 21.931525, 21.959276, 21.986850, 22.014251, 22.041480,
	22.068539, 22.095431, 22.122157, 22.148720, 22.175121, 22.201363, 22.227447,
	22.253376, 22.279150, 22.304773, 22.330245, 22.355569, 22.380745, 22.405777,
	22.430665, 22.455412, 22.480018, 22.504486, 22.528817, 22.553012, 22.577073,
	22.601001, 22.624798, 22.648466, 22.672005, 22.695418, 22.718705, 22.741867,
	22.764907, 22.787826, 22.810623, 22.833302, 22.855863, 22.878308, 22.900637,
	22.922852, 22.944954, 22.966944, 22.988823, 23.010592, 23.032253, 23.053807,
	23.075254, 23.096595, 23.117832, 23.138966, 23.159998, 23.180928, 23.201758,
	23.222488, 23.243120, 23.263654, 23.284092, 23.304434, 23.324681, 23.344834,
	23.364894, 23.384862, 23.404739, 23.424524, 23.444221, 23.463828, 23.483347,
	23.502779, 23.522124, 23.541384, 23.560558, 23.579648, 23.598655, 23.617578,
	23.636420, 23.655180, 23.673860, 23.692460, 23.710980, 23.729422, 23.747785,
	23.766072, 23.784281, 23.802415, 23.820473, 23.838457, 23.856366, 23.874202,
	23.891964, 23.909655, 23.927274, 23.944821, 23.962298, 23.979705, 23.997042,
	24.014311, 24.031511, 24.048643, 24.065708, 24.082706};

/* TODO */
int32_t AlignedEndPrint(AlignedEnd *a,
		FILE *outputFP,
		int32_t space,
		int32_t binaryOutput)
{
	int32_t i;
	assert(NULL != a->read);
	if(binaryOutput == TextOutput) {
		if(fprintf(outputFP, "%s\t%s\t%d\n",
					a->read,
					a->qual,
					a->numEntries) < 0) {
			return EOF;
		}
	}
	else {
		if(fwrite(&a->readLength, sizeof(int32_t), 1, outputFP) != 1 ||
				fwrite(&a->qualLength, sizeof(int32_t), 1, outputFP) != 1 ||
				fwrite(a->read, sizeof(char), a->readLength, outputFP) != a->readLength ||
				fwrite(a->qual, sizeof(char), a->qualLength, outputFP) != a->qualLength ||
				fwrite(&a->numEntries, sizeof(int32_t), 1, outputFP) != 1) {
			return EOF;
		}
	}

	for(i=0;i<a->numEntries;i++) {
		if(EOF == AlignedEntryPrint(&a->entries[i],
					outputFP,
					space,
					binaryOutput)) {
			return EOF;
		}
	}

	return 1;
}

/* TODO */
int32_t AlignedEndRead(AlignedEnd *a,
		FILE *inputFP,
		int32_t space,
		int32_t binaryInput)
{
	char *FnName = "AlignedEndRead";
	char read[SEQUENCE_LENGTH]="\0";
	char qual[SEQUENCE_LENGTH]="\0";
	int32_t i;

	if(binaryInput == TextInput) {
		if(fscanf(inputFP, "%s\t%s\t%d",
					read,
					qual,
					&a->numEntries) < 0) {
			return EOF;
		}

		a->readLength = strlen(read);
		a->qualLength = strlen(qual);

		/* Allocate memory for the alignment */
		if(a->read == NULL) {
			a->read = malloc(sizeof(char)*(1+a->readLength));
			if(NULL == a->read) {
				PrintError(FnName,
						"a->read",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}
		if(a->qual == NULL) {
			a->qual = malloc(sizeof(char)*(1+a->qualLength));
			if(NULL == a->qual) {
				PrintError(FnName,
						"a->qual",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}
		/* Copy over */
		strcpy(a->read, read);
		strcpy(a->qual, qual);
	}
	else {
		if(fread(&a->readLength, sizeof(int32_t), 1, inputFP) != 1 ||
				fread(&a->qualLength, sizeof(int32_t), 1, inputFP) != 1) {
			return EOF;
		}
		/* Allocate memory for the alignment */
		if(a->read == NULL) {
			a->read = malloc(sizeof(char)*(1+a->readLength));
			if(NULL == a->read) {
				PrintError(FnName,
						"a->read",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}
		if(a->qual == NULL) {
			a->qual = malloc(sizeof(char)*(1+a->qualLength));
			if(NULL == a->qual) {
				PrintError(FnName,
						"a->qual",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
		}

		if(fread(a->read, sizeof(char), a->readLength, inputFP) != a->readLength ||
				fread(a->qual, sizeof(char), a->qualLength, inputFP) != a->qualLength ||
				fread(&a->numEntries, sizeof(int32_t), 1, inputFP) != 1) {
			PrintError(FnName,
					"a->reads, a->qual, and a->numEntries",
					"Could not read from file",
					Exit,
					ReadFileError);
		}
		/* Add the null terminator to strings */
		a->read[a->readLength]='\0';
		a->qual[a->qualLength]='\0';
	}

	/* Allocate room for the the entries */
	AlignedEndAllocate(a,
			a->numEntries);

	for(i=0;i<a->numEntries;i++) {
		AlignedEntryInitialize(&a->entries[i]);
		if(EOF == AlignedEntryRead(&a->entries[i],
					inputFP,
					space,
					binaryInput)) {
			PrintError(FnName,
					"a->entries[i]",
					"Could not read from file",
					Exit,
					ReadFileError);
		}
	}

	return 1;
}

/* TODO */
int32_t AlignedEndRemoveDuplicates(AlignedEnd *end,
		int32_t sortOrder)
{
	/*
	char *FnName="AlignedEndRemoveDuplicates";
	*/
	int32_t i, prevIndex;

	if(end->numEntries > 1) {
		/* Sort the entries */
		AlignedEndQuickSort(end, sortOrder, 0);
		/*
		AlignedEndMergeSort(end, sortOrder, 0);
		 */

		/* Remove duplicates */
		prevIndex=0;
		for(i=1;i<end->numEntries;i++) {
			if(AlignedEntryCompareAtIndex(end->entries, prevIndex, end->entries, i, sortOrder)==0) {
				/* Do nothing */
			}
			else {
				/* Increment prevIndex */
				prevIndex++;
				/* Copy to prevIndex (incremented) */
				AlignedEntryCopyAtIndex(end->entries, prevIndex, end->entries, i);
			}
		}

		/* Reallocate */
		AlignedEndReallocate(end, prevIndex+1);
	}
	return end->numEntries;
}

/* TODO */
/* Log-n space */
/* Do not use, since it is buggy and has not been updated lately */  
void AlignedEndQuickSort(AlignedEnd *a,
		int32_t sortOrder,
		int32_t showPercentComplete)
{

	double curPercent = 0.0;
	AlignedEntryQuickSort(&a->entries,
			0,
			a->numEntries-1,
			sortOrder,
			showPercentComplete,
			&curPercent,
			a->numEntries);
}

/* TODO */
/* O(n) space, but really double */
void AlignedEndMergeSort(AlignedEnd *a,
		int32_t sortOrder,
		int32_t showPercentComplete)
{
	double curPercent = 0.0;
	AlignedEntryMergeSort(&a->entries,
			0,
			a->numEntries-1,
			sortOrder,
			showPercentComplete,
			&curPercent,
			a->numEntries);
}

/* TODO */
int32_t AlignedEndCompare(AlignedEnd *a,
		AlignedEnd *b, 
		int32_t sortOrder)
{
	assert(1 == a->numEntries);
	assert(1 == b->numEntries);

	return AlignedEntryCompare(&a->entries[0],
			&b->entries[0],
			sortOrder);
}

/* TODO */
void AlignedEndCopyAtIndex(AlignedEnd *dest, int32_t destIndex, AlignedEnd *src, int32_t srcIndex)
{
	if(dest != src || srcIndex != destIndex) {
		AlignedEndCopy(&(dest[destIndex]), &(src[srcIndex]));
	}
}

/* TODO */
void AlignedEndCopy(AlignedEnd *dest, AlignedEnd *src)
{
	char *FnName = "AlignedEndCopy";
	int32_t i;
	if(src != dest) {
		assert(src->read != NULL);
		/* read */
		dest->readLength = src->readLength;
		dest->read = realloc(dest->read, sizeof(char)*(src->readLength+1));
		if(NULL == dest->read) {
			PrintError(FnName,
					"dest->read",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		assert(src->read != NULL);
		strcpy(dest->read, src->read);
		/* qual */
		dest->qualLength = src->qualLength;
		dest->qual = realloc(dest->qual, sizeof(char)*(src->qualLength+1));
		if(NULL == dest->qual) {
			PrintError(FnName,
					"dest->qual",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		assert(src->qual != NULL);
		strcpy(dest->qual, src->qual);
		/* Reallocate */
		AlignedEndReallocate(dest,
				src->numEntries);
		/* Copy entries */
		for(i=0;i<dest->numEntries;i++) {
			AlignedEntryCopy(&dest->entries[i], 
					&src->entries[i]);
		}
	}
}

void AlignedEndAllocate(AlignedEnd *a,
		int32_t numEntries)
{
	char *FnName="AlignedEndAllocate";
	int32_t i;

	/* Allocate */
	a->entries = malloc(sizeof(AlignedEntry)*numEntries);
	if(NULL == a->entries && 0 < numEntries) {
		PrintError(FnName,
				"a->entries",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}

	a->numEntries = numEntries;
	for(i=0;i<a->numEntries;i++) {
		AlignedEntryInitialize(&a->entries[i]);
	}

}

void AlignedEndReallocate(AlignedEnd *a,
		int32_t numEntries)
{
	char *FnName="AlignedEndReallocate";
	int32_t i;

	if(numEntries < a->numEntries) {
		for(i=numEntries;i<a->numEntries;i++) {
			AlignedEntryFree(&a->entries[i]);
		}
	}

	/* Reallocate */
	a->entries = realloc(a->entries, sizeof(AlignedEntry)*numEntries);
	if(NULL == a->entries && 0 < numEntries) {
		PrintError(FnName,
				"a->entries",
				"Could not reallocate memory",
				Exit,
				ReallocMemory);
	}

	for(i=a->numEntries;i<numEntries;i++) {
		AlignedEntryInitialize(&a->entries[i]);
	}

	a->numEntries = numEntries;

}

void AlignedEndFree(AlignedEnd *a)
{
	int32_t i;

	free(a->read);
	free(a->qual);
	for(i=0;i<a->numEntries;i++) {
		AlignedEntryFree(&a->entries[i]);
	}
	AlignedEndInitialize(a);
}

void AlignedEndInitialize(AlignedEnd *a) 
{
	a->read=NULL;
	a->readLength=0;
	a->qual=NULL;
	a->qualLength=0;
	a->numEntries=0;
	a->entries=NULL;
}

void AlignedEndAssignMappingQualities(AlignedEnd *end,
		ScoringMatrix *sm,
		QS *qs,
		RGIndexAccuracyMismatchProfile *profile,
		int32_t space)
{
	char *FnName="AlignedEndAssignMappingQualities";
	/* Imitation is the sincerest form of flattery... 
	 * This follows the maq-0.7.1 mapping qualities approximations.
	 *
	 * Assumptions:
	 * 		- Given an error probabilitiy e, (1-e) is approx 0.
	 * 		- Illumina: mismatches are assumed to be errors etc.
	 * 		- ABI: color errors are errors, mismatches are 
	 * 		(possible) variants etc.
	 *
	 * */
	int32_t i, j, ctr;
	int32_t scores[2] = {INT_MIN, INT_MIN}; /* Scores */ 
	int32_t num[2] = {0, 0}; /* # of each */
	int32_t maxQ[2] = {INT_MIN, INT_MIN}; /* max quality */
	int32_t curScore=0;
	int32_t *quals=NULL;
	int32_t mappingQ, altMappingQ;
	int32_t maxMM, curMM, diffMM;
	double minScore = 0.0;

	if(end->numEntries <= 0) {
		return;
	}

	quals = calloc(sizeof(int32_t), end->numEntries);
	if(NULL == quals) {
		PrintError(FnName,
				"quals",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Get best and second best */
	for(i=0;i<end->numEntries;i++) {
		quals[i] = 0;
		curScore = (int32_t)end->entries[i].score;
		if(scores[0] <= curScore ||
				(curScore < scores[0] && scores[1] <= curScore)) {
			/* Get sum of mismatch/color error qualities */
			for(j=ctr=0;j<end->entries[i].length;j++) {
				if(GAP != end->entries[i].reference[j] &&
						GAP != end->entries[i].read[j]) { /* Not an indel */
					if((NTSpace == space && end->entries[i].reference[j] != end->entries[i].read[j]) || /* Mismatch */
							(ColorSpace == space && end->entries[i].colorError[j] != GAP)) { /* Color error */
						/* Add quality score */
						quals[i] += CHAR2QUAL(end->qual[ctr]);
					}
					ctr++;
				}
			}

			/* Bin appropriately */
			if(curScore == scores[0]) {
				num[0]++;
				if(maxQ[0] < quals[i] ) {
					maxQ[0] = quals[i];
				}
			}
			else if(scores[0] < curScore) {
				num[0]=1;
				scores[0] = curScore;
				maxQ[0] = quals[i];
			}
			else if(curScore == scores[1]) {
				num[1]++;
				if(maxQ[1] < quals[i] ) {
					maxQ[1] = quals[i];
				}
			}
			else if(scores[1] < curScore) {
				num[1]++;
				scores[1] = curScore;
				maxQ[1] = quals[i];
			}
		}
	}

	assert(1 <= num[0]);

	for(i=0;i<end->numEntries;i++) {
		curScore = (int32_t)end->entries[i].score;
		mappingQ = altMappingQ = 0;
		if(curScore == scores[0] &&
				1 == num[0]) { /* Uniquely the best */

			/* Get primary mapping quality */
			if(0 == num[1]) { /* No second best exists */
				mappingQ = 99;
			}
			else {
				end->entries[i].mappingQuality = 99;
				/* Second best exists */
				/* quality_2 - quality_cur - int(4.343 * log(num_quality_2) + 0.5) */
				if(num[1] <= 256) {
					mappingQ = maxQ[1] - quals[i] - (int)(logN[num[1]] + 0.5);
				}
				else {
					mappingQ = maxQ[1] - quals[i] - (int)(4.343*log(num[1]) + 0.5);
				}
			}
			/* Get alternative mapping quality */
			minScore = DBL_MIN;
			assert(end->readLength <= profile->maxReadLength);
			maxMM = profile->maxMismatches[end->readLength];
			assert(0 <= maxMM);
			/* Get number of mismatches (NT space) or color errors (Color space)
			 * in the alignment */
			curMM=0;
			for(j=0;j<end->entries[i].length;j++) {
				if(NTSpace == space) {
					if(GAP != end->entries[i].read[j] &&
							GAP != end->entries[i].reference[j] &&
							end->entries[i].read[j] != end->entries[i].reference[j]) {
						curMM++;
					}
				}
				else {
					if(GAP != end->entries[i].colorError[j]) {
						curMM++;
					}
				}
			}
			/* Get minimum possible score with only mismatches or color errors given curMM */
			assert(0 <= curMM);
			if(NTSpace == space) {
				minScore = (end->readLength - maxMM)*sm->maxNTScore + curMM*sm->minNTScore;
			}
			else {
				/* Minus one due to the adaptor */
				minScore = (end->readLength - maxMM - 1)*sm->maxColorScore + curMM*sm->minColorScore;
			}
			/* Get difference between the number of mismatches found and the the number of mismatches
			 * BFAST will confidently allow */
			diffMM = maxMM + 1 - curMM; 
			if(0 <= diffMM) {
				altMappingQ = 2*maxMM + (diffMM*quals[i]/end->qualLength - 13) + QSGet(qs, (int32_t)(end->entries[i].score - minScore));
			}
			else {
				altMappingQ = 0;
			}

			/* Compare primary and alternative mapping qualities */
			if(0 < altMappingQ && altMappingQ < mappingQ) {
				mappingQ = altMappingQ;
			}
			if(mappingQ < 0) {
				mappingQ = 0;
			}
			else if(99 < mappingQ) {
				mappingQ = 99;
			}
		}
		/* Store mapping quality */
		end->entries[i].mappingQuality = mappingQ;
	}

	free(quals);
}
