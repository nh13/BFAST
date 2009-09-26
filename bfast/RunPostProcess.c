#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <zlib.h>
#include <limits.h>
#include "BLibDefinitions.h"
#include "BError.h"
#include "AlignedRead.h"
#include "AlignedEnd.h"
#include "AlignedReadConvert.h"
#include "ScoringMatrix.h"
#include "RunPostProcess.h"

/* TODO */
void ReadInputFilterAndOutput(RGBinary *rg,
		char *inputFileName,
		int algorithm,
		int queueLength,
		int outputFormat,
		char *outputID,
		char *unmappedFileName)
{
	char *FnName="ReadInputFilterAndOutput";
	gzFile fp=NULL;
	int32_t i, j;
	int64_t counter, foundType, numUnmapped, numReported;
	AlignedRead *aBuffer;
	int32_t aBufferLength=0;
	int32_t numRead, aBufferIndex;
	gzFile fpReportedGZ=NULL;
	FILE *fpReported=NULL;
	gzFile fpUnmapped=NULL;

	/* Open the input file */
	if(NULL == inputFileName) {
		if(!(fp=gzdopen(fileno(stdin), "rb"))) {
			PrintError(FnName, "stdin", "Could not open stdin for reading", Exit, OpenFileError);
		}
	}
	else {
		if(!(fp=gzopen(inputFileName, "rb"))) {
			PrintError(FnName, inputFileName, "Could not open inputFileName for reading", Exit, OpenFileError);
		}
	}

	/* Open output files, if necessary */
	if(BAF == outputFormat) {
		if(!(fpReportedGZ=gzdopen(fileno(stdout), "wb"))) {
			PrintError(FnName, "stdout", "Could not open stdout for writing", Exit, OpenFileError);
		}
	}
	else {
		if(!(fpReported=fdopen(fileno(stdout), "wb"))) {
			PrintError(FnName, "stdout", "Could not open stdout for writing", Exit, OpenFileError);
		}
	}
	if(NULL != unmappedFileName) {
		if(!(fpUnmapped=gzopen(unmappedFileName, "wb"))) {
			PrintError(FnName, unmappedFileName, "Could not open unmappedFileName for writing", Exit, OpenFileError);
		}
	}

	AlignedReadConvertPrintHeader(fpReported, rg, outputFormat);

	aBufferLength=queueLength;
	aBuffer=malloc(sizeof(AlignedRead)*aBufferLength);
	if(NULL == aBuffer) {
		PrintError(FnName, "aBuffer", "Could not allocate memory", Exit, MallocMemory);
	}

	/* Go through each read */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Processing reads, currently on:\n0");
	}
	counter = numReported = numUnmapped = numRead = 0;

	while(0 != (numRead = GetAlignedReads(fp, aBuffer, aBufferLength))) {

		for(aBufferIndex=0;aBufferIndex<numRead;aBufferIndex++) {

			if(VERBOSE >= 0 && counter%ALIGNENTRIES_READ_ROTATE_NUM==0) {
				fprintf(stderr, "\r%lld",
						(long long int)counter);
			}

			/* Filter */
			foundType=FilterAlignedRead(&aBuffer[aBufferIndex],
					algorithm);

			if(NoneFound == foundType) {
				/* Print to Not Reported file */
				if(NULL != unmappedFileName) {
					AlignedReadPrint(&aBuffer[aBufferIndex], fpUnmapped);
				}
				numUnmapped++;

				/* Free the alignments for output */
				for(i=0;i<aBuffer[aBufferIndex].numEnds;i++) {
					for(j=0;j<aBuffer[aBufferIndex].ends[i].numEntries;j++) {
						AlignedEntryFree(&aBuffer[aBufferIndex].ends[i].entries[j]);
					}
					aBuffer[aBufferIndex].ends[i].numEntries=0;
				}
			}
			else {
				numReported++;
			}

			/* Increment counter */
			counter++;
		}

		if(VERBOSE >= 0) {
			fprintf(stderr, "\r%lld",
					(long long int)counter);
		}

		/* Print to Output file */
		for(aBufferIndex=0;aBufferIndex<numRead;aBufferIndex++) {
			AlignedReadConvertPrintOutputFormat(&aBuffer[aBufferIndex], rg, fpReported, fpReportedGZ, (NULL == outputID) ? "" : outputID, outputFormat, BinaryOutput);

			/* Free memory */
			AlignedReadFree(&aBuffer[aBufferIndex]);
		}
	}
	if(VERBOSE>=0) {
		fprintf(stderr, "\r%lld\n",
				(long long int)counter);
	}

	/* Close output files, if necessary */
	if(BAF == outputFormat) {
		gzclose(fpReportedGZ);
	}
	else {
		fclose(fpReported);
	}
	if(NULL != unmappedFileName) {
		gzclose(fpUnmapped);
	}
	/* Close the input file */
	gzclose(fp);
	free(aBuffer);

	if(VERBOSE>=0) {
		fprintf(stderr, "Found %lld mappings after filtering.\n%lld reads were left unmapped\n",
				(long long int)numReported,
				(long long int)numUnmapped);
	}
}

int32_t GetAlignedReads(gzFile fp, AlignedRead *aBuffer, int32_t maxToRead)
{
	int32_t numRead=0;
	while(numRead < maxToRead) {
		AlignedReadInitialize(&aBuffer[numRead]);
		if(EOF == AlignedReadRead(&aBuffer[numRead], fp)) {
			return numRead;
		}
		numRead++;
	}
	return numRead;
}

/* TODO */
int FilterAlignedRead(AlignedRead *a,
		int algorithm) 
{
	char *FnName="FilterAlignedRead";
	int foundType;
	int32_t *foundTypes=NULL;
	AlignedRead tmpA;
	int32_t i, j, ctr;
	int32_t best, bestIndex, numBest;

	AlignedReadInitialize(&tmpA);

	/* We should only modify "a" if it is going to be reported */ 
	/* Copy in case we do not find anything to report */
	AlignedReadCopy(&tmpA, a);


	foundType=NoneFound;
	foundTypes=malloc(sizeof(int32_t)*tmpA.numEnds);
	if(NULL == foundTypes) {
		PrintError(FnName, "foundTypes", "Could not allocate memory", Exit, MallocMemory);
	}
	for(i=0;i<tmpA.numEnds;i++) {
		foundTypes[i]=NoneFound;
	}

	/* Pick alignment for each end individually (is this a good idea?) */
	for(i=0;i<tmpA.numEnds;i++) {
		/* Choose each end */
		switch(algorithm) {
			case NoFiltering:
			case AllNotFiltered:
				foundTypes[i] = (0<tmpA.ends[i].numEntries)?Found:NoneFound;
				break;
			case Unique:
				foundTypes[i]=(1==tmpA.ends[i].numEntries)?Found:NoneFound;
				break;
			case BestScore:
			case BestScoreAll:
				best = INT_MIN;
				bestIndex = -1;
				numBest = 0;
				for(j=0;j<tmpA.ends[i].numEntries;j++) {
					if(best < tmpA.ends[i].entries[j].score) {
						best = tmpA.ends[i].entries[j].score;
						bestIndex = j;
						numBest = 1;
					}
					else if(best == tmpA.ends[i].entries[j].score) {
						numBest++;
					}
				}
				if(BestScore == algorithm &&
						1 == numBest) {
					foundTypes[i] = Found;
					/* Copy to front */
					AlignedEntryCopy(&tmpA.ends[i].entries[0], 
							&tmpA.ends[i].entries[bestIndex]);
					AlignedEndReallocate(&tmpA.ends[i], 1);
				}
				else if(BestScoreAll == algorithm) {
					foundTypes[i] = Found;
					ctr=0;
					for(j=0;j<tmpA.ends[i].numEntries;j++) {
						if(tmpA.ends[i].entries[j].score == best) {
							if(ctr != j) {
								AlignedEntryCopy(&tmpA.ends[i].entries[ctr], 
										&tmpA.ends[i].entries[j]);
							}
							ctr++;
						}
					}
					assert(ctr == numBest);
					AlignedEndReallocate(&tmpA.ends[i], numBest);
				}
				break;
			default:
				PrintError(FnName, "algorithm", "Could not understand algorithm", Exit, OutOfRange);
				break;
		}
		/* Free if not found */
		if(NoneFound == foundTypes[i]) {
			AlignedEndReallocate(&tmpA.ends[i],
					0);
		}
	}

	if(1 == tmpA.numEnds) {
		foundType=foundTypes[0];
	}
	else {
		/* Call found if at least one has been found */
		foundType=NoneFound;
		for(i=0;NoneFound==foundType && i<tmpA.numEnds;i++) {
			if(Found == foundTypes[i]) {
				foundType=Found;
				break;
			}
		}
	}

	/* If we found, then copy back */
	if(NoneFound != foundType) {
		AlignedReadFree(a);
		AlignedReadCopy(a, &tmpA);
	}
	AlignedReadFree(&tmpA);
	free(foundTypes);

	return foundType;
}
