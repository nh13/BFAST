#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <zlib.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include "BLibDefinitions.h"
#include "BLib.h"
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
		int pairedEndInfer,
		int numThreads,
		int queueLength,
		int outputFormat,
		char *outputID,
		char *readGroup,
		char *unmappedFileName)
{
	char *FnName="ReadInputFilterAndOutput";
	gzFile fp=NULL;
	int32_t i;
	int32_t numUnmapped=0, numReported=0;
	gzFile fpReportedGZ=NULL;
	FILE *fpReported=NULL;
	gzFile fpUnmapped=NULL;
	PEDBins bins;
	int32_t *mappedEndCounts=NULL;
	int32_t mappedEndCountsNumEnds=-1;
	int32_t pairedEndInferRescued=0;
	char *readGroupString=NULL;
	pthread_mutex_t inputFP_mutex, outputFP_mutex;
	pthread_t *threads=NULL;
	int errCode;
	void *status=NULL;
	PostProcessThreadData *data=NULL;

	if(NULL != readGroup) {
		readGroupString=ParseReadGroup(readGroup);
	}

	/* Get the PEDBins if necessary */
	if(1 == pairedEndInfer) {
		PEDBinsInitialize(&bins);
		pairedEndInfer = GetPEDBins(inputFileName, algorithm, queueLength, &bins);
		PEDBinsMakeIntoProbability(&bins);
	}

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

	AlignedReadConvertPrintHeader(fpReported, rg, outputFormat, readGroup);

	// Initialize mutex
	pthread_mutex_init(&outputFP_mutex, NULL);
	pthread_mutex_init(&inputFP_mutex, NULL);

	/* Allocate memory for threads */
	threads=malloc(sizeof(pthread_t)*numThreads);
	if(NULL==threads) {
		PrintError(FnName, "threads", "Could not allocate memory", Exit, MallocMemory);
	}
	/* Allocate memory to pass data to threads */
	data=malloc(sizeof(PostProcessThreadData)*numThreads);
	if(NULL==data) {
		PrintError(FnName, "data", "Could not allocate memory", Exit, MallocMemory);
	}

	/* Initialize thread data */
	for(i=0;i<numThreads;i++) {
		data[i].rg = rg;
		data[i].bins = &bins;
		data[i].algorithm = algorithm;
		data[i].pairedEndInfer = pairedEndInfer;
		data[i].queueLength = queueLength;
		data[i].outputFormat = outputFormat;
		data[i].readGroupString = readGroupString;
		data[i].outputID = outputID;
		data[i].mappedEndCounts = &mappedEndCounts;
		data[i].mappedEndCountsNumEnds = &mappedEndCountsNumEnds;
		data[i].numReported = &numReported;
		data[i].numUnmapped = &numUnmapped;
		data[i].pairedEndInferRescued = &pairedEndInferRescued;

		data[i].inputFP_mutex = &inputFP_mutex;
		data[i].outputFP_mutex = &inputFP_mutex;
		data[i].fpIn = fp;
		data[i].fpReported = fpReported;
		data[i].fpReportedGZ = fpReportedGZ;
		data[i].fpUnmapped = fpUnmapped;
		data[i].threadID = i+1;
	}

	/* Go through each read */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Processing reads, currently on:\n");
	}

	/* Open threads */
	for(i=0;i<numThreads;i++) {
		/* Start thread */
		errCode = pthread_create(&threads[i], /* thread struct */
				NULL, /* default thread attributes */
				ReadInputFilterAndOutputThread, /* start routine */
				&data[i]); /* data to routine */
		if(0!=errCode) {
			PrintError(FnName, "pthread_create: errCode", "Could not start thread", Exit, ThreadError);
		}
	}
	/* Wait for threads to return */
	for(i=0;i<numThreads;i++) {
		/* Wait for the given thread to return */
		errCode = pthread_join(threads[i],
				&status);
		/* Check the return code of the thread */
		if(0!=errCode) {
			PrintError(FnName, "pthread_join: errCode", "Thread returned an error", Exit, ThreadError);
		}
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\n");
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

	if(VERBOSE>=0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Found %10lld reads with no ends mapped.\n", 
				(long long int)numUnmapped);
		if(!(mappedEndCountsNumEnds < 1 || numUnmapped == mappedEndCounts[0])) {
			fprintf(stderr, "%d < 1 || %d == %d\n",
					mappedEndCountsNumEnds,
					(int)numUnmapped,
					mappedEndCounts[0]);
		}
		assert(mappedEndCountsNumEnds < 1 || numUnmapped == mappedEndCounts[0]);
		for(i=1;i<=mappedEndCountsNumEnds;i++) {
			if(1 == i) fprintf(stderr, "Found %10d reads with %2d end mapped.\n", mappedEndCounts[i], i);
			else fprintf(stderr, "Found %10d reads with %2d ends mapped.\n", mappedEndCounts[i], i);
		}
		fprintf(stderr, "Found %10lld reads with at least one end mapping.\n",
				(long long int)numReported);
		if(1 == pairedEndInfer) {
			fprintf(stderr, "Rescued %d ends using the paired end insert distribution.\n",
					pairedEndInferRescued);
		}
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Free */
	if(1 == pairedEndInfer) {
		PEDBinsFree(&bins);
	}
	free(mappedEndCounts);
	free(readGroupString);
	free(threads);
	free(data);
}

void *ReadInputFilterAndOutputThread(void *arg)
{
	char *FnName="ReadInputFilterAndOutputThread";
	PostProcessThreadData *data = (PostProcessThreadData*)arg;
	PEDBins *bins = data->bins;
	RGBinary *rg = data->rg;
	int algorithm = data->algorithm;
	int pairedEndInfer = data->pairedEndInfer;
	int queueLength = data->queueLength;
	int outputFormat = data->outputFormat;
	char *outputID = data->outputID;
	char *readGroupString = data->readGroupString;
	int32_t **mappedEndCounts = data->mappedEndCounts;
	int32_t *mappedEndCountsNumEnds = data->mappedEndCountsNumEnds;
	int32_t *numReported = data->numReported;
	int32_t *numUnmapped = data->numUnmapped;
	int32_t *pairedEndInferRescued = data->pairedEndInferRescued;
	pthread_mutex_t *inputFP_mutex = data->inputFP_mutex;
	gzFile fpIn = data->fpIn;
	FILE *fpReported = data->fpReported;
	gzFile fpReportedGZ = data->fpReportedGZ;
	gzFile fpUnmapped = data->fpUnmapped;
	pthread_mutex_t *outputFP_mutex = data->outputFP_mutex;
	int32_t threadID = data->threadID;

	int32_t i, j;
	int32_t numEnds=0;
	int64_t counter, foundType;
	AlignedRead *aBuffer=NULL;
	int32_t aBufferLength=0;
	int32_t numRead, aBufferIndex;
	int32_t **numEntries=NULL;
	int32_t *numEntriesN=NULL;

	aBufferLength=queueLength;
	aBuffer=malloc(sizeof(AlignedRead)*aBufferLength);
	if(NULL == aBuffer) {
		PrintError(FnName, "aBuffer", "Could not allocate memory", Exit, MallocMemory);
	}
	numEntries=malloc(sizeof(int32_t*)*aBufferLength);
	if(NULL == numEntries) {
		PrintError(FnName, "numEntries", "Could not allocate memory", Exit, MallocMemory);
	}
	numEntriesN=malloc(sizeof(int32_t)*aBufferLength);
	if(NULL == numEntriesN) {
		PrintError(FnName, "numEntriesN", "Could not allocate memory", Exit, MallocMemory);
	}
	for(i=0;i<aBufferLength;i++) {
		numEntries[i] = NULL;
		numEntries[i] = 0;
	}

	counter = numRead = 0;

	while(0 != (numRead = GetAlignedReads(fpIn, aBuffer, aBufferLength, inputFP_mutex))) {

		for(aBufferIndex=0;aBufferIndex<numRead;aBufferIndex++) {

			if(VERBOSE >= 0 && counter%ALIGNENTRIES_READ_ROTATE_NUM==0) {
				fprintf(stderr, "\rthreadID:%d\t[%lld]",
						threadID,
						(long long int)counter);
			}

			// Store the original # of entries for SAM output
			if(aBuffer[aBufferIndex].numEnds < numEntriesN[aBufferIndex]) {
				numEntriesN[aBufferIndex] = aBuffer[aBufferIndex].numEnds;
				numEntries[aBufferIndex]=realloc(numEntries[aBufferIndex], sizeof(int32_t)*numEntriesN[aBufferIndex]);
				if(NULL == numEntries[aBufferIndex]) {
					PrintError(FnName, "numEntries[aBufferIndex]", "Could not reallocate memory", Exit, ReallocMemory);
				}
				for(i=0;i<aBuffer[aBufferIndex].numEnds;i++) {
					numEntries[aBufferIndex][i] = aBuffer[aBufferIndex].ends[i].numEntries;
				}
			}

			/* Filter */
			foundType=FilterAlignedRead(&aBuffer[aBufferIndex],
					algorithm,
					pairedEndInfer,
					bins,
					pairedEndInferRescued);

			numEnds=0;
			if(NoneFound == foundType) {
				/* Print to Not Reported file */
				if(NULL != fpUnmapped) {
					AlignedReadPrint(&aBuffer[aBufferIndex], fpUnmapped);
				}

				/* Free the alignments for output */
				for(i=0;i<aBuffer[aBufferIndex].numEnds;i++) {
					for(j=0;j<aBuffer[aBufferIndex].ends[i].numEntries;j++) {
						AlignedEntryFree(&aBuffer[aBufferIndex].ends[i].entries[j]);
					}
					aBuffer[aBufferIndex].ends[i].numEntries=0;
				}
			}
			else {
				(*numReported)++;
			}

			/* Increment counter */
			counter++;
		}

		/* Print to Output file */
		pthread_mutex_lock(outputFP_mutex);
		for(aBufferIndex=0;aBufferIndex<numRead;aBufferIndex++) {
			// Get the # of ends
			numEnds = 0;
			for(i=0;i<aBuffer[aBufferIndex].numEnds;i++) {
				if(0 < aBuffer[aBufferIndex].ends[i].numEntries) {
					numEnds++;
				}
			}
			if(0 == numEnds) {
				(*numUnmapped)++;
			}
			// Clean up
			if((*mappedEndCountsNumEnds) < numEnds) {
				// Reallocate
				(*mappedEndCounts) = realloc((*mappedEndCounts), sizeof(int32_t)*(1+numEnds));
				if(NULL == (*mappedEndCounts)) {
					PrintError(FnName, "mappedEndCounts", "Could not reallocate memory", Exit, ReallocMemory);
				}
				// Initialize
				for(i=1+(*mappedEndCountsNumEnds);i<=numEnds;i++) {
					(*mappedEndCounts)[i] = 0;
				}
				(*mappedEndCountsNumEnds) = numEnds;
			}
			(*mappedEndCounts)[numEnds]++;

			AlignedReadConvertPrintOutputFormat(&aBuffer[aBufferIndex], rg, fpReported, fpReportedGZ, (NULL == outputID) ? "" : outputID, readGroupString, algorithm, numEntries[aBufferIndex], outputFormat, BinaryOutput);

			/* Free memory */
			AlignedReadFree(&aBuffer[aBufferIndex]);
		}
		pthread_mutex_unlock(outputFP_mutex);
		if(VERBOSE >= 0 && counter%ALIGNENTRIES_READ_ROTATE_NUM==0) {
			fprintf(stderr, "\rthreadID:%d\t[%lld]",
					threadID,
					(long long int)counter);
		}
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\rthreadID:%d\t[%lld]",
				threadID,
				(long long int)counter);
	}

	free(aBuffer);
	for(i=0;i<aBufferLength;i++) {
		free(numEntries[i]);
	}
	free(numEntries);
	free(numEntriesN);
	return arg;
}

int32_t GetAlignedReads(gzFile fp, AlignedRead *aBuffer, int32_t maxToRead, pthread_mutex_t *inputFP_mutex)
{
	if(NULL != inputFP_mutex) pthread_mutex_lock(inputFP_mutex);
	int32_t numRead=0;
	while(numRead < maxToRead) {
		AlignedReadInitialize(&aBuffer[numRead]);
		if(EOF == AlignedReadRead(&aBuffer[numRead], fp)) {
			break;
		}
		numRead++;
	}
	if(NULL != inputFP_mutex) pthread_mutex_unlock(inputFP_mutex);
	return numRead;
}

/* TODO */
int FilterAlignedRead(AlignedRead *a,
		int algorithm,
		int  pairedEndInfer,
		PEDBins *b,
		int32_t *pairedEndInferRescued)
{
	char *FnName="FilterAlignedRead";
	int foundType;
	int32_t *foundTypes=NULL;
	AlignedRead tmpA;
	int32_t i, j, ctr;
	int32_t best, bestIndex, numBest;
	int32_t prob, bestProb;

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

	/* If we are going to use the paired end to infer the other end, then
	 * use the BestScoreAll algorithm and try to infer. */
	if(1 == pairedEndInfer && 2 == tmpA.numEnds) {
		assert(BestScore == algorithm);
		algorithm = BestScoreAll;
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
				else if(BestScoreAll == algorithm &&
						1 <= numBest) {
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

	/* Now see if need to infer one end using the other */
	if(1 == pairedEndInfer && 2 == tmpA.numEnds) {
		if(1 == tmpA.ends[0].numEntries &&
				1 == tmpA.ends[1].numEntries) {
			foundType = Found;
		}
		else if((1 == tmpA.ends[0].numEntries && 1 < tmpA.ends[1].numEntries) ||
				(1 < tmpA.ends[0].numEntries && 1 == tmpA.ends[1].numEntries)) {
			int32_t endOne, endTwo;
			if(1 == tmpA.ends[0].numEntries &&
					1 < tmpA.ends[1].numEntries) {
				endOne = 0;
				endTwo = 1;
			}
			else {
				endOne = 1;
				endTwo = 0;
			}
			// Get entry with the highest probability
			bestIndex = -1;
			bestProb = -1.0;
			numBest = 0;
			for(i=0;i<tmpA.ends[endTwo].numEntries;i++) {
				prob = PEDBinsGetProbability(b,
						tmpA.ends[endOne].entries[0].contig,
						tmpA.ends[endOne].entries[0].position,
						tmpA.ends[endTwo].entries[i].contig,
						tmpA.ends[endTwo].entries[i].position);
				if(bestProb < prob) {
					bestProb = prob;
					bestIndex = i;
					numBest = 1;
				}
				else if(prob < bestProb) {
					// Ignore
				}
				else {
					numBest++;
				}
			}
			// Check if this helped
			if(1 == numBest) { // Keep
				/* Copy to front */
				foundTypes[endTwo] = Found;
				AlignedEntryCopy(&tmpA.ends[endTwo].entries[0],
						&tmpA.ends[endTwo].entries[bestIndex]);
				AlignedEndReallocate(&tmpA.ends[endTwo], 1);
				tmpA.ends[endTwo].entries[0].mappingQuality = bestProb * tmpA.ends[endOne].entries[0].mappingQuality / b->numDistances;
				// Updating the mapping quality (how?)
				(*pairedEndInferRescued)++;
			}
			else { // Remove
				AlignedEndReallocate(&tmpA.ends[endTwo], 0);
				foundTypes[endTwo] = NoneFound;
			}
			foundType = Found; // We found at least one end!
		}
		else {
			/* Check each end individually */
			foundType=NoneFound;
			for(i=0;i<tmpA.numEnds;i++) {
				if(1 == tmpA.ends[i].numEntries) { // Keep
					foundTypes[i] = foundType = Found;
				}
				else { // Remove
					AlignedEndReallocate(&tmpA.ends[i], 0);
					foundTypes[i] = NoneFound;
				}
			}
		}
	}
	else if(1 == tmpA.numEnds) {
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

int32_t GetPEDBins(char *inputFileName, 
		int algorithm,
		int queueLength,
		PEDBins *b)
{
	char *FnName="GetPEDBins";
	gzFile fp=NULL;
	int64_t counter, foundType;
	AlignedRead *aBuffer;
	int32_t aBufferLength=0;
	int32_t numRead, aBufferIndex;

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

	aBufferLength=queueLength;
	aBuffer=malloc(sizeof(AlignedRead)*aBufferLength);
	if(NULL == aBuffer) {
		PrintError(FnName, "aBuffer", "Could not allocate memory", Exit, MallocMemory);
	}

	/* Go through each read */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Estimating paired end distance, currently on:\n0");
	}
	counter = numRead = 0;

	// TODO: make this multi-threaded

	while(0 != (numRead = GetAlignedReads(fp, aBuffer, aBufferLength, NULL))) {

		for(aBufferIndex=0;aBufferIndex<numRead;aBufferIndex++) {

			if(2 == aBuffer[aBufferIndex].numEnds) { // Only paired end data

				if(VERBOSE >= 0 && counter%ALIGNENTRIES_READ_ROTATE_NUM==0) {
					fprintf(stderr, "\r%lld",
							(long long int)counter);
				}

				/* Filter */
				foundType=FilterAlignedRead(&aBuffer[aBufferIndex],
						algorithm,
						0,
						NULL,
						NULL);

				if(Found == foundType) {
					assert(2 == aBuffer[aBufferIndex].numEnds);
					/* Must only have one alignment per end and on the same contig.
					 * There is a potential this will be inferred incorrectly under
					 * many scenarios.  Be careful! */
					if(1 == aBuffer[aBufferIndex].ends[0].numEntries &&
							1 == aBuffer[aBufferIndex].ends[1].numEntries &&
							aBuffer[aBufferIndex].ends[0].entries[0].contig == aBuffer[aBufferIndex].ends[1].entries[0].contig) {
						PEDBinsInsert(b, 
								fabs(aBuffer[aBufferIndex].ends[1].entries[0].position - aBuffer[aBufferIndex].ends[0].entries[0].position));
					}
				}
			}

			/* Increment counter */
			counter++;
		}

		if(VERBOSE >= 0) {
			fprintf(stderr, "\r%lld",
					(long long int)counter);
		}

		/* Free buffer */
		for(aBufferIndex=0;aBufferIndex<numRead;aBufferIndex++) {
			/* Free memory */
			AlignedReadFree(&aBuffer[aBufferIndex]);
		}

		/* Do we really need this many? */
		if(MAX_PEDBINS_SIZE < counter) {
			break;
		}
	}
	if(VERBOSE>=0) {
		fprintf(stderr, "\r%lld\n",
				(long long int)counter);
	}

	/* Close the input file */
	gzclose(fp);
	free(aBuffer);

	if(b->numDistances < MIN_PEDBINS_SIZE) {
		fprintf(stderr, "Found only %d distances to infer the insert size distribution\n", b->numDistances);
		PrintError(FnName, "b->numDistances", "Not enough distances to infer insert size distribution", Warn, OutOfRange);
		PEDBinsFree(b);
		return 0;
	}

	if(VERBOSE>=0) {
		// Print Statistics
		PEDBinsPrintStatistics(b, stderr);
	}
	return 1;
}

void PEDBinsInitialize(PEDBins *b)
{
	b->minDistance = INT_MAX;
	b->maxDistance = INT_MIN;
	b->bins = NULL;
	b->numDistances = 0;
}

void PEDBinsFree(PEDBins *b) {
	free(b->bins);
	PEDBinsInitialize(b);
}

void PEDBinsInsert(PEDBins *b,
		int32_t distance)
{
	char *FnName="PEDBinsInsert";
	int32_t prevMinDistance, prevMaxDistance;
	int32_t i;

	// TODO: there is a bug here with an off-by-one error
	// such that b->numDistances does not equal the sum of
	// all the bins.

	if(distance < MIN_PEDBINS_DISTANCE ||
			MAX_PEDBINS_DISTANCE < distance) {
		return;
	}

	if(distance < b->minDistance) {
		if(INT_MIN < b->maxDistance) { // Bins exist
			prevMinDistance = b->minDistance;
			b->minDistance = distance;
			b->bins = realloc(b->bins, sizeof(int32_t)*(b->maxDistance - b->minDistance + 1));
			if(NULL == b->bins) {
				PrintError(FnName, "b->bins", "Could not reallocate memory", Exit, ReallocMemory);
			}
			// Move over old bins
			for(i=b->maxDistance - prevMinDistance;0<=i;i--) {
				b->bins[i + (prevMinDistance - b->minDistance)] = b->bins[i];
			}
			// Initialize new bins
			for(i=0;i < prevMinDistance - b->minDistance;i++) {
				b->bins[i] = 0;
			}
		}
		else { // No bins
			b->minDistance = b->maxDistance = distance;
			b->bins = malloc(sizeof(int32_t));
			if(NULL == b->bins) {
				PrintError(FnName, "b->bins", "Could not allocate memory", Exit, MallocMemory);
			}
		}
	}
	else if(b->maxDistance < distance) {
		if(b->minDistance < INT_MAX) { // Bins exist
			prevMaxDistance = b->maxDistance;
			b->maxDistance = distance;
			b->bins = realloc(b->bins, sizeof(int32_t)*(b->maxDistance - b->minDistance + 1));
			if(NULL == b->bins) {
				PrintError(FnName, "b->bins", "Could not reallocate memory", Exit, ReallocMemory);
			}
			// Initialize new bins
			for(i=prevMaxDistance-b->minDistance+1;i<b->maxDistance-b->minDistance+1;i++) {
				b->bins[i] = 0;
			}
		}
		else { // No bins
			b->minDistance = b->maxDistance = distance;
			b->bins = malloc(sizeof(int32_t));
			if(NULL == b->bins) {
				PrintError(FnName, "b->bins", "Could not allocate memory", Exit, MallocMemory);
			}
		}
	}

	// Add to bin
	b->bins[distance - b->minDistance]++;
	b->numDistances++;
}

void PEDBinsPrintStatistics(PEDBins *b, FILE *fp)
{
	// Mean, Range, and SD
	int32_t i;
	double mean, sd;

	// Mean
	mean = 0.0;
	for(i=0;i<b->maxDistance-b->minDistance+1;i++) {
		mean += (b->minDistance + i)*b->bins[i];
	}
	mean /= b->numDistances;

	// SD
	sd = 0.0;
	for(i=0;i<b->maxDistance-b->minDistance+1;i++) {
		sd += b->bins[i]*((b->minDistance + i) - mean)*((b->minDistance + i) - mean);
	}
	sd /= b->numDistances-1;
	sd = sqrt(sd);

	if(0<=VERBOSE) {
		fprintf(stderr, "Used %d paired end distances to infer the insert size distribution.\n",
				b->numDistances);
		fprintf(stderr, "The paired end distance range was from %d to %d.\n",
				b->minDistance, b->maxDistance);
		fprintf(stderr, "The paired end distance mean and standard deviation was %.2lf and %.2lf.\n",
				mean, sd);
		fprintf(stderr, "%s", BREAK_LINE);
	}
}

void PEDBinsMakeIntoProbability(PEDBins *b)
{
	// Mean, Range, and SD
	int32_t i;
	double mean;
	int32_t distance, mid, sum;

	// Mean
	mean = 0.0;
	for(i=0;i<b->maxDistance-b->minDistance+1;i++) {
		mean += (b->minDistance + i)*b->bins[i];
	}
	mean /= b->numDistances;

	mid = (int)mean;
	if(b->maxDistance - mid < mid - b->minDistance) {
		distance = mid - b->minDistance + 1;
	}
	else {
		distance = b->maxDistance - mid + 1;
	}

	sum=0;
	while(0 != distance) {
		// Compute sum
		if(distance <= b->maxDistance - mid) {
			assert(mid + distance - b->minDistance < b->maxDistance - b->minDistance + 1);
			sum += b->bins[mid + distance - b->minDistance];
		}
		if(distance <= mid - b->minDistance) {
			assert(0 <= mid - distance - b->minDistance);
			sum += b->bins[mid - distance - b->minDistance];
		}
		// Put in sum
		if(distance <= b->maxDistance - mid) {
			assert(mid + distance - b->minDistance < b->maxDistance - b->minDistance + 1);
			b->bins[mid + distance - b->minDistance] = sum;
		}
		if(distance <= mid - b->minDistance) {
			assert(0 <= mid - distance - b->minDistance);
			b->bins[mid - distance - b->minDistance] = sum;
		}
		/*
		   fprintf(stderr, "distance=%d\tsum=%d\n",
		   distance, sum);
		   */
		// Decrement the distance
		distance--;
	}
	// Take care of the final case
	sum += b->bins[mid - b->minDistance];
	b->bins[mid - b->minDistance] = sum;
}

// Not normalized
int32_t PEDBinsGetProbability(PEDBins *b, int32_t contigA, int32_t positionA, int32_t contigB, int32_t positionB)
{
	int32_t distance;

	// Must be on the same contig
	if(contigA != contigB) {
		return 0.0;
	}

	// Get the distance
	distance = (int)fabs(positionA - positionB);

	// Distance is out of bounds
	if(distance < b->minDistance ||
			b->maxDistance < distance) {
		return 0.0;
	}

	return b->bins[distance - b->minDistance];
}
