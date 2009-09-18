#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#include <zlib.h>

#endif
#include "AlignedEntry.h"
#include "AlignedRead.h"
#include "AlignedEnd.h"
#include "../BLibDefinitions.h"
#include "../BLib.h"
#include "../BError.h"
#include "bsort.h"

#define Name "bsort"
#define BSORT_ROTATE_NUM 100000
#define BSORT_MAX_LINE_LENGTH 100

/* Sorts a bfast report file.
 * */

void TmpGZFileOpen(TmpGZFile *tmpFile,
		char *tmpDir)
{
	TmpGZFileInitialize(tmpFile);
	tmpFile->FP = OpenTmpGZFile(tmpDir, &tmpFile->FileName);
}

void TmpGZFileClose(TmpGZFile *tmpFile) 
{
	CloseTmpGZFile(&tmpFile->FP, &tmpFile->FileName, 1);
	TmpGZFileInitialize(tmpFile);
}

void TmpGZFileInitialize(TmpGZFile *tmpFile)
{
	tmpFile->FP = NULL;
	tmpFile->FileName = NULL;
	tmpFile->startContig = INT_MAX;
	tmpFile->startPos = INT_MAX;
	tmpFile->endContig= 0;
	tmpFile->endPos = 0;
	tmpFile->numEntries = 0;
}

void TmpGZFileUpdateMetaData(TmpGZFile *tmpFile,
		AlignedRead *a)
{
	int32_t minIndex=INT_MIN;
	int32_t i;

	for(i=0;i<a->numEnds;i++) {
		if(0 < a->ends[i].numEntries) {
			assert(1 == a->ends[i].numEntries);
			if(minIndex < 0 ||
					AlignedEntryCompare(&a->ends[i].entries[0], 
						&a->ends[minIndex].entries[0], 
						AlignedEntrySortByContigPos) < 0) {
				minIndex = i;
			}
		}
	}
	assert(0 <= minIndex);

	TmpGZFileUpdateMetaDataHelper(tmpFile, &a->ends[minIndex].entries[0]);
	tmpFile->numEntries++;
}

void TmpGZFileUpdateMetaDataHelper(TmpGZFile *tmpFile,
		AlignedEntry *a)
{
	if(a->contig < tmpFile->startContig || 
			(a->contig == tmpFile->startContig && a->position < tmpFile->startPos)) {
		tmpFile->startContig = a->contig;
		tmpFile->startPos = a->position;
	}
	if(tmpFile->endContig < a->contig ||
			(tmpFile->endContig == a->contig && tmpFile->endPos < a->position)) {
		tmpFile->endContig = a->contig;
		tmpFile->endPos = a->position;
	}
}

void MoveAllIntoTmpGZFile(char *inputFileName, 
		TmpGZFile *tmpFile,
		char *tmpDir)
{
	char *FnName="MoveAllIntoTmpGZFile";
	int32_t i;
	int64_t counter=0;
	gzFile fpIn;
	AlignedRead a;

	/* Open tmp file */
	TmpGZFileOpen(tmpFile, tmpDir);

	/* Open the input file */
	if(!(fpIn=gzopen(inputFileName, "rb"))) {
		PrintError(Name,
				inputFileName,
				"Could not open file for reading",
				Exit,
				OpenFileError);
	}

	/* Move all entries into the tmp file */
	fprintf(stderr, "Moving all entries into a tmp file.  Currently on read:\n0");
	AlignedReadInitialize(&a);
	while(EOF != AlignedReadRead(&a, fpIn)) {
		if(counter%BSORT_ROTATE_NUM==0) {
			fprintf(stderr, "\r%lld",
					(long long int)counter);
		}
		counter++;

		/* Store AlignedRead */
		for(i=0;i<a.numEnds;i++) {
			if(0 != a.ends[i].numEntries &&
					1 != a.ends[i].numEntries) {
				PrintError(FnName,
						a.readName,
						"Read was not uniquely aligned",
						Exit,
						OutOfRange);
			}
		}
		AlignedReadPrint(&a, 
				tmpFile->FP);
		TmpGZFileUpdateMetaData(tmpFile, 
				&a);
		AlignedReadFree(&a);
	}
	fprintf(stderr, "\r%lld\n",
			(long long int)counter);
	fprintf(stderr, "Sorting range is contig%lld:%lld to contig%lld:%lld.\n",
			(long long int)tmpFile->startContig,
			(long long int)tmpFile->startPos,
			(long long int)tmpFile->endContig,
			(long long int)tmpFile->endPos
		   );

	/* Close the input file */
	gzclose(fpIn);
}

void SplitEntriesAndPrint(gzFile outputFP,
		TmpGZFile *tmpFile, 
		char *tmpDir,
		int32_t maxNumEntries,
		int32_t numThreads)
{
	char *FnName="SplitEntriesAndPrint";
	int64_t meanPos;
	int32_t meanContig;
	AlignedRead *entries=NULL;
	AlignedRead **entriesPtr=NULL;
	AlignedRead a;
	int64_t numEntries=0;
	TmpGZFile belowTmpGZFile, aboveTmpGZFile;
	int32_t endContig, endPos;
	int32_t belowMinPos, belowMinContig, belowMaxPos, belowMaxContig;
	int32_t aboveMinPos, aboveMinContig, aboveMaxPos, aboveMaxContig;
	AlignedEntry *tmpAlignedEntry=NULL;
	int64_t i, j, numPrinted;
	ThreadSortData *sortData=NULL;
	ThreadSortData *mergeData=NULL;
	pthread_t *threads=NULL;
	int errCode;
	void *status;
	int32_t curNumThreads = numThreads;
	int32_t curMergeIteration, curThread;
	int32_t minIndex=INT_MIN;

	assert(tmpFile->startContig < tmpFile->endContig ||
			(tmpFile->startContig == tmpFile->endContig && tmpFile->startPos <= tmpFile->endPos));
	assert(0 < tmpFile->numEntries);
	/*
	   if(tmpFile->numEntries <= 0) {
	   return;
	   }
	   */

	/* Move to the beginning of the tmp file */
	CloseTmpGZFile(&tmpFile->FP, &tmpFile->FileName, 0);
	if(!(tmpFile->FP=gzopen(tmpFile->FileName, "rb"))) {
		PrintError(FnName,
				tmpFile->FileName,
				"Could not re-open file for reading",
				Exit,
				OpenFileError);
	}

	/* Check if we should print or split */
	if(tmpFile->numEntries <= maxNumEntries) {
		/* Sort and print */
		numPrinted=PrintContigPos(stderr,
				tmpFile->startContig,
				tmpFile->startPos);
		numPrinted+=fprintf(stderr, ": sorting");
		assert(numPrinted<=BSORT_MAX_LINE_LENGTH);
		while(numPrinted<=BSORT_MAX_LINE_LENGTH) {
			fprintf(stderr, " ");
			numPrinted++;
		}
		endContig = tmpFile->endContig;
		endPos = tmpFile->endPos;
		assert(0 < tmpFile->numEntries);

		/* Allocate memory for the entries */
		entriesPtr = malloc(sizeof(AlignedRead*)*tmpFile->numEntries);
		if(NULL == entriesPtr) {
			PrintError(FnName,
					"entriesPtr",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		entries = malloc(sizeof(AlignedRead)*tmpFile->numEntries);
		if(NULL == entries) {
			PrintError(FnName,
					"entries",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Initialize */
		for(i=0;i<tmpFile->numEntries;i++) {
			entriesPtr[i] = &entries[i];
			AlignedReadInitialize(entriesPtr[i]);
		}

		/* Read in, sort, and print */
		numEntries = 0;
		while(numEntries < tmpFile->numEntries &&
				EOF != AlignedReadRead(entriesPtr[numEntries],
					tmpFile->FP)) {
			assert(numEntries < tmpFile->numEntries);
			numEntries++;
		}
		assert(numEntries == tmpFile->numEntries);

		/* Close the file */
		TmpGZFileClose(tmpFile);
		/* Sort */
		if(numEntries < BSORT_THREADED_SORT_MIN ||
				numThreads <= 1) {
			/* Ignore threading */
			AlignedReadMergeSortAll(entriesPtr, 
					0,
					numEntries-1);
		}
		else {
			/* Should check that the number of threads is a power of 4 since we split
			 * in half in both sorts. */
			assert(IsAPowerOfTwo(numThreads)==1);

			/* Allocate memory for the thread arguments */
			sortData = malloc(sizeof(ThreadSortData)*numThreads);
			if(NULL==sortData) {
				PrintError(FnName,
						"sortData",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
			/* Allocate memory for the thread point32_ters */
			threads = malloc(sizeof(pthread_t)*numThreads);
			if(NULL==threads) {
				PrintError(FnName,
						"threads",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}

			/* Initialize sortData */
			for(i=0;i<numThreads;i++) {
				sortData[i].entriesPtr=entriesPtr;
				sortData[i].threadID = i;
				sortData[i].low = i*(numEntries/numThreads);
				sortData[i].high = (i+1)*(numEntries/numThreads)-1;
				assert(sortData[i].low >= 0 && sortData[i].high < numEntries);
			}
			sortData[0].low = 0;
			sortData[numThreads-1].high = numEntries-1;

			/* Check that we split correctly */
			for(i=1;i<numThreads;i++) {
				assert(sortData[i-1].high < sortData[i].low);
			}
			/* Create threads */
			for(i=0;i<numThreads;i++) {
				/* Start thread */
				errCode = pthread_create(&threads[i], /* thread struct */
						NULL, /* default thread attributes */
						SortAlignedReadHelper, /* start routine */
						(void*)(&sortData[i])); /* sortData to routine */
				if(0!=errCode) {
					PrintError(FnName,
							"pthread_create: errCode",
							"Could not start thread",
							Exit,
							ThreadError);
				}
			}

			/* Wait for the threads to finish */
			for(i=0;i<numThreads;i++) {
				/* Wait for the given thread to return */
				errCode = pthread_join(threads[i],
						&status);
				/* Check the return code of the thread */
				if(0!=errCode) {
					PrintError(FnName,
							"pthread_join: errCode",
							"Thread returned an error",
							Exit,
							ThreadError);
				}
			}
			/* Free memory for the threads */
			free(threads);
			threads=NULL;

			/* Now we must merge the results from the threads */
			/* Merge intelligently i.e. merge recursively so 
			 * there are only nlogn merges where n is the 
			 * number of threads. */
			curNumThreads = numThreads;
			for(i=1, curMergeIteration=1;i<numThreads;i=i*2, curMergeIteration++) { /* The number of merge iterations */
				curNumThreads /= 2; /* The number of threads to spawn */
				/* Allocate memory for the thread arguments */
				mergeData = malloc(sizeof(ThreadRGIndexMergeData)*curNumThreads);
				if(NULL==mergeData) {
					PrintError(FnName,
							"mergeData",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}
				/* Allocate memory for the thread point32_ters */
				threads = malloc(sizeof(pthread_t)*curNumThreads);
				if(NULL==threads) {
					PrintError(FnName,
							"threads",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}
				/* Initialize data for threads */
				for(j=0,curThread=0;j<numThreads;j+=2*i,curThread++) {
					mergeData[curThread].entriesPtr = entriesPtr;
					mergeData[curThread].threadID = curThread;
					/* Use the same bounds as was used in the sort */
					mergeData[curThread].low = sortData[j].low;
					mergeData[curThread].mid = sortData[i+j].low-1;
					mergeData[curThread].high = sortData[j+2*i-1].high;
				}
				/* Check that we split correctly */
				for(j=1;j<curNumThreads;j++) {
					if(mergeData[j-1].high >= mergeData[j].low) {
						PrintError(FnName,
								NULL,
								"mergeData[j-1].high >= mergeData[j].low",
								Exit,
								OutOfRange);
					}
				}
				/* Create threads */
				for(j=0;j<curNumThreads;j++) {
					/* Start thread */
					errCode = pthread_create(&threads[j], /* thread struct */
							NULL, /* default thread attributes */
							MergeAlignedReadHelper, /* start routine */
							(void*)(&mergeData[j])); /* sortData to routine */
					if(0!=errCode) {
						PrintError(FnName,
								"pthread_create: errCode",
								"Could not start thread",
								Exit,
								ThreadError);
					}
				}

				/* Wait for the threads to finish */
				for(j=0;j<curNumThreads;j++) {
					/* Wait for the given thread to return */
					errCode = pthread_join(threads[j],
							&status);
					/* Check the return code of the thread */
					if(0!=errCode) {
						PrintError(FnName,
								"pthread_join: errCode",
								"Thread returned an error",
								Exit,
								ThreadError);
					}
				}

				/* Free memory for the merge data */
				free(mergeData);
				mergeData=NULL;
				/* Free memory for the threads */
				free(threads);
				threads=NULL;
			}

			/* Free memory for sort data */
			free(sortData);
			sortData=NULL;

		}

		/* Print and Free memory */
		for(i=0;i<numEntries;i++) {
			/* Print */
			AlignedReadPrint(entriesPtr[i],
					outputFP);
			/* Free memory */
			AlignedReadFree(entriesPtr[i]);
			entriesPtr[i]=NULL;
		}
		free(entriesPtr);
		entriesPtr=NULL;
		free(entries);
		entries=NULL;
	}
	else if(tmpFile->startContig == tmpFile->endContig && 
			tmpFile->startPos == tmpFile->endPos) {
		PrintError(FnName,
				NULL,
				"Could not split the file any further.  Try increasing your the maximum number of entries.",
				Exit,
				OutOfRange);
	}
	else {
		/* Split and recurse */
		/*
		   fprintf(stderr, "\rSplitting file of size %lld.",
		   (long long int)tmpFile->numEntries);
		   */

		/* Initialize */
		AlignedReadInitialize(&a);
		TmpGZFileOpen(&belowTmpGZFile, tmpDir);
		TmpGZFileOpen(&aboveTmpGZFile, tmpDir);
		if(tmpFile->startContig == tmpFile->endContig) {
			meanPos = (tmpFile->startPos + tmpFile->endPos)/2;
			numPrinted=PrintContigPos(stderr,
					tmpFile->startContig,
					tmpFile->startPos);
			numPrinted+=fprintf(stderr, ": splitting file of size %lld by position %lld on contig %lld",
					(long long int)tmpFile->numEntries,
					(long long int)meanPos,
					(long long int)tmpFile->startContig
					);
			belowTmpGZFile.startContig = tmpFile->startContig;
			belowTmpGZFile.startPos = tmpFile->startPos;
			belowTmpGZFile.endContig = tmpFile->endContig;
			belowTmpGZFile.endPos = meanPos;
			aboveTmpGZFile.startContig = tmpFile->startContig;
			aboveTmpGZFile.startPos = meanPos + 1;
			aboveTmpGZFile.endContig = tmpFile->endContig;
			aboveTmpGZFile.endPos = tmpFile->endPos;
		}
		else {
			meanContig = (tmpFile->startContig + tmpFile->endContig)/2;
			numPrinted=PrintContigPos(stderr,
					tmpFile->startContig,
					tmpFile->startPos);
			numPrinted+=fprintf(stderr, ": splitting file of size %lld by contig %d",
					(long long int)tmpFile->numEntries,
					meanContig
					);
			belowTmpGZFile.startContig = tmpFile->startContig;
			belowTmpGZFile.startPos = tmpFile->startPos;
			belowTmpGZFile.endContig = meanContig;
			belowTmpGZFile.endPos = INT_MAX-1;
			aboveTmpGZFile.startContig = meanContig + 1;
			aboveTmpGZFile.startPos = 1;
			aboveTmpGZFile.endContig = tmpFile->endContig;
			aboveTmpGZFile.endPos = tmpFile->endPos;
		}
		assert(numPrinted<=BSORT_MAX_LINE_LENGTH);
		while(numPrinted<=BSORT_MAX_LINE_LENGTH) {
			fprintf(stderr, " ");
			numPrinted++;
		}

		belowMinPos = INT_MAX;
		belowMinContig = INT_MAX;
		belowMaxPos = 0;
		belowMaxContig = 0;
		aboveMinPos = INT_MAX;
		aboveMinContig = INT_MAX;
		aboveMaxPos = 0;
		aboveMaxContig = 0;

		/* Split */
		while(EOF != AlignedReadRead(&a,
					tmpFile->FP)) {

			for(i=0;i<a.numEnds;i++) {
				if(0 < a.ends[i].numEntries) {
					assert(1 == a.ends[i].numEntries);
					if(minIndex < 0 ||
							AlignedEntryCompare(&a.ends[i].entries[0], 
								&a.ends[minIndex].entries[0], 
								AlignedEntrySortByContigPos) < 0) {
						minIndex = i;
					}
				}
			}
			assert(0 <= minIndex);
			tmpAlignedEntry = &a.ends[minIndex].entries[0];

			/* Print to the appropriate file */
			if(tmpAlignedEntry->contig < belowTmpGZFile.endContig ||
					(tmpAlignedEntry->contig == belowTmpGZFile.endContig && tmpAlignedEntry->position < belowTmpGZFile.endPos)) {
				/* Print */
				AlignedReadPrint(&a,
						belowTmpGZFile.FP);
				belowTmpGZFile.numEntries++;

				/* Update bounds */
				if(tmpAlignedEntry->contig < belowMinContig ||
						(tmpAlignedEntry->contig == belowMinContig && tmpAlignedEntry->position < belowMinPos)) {
					belowMinContig = tmpAlignedEntry->contig;
					belowMinPos = tmpAlignedEntry->position;
				}
				if(belowMaxContig < tmpAlignedEntry->contig || 
						(belowMaxContig == tmpAlignedEntry->contig && belowMaxPos < tmpAlignedEntry->position)) {
					belowMaxContig = tmpAlignedEntry->contig;
					belowMaxPos = tmpAlignedEntry->position;
				}
			}
			else {
				/* Print */
				AlignedReadPrint(&a,
						aboveTmpGZFile.FP);
				aboveTmpGZFile.numEntries++;

				/* Update bounds */
				if(tmpAlignedEntry->contig < aboveMinContig ||
						(tmpAlignedEntry->contig == aboveMinContig && tmpAlignedEntry->position < aboveMinPos)) {
					aboveMinContig = tmpAlignedEntry->contig;
					aboveMinPos = tmpAlignedEntry->position;
				}
				if(aboveMaxContig < tmpAlignedEntry->contig || 
						(aboveMaxContig == tmpAlignedEntry->contig && aboveMaxPos < tmpAlignedEntry->position)) {
					aboveMaxContig = tmpAlignedEntry->contig;
					aboveMaxPos = tmpAlignedEntry->position;
				}
			}
			AlignedReadFree(&a);
		}

		/* Close tmp file */
		TmpGZFileClose(tmpFile);

		/* Close tmp files if we are not going to use them */
		if(belowTmpGZFile.numEntries <= 0) {
			assert(belowTmpGZFile.startContig < belowMinContig ||
					(belowTmpGZFile.startContig == belowMinContig && belowTmpGZFile.startPos <= belowMinPos));
			assert(belowMaxContig < belowTmpGZFile.endContig ||
					(belowMaxContig == belowTmpGZFile.endContig && belowMaxPos <= belowTmpGZFile.endPos));
			belowTmpGZFile.startContig = belowMinContig;
			belowTmpGZFile.startPos = belowMinPos;
			belowTmpGZFile.endContig = belowMaxContig;
			belowTmpGZFile.endPos = belowMaxPos;
			TmpGZFileClose(&belowTmpGZFile);
		}
		if(aboveTmpGZFile.numEntries <= 0) {
			assert(aboveTmpGZFile.startContig < belowMinContig ||
					(aboveTmpGZFile.startContig == belowMinContig && belowTmpGZFile.startPos <= belowMinPos));
			assert(aboveMaxContig < belowTmpGZFile.endContig ||
					(aboveMaxContig == belowTmpGZFile.endContig && belowMaxPos <= belowTmpGZFile.endPos));
			aboveTmpGZFile.startContig = aboveMinContig;
			aboveTmpGZFile.startPos = aboveMinPos;
			aboveTmpGZFile.endContig = aboveMaxContig;
			aboveTmpGZFile.endPos = aboveMaxPos;
			TmpGZFileClose(&aboveTmpGZFile);
		}

		/* Recurse on the two */
		if(0 < belowTmpGZFile.numEntries) {
			/* First update based on bounds */
			SplitEntriesAndPrint(outputFP,
					&belowTmpGZFile,
					tmpDir,
					maxNumEntries,
					numThreads);
		}
		if(0 < aboveTmpGZFile.numEntries) {
			SplitEntriesAndPrint(outputFP,
					&aboveTmpGZFile,
					tmpDir,
					maxNumEntries,
					numThreads);
		}
	}
}

void *SortAlignedReadHelper(void *arg)
{
	ThreadSortData *data = (ThreadSortData*)arg;
	AlignedRead **entriesPtr = data->entriesPtr;
	int64_t low = data->low;
	int64_t high = data->high;

	AlignedReadMergeSortAll(entriesPtr, 
			low,
			high);
	return arg;
}

void *MergeAlignedReadHelper(void *arg)
{
	ThreadSortData *data = (ThreadSortData*)arg;
	AlignedRead **entriesPtr = data->entriesPtr;
	int64_t low = data->low;
	int64_t mid = data->mid;
	int64_t high = data->high;

	AlignedReadMergeAll(entriesPtr, 
			low,
			mid,
			high);
	return arg;
}

int main(int argc, char *argv[])
{
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	int32_t maxNumEntries=0;
	int32_t numThreads=0;
	char tmpDir[MAX_FILENAME_LENGTH]="\0";
	TmpGZFile tmpFile;
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	gzFile outputFP=NULL;
	char *last=NULL;

	if(argc == 5) {
		strcpy(inputFileName, argv[1]);
		maxNumEntries = atoi(argv[2]);
		numThreads=atoi(argv[3]);
		strcpy(tmpDir, argv[4]);

		assert(1 < maxNumEntries);
		assert(1 <= numThreads);
		if(1 < numThreads) {
			/* Should check that the number of threads is a power of 4 since we split
			 * in half in both sorts. */
			assert(IsAPowerOfTwo(numThreads)==1);
		}

		/* Move all into a tmp file */
		fprintf(stderr, "%s", BREAK_LINE);
		MoveAllIntoTmpGZFile(inputFileName, &tmpFile, tmpDir);
		fprintf(stderr, "%s", BREAK_LINE);

		/* Create output file name */
		last = StrStrGetLast(inputFileName,
				BFAST_ALIGNED_FILE_EXTENSION);
		if(NULL == last) {
			PrintError(Name,
					inputFileName,
					"Could not recognize file extension",
					Exit,
					OutOfRange);
		}
		strncpy(outputFileName, inputFileName, (last - inputFileName));
		strcat(outputFileName, "sorted.");
		strcat(outputFileName, BFAST_ALIGNED_FILE_EXTENSION);

		if(!(outputFP = gzopen(outputFileName, "wb"))) {
			PrintError(Name,
					outputFileName,
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}

		/* Split entries and print */
		fprintf(stderr, "Performing initial split, this could take some time.\n");
		SplitEntriesAndPrint(outputFP,
				&tmpFile,
				tmpDir,
				maxNumEntries,
				numThreads);
		fprintf(stderr, "\n");

		/* Close files */
		gzclose(outputFP);

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully!\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast report file name>\n");
		fprintf(stderr, "\t<maximum number of entries when sorting>\n");
		fprintf(stderr, "\t<number of threads>\n");
		fprintf(stderr, "\t<tmp directory>\n");
	}
	return 0;
}
