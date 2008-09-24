#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <pthread.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BLib.h"
#include "../blib/BError.h"
#include "../blib/RGIndex.h"
#include "../blib/RGRanges.h"
#include "../blib/RGMatch.h"
#include "../blib/RGReads.h"
#include "bindexdist.h"

#define Name "bindexdist"
#define BINDEXDIST_ROTATE_NUM 10000
#define BINDEXDIST_SORT_ROTATE_INC 0.01

/* Prints each unique read from the genome and the number 
 * of times it occurs, where the genome is contained in 
 * the bfast index file.
 * */

int main(int argc, char *argv[]) 
{
	FILE *fp=NULL;
	char indexFileName[MAX_FILENAME_LENGTH]="\0";
	char distributionFileName[MAX_FILENAME_LENGTH]="\0";
	char rgFileName[MAX_FILENAME_LENGTH]="\0";
	char outputDir[MAX_FILENAME_LENGTH]="\0";
	char outputID[MAX_FILENAME_LENGTH]="\0";
	char tmpDir[MAX_FILENAME_LENGTH]="\0";
	int numMismatches = 0;
	int numThreads = 0;

	if(argc == 8) {
		RGBinary rg;
		RGIndex index;

		strcpy(rgFileName, argv[1]);
		strcpy(indexFileName, argv[2]);
		numMismatches = atoi(argv[3]);
		strcpy(outputDir, argv[4]);
		strcpy(outputID, argv[5]);
		strcpy(tmpDir, argv[6]);
		numThreads = atoi(argv[7]);

		/* Create the distribution file name */
		sprintf(distributionFileName, "%s%s.dist.%d",
				outputDir,
				outputID,
				numMismatches);

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Starting %s.\n", Name);

		/* Read in the rg binary file */
		RGBinaryReadBinary(&rg, rgFileName);

		/* Read the index */
		fprintf(stderr, "Reading in index from %s.\n",
				indexFileName);
		if(!(fp=fopen(indexFileName, "rb"))) {
			PrintError(Name,
					indexFileName,
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}
		RGIndexRead(fp, &index, 1);
		fclose(fp);

		fprintf(stderr, "%s", BREAK_LINE);
		PrintDistribution(&index, 
				&rg, 
				distributionFileName,
				numMismatches,
				tmpDir,
				numThreads); 
		fprintf(stderr, "%s", BREAK_LINE);

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Cleaning up.\n");
		/* Delete the index */
		RGIndexDelete(&index);
		/* Delete the rg */
		RGBinaryDelete(&rg);
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully!\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: bindexdist [OPTIONS]\n");
		fprintf(stderr, "\t\t<bfast reference genome file name>\n");
		fprintf(stderr, "\t\t<bfast index file name>\n");
		fprintf(stderr, "\t\t<number of mismatches>\n");
		fprintf(stderr, "\t\t<output directory>\n");
		fprintf(stderr, "\t\t<output id>\n");
		fprintf(stderr, "\t\t<tmp file directory>\n");
		fprintf(stderr, "\t\t<number of threads>\n");
	}

	return 0;
}

void PrintDistribution(RGIndex *index, 
		RGBinary *rg,
		char *distributionFileName,
		int numMismatches,
		char *tmpDir,
		int numThreads)
{
	char *FnName = "PrintDistribution";
	FILE *fp;
	int64_t startIndex = 0;
	int64_t endIndex = index->length-1;
	int64_t curIndex=0, nextIndex=0;
	int64_t counter=0;
	int64_t numDifferent = 0;
	int64_t numForward, numReverse;
	char *read=NULL;
	char *reverseRead=NULL;
	int64_t i, j;
	char **reads=NULL;
	int64_t *readCounts=NULL;
	int64_t numReads=0;
	pthread_t *threads=NULL;
	ThreadData *data=NULL;
	int errCode;
	void *status;

	fprintf(stderr, "Out of %lld, currently on:\n0",
			(long long int)(endIndex - startIndex + 1));
	/* Go through every possible read in the genome using the index */
	for(curIndex=startIndex, nextIndex=startIndex, counter=0, numDifferent=0;
			curIndex <= endIndex;
			curIndex = nextIndex) {
		if(counter >= BINDEXDIST_ROTATE_NUM) {
			fprintf(stderr, "\r%10lld", 
					(long long int)(curIndex-startIndex));
			counter -= BINDEXDIST_ROTATE_NUM;
		}
		/* Get the matches for the contig/pos */
		GetMatchesFromContigPos(index,
				rg,
				(index->contigType==Contig_8)?(index->contigs_8[curIndex]):(index->contigs_32[curIndex]),
				index->positions[curIndex],
				numMismatches,
				&numForward, 
				&numReverse,
				&read,
				&reverseRead);
		assert(numForward + numReverse> 0);

		nextIndex += numForward;
		counter += numForward;
		/* In case reverse is zero */
		if(numForward <= 0) {
			nextIndex++;
			counter++;
		}

		/* Reallocate memory */
		numReads+=2; /* One for both strands */
		reads = realloc(reads, sizeof(char*)*numReads);
		if(NULL==reads) {
			PrintError(FnName,
					"reads",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		reads[numReads-2] = malloc(sizeof(char)*(index->width+1));
		if(NULL==reads[numReads-2]) {
			PrintError(FnName,
					"reads[numReads-1]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		reads[numReads-1] = malloc(sizeof(char)*(index->width+1));
		if(NULL==reads[numReads-1]) {
			PrintError(FnName,
					"reads[numReads-1]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		readCounts = realloc(readCounts, sizeof(int64_t)*numReads);
		if(NULL==readCounts) {
			PrintError(FnName,
					"readCounts",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Copy over */
		assert(strlen(read) < SEQUENCE_LENGTH);
		assert(strlen(reverseRead) < SEQUENCE_LENGTH);
		strcpy(reads[numReads-2], read);
		strcpy(reads[numReads-1], reverseRead);
		readCounts[numReads-1] = numForward+numReverse;
		readCounts[numReads-2] = numForward+numReverse;

		/* Free memory */
		free(read);
		read = NULL;
		free(reverseRead);
		reverseRead = NULL;
	}
	fprintf(stderr, "\r%10lld\n", 
			(long long int)(curIndex-startIndex+1));

	/* Allocate memory for threads */
	threads=malloc(sizeof(pthread_t)*numThreads);
	if(NULL==threads) {
		PrintError(FnName,
				"threads",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory to pass data to threads */
	data=malloc(sizeof(ThreadData)*numThreads);
	if(NULL==data) {
		PrintError(FnName,
				"data",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	for(i=0;i<numThreads;i++) {
		data[i].reads = reads;
		data[i].readCounts = readCounts;
		data[i].low = i*(numReads/numThreads);
		data[i].high = (i+1)*(numReads/numThreads)-1;
		data[i].tmpDir = tmpDir;
		data[i].readLength = index->width;
		data[i].showPercentComplete = 0;
	}
	data[0].low = 0;
	data[numThreads-1].high = numReads-1;
	data[numThreads-1].showPercentComplete = 1;

	/* Open threads */
	for(i=0;i<numThreads;i++) {
		/* Start thread */
		errCode = pthread_create(&threads[i], /* thread struct */
				NULL, /* default thread attributes */
				MergeSortReads, /* start routine */
				&data[i]); /* data to routine */
		if(0!=errCode) {
			PrintError(FnName,
					"pthread_create: errCode",
					"Could not start thread",
					Exit,
					ThreadError);
		}
	}
	/* Wait for threads to return */
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

	/* Merge results from the sorts */
	fprintf(stderr, "\rMerging sorts from threads...    \n");
	for(j=1;j<numThreads;j=j*2) {
		for(i=0;i<numThreads;i+=2*j) {
			MergeHelper(reads,
					readCounts,
					data[i].low,
					data[i+j].low-1,
					data[i+2*j-1].high,
					tmpDir,
					index->width);
		}
	}
	fprintf(stderr, "Sorting complete.\n");
	fprintf(stderr, "%s", BREAK_LINE);

	/* Open the output file */
	fprintf(stderr, "Outputting to %s.\n",
			distributionFileName);
	if(!(fp = fopen(distributionFileName, "w"))) {
		PrintError(FnName,
				distributionFileName,
				"Could not open file for writing",
				Exit,
				OpenFileError);
	}

	/* Print */
	for(i=0;i<numReads;i++) {
		fprintf(fp, "%s\t%lld\n",
				reads[i],
				(long long int)readCounts[i]);
	}

	/* Close the file */
	fclose(fp);

	/* Free memory */
	for(i=0;i<numReads;i++) {
		free(reads[i]);
		reads[i]=NULL;
	}
	free(reads);
	reads=NULL;
	free(readCounts);
	readCounts=NULL;
}

/* Get the matches for the contig/pos */
void GetMatchesFromContigPos(RGIndex *index,
		RGBinary *rg,
		uint32_t curContig,
		uint32_t curPos,
		int numMismatches,
		int64_t *numForward,
		int64_t *numReverse, 
		char **read,
		char **reverseRead)
{
	char *FnName = "GetMatchesFromContigPos";
	int readLength = index->width;
	int returnLength, returnPosition;
	int i;
	RGReads reads;
	RGRanges ranges;

	/* Initialiez reads */
	RGReadsInitialize(&reads);
	RGRangesInitialize(&ranges);

	/* Get the read */
	RGBinaryGetReference(rg,
			curContig,
			curPos,
			FORWARD,
			0,
			read,
			readLength,
			&returnLength,
			&returnPosition);
	assert(returnLength == readLength);
	assert(returnPosition == curPos);

	/* First generate the perfect match for the forward and
	 * reverse strand */

	reverseRead = malloc(sizeof(char)*(returnLength+1));
	if(NULL==reverseRead) {
		PrintError(FnName,
				"reverseRead",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	GetReverseComplimentAnyCase((*read),
			(*reverseRead),
			readLength);

	RGReadsGeneratePerfectMatch((*read),
			readLength,
			FORWARD,
			0,
			index,
			&reads);
	RGReadsGeneratePerfectMatch((*reverseRead),
			readLength,
			REVERSE,
			0,
			index,
			&reads);

	if(numMismatches > 0) {
		/* Generate reads with the necessary mismatches for 
		 *          * both the forward and reverse strands */
		RGReadsGenerateMismatches((*reverseRead),
				readLength,
				REVERSE,
				0,
				numMismatches,
				index,
				&reads);
		RGReadsGenerateMismatches((*read),
				readLength,
				FORWARD,
				0,
				numMismatches,
				index,
				&reads);
	}

	for(i=0;i<reads.numReads;i++) {
		/* Get the matches for the read */
		RGIndexGetRanges(index,
				rg,
				reads.reads[i],
				reads.readLength[i],
				reads.strand[i],
				reads.offset[i],
				&ranges);
	}

	/* Remove duplicates */
	RGRangesRemoveDuplicates(&ranges);

	/* Error check */
	if(ranges.numEntries <= 0) {
		PrintError(FnName,
				"ranges",
				"Returned zero ranges",
				Exit,
				OutOfRange);
	}

	/* Return the number of FORWARD strand matches so that we can skip over */
	(*numForward) = (*numReverse) = 0;
	for(i=0;i<ranges.numEntries;i++) {
		switch(ranges.strand[i]) {
			case FORWARD:
				(*numForward) += ranges.endIndex[i] - ranges.startIndex[i] + 1;
				break;
			case REVERSE:
				(*numReverse) += ranges.endIndex[i] - ranges.startIndex[i] + 1;
				break;
			default:
				PrintError(FnName,
						"m->strand[i]",
						"Could not understand strand",
						Exit,
						OutOfRange);
				break;
		}
	}

	/* Null out the gaps */
	for(i=0;i<index->width;i++) {
		switch(index->mask[i]) {
			case 0:
				(*read)[i] = 'N';
				break;
			case 1:
				break;
			default:
				PrintError(FnName,
						"index->masks[i]",
						"Could not understand mask",
						Exit,
						OutOfRange);
		}
	}

	/* Free memory */
	RGReadsFree(&reads);
	RGRangesFree(&ranges);
}

void *MergeSortReads(void *arg)
{
	ThreadData *data = (ThreadData*)arg;
	char **reads = data->reads;
	int64_t *readCounts = data->readCounts;
	int64_t low = data->low;
	int64_t high = data->high;
	char *tmpDir = data->tmpDir;
	int readLength = data->readLength;

	double curPercentComplete = 0.0;
	/* Call helper */
	if(data->showPercentComplete == 1) {
		fprintf(stderr, "\r%3.2lf percent complete", 0.0);
	}

	MergeSortReadsHelper(reads,
			readCounts,
			low,
			high,
			low,
			high - low,
			data->showPercentComplete,
			&curPercentComplete,
			tmpDir,
			readLength);
	if(data->showPercentComplete == 1) {
		fprintf(stderr, "\r");
		fprintf(stderr, "thread %3.2lf percent complete", 100.0);
	}

	return NULL;
}

void MergeSortReadsHelper(char **reads,
		int64_t *readCounts,
		int64_t low,
		int64_t high,
		int64_t startLow,
		int64_t total,
		int showPercentComplete,
		double *curPercentComplete,
		char *tmpDir,
		int readLength)
{
	int64_t mid = (low + high)/2;
	if(low >= high) {

		if(showPercentComplete == 1) {
			assert(NULL!=curPercentComplete);
			if((*curPercentComplete) < 100.0*((double)(low - startLow))/total) {
				while((*curPercentComplete) < 100.0*((double)(low - startLow))/total) {
					(*curPercentComplete) += BINDEXDIST_SORT_ROTATE_INC;
				}
				PrintPercentCompleteShort((*curPercentComplete));
			}
		}
		return;
	}

	/* Sort recursively */
	MergeSortReadsHelper(reads,
			readCounts,
			low,
			mid,
			startLow, 
			total,
			showPercentComplete,
			curPercentComplete,
			tmpDir,
			readLength);
	MergeSortReadsHelper(reads,
			readCounts,
			mid+1,
			high,
			startLow, 
			total,
			showPercentComplete,
			curPercentComplete,
			tmpDir,
			readLength);

	/* Merge the two lists */
	MergeHelper(reads,
			readCounts,
			low,
			mid,
			high,
			tmpDir,
			readLength);
}

void MergeHelper(char **reads,
		int64_t *readCounts,
		int64_t low,
		int64_t mid,
		int64_t high,
		char *tmpDir,
		int readLength)
{
	char *FnName = "MergeHelper";
	int64_t i=0;
	char **tmpReads=NULL;
	int64_t *tmpReadCounts=NULL;
	int64_t startUpper, startLower, endUpper, endLower;
	int64_t ctr=0;
	FILE *tmpLowerFP=NULL;
	FILE *tmpUpperFP=NULL;
	char *tmpLowerFileName=NULL;
	char *tmpUpperFileName=NULL;
	char tmpLowerRead[SEQUENCE_LENGTH]="\0";
	char tmpUpperRead[SEQUENCE_LENGTH]="\0";
	long long int tmpLowerReadCount=0;
	long long int tmpUpperReadCount=0;
	int eofLower, eofUpper;

	/* Merge the two lists */
	/* Since we want to keep space requirement small, use an upper bound on memory,
	 * so that we use tmp files when memory requirements become to large */
	if( (high-low+1)*(sizeof(int64_t) + sizeof(char*)) <= ONE_GIGABYTE) {

		/* Use memory */
		tmpReads = malloc(sizeof(char*)*(high-low+1));
		if(NULL == tmpReads) {
			PrintError(FnName,
					"tmpReads",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		tmpReadCounts = malloc(sizeof(int64_t)*(high-low+1));
		if(NULL == tmpReadCounts) {
			PrintError(FnName,
					"tmpReadCounts",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}

		/* Merge */
		startLower = low;
		endLower = mid;
		startUpper = mid+1;
		endUpper = high;
		ctr=0;
		while( (startLower <= endLower) && (startUpper <= endUpper) ) {
			if(strcmp(reads[startLower], reads[startUpper]) <= 0) {
				tmpReads[ctr] = reads[startLower];
				tmpReadCounts[ctr] = readCounts[startLower];
				startLower++;
			}
			else {
				tmpReads[ctr] = reads[startUpper];
				tmpReadCounts[ctr] = readCounts[startUpper];
				startUpper++;
			}
			ctr++;
		}
		while(startLower <= endLower) {
			tmpReads[ctr] = reads[startLower];
			tmpReadCounts[ctr] = readCounts[startLower];
			startLower++;
			ctr++;
		}
		while(startUpper <= endUpper) {
			tmpReads[ctr] = reads[startUpper];
			tmpReadCounts[ctr] = readCounts[startUpper];
			startUpper++;
			ctr++;
		}
		/* Copy back */
		for(i=low, ctr=0;
				i<=high;
				i++, ctr++) {
			reads[i] = tmpReads[ctr];
			readCounts[i] = tmpReadCounts[ctr];
		}

		/* Free memory */
		free(tmpReads);
		tmpReads=NULL;
		free(tmpReadCounts);
		tmpReadCounts=NULL;
	}
	else {
		/* Use tmp files */
		assert(sizeof(int64_t) == sizeof(long long int));

		/* Open tmp files */
		tmpLowerFP = OpenTmpFile(tmpDir, &tmpLowerFileName);
		tmpUpperFP = OpenTmpFile(tmpDir, &tmpUpperFileName);

		/* Print to tmp files */
		for(i=low;i<=mid;i++) {
			if(0 > fprintf(tmpLowerFP, "%s\t%lld\n", 
						reads[i],
						(long long int)readCounts[i])) { 
				PrintError(FnName,
						NULL,
						"Could not write to tmp lower file",
						Exit,
						WriteFileError);
			}
		}
		for(i=mid+1;i<=high;i++) {
			if(0 > fprintf(tmpUpperFP, "%s\t%lld\n",
						reads[i],
						(long long int)readCounts[i])) {
				PrintError(FnName,
						NULL,
						"Could not write to tmp upper file",
						Exit,
						WriteFileError);
			}
		}

		/* Move to beginning of the files */
		fseek(tmpLowerFP, 0 , SEEK_SET);
		fseek(tmpUpperFP, 0 , SEEK_SET);

		/* Merge tmp files back into index */
		/* Get first contig/pos */
		if(0 > fscanf(tmpLowerFP, "%s %lld\n",
					tmpLowerRead,
					&tmpLowerReadCount)) {
			PrintError(FnName,
					NULL,
					"Could not read in tmp lower",
					Exit,
					ReadFileError);
		}
		if(0 > fscanf(tmpUpperFP, "%s %lld\n",
					tmpUpperRead,
					&tmpUpperReadCount)) {
			PrintError(FnName,
					NULL,
					"Could not read in tmp upper",
					Exit,
					ReadFileError);
		}
		for(i=low, ctr=0, eofLower = 0, eofUpper = 0;
				i<=high &&
				eofLower == 0 &&
				eofUpper == 0;
				i++, ctr++) {
			if(strcmp(tmpLowerRead, tmpUpperRead) <= 0) {
				/* Copy lower */
				strcpy(reads[i], tmpLowerRead);
				readCounts[i] = tmpLowerReadCount;
				/* Get new tmpLower */
				if(0 > fscanf(tmpLowerFP, "%s %lld\n",
							tmpLowerRead,
							&tmpLowerReadCount)) {
					eofLower = 1;
				}
			}
			else {
				/* Copy upper */
				strcpy(reads[i], tmpUpperRead);
				readCounts[i] = tmpUpperReadCount;
				/* Get new tmpUpper */
				if(0 > fscanf(tmpUpperFP, "%s %lld\n",
							tmpUpperRead,
							&tmpUpperReadCount)) {
					eofUpper = 1;
				}
			}
		}
		while(eofLower != 1) {
			/* Copy lower */
			strcpy(reads[i], tmpLowerRead);
			readCounts[i] = tmpLowerReadCount;
			/* Get new tmpLower */
			if(0 > fscanf(tmpLowerFP, "%s %lld\n",
						tmpLowerRead,
						&tmpLowerReadCount)) {
				eofLower = 1;
			}
			i++;
			ctr++;
		}
		while(eofUpper != 1) {
			/* Copy upper */
			strcpy(reads[i], tmpUpperRead);
			readCounts[i] = tmpUpperReadCount;
			/* Get new tmpUpper */
			if(0 > fscanf(tmpUpperFP, "%s %lld\n",
						tmpUpperRead,
						&tmpUpperReadCount)) {
				eofUpper = 1;
			}
			i++;
			ctr++;
		}
		assert(ctr == (high - low + 1));
		assert(i == high + 1);

		/* Close tmp files */
		CloseTmpFile(&tmpLowerFP, &tmpLowerFileName);
		CloseTmpFile(&tmpUpperFP, &tmpUpperFileName);
	}
}
