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
#define BINDEXDIST_ROTATE_NUM 1000000

int main(int argc, char *argv[]) 
{
	FILE *fp=NULL;
	char indexFileName[MAX_FILENAME_LENGTH]="\0";
	char distributionFileName[MAX_FILENAME_LENGTH]="\0";
	char rgFileName[MAX_FILENAME_LENGTH]="\0";
	char tmpDir[MAX_FILENAME_LENGTH]="\0";
	int numMismatches = 0;

	if(argc == 5) {
		RGBinary rg;
		RGIndex index;

		strcpy(rgFileName, argv[1]);
		strcpy(indexFileName, argv[2]);
		numMismatches = atoi(argv[3]);
		strcpy(tmpDir, argv[4]);

		/* Create the distribution file name */
		strcpy(distributionFileName, indexFileName);
		strcat(distributionFileName, ".dist");

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
				tmpDir);
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
		fprintf(stderr, "\t\t<reference genome file name>\n");
		fprintf(stderr, "\t\t<index file name>\n");
		fprintf(stderr, "\t\t<number of mismatches>\n");
		fprintf(stderr, "\t\t<tmp file directory>\n");
	}

	return 0;
}

void PrintDistribution(RGIndex *index, 
		RGBinary *rg,
		char *distributionFileName,
		int numMismatches,
		char *tmpDir)
{
	char *FnName = "PrintDistribution";
	FILE *fp;
	int64_t startIndex = 0;
	int64_t endIndex = index->length-1;
	int64_t curIndex=0, nextIndex=0;
	int64_t counter=0;
	int64_t numDifferent = 0;
	int64_t numForward, numReverse;
	char read[SEQUENCE_LENGTH] = "\0";
	char reverseRead[SEQUENCE_LENGTH] = "\0";
	int numTmpFiles=2;
	int i;
	char **tmpFileNames;
	FILE **tmpFPs;

	/* Allocate memory for temporary file pointer and names */
	tmpFileNames = malloc(sizeof(char*)*numTmpFiles);
	if(NULL==tmpFileNames) {
		PrintError(FnName,
				"tmpFileNames",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	tmpFPs = malloc(sizeof(FILE*)*numTmpFiles);
	if(NULL==tmpFPs) {
		PrintError(FnName,
				"tmpFPs",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Open temporary files */
	for(i=0;i<numTmpFiles;i++) {
		tmpFPs[i] = OpenTmpFile(tmpDir,
				&tmpFileNames[i]);
	}

	/* Go through every possible read in the genome using the index */
	for(curIndex=startIndex, nextIndex=startIndex, counter=0, numDifferent=0;
			curIndex <= endIndex;
			curIndex = nextIndex) {
		if(counter >= BINDEXDIST_ROTATE_NUM) {
			fprintf(stderr, "\r%10lld", 
					(long long int)(curIndex-startIndex));
			counter -= BINDEXDIST_ROTATE_NUM;
		}
		/* Get the matches for the chr/pos */
		GetMatchesFromChrPos(index,
				rg,
				index->chromosomes[curIndex],
				index->positions[curIndex],
				numMismatches,
				&numForward, 
				&numReverse,
				read,
				reverseRead);
		assert(numForward > 0);

		nextIndex += numForward;
		counter += numForward;

		/* Print the number of matches to file */
		assert(numTmpFiles == 2); /* One for + strand, one for - strand */
		fprintf(tmpFPs[0], "%s\t%lld\n",
				read,
				(long long int)(numForward+numReverse));
		fprintf(tmpFPs[1], "%s\t%lld\n",
				reverseRead,
				(long long int)(numForward+numReverse));
	}
	fprintf(stderr, "\r%10lld\n", 
			(long long int)(curIndex-startIndex));

	/* Reverse the "reverse" strand file */
	ReverseFile(&tmpFPs[1], &tmpFileNames[1], tmpDir);

	/* Open the output file */
	if(!(fp = fopen(distributionFileName, "w"))) {
		PrintError(FnName,
				distributionFileName,
				"Could not open file for writing",
				Exit,
				OpenFileError);
	}

	/* Merge the files in to the output */
	MergeFiles(tmpFPs[0],
			tmpFPs[1],
			fp);

	/* Close the file */
	fclose(fp);

	/* Close temporary files */
	for(i=0;i<numTmpFiles;i++) {
		CloseTmpFile(&tmpFPs[i],
				&tmpFileNames[i]);
	}
}

/* Get the matches for the chr/pos */
void GetMatchesFromChrPos(RGIndex *index,
		RGBinary *rg,
		uint32_t curChr,
		uint32_t curPos,
		int numMismatches,
		int64_t *numForward,
		int64_t *numReverse, 
		char *read,
		char *reverseRead)
{
	char *FnName = "GetMatchesFromChrPos";
	int readLength = index->totalLength;
	int returnLength, returnPosition;
	int i;
	RGReads reads;
	RGRanges ranges;

	/* Initialiez reads */
	RGReadsInitialize(&reads);
	RGRangesInitialize(&ranges);

	/* Get the read */
	RGBinaryGetSequence(rg,
			curChr,
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
	GetReverseComplimentAnyCase(read,
			reverseRead,
			readLength);
	RGReadsGeneratePerfectMatch(read,
			readLength,
			FORWARD,
			0,
			index->numTiles,
			index->tileLengths,
			index->gaps,
			index->totalLength,
			&reads);
	RGReadsGeneratePerfectMatch(reverseRead,
			readLength,
			REVERSE,
			0,
			index->numTiles,
			index->tileLengths,
			index->gaps,
			index->totalLength,
			&reads);

	if(numMismatches > 0) {
		/* Generate reads with the necessary mismatches for 
		 *          * both the forward and reverse strands */
		RGReadsGenerateMismatches(reverseRead,
				readLength,
				REVERSE,
				0,
				index->numTiles,
				index->tileLengths,
				index->gaps,
				index->totalLength,
				numMismatches,
				&reads);
		RGReadsGenerateMismatches(read,
				readLength,
				FORWARD,
				0,
				index->numTiles,
				index->tileLengths,
				index->gaps,
				index->totalLength,
				numMismatches,
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
	assert((*numForward)>0);

	/* Free memory */
	RGReadsFree(&reads);
	RGRangesFree(&ranges);
}

void ReverseFile(FILE **fp,
		char **fileName,
		char *tmpDir)
{
	char *FnName="ReverseFile";
	FILE *tmpFP=NULL;
	char *tmpFileName=NULL;
	char read[SEQUENCE_LENGTH]="\0";
	long long int count=0;

	/* Move to the beginning of the file to be read */
	fseek((*fp), 0, SEEK_SET);

	/* Open a tmp file */
	tmpFP=OpenTmpFile(tmpDir,
			&tmpFileName);

	/* Write to tmp file */
	while(0 < fscanf((*fp), 
				"%s %lld",
				read,
				&count)) {
		if(0 > fprintf(tmpFP, "%s\t%lld\n",
					read,
					count)) {
			PrintError(FnName,
					NULL,
					"Could not write read or count",
					Exit,
					WriteFileError);
		}
	}

	/* Close the tmp file that we read in from */
	CloseTmpFile(fp, fileName); 

	/* Now adjust memory pointers */
	(*fp) = tmpFP;
	(*fileName) = tmpFileName;
	tmpFP = NULL;
	tmpFileName = NULL;
}

void MergeFiles(FILE *fpIn1,
		FILE *fpIn2,
		FILE *fpOut)
{
	char read1[SEQUENCE_LENGTH]="\0";
	char read2[SEQUENCE_LENGTH]="\0";
	long long int count1, count2;
	int eof1=0, eof2=0;
	int getCount1=1, getCount2=1;

	/* Move to the beginning of the files */
	fseek(fpIn1, 0, SEEK_SET);
	fseek(fpIn2, 0, SEEK_SET);

	while(eof1 != 1 || 
			eof2 != 1) {
		/* Get data */
		if(getCount1 == 1 && eof1 != 1) {
			if(EOF == fscanf(fpIn1, 
						"%s %lld",
						read1,
						&count1)) {
				eof1 = 1;
			}
		}
		if(getCount2 == 1 && eof2 != 1) {
			if(EOF == fscanf(fpIn2, 
						"%s %lld",
						read2,
						&count2)) {
				eof2 = 1;
			}
		}
		if(eof1 != 1 && eof2 != 1) {
			/* Compare */
			int cmp = strcmp(read1, read2);
			if(cmp == 0) {
				fprintf(fpOut, "%s\t%lld\n",
						read1,
						(count1 + count2));
				getCount1 = getCount2 = 1;
			}
			else if(cmp < 0) {
				fprintf(fpOut, "%s\t%lld\n",
						read1,
						count1); 
				getCount1 = 1;
				getCount2 = 0;
			}
			else {
				fprintf(fpOut, "%s\t%lld\n",
						read2,
						count2); 
				getCount1 = 0;
				getCount2 = 1;
			}
		}
		else if(eof1 != 1 && eof2 != 1) {
			fprintf(fpOut, "%s\t%lld\n",
					read1,
					count1); 
			getCount1 = 1;
			getCount2 = 0;
		}
		else if(eof1 == 1 && eof2 != 1) {
			fprintf(fpOut, "%s\t%lld\n",
					read2,
					count2); 
			getCount1 = 0;
			getCount2 = 1;
		}
	}
}
