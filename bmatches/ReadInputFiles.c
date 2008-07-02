#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "../blib/RGMatch.h"
#include "Definitions.h"
#include "ReadInputFiles.h"


/* TODO */
/* Read the next read from the stream */
int GetNextRead(FILE *fp, char **read, int *readLength, char **pairedRead, int *pairedReadLength,  char**readName, int pairedEnd)
{
	if(pairedEnd == 1) {
		if(EOF==fscanf(fp, "%s", (*readName)) || EOF==fscanf(fp, "%s", (*read)) || EOF==fscanf(fp, "%s", (*pairedRead))) {
			return EOF;
		}
		(*readLength) = strlen((*read));
		(*pairedReadLength) = strlen((*pairedRead));
	}
	else {
		if(EOF==fscanf(fp, "%s", (*readName)) || EOF==fscanf(fp, "%s", (*read))) {
			return EOF;
		}
		(*readLength) = strlen((*read));
		(*pairedRead)[0]='\0';
		(*pairedReadLength)=0;
	}
	return 1;
}

/* TODO */
int WriteRead(FILE *fp, char *readName, char *read, char *pairedRead, int pairedEnd)
{
	/* Print read */
	if(fprintf(fp, "%s\n", readName) < 0) {
		return EOF;
	}
	/* Print read */
	if(fprintf(fp, "%s\n", read) < 0) {
		return EOF;
	}
	/* Print paired read */
	if(pairedEnd==1) {
		if(fprintf(fp, "%s\n", pairedRead) < 0) {
			return EOF;
		}
	}
	return 1;
}

/* TODO */
void WriteReadsToTempFile(FILE *seqFP,
		FILE ***tempSeqFPs, /* One for each thread */ 
		char ***tempSeqFileNames,
		int startReadNum, 
		int endReadNum, 
		int pairedEnd,
		int numThreads,
		char *tmpDir,
		int timing,
		int *numWritten,
		int *numFiltered)
{
	char *FnName = "WriteReadsToTempFile";
	int i;
	int curSeqFPIndex=0;
	char *readName=NULL;
	char *read=NULL;
	char *pairedRead=NULL;
	int readLength;
	int pairedReadLength;
	int curReadNum = 1;
	time_t startTime = time(NULL);
	time_t endTime;

	(*numFiltered)=0;
	(*numWritten)=0;

	/* Allocate memory for the read name, read and paired read */
	readName = malloc(sizeof(char)*SEQUENCE_NAME_LENGTH);
	if(NULL == readName) {
		PrintError(FnName,
				"readName",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	read = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL == read) {
		PrintError(FnName,
				"read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	pairedRead = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL == pairedRead) {
		PrintError(FnName,
				"pairedRead",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
			

	/* open one temporary file, one for each thread */
	for(i=0;i<numThreads;i++) {
		assert((*tempSeqFPs)[i] == NULL);
		(*tempSeqFPs)[i] = OpenTmpFile(tmpDir, &(*tempSeqFileNames)[i]);
	}

	while((endReadNum<=0 || endReadNum >= curReadNum) && 
			EOF != GetNextRead(seqFP, &read, &readLength, &pairedRead, &pairedReadLength, &readName, pairedEnd)) {
		/* Get which tmp file to put the read in */
		curSeqFPIndex = (curReadNum-1)%numThreads;

		/* Print only if we are within the desired limit and the read checks out */
		if( (startReadNum<=0 || curReadNum >= startReadNum)  /* Only if we are within the bounds for the reads */
				&& (1 == UpdateRead(read, readLength)) /* The read is valid */
				&& (0 == pairedEnd || 1 == UpdateRead(pairedRead, pairedReadLength))) { /* The paired read is valid */
			if(EOF == WriteRead((*tempSeqFPs)[curSeqFPIndex], readName, read, pairedRead, pairedEnd)) {
				PrintError(FnName,
						NULL,
						"Could not write read",
						Exit,
						WriteFileError);
			}
			(*numWritten)++;
		}
		else {
			(*numFiltered)++;
		}
		/* Increment read number */
		curReadNum++;
	}

	/* reset pointer to temp files to the beginning of the file */
	for(i=0;i<numThreads;i++) {
		fseek((*tempSeqFPs)[i], 0, SEEK_SET);
	}

	/* Free memory */
	free(readName);
	readName=NULL;
	free(read);
	readName=NULL;
	free(pairedRead);
	pairedRead=NULL;

	endTime = time(NULL);
	int seconds = endTime - startTime;
	int hours = seconds/3600;
	seconds -= hours*3600;
	int minutes = seconds/60;
	seconds -= minutes*60;

	if(timing == 1) {
		fprintf(stderr, "Reading to temp file took: %d hours, %d minutes and %d seconds.\n",
				hours,
				minutes,
				seconds
			   );
	}
}

/* TODO */
/* Go through the temporary output file and output those reads that have 
 * at least one match to the final output file.  For those reads that have
 * zero matches, output them to the temporary read file *
 * */
int ReadTempReadsAndOutput(FILE *tempOutputFP,
		FILE *outputFP,
		FILE *tempSeqFP,
		int pairedEnd,
		int binaryOutput)
{
	char *FnName = "ReadTempReadsAndOutput";
	char readName[SEQUENCE_NAME_LENGTH];
	char read[SEQUENCE_LENGTH];
	char pairedRead[SEQUENCE_LENGTH];
	RGMatch readMatch;
	RGMatch pairedReadMatch;
	int numReads = 0;
	int numOutputted=0;

	/* Initialize match structures */
	RGMatchInitialize(&readMatch);
	RGMatchInitialize(&pairedReadMatch);

	/* Go to the beginning of the temporary output file */
	fseek(tempOutputFP, 0, SEEK_SET);

	while(RGMatchRead(tempOutputFP, 
				readName, 
				read, 
				pairedRead, 
				&readMatch, 
				&pairedReadMatch,
				pairedEnd,
				binaryOutput)!=EOF) {
		/* Output if either pair has more than one entry */
		if(readMatch.numEntries > 0 ||
				pairedReadMatch.numEntries > 0) {
			/* Output to final output file */
			RGMatchPrint(outputFP,
					readName,
					read,
					pairedRead,
					&readMatch,
					&pairedReadMatch,
					pairedEnd,
					binaryOutput);
			numOutputted++;
		}
		else {
			/* Put back in the read file */
			if(EOF == WriteRead(tempSeqFP, readName, read, pairedRead, pairedEnd)) {
				PrintError(FnName,
						NULL,
						"Could not write read.",
						Exit,
						WriteFileError);
			}
			numReads++;
		}

		/* Free match memory and reinitialize match structures */
		RGMatchFree(&readMatch);
		RGMatchFree(&pairedReadMatch);
	}
	return numReads;
}

/* TODO */
/* Reads in a RGIndexfrom file */
void ReadRGIndex(char *rgIndexFileName, RGIndex *index, int binaryInput)
{
	FILE *fp;

	if(VERBOSE >= 0) {
		fprintf(stderr, "Reading index from %s.\n",
				rgIndexFileName);
	}

	/* open file */
	if((fp=fopen(rgIndexFileName, "r"))==0) {
		PrintError("ReadRGIndex",
				rgIndexFileName,
				"Could not open rgIndexFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Read from file */
	RGIndexRead(fp, index, binaryInput);

	/* close file */
	fclose(fp);

	if(VERBOSE >= 0) {
		fprintf(stderr, "Read index from %s.\n",
				rgIndexFileName);
	}
}

/* TODO */
int ReadFileNames(char *listFileName, char ***fileNames)
{
	char tempFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fp;
	int numFileNames=0;

	if(VERBOSE>0) {
		fprintf(stderr, "Reading in file names from %s.\n",
				listFileName);
	}

	/* open file */
	if((fp=fopen(listFileName, "r"))==0) {
		PrintError("ReadFileNames",
				listFileName,
				"Could not open listFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Read in the file names */
	while(fscanf(fp, "%s", tempFileName)!=EOF) {
		numFileNames++;
		(*fileNames) = (char**)realloc((*fileNames), sizeof(char*)*numFileNames);
		(*fileNames)[numFileNames-1] = (char*)malloc(sizeof(char)*MAX_FILENAME_LENGTH);
		strcpy((*fileNames)[numFileNames-1], tempFileName);
		if(VERBOSE>0) {
			fprintf(stderr, "file name %d:%s\n", 
					numFileNames,
					(*fileNames)[numFileNames-1]);
		}
	}

	/* close file */
	fclose(fp);

	if(VERBOSE>0) {
		fprintf(stderr, "Read %d file names from %s.\n", numFileNames, listFileName);
	}

	return numFileNames;
}

/* TODO */
int ReadOffsets(char *offsetsFileName, int **offsets) 
{
	FILE *fp;
	int numOffsets=0;
	int tempInt;

	if(VERBOSE>0) {
		fprintf(stderr, "Reading in offsets from %s.\n",
				offsetsFileName);
	}

	/* open file */
	if((fp=fopen(offsetsFileName, "r"))==0) {
		PrintError("ReadOffsets",
				offsetsFileName,
				"Could not open offsetsFileName for reading",
				Exit,
				OpenFileError);
	}

	numOffsets=0;
	while(fscanf(fp, "%d", &tempInt)!=EOF) {
		numOffsets++;
		(*offsets)=(int*)realloc((*offsets), sizeof(int)*numOffsets);
		(*offsets)[numOffsets-1] = tempInt;
		if(VERBOSE>0) {
			fprintf(stderr, "Offset %d:%d\n", numOffsets, (*offsets)[numOffsets-1]);
		}
	}
	assert(numOffsets>0);

	/* close file */
	fclose(fp);

	if(VERBOSE>0) {
		fprintf(stderr, "Read %d offsets from %s.\n", numOffsets, offsetsFileName);
	}

	return numOffsets;
}
