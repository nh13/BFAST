#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "../blib/RGMatch.h"
#include "../blib/RGMatches.h"
#include "Definitions.h"
#include "ReadInputFiles.h"

/* TODO */
/* Read the next read from the stream */
int GetNextRead(FILE *fp, 
		RGMatches *m,
		int pairedEnd)
{
	char *FnName = "GetNextRead";
	char readName[SEQUENCE_NAME_LENGTH]="\0";
	char readOne[SEQUENCE_LENGTH]="\0";
	char readTwo[SEQUENCE_LENGTH]="\0";

	m->pairedEnd = pairedEnd;

	/* Read in read name and read */
	if(NULL==fgets(readName, SEQUENCE_NAME_LENGTH-1, fp)) {
		return EOF;
	}
	if(NULL==fgets(readOne, SEQUENCE_LENGTH-1, fp)) {
		PrintError(FnName,
				"readOne",
				"Could not read from file",
				Exit,
				ReadFileError);
	}

	/* Trim leading and following whitespaces */
	m->readNameLength = StringTrimWhiteSpace(readName);
	m->matchOne.readLength = StringTrimWhiteSpace(readOne);

	/* Allocate memory */
	m->readName = malloc(sizeof(int8_t)*(m->readNameLength+1));
	if(NULL == m->readName) {
		PrintError(FnName,
				"m->readName",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	m->matchOne.read = malloc(sizeof(int8_t)*(m->matchOne.readLength+1));
	if(NULL == m->matchOne.read) {
		PrintError(FnName,
				"m->matchOne.read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Copy over */
	strcpy((char*)m->readName, readName);
	strcpy((char*)m->matchOne.read, readOne);

	/* Read in paired end if necessary */
	if(pairedEnd == 1) {
		/* Read in read number two */
		if(NULL==fgets(readTwo, SEQUENCE_LENGTH-1, fp)) {
			PrintError(FnName,
					"readTwo",
					"Could not read from file",
					Exit,
					ReadFileError);
		}
		/* Copy over read two length as well as adjusting the 
		 * read to get rid of the "\n" character. */
		m->matchTwo.readLength = strlen(readTwo)-1;
		readTwo[m->matchTwo.readLength]='\0';
		/* Allocate memory */
		m->matchTwo.read = malloc(sizeof(int8_t)*(m->matchTwo.readLength+1));
		if(NULL == m->matchTwo.read) {
			PrintError(FnName,
					"m->matchTwo.read",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Copy over */
		strcpy((char*)m->matchTwo.read, readTwo);
	}
	return 1;
}

/* TODO */
int WriteRead(FILE *fp, RGMatches *m, int pairedEnd) 
{
	assert(m->pairedEnd == pairedEnd);

	/* Print read */
	if(fprintf(fp, "%s\n", m->readName) < 0) {
		return EOF;
	}
	/* Print read */
	if(fprintf(fp, "%s\n", m->matchOne.read) < 0) {
		return EOF;
	}
	/* Print paired read */
	if(pairedEnd==1) {
		if(fprintf(fp, "%s\n", m->matchTwo.read) < 0) {
			return EOF;
		}
	}
	return 1;
}

/* TODO */
void WriteReadsToTempFile(FILE *seqFP,
		FILE *seqFilteredFP,
		FILE ***tempSeqFPs, /* One for each thread */ 
		char ***tempSeqFileNames,
		int startReadNum, 
		int endReadNum, 
		int pairedEnd,
		int numThreads,
		char *tmpDir,
		int *numWritten,
		int *numFiltered)
{
	char *FnName = "WriteReadsToTempFile";
	int i;
	int curSeqFPIndex=0;
	int curReadNum = 1;
	RGMatches m;
	char curLine[MAX_HEADER_LENGTH]="\0";
	fpos_t curFilePos;


	/* Read in any lines that begin with # */
	do {
		/* Get current file position */
		if(0!=fgetpos(seqFP, &curFilePos)) {
			PrintError(FnName,
					"fgetpos",
					"Could not get position in file",
					Exit,
					OutOfRange);
		}
	} while(NULL!=fgets(curLine, MAX_HEADER_LENGTH, seqFP) && curLine[0]=='#');
	/* Restore position */
	if(0!=fsetpos(seqFP, &curFilePos)) {
		PrintError(FnName,
				"fsetpos",
				"Could not set position in the file",
				Exit,
				OutOfRange);
	}

	(*numFiltered)=0;
	(*numWritten)=0;
	RGMatchesInitialize(&m);

	/* Open one temporary file, one for each thread */
	for(i=0;i<numThreads;i++) {
		assert((*tempSeqFPs)[i] == NULL);
		(*tempSeqFPs)[i] = OpenTmpFile(tmpDir, &(*tempSeqFileNames)[i]);
		assert((*tempSeqFPs)[i] != NULL);
	}

	while((endReadNum<=0 || endReadNum >= curReadNum) && 
			EOF != GetNextRead(seqFP, &m, pairedEnd)) {
		/* Get which tmp file to put the read in */
		curSeqFPIndex = (curReadNum-1)%numThreads;
			
		/* Print only if we are within the desired limit and the read checks out */
		if( (startReadNum<=0 || curReadNum >= startReadNum)  /* Only if we are within the bounds for the reads */
				&& (1 == UpdateRead((char*)m.matchOne.read, m.matchOne.readLength)) /* The first read is valid */
				&& (0 == pairedEnd || 1 == UpdateRead((char*)m.matchTwo.read, m.matchTwo.readLength))) { /* The second read is valid */
			/* Print */
			if(EOF == WriteRead((*tempSeqFPs)[curSeqFPIndex], &m, pairedEnd)) {
				PrintError(FnName,
						NULL,
						"Could not write read",
						Exit,
						WriteFileError);
			}
			(*numWritten)++;
		}
		else {
			assert(seqFilteredFP != NULL);
			/* Write to filtered read file */
			if(EOF == WriteRead(seqFilteredFP, &m, pairedEnd)) {
				PrintError(FnName,
						NULL,
						"Could not write read",
						Exit,
						WriteFileError);
			}
			(*numFiltered)++;
		}
		/* Increment read number */
		curReadNum++;
		/* Free memory */
		RGMatchesFree(&m);
	}

	/* reset pointer to temp files to the beginning of the file */
	for(i=0;i<numThreads;i++) {
		fseek((*tempSeqFPs)[i], 0, SEEK_SET);
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
	RGMatches m;
	int numReads = 0;
	int numOutputted=0;

	/* Initialize */
	RGMatchesInitialize(&m);

	/* Go to the beginning of the temporary output file */
	fseek(tempOutputFP, 0, SEEK_SET);

	while(RGMatchesRead(tempOutputFP, 
				&m,
				pairedEnd,
				binaryOutput)!=EOF) {
		/* Output if either pair has more than one entry */
		if(m.matchOne.numEntries > 0 ||
				m.matchTwo.numEntries > 0) {
			/* Output to final output file */
			RGMatchesPrint(outputFP,
					&m,
					pairedEnd,
					binaryOutput);
			numOutputted++;
		}
		else {
			/* Put back in the read file */
			if(EOF == WriteRead(tempSeqFP, &m, pairedEnd)) {
				PrintError(FnName,
						NULL,
						"Could not write read.",
						Exit,
						WriteFileError);
			}
			numReads++;
		}

		/* Free memory */
		RGMatchesFree(&m);
	}
	return numReads;
}

/* TODO */
/* Reads in a RGIndexfrom file */
void ReadRGIndex(char *rgIndexFileName, RGIndex *index, int binaryInput, int space)
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

	if(index->space != space) {
		PrintError("space",
				rgIndexFileName,
				"The index has a different space parity than specified",
				Exit,
				OutOfRange);
	}

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
