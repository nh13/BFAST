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
/* Read the next read end from the stream */
int GetNextRead(FILE *fp, 
		char *readName,
		char *read,
		char *qual)
{
	char *FnName = "GetNextRead";
	char comment[SEQUENCE_NAME_LENGTH]="\0";

	/* Move to the beginning of the read name */
	assert(0 == feof(fp));
	while(0 == feof(fp)
			&& '@' != fgetc(fp)) {
	}
	/* Read in read name */
	if(0 != feof(fp) ||
			NULL==fgets(readName, SEQUENCE_NAME_LENGTH-1, fp)) {
		return EOF;
	}
	StringTrimWhiteSpace(readName);
	/* Read in sequence */
	if(NULL==fgets(read, SEQUENCE_LENGTH-1, fp)) {
		PrintError(FnName, 
				"read",
				"Could not read in read",
				Exit,
				ReadFileError);
	}
	StringTrimWhiteSpace(read);
	/* Read in comment line */
	if(NULL==fgets(comment, SEQUENCE_LENGTH-1, fp)) {
		/* Ignore */
		PrintError(FnName, 
				"comment",
				"Could not read in comment",
				Exit,
				ReadFileError);
	}
	/* Read in qualities */
	if(NULL==fgets(qual, SEQUENCE_LENGTH-1, fp)) {
		PrintError(FnName, 
				"qual",
				"Could not read in qualities",
				Exit,
				ReadFileError);
	}
	StringTrimWhiteSpace(qual);
	return 1;
}

/* TODO */
/* Read all the ends for the next read from the stream */
int GetRead(FILE *fp, 
		RGMatches *m)
{
	char *FnName = "GetRead";
	fpos_t curPos;
	char readName[SEQUENCE_NAME_LENGTH]="\0";
	char read[SEQUENCE_LENGTH]="\0";
	char qual[SEQUENCE_LENGTH]="\0";
	int32_t foundAllEnds;

	if(0 != feof(fp)) {
		return EOF;
	}
	assert(0 == fgetpos(fp, &curPos));
	foundAllEnds=0;
	while(0 == foundAllEnds &&
			EOF != GetNextRead(fp,
				readName,
				read,
				qual)) {
		/* Inset read name if this is the first end */
		if(0 == m->numEnds) {
			/* Allocate memory */
			m->readNameLength = strlen(readName);
			m->readName = malloc(sizeof(char)*(m->readNameLength+1));
			if(NULL == m->readName) {
				PrintError(FnName,
						"m->readName",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
			strcpy(m->readName, readName);
		}
		/* Add end if necessary */
		if(0 == strcmp(m->readName, readName)) {
			/* Save position */
			assert(0 == fgetpos(fp, &curPos));
			/* Reallocate */
			m->numEnds++;
			m->ends = realloc(m->ends, sizeof(RGMatch)*m->numEnds);
			if(NULL == m->ends) {
				PrintError(FnName,
						"m->ends",
						"Could not reallocate memory",
						Exit,
						ReallocMemory);
			}
			RGMatchInitialize(&m->ends[m->numEnds-1]);
			m->ends[m->numEnds-1].readLength = strlen(read);
			m->ends[m->numEnds-1].qualLength = strlen(qual);
			m->ends[m->numEnds-1].read = malloc(sizeof(char)*(m->ends[m->numEnds-1].readLength+1));
			if(NULL == m->ends[m->numEnds-1].read) {
				PrintError(FnName,
						"m->ends[m->numEnds-1].read",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
			strcpy(m->ends[m->numEnds-1].read, read);
			m->ends[m->numEnds-1].qual = malloc(sizeof(char)*(m->ends[m->numEnds-1].qualLength+1));
			if(NULL == m->ends[m->numEnds-1].qual) {
				PrintError(FnName,
						"m->ends[m->numEnds-1].qual",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
			strcpy(m->ends[m->numEnds-1].qual, qual);
		}
		else {
			/* Found all ends */
			foundAllEnds = 1;
			/* Restore position */
			assert(0 == fsetpos(fp, &curPos));
		}
	}

	if(0 == m->numEnds) {
		PrintError(FnName,
				"0 == m->numEnds",
				"Did not find any reads and qualities",
				Exit,
				ReadFileError);
	}

	return 1;
}

/* TODO */
int WriteRead(FILE *fp, RGMatches *m)
{
	int32_t i;

	/* Print read */
	for(i=0;i<m->numEnds;i++) {
		if(fprintf(fp, "@%s\n%s\n+\n%s\n", 
					m->readName,
					m->ends[i].read,
					m->ends[i].qual) < 0) { 
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
		int numThreads,
		char *tmpDir,
		int *numWritten,
		int *numFiltered,
		int32_t space)
{
	char *FnName = "WriteReadsToTempFile";
	int i;
	int curSeqFPIndex=0;
	int curReadNum = 1;
	int32_t isValidRead = 0;
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
			EOF != GetRead(seqFP, &m)) {
		/* Get which tmp file to put the read in */
		curSeqFPIndex = (curReadNum-1)%numThreads;

		/* Print only if we are within the desired limit and the read checks out */
		if(startReadNum<=0 || curReadNum >= startReadNum) {/* Only if we are within the bounds for the reads */

			for(i=0,isValidRead=1;1==isValidRead && i<m.numEnds;i++) {
				if(0 == UpdateRead(m.ends[i].read,
							m.ends[i].readLength)) {
					isValidRead= 0;
				}
			}

			if(1 == isValidRead) {
				/* Print */
				if(EOF == WriteRead((*tempSeqFPs)[curSeqFPIndex], &m)) {
					PrintError(FnName,
							NULL,
							"Could not write read",
							Exit,
							WriteFileError);
				}
				(*numWritten)++;
			}
		}
		else {
			assert(seqFilteredFP != NULL);
			/* Write to filtered read file */
			if(EOF == WriteRead(seqFilteredFP, &m)) {
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
		int binaryOutput)
{
	char *FnName = "ReadTempReadsAndOutput";
	RGMatches m;
	int32_t i;
	int numReads = 0;
	int numOutputted=0;
	int hasEntries=0;

	/* Initialize */
	RGMatchesInitialize(&m);

	/* Go to the beginning of the temporary output file */
	fseek(tempOutputFP, 0, SEEK_SET);

	while(RGMatchesRead(tempOutputFP, 
				&m,
				binaryOutput)!=EOF) {
		/* Output if any end has more than one entry */
		for(i=hasEntries=0;0==hasEntries && i<m.numEnds;i++) {
			if(0 < m.ends[i].numEntries) {
				hasEntries=1;
			}
		}
		if(1 == hasEntries) {
			/* Output to final output file */
			RGMatchesPrint(outputFP,
					&m,
					binaryOutput);
			numOutputted++;
		}
		else {
			/* Put back in the read file */
			if(EOF == WriteRead(tempSeqFP, &m)) {
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

void ReadRGIndex(char *rgIndexFileName, RGIndex *index, int space)
{

	/* Read from file */
	RGIndexRead(index, rgIndexFileName);

	if(index->space != space) {
		PrintError("space",
				rgIndexFileName,
				"The index has a different space parity than specified",
				Exit,
				OutOfRange);
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
