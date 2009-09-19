#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <zlib.h>
#include "BError.h"
#include "BLib.h"
#include "RGMatch.h"
#include "RGMatches.h"
#include "MatchesReadInputFiles.h"

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
		PrintError(FnName, "read", "Could not read in read", Exit, ReadFileError);
	}
	StringTrimWhiteSpace(read);
	/* Read in comment line */
	if(NULL==fgets(comment, SEQUENCE_LENGTH-1, fp)) {
		/* Ignore */
		PrintError(FnName, "comment", "Could not read in comment", Exit, ReadFileError);
	}
	/* Read in qualities */
	if(NULL==fgets(qual, SEQUENCE_LENGTH-1, fp)) {
		PrintError(FnName, "qual", "Could not read in qualities", Exit, ReadFileError);
	}
	StringTrimWhiteSpace(qual);

	return 1;
}

/* TODO */
/* Read all the ends for the next read from the stream */
int GetRead(FILE *fp, 
		RGMatches *m,
		int space)
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
		/* Insert read name if this is the first end */
		if(0 == m->numEnds) {
			/* Allocate memory */
			m->readNameLength = strlen(readName);
			m->readName = malloc(sizeof(char)*(m->readNameLength+1));
			if(NULL == m->readName) {
				PrintError(FnName, "m->readName", "Could not allocate memory", Exit, MallocMemory);
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
				PrintError(FnName, "m->ends", "Could not reallocate memory", Exit, ReallocMemory);
			}
			RGMatchInitialize(&m->ends[m->numEnds-1]);
			m->ends[m->numEnds-1].readLength = strlen(read);
			m->ends[m->numEnds-1].qualLength = strlen(qual);
			m->ends[m->numEnds-1].read = malloc(sizeof(char)*(m->ends[m->numEnds-1].readLength+1));
			if(NULL == m->ends[m->numEnds-1].read) {
				PrintError(FnName, "m->ends[m->numEnds-1].read", "Could not allocate memory", Exit, MallocMemory);
			}
			strcpy(m->ends[m->numEnds-1].read, read);
			m->ends[m->numEnds-1].qual = malloc(sizeof(char)*(m->ends[m->numEnds-1].qualLength+1));
			if(NULL == m->ends[m->numEnds-1].qual) {
				PrintError(FnName, "m->ends[m->numEnds-1].qual", "Could not allocate memory", Exit, MallocMemory);
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
		if(0 != feof(fp)) {
			return EOF;
		}
		PrintError(FnName, "0 == m->numEnds", "Did not find any reads and qualities", Exit, ReadFileError);
	}

	if(0 == IsValidRead(m, space)) {
		PrintError(FnName, NULL, "The input read was not in the proper format", Exit, OutOfRange);
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
		FILE ***tempSeqFPs, /* One for each thread */ 
		char ***tempSeqFileNames,
		int startReadNum, 
		int endReadNum, 
		int numThreads,
		char *tmpDir,
		int *numWritten,
		int32_t space)
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
			PrintError(FnName, "fgetpos", "Could not get position in file", Exit, OutOfRange);
		}
	} while(NULL!=fgets(curLine, MAX_HEADER_LENGTH, seqFP) && curLine[0]=='#');
	/* Restore position */
	if(0!=fsetpos(seqFP, &curFilePos)) {
		PrintError(FnName, "fsetpos", "Could not set position in the file", Exit, OutOfRange);
	}

	(*numWritten)=0;
	RGMatchesInitialize(&m);

	/* Open one temporary file, one for each thread */
	for(i=0;i<numThreads;i++) {
		assert((*tempSeqFPs)[i] == NULL);
		(*tempSeqFPs)[i] = OpenTmpFile(tmpDir, &(*tempSeqFileNames)[i]);
		assert((*tempSeqFPs)[i] != NULL);
	}

	while((endReadNum<=0 || endReadNum >= curReadNum) && 
			EOF != GetRead(seqFP, &m, space)) {
		/* Get which tmp file to put the read in */
		curSeqFPIndex = (curReadNum-1)%numThreads;

		/* Print only if we are within the desired limit and the read checks out */
		if(startReadNum<=0 || curReadNum >= startReadNum) {/* Only if we are within the bounds for the reads */

			/* Print */
			if(EOF == WriteRead((*tempSeqFPs)[curSeqFPIndex], &m)) {
				PrintError(FnName, NULL, "Could not write read", Exit, WriteFileError);
			}
			(*numWritten)++;

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
int ReadTempReadsAndOutput(gzFile tempOutputFP,
		char *tempOutputFileName,
		gzFile outputFP,
		FILE *tempSeqFP)
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
	CloseTmpGZFile(&tempOutputFP,
			&tempOutputFileName,
			0);
	if(!(tempOutputFP=gzopen(tempOutputFileName, "rb"))) {
		PrintError(FnName, tempOutputFileName, "Could not re-open file for reading", Exit, OpenFileError);
	}

	while(RGMatchesRead(tempOutputFP, 
				&m)!=EOF) {
		/* Output if any end has more than one entry */
		for(i=hasEntries=0;0==hasEntries && i<m.numEnds;i++) {
			if(0 < m.ends[i].numEntries) {
				hasEntries=1;
			}
		}
		/* Output to final output file */
		if(1 == hasEntries) {
			RGMatchesPrint(outputFP,
					&m);
			numOutputted++;
		}
		else {
			/* Put back in the read file */
			if(EOF == WriteRead(tempSeqFP, &m)) {
				PrintError(FnName, NULL, "Could not write read.", Exit, WriteFileError);
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
		PrintError("space", rgIndexFileName, "The index has a different space parity than specified", Exit, OutOfRange);
	}
}

/* TODO */
int GetIndexFileNames(char *fastaFileName, int32_t space, char *indexes, char ***fileNames)
{
	char *FnName="GetIndexFileNames";
	char prefix[MAX_FILENAME_LENGTH]="\0";
	int32_t i, j, numFiles=0;
	int32_t *indexNumbers=NULL;
	int32_t numIndexNumbers=0;
	int32_t maxBin;
	char *pch=NULL;

	assert(NULL != fastaFileName);

	/* Build the index file names */
	strcpy(prefix, fastaFileName);
	strcat(prefix, ".");
	strcat(prefix, SPACENAME(space));

	if(NULL != indexes) { // Tokenize
		pch = strtok(indexes, ",");
		while(pch != NULL) {
			numIndexNumbers++;
			indexNumbers=realloc(indexNumbers, sizeof(int32_t)*numIndexNumbers);
			if(NULL == indexNumbers) {
				PrintError(FnName, "indexNumbers", "Could not reallocate memory", Exit, ReallocMemory);			
			}			
			indexNumbers[numIndexNumbers-1]=atoi(pch);			
			if(0 == indexNumbers[numIndexNumbers-1]) {				
				PrintError(FnName, pch, "Could not understand index number", Exit, OutOfRange);
			}
			pch = strtok(NULL, ",");
		}

		for(i=0;i<numIndexNumbers;i++) {
			maxBin = GetBIFMaximumBin(prefix, indexNumbers[i]);
			if(0 == maxBin) {
				fprintf(stderr, "Index number: %d\n", indexNumbers[i]);
				PrintError(FnName, NULL, "The index does not exist", Exit, OutOfRange);			
			}			
			(*fileNames) = realloc((*fileNames), sizeof(char*)*(numFiles+maxBin));			
			if(NULL == (*fileNames)) {				
				PrintError(FnName, "fileNames", "Could not reallocate memory", Exit, ReallocMemory);
			}
			/* Insert file names */
			for(j=1;j<=maxBin;j++) {
				(*fileNames)[numFiles] = malloc(sizeof(char)*MAX_FILENAME_LENGTH);
				if(NULL == (*fileNames)[numFiles]) {
					PrintError(FnName, "fileNames[j]", "Could not allocate memory", Exit, MallocMemory);				
				}				
				sprintf((*fileNames)[numFiles], "%s.%d.%d.%s", prefix, indexNumbers[i], 
						j,
						BFAST_INDEX_FILE_EXTENSION);
				if(0 == FileExists((*fileNames)[numFiles])) {
					PrintError(FnName, (*fileNames)[numFiles], "The index does not exist", Exit, OutOfRange);				
				}				
				numFiles++;			
			}		
		}
		free(indexNumbers);
	}
	else {
		i=1;
		numIndexNumbers=0;
		while(1) { // ^^
			maxBin = GetBIFMaximumBin(prefix, i);
			if(0 == maxBin) {
				break;
			}
			(*fileNames) = realloc((*fileNames), sizeof(char*)*(numFiles+maxBin));
			if(NULL == (*fileNames)) {
				PrintError(FnName, "fileNames", "Could not reallocate memory", Exit, ReallocMemory);			
			}			
			/* Insert file names */			
			for(j=1;j<=maxBin;j++) {				
				(*fileNames)[numFiles] = malloc(sizeof(char)*MAX_FILENAME_LENGTH);
				if(NULL == (*fileNames)[numFiles]) {
					PrintError(FnName, "fileNames[j]", "Could not allocate memory", Exit, MallocMemory);				
				}				
				sprintf((*fileNames)[numFiles], "%s.%d.%d.%s", prefix, i, 
						j,
						BFAST_INDEX_FILE_EXTENSION);
				if(0 == FileExists((*fileNames)[numFiles])) {
					PrintError(FnName, (*fileNames)[numFiles], "Missing Bin: The index does not exist", Exit, OutOfRange);				
				}				
				numFiles++;			
			}			
			i++;
			numIndexNumbers++;
		}
		free(indexNumbers);
	}

	if(0 == numIndexNumbers) {
		PrintError(FnName, prefix, "Could not find any indexes with the given prefix", Exit, OutOfRange);	
	}	
	//for(i=0;i<numFiles;i++) fprintf(stderr, "f[%d]=[%s]\n", i, (*fileNames)[i]);
	if(VERBOSE>=0) {
		if(1 == numIndexNumbers && 1 == numFiles) fprintf(stderr, "Found %d index (%d file).\n", numIndexNumbers, numFiles);
		else if(1 == numIndexNumbers && 1 != numFiles) fprintf(stderr, "Found %d index (%d total files).\n", numIndexNumbers, numFiles);
		else fprintf(stderr, "Found %d index (%d total files).\n", numIndexNumbers, numFiles);
	}

	return numFiles;
}

/* TODO */
int ReadOffsets(char *offsetsInput, int **offsets) 
{
	char *FnName="ReadOffsets";
	int numOffsets=0;
	char *pch=NULL;
	int start, end, offset;

	if(NULL != strchr(offsetsInput, ',')) { // use comma-separated offsets
		pch = strtok(offsetsInput, ",");
		while(pch != NULL) {
			numOffsets++;
			(*offsets)=(int*)realloc((*offsets), sizeof(int)*numOffsets);
			(*offsets)[numOffsets-1] = atoi(pch);
			if(0 == (*offsets)[numOffsets-1] && 0 == strcmp("0", pch)) {
				PrintError(FnName, pch, "Could not understand offset", Exit, OutOfRange);		
			}		
			else if((*offsets)[numOffsets-1] < 0) {			
				PrintError(FnName, pch, "Offset was negative", Exit, OutOfRange);		
			}
			pch = strtok(NULL, ",");
		}
	}
	else if(NULL != strchr(offsetsInput, '-')) { // use %d-%d
		pch = strtok(offsetsInput, "-");
		start = end = -1;
		while(pch != NULL) {
			offset = atoi(pch);
			if(0 == offset && strcmp("0", pch)) {
				PrintError(FnName, pch, "Could not understand offset", Exit, OutOfRange);		
			}		
			else if(offset < 0) {
				PrintError(FnName, pch, "Offset was negative", Exit, OutOfRange);		
			}
			if(start < 0) {
				start = offset;
			}
			else if(end < 0) {
				end = offset;
			}
			else {
				PrintError(FnName, NULL, "Offsets were not given in the correct format", Exit, OutOfRange);
			}
			pch = strtok(NULL, "-");
		}
		numOffsets=end-start+1;
		(*offsets)=(int*)realloc((*offsets), sizeof(int)*numOffsets);
		while(start<=end) {
			(*offsets)[end-start] = start;
			start++;
		}
	}
	if(numOffsets <= 0) {
		PrintError(FnName, NULL, "Could not find any offsets", Exit, OutOfRange);
	}

	if(VERBOSE>=0) {
		fprintf(stderr, "Read %d offsets.\n", numOffsets);
	}

	return numOffsets;
}

int32_t GetReads(FILE *fp, RGMatches *m, int32_t maxToRead, int32_t space) 
{
	int32_t numRead = 0;

	while(numRead < maxToRead) {
		RGMatchesInitialize(&(m[numRead]));
		if(EOF == GetRead(fp, &(m[numRead]), space)) {
			return numRead;
		}
		numRead++;
	}
	return numRead;
}
