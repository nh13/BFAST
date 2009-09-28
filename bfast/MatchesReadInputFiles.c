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
/* Read the next read end from the stream */
int32_t GetNextReadBuffered(char **buffer,
		int32_t numLines,
		char *readName,
		char *read,
		char *qual)
{
	char *FnName = "GetNextReadBuffered";
	char comment[SEQUENCE_NAME_LENGTH]="\0";

	if('>' == buffer[0][0]) { // FASTA format
		// We could support FASTA by just putting in dummy values for the quals (?)
		PrintError(FnName, "read", "Read was in FASTA format.  Currently not supported", Exit, OutOfRange);
	}

	/* Move to the beginning of the read name */
	if(numLines < 4) {
		return -1;
	}

	/* Read in read name */
	strcpy(readName, buffer[0] + 1); // Ignore the '@'
	StringTrimWhiteSpace(readName);
	/* Read in sequence */
	strcpy(read, buffer[1]);
	StringTrimWhiteSpace(read);
	/* Read in comment line */
	strcpy(comment, buffer[2]);
	/* Read in qualities */
	strcpy(qual, buffer[3]);
	StringTrimWhiteSpace(qual);

	return 4;
}

/* TODO */
/* Read all the ends for the next read from the buffer */
int32_t GetReadBuffered(char **buffer,
		int32_t numLines,
		RGMatches *m,
		int space,
		int32_t abort)
{
	char *FnName = "GetReadBuffered";
	char readName[SEQUENCE_NAME_LENGTH]="\0";
	char read[SEQUENCE_LENGTH]="\0";
	char qual[SEQUENCE_LENGTH]="\0";
	int32_t foundAllEnds, curLine, numLinesRead;

	if(0 == numLines) {
		return -1;
	}
	foundAllEnds=0;

	curLine = 0;
	while(0 == foundAllEnds && curLine < numLines) {
		/* Get some reads */
		numLinesRead = GetNextReadBuffered(buffer + curLine,
				(numLines - curLine),
				readName,
				read,
				qual);
		if(numLinesRead < 0) {
			break;
		}

		// Add the number of lines read for now
		curLine += numLinesRead;

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
			// Subtract the number of ends reached since we don't want the last read
			curLine -= numLinesRead;
		}
	}
	
	if(0 == m->numEnds) { // No ends found
		return 0;
	}

	// Did not find all the ends since there wasn't enough buffer. Signal more buffer.
	if(0 == foundAllEnds && 0 == abort) { 
		return curLine;
	}

	// Check read
	if(0 == IsValidRead(m, space)) {
		PrintError(FnName, NULL, "The input read was not in the proper format", Exit, OutOfRange);
	}

	return curLine;
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
	char **lineBuffer=NULL;
	int32_t numLines = 0;
	int32_t curNumLines, numLinesProcessed, continueReading;

	/* Allocate the buffer */
	lineBuffer = malloc(sizeof(char*)*READS_BUFFER_LENGTH);
	if(NULL == lineBuffer) {
		PrintError(FnName, "lineBuffer", "Could not allocate memory", Exit, MallocMemory);
	}
	for(i=0;i<READS_BUFFER_LENGTH;i++) {
		lineBuffer[i] = malloc(sizeof(char)*(1+SEQUENCE_LENGTH));
		if(NULL == lineBuffer[i]) {
			PrintError(FnName, "lineBuffer[i]", "Could not allocate memory", Exit, MallocMemory);
		}
	}

	/* Open one temporary file, one for each thread */
	for(i=0;i<numThreads;i++) {
		assert((*tempSeqFPs)[i] == NULL);
		(*tempSeqFPs)[i] = OpenTmpFile(tmpDir, &(*tempSeqFileNames)[i]);
		assert((*tempSeqFPs)[i] != NULL);
	}

	/* Get the reads */
	continueReading = 1;
	(*numWritten)=0;
	RGMatchesInitialize(&m);
	do {
		/* Read in the buffer */
		for(i=numLines;i<READS_BUFFER_LENGTH;i++) {
			if(EOF == fscanf(seqFP, "%s", lineBuffer[i])) {
				continueReading = 0; break;
			}
			numLines++;
		}

		/* Process the buffer */
		curNumLines=0;
		while(curNumLines < numLines) {
			if(lineBuffer[curNumLines][0] != '#') { // skip comments
				RGMatchesInitialize(&m);

				/* Read in */
				numLinesProcessed=GetReadBuffered(lineBuffer + curNumLines, 
						(numLines - curNumLines),
						&m,
						space,
						continueReading);
				if(0 == numLinesProcessed) {
					PrintError(FnName, "numLinesProcessed == 0", "Could not read input file", Exit, OutOfRange);
				}
				
				/* Print only if we are within the desired limit and the read checks out */
				if((endReadNum<=0 || endReadNum >= curReadNum)) {
					/* Get which tmp file to put the read in */
					curSeqFPIndex = (curReadNum-1)%numThreads;

					if(EOF == WriteRead((*tempSeqFPs)[curSeqFPIndex], &m)) {
						PrintError(FnName, NULL, "Could not write read", Exit, WriteFileError);
					}
					(*numWritten)++;
				}

				curNumLines += numLinesProcessed;

				/* Increment read number */
				curReadNum++;
				/* Free memory */
				RGMatchesFree(&m);
				
			}
			else {
				curNumLines++;
			}
		}

		// DEBUGGING
		if(0 == continueReading) {
			if(curNumLines != numLines) {
				fprintf(stderr, "\n%d != %d\n", curNumLines, numLines);
				PrintError(FnName, "curNumLines != numLines", "Please report", Exit, OutOfRange);
			}
		}

		/* Copy over unused */
		for(i=0;i<(numLines - curNumLines);i++) {
			strcpy(lineBuffer[i], lineBuffer[curNumLines]);
		}
		numLines -= curNumLines;
	} while(continueReading == 1);

	/* Free the buffer */
	for(i=0;i<READS_BUFFER_LENGTH;i++) {
		free(lineBuffer[i]);
	}
	free(lineBuffer);

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
int GetIndexFileNames(char *fastaFileName, 
		int32_t space, 
		char *indexes, 
		char ***fileNames,
		int32_t ***indexIDs)
{
	char *FnName="GetIndexFileNames";
	char prefix[MAX_FILENAME_LENGTH]="\0";
	int32_t i, j, numFiles=0;
	int32_t *indexNumbers=NULL;
	int32_t numIndexNumbers=0;
	int32_t maxBin;

	assert(NULL != fastaFileName);

	/* Build the index file names */
	strcpy(prefix, fastaFileName);
	strcat(prefix, ".");
	strcat(prefix, SPACENAME(space));

	if(NULL != indexes) { // Tokenize
		indexNumbers = GetNumbersFromString(indexes, &numIndexNumbers);
		for(i=0;i<numIndexNumbers;i++) {
			if(indexNumbers[i] <= 0) {
				PrintError(FnName, indexes, "Could not understand index number", Exit, OutOfRange);
			}
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
				(*indexIDs) = realloc((*indexIDs), sizeof(int32_t*)*(1+numFiles));
				if(NULL == (*indexIDs)) {
					PrintError(FnName, "(*indexIDs)", "Could not reallocate memory", Exit, ReallocMemory);
				}
				(*indexIDs)[numFiles] = malloc(sizeof(int32_t)*2);
				if(NULL == (*indexIDs)[numFiles]) {
					PrintError(FnName, "(*indexIDs)[numFiles]", "Could not allocate memory", Exit, MallocMemory);
				}
				(*indexIDs)[numFiles][0] = indexNumbers[i];
				(*indexIDs)[numFiles][1] = j;
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

				(*indexIDs) = realloc((*indexIDs), sizeof(int32_t*)*(1+numFiles));
				if(NULL == (*indexIDs)) {
					PrintError(FnName, "(*indexIDs)", "Could not reallocate memory", Exit, ReallocMemory);
				}
				(*indexIDs)[numFiles] = malloc(sizeof(int32_t)*2);
				if(NULL == (*indexIDs)[numFiles]) {
					PrintError(FnName, "(*indexIDs)[numFiles]", "Could not allocate memory", Exit, MallocMemory);
				}
				(*indexIDs)[numFiles][0] = i;
				(*indexIDs)[numFiles][1] = j;

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
int32_t ReadOffsets(char *offsetsInput, int32_t **offsets) 
{
	char *FnName="ReadOffsets";
	int numOffsets=0;
	int32_t i;

	(*offsets) = GetNumbersFromString(offsetsInput, &numOffsets);

	if(NULL == (*offsets)) {
		PrintError(FnName, offsetsInput, "Could not parse the offsets", Exit, OutOfRange);
	}
	else if(numOffsets <= 0) {
		PrintError(FnName, NULL, "Could not find any offsets", Exit, OutOfRange);
	}

	// Check input
	for(i=0;i<numOffsets;i++) {
		if((*offsets)[i] < 0) {			
			PrintError(FnName, offsetsInput, "Offset was negative", Exit, OutOfRange);		
		}
		else if(0 < i && (*offsets)[i] <= (*offsets)[i-1]) {
			PrintError(FnName, offsetsInput, "Offset were not in increasing order", Exit, OutOfRange);		
		}
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
