#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include "BLibDefinitions.h"
#include "BError.h"
#include "RGMatch.h"
#include "RGMatches.h"

#define RGMATCHES_CHECK 0

/* TODO */
int32_t RGMatchesRead(FILE *fp,
		RGMatches *m,
		int32_t pairedEnd,
		int32_t binaryInput)
{
	char *FnName = "RGMatchesRead";

	/* Read the matches from the input file */
	if(binaryInput == TextInput) {
		/* Read paired end */
		if(fscanf(fp, "%d", &m->pairedEnd)==EOF) {
			return EOF;
		}
		/* Check paired end */
		if(pairedEnd != PairedEndDoesNotMatter && 
				m->pairedEnd != pairedEnd) {
			PrintError(FnName,
					"pairedEnd",
					"Error.  Paired end did not match",
					Exit,
					OutOfRange);
		}
		/* Read read name length */
		if(fscanf(fp, "%d", &m->readNameLength)==EOF) {
			PrintError(FnName,
					"m->readNameLength",
					"Could not read name length",
					Exit,
					ReadFileError);
		}
		assert(m->readNameLength < SEQUENCE_NAME_LENGTH);
		assert(m->readNameLength > 0);
		/* Allocate memory for the read name */
		m->readName = malloc(sizeof(int8_t)*(m->readNameLength + 1));
		if(NULL == m->readName) {
			PrintError(FnName,
					"m->readName",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}

		/* Read read name */
		if(EOF==fscanf(fp, "%s", m->readName)) {
			PrintError(FnName,
					"m->readName",
					"Could not read in read name",
					Exit,
					ReadFileError);
		}
	}
	else {
		/* Read paired end */
		if(fread(&m->pairedEnd, sizeof(int32_t), 1, fp)!=1) {
			if(feof(fp) != 0) {
				return EOF;
			}
			else {
				PrintError(FnName,
						"paired end",
						"Could not read in paired end",
						Exit,
						ReadFileError);
			}
		}
		/* Check paired end */
		if(pairedEnd != PairedEndDoesNotMatter && 
		m->pairedEnd != pairedEnd) {
			PrintError(FnName,
					"pairedEnd",
					"Error.  Paired end did not match",
					Exit,
					OutOfRange);
		}

		/* Read read name length */
		if(fread(&m->readNameLength, sizeof(int32_t), 1, fp)!=1) {
			PrintError(FnName,
					"m->readNameLength",
					"Could not read in read name length",
					Exit,
					ReadFileError);
		}
		assert(m->readNameLength < SEQUENCE_NAME_LENGTH);
		assert(m->readNameLength > 0);

		/* Allocate memory for the read name */
		m->readName = malloc(sizeof(int8_t)*(m->readNameLength + 1));
		if(NULL == m->readName) {
			PrintError(FnName,
					"m->readName",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}

		/* Read in read name */
		if(fread(m->readName, sizeof(int8_t), m->readNameLength, fp)!=m->readNameLength) {
			PrintError(FnName,
					"m->readName",
					"Could not read in read name",
					Exit,
					ReadFileError);
		}
		m->readName[m->readNameLength]='\0';
	}

	/* Read match one */
	RGMatchRead(fp,
			&m->matchOne,
			binaryInput);
	if(1==m->pairedEnd) {
		/* Read match two if necessary */
		RGMatchRead(fp,
				&m->matchTwo,
				binaryInput);
	}

	/* Check m */
	if(1==RGMATCHES_CHECK) {
		RGMatchesCheck(m);
	}

	return 1;
}

/* TODO */
void RGMatchesPrint(FILE *fp,
		RGMatches *m,
		int32_t binaryOutput)
{
	char *FnName = "RGMatchesPrint";
	assert(fp!=NULL);

	/* Check m */
	/*
	   if(1==RGMATCHES_CHECK) {
	   RGMatchesCheck(m);
	   }
	   */

	/* Print the matches to the output file */
	if(binaryOutput == TextInput) {
		/* Print paired end, read name length, and read name */
		if(0 > fprintf(fp, "%d %d %s\n",
					m->pairedEnd,
					m->readNameLength,
					m->readName)) {
			PrintError(FnName,
					NULL,
					"Could not write m->pairedEnd, m->readNameLength, and m->readName",
					Exit,
					WriteFileError);
		}
	}
	else {
		/* Print paired end, read name length, and read name */
		if(fwrite(&m->pairedEnd, sizeof(int32_t), 1, fp) != 1 ||
				fwrite(&m->readNameLength, sizeof(int32_t), 1, fp) != 1 ||
				fwrite(m->readName, sizeof(int8_t), m->readNameLength, fp) != m->readNameLength) {
			PrintError(FnName,
					NULL,
					"Could not write m->pairedEnd, m->readNameLength, and m->readName",
					Exit,
					WriteFileError);
		}
	}

	/* Print match one */
	RGMatchPrint(fp,
			&m->matchOne,
			binaryOutput);
	/* Print match two if necessary */
	if(m->pairedEnd == 1) {
		RGMatchPrint(fp,
				&m->matchTwo,
				binaryOutput);
	}
}

/* TODO */
void RGMatchesRemoveDuplicates(RGMatches *m,
		int32_t maxNumMatches)
{
	RGMatchRemoveDuplicates(&m->matchOne, maxNumMatches);
	assert(m->matchOne.numEntries <= maxNumMatches);
	if(1 == m->pairedEnd) {
		RGMatchRemoveDuplicates(&m->matchTwo, maxNumMatches);
		assert(m->matchTwo.numEntries <= maxNumMatches);
	}

	/* Check m */
	/*
	   if(1==RGMATCHES_CHECK) {
	   RGMatchesCheck(m);
	   }
	   */
}

/* TODO */
/* Merges matches from the same read */
int32_t RGMatchesMergeFilesAndOutput(FILE **tempFPs,
		int32_t numFiles,
		FILE *outputFP,
		int32_t pairedEnd,
		int32_t binaryOutput,
		int32_t maxNumMatches)
{
	char *FnName="RGMatchesMergeFilesAndOutput";
	int32_t i;
	int32_t counter;
	RGMatches matches;
	RGMatches tempMatches;
	int32_t numMatches=0;
	int32_t numFinished = 0;

	/* Initialize matches */
	RGMatchesInitialize(&matches);
	RGMatchesInitialize(&tempMatches);

	/* Seek to the beginning of the files */
	for(i=0;i<numFiles;i++) {
		fseek(tempFPs[i], 0, SEEK_SET);
	}

	/* Read in each sequence/match one at a time */
	counter = 0;
	if(VERBOSE >=0) {
		fprintf(stderr, "\r[%d]", 0);
	}
	while(numFinished == 0) {
		if(VERBOSE >=0 && counter%RGMATCH_MERGE_ROTATE_NUM == 0) {
			fprintf(stderr, "\r[%d]", counter);
		}
		counter++;

		/* Read matches for one read from each file */ 
		for(i=0;i<numFiles;i++) {
			if(RGMatchesRead(tempFPs[i],
						&tempMatches,
						pairedEnd,
						binaryOutput)==EOF) {
				numFinished++;
			}
			else {
				if(matches.readName != NULL &&
						strcmp((char*)matches.readName, (char*)tempMatches.readName)!=0) {
					PrintError(FnName,
							NULL,
							"Read names do not match",
							Exit,
							OutOfRange);
				}
				/* Append temp matches to matches */
				RGMatchesAppend(&tempMatches, &matches);
			}

			/* Free temp matches */
			RGMatchesFree(&tempMatches);
		}
		/* We must finish all at the same time */
		assert(numFinished == 0 || numFinished == numFiles);

		if(numFinished == 0) {
			/* Remove duplicates */
			RGMatchesRemoveDuplicates(&matches, maxNumMatches);

			/* Print to output file */
			if(matches.matchOne.numEntries > 0 ||
					(1==pairedEnd && matches.matchTwo.numEntries > 0)) {
				numMatches++;
			}
			RGMatchesPrint(outputFP,
					&matches,
					binaryOutput);
		}
		/* Free memory */
		RGMatchesFree(&matches);
	}

	if(VERBOSE >=0) {
		fprintf(stderr, "\r[%d]... completed.\n", counter-1);
	}

	return numMatches;
}

/* TODO */
/* Merges matches from different reads */
int32_t RGMatchesMergeThreadTempFilesIntoOutputTempFile(FILE **threadFPs,
		int32_t numThreads,
		FILE *outputFP,
		int32_t pairedEnd,
		int32_t binaryOutput)
{
	char *FnName = "RGMatchesMergeThreadTempFilesIntoOutputTempFile";
	int32_t counter;
	int32_t i;
	RGMatches matches;
	int32_t numFinished;
	int *finished=NULL;

	/* Initialize matches */
	RGMatchesInitialize(&matches);

	/* Initialize thread file pointers */
	for(i=0;i<numThreads;i++) {
		fseek(threadFPs[i], 0, SEEK_SET);
	}

	/* Allocate memory for the finished array */
	finished = malloc(sizeof(int)*numThreads);
	if(NULL == finished) {
		PrintError(FnName,
				"finished",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Initialize finished array */
	for(i=0;i<numThreads;i++) {
		finished[i] = 0;
	}

	counter = 0;
	numFinished = 0;
	if(VERBOSE >= 0) {
		fprintf(stderr, "\r[%d]", counter);
	}
	while(numFinished < numThreads) {
		/* For each thread */
		for(i=0;i<numThreads;i++) {
			/* Only try reading from those that are not finished */
			if(0 == finished[i]) {
				if(RGMatchesRead(threadFPs[i],
							&matches,
							pairedEnd,
							binaryOutput)==EOF) {
					finished[i] = 1;
					numFinished++;
				}
				else {
					assert(numFinished < numThreads);

					/*
					   if(match.numEntries > 0) {
					   counter++;
					   }
					   */
					if(VERBOSE >=0 && counter%RGMATCH_MERGE_ROTATE_NUM == 0) {
						fprintf(stderr, "\r[%d]", counter);
					}
					counter++;

					RGMatchesPrint(outputFP,
							&matches,
							binaryOutput);

				}
				/* Free memory */
				RGMatchesFree(&matches);
			}
		}
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\r[%d]\n", counter);
	}
	assert(numFinished == numThreads);
	for(i=0;i<numThreads;i++) {
		assert(1==finished[i]);
	}

	/* Free memory */
	free(finished);

	return counter;
}

/* TODO */
void RGMatchesAppend(RGMatches *src, RGMatches *dest)
{
	char *FnName = "RGMatchesAppend";
	/* Check that we are not appending to ourselves */
	assert(src != dest);

	/* Check to see if we need to add in the read name */
	if(dest->readNameLength <= 0) {
		assert(dest->readName == NULL);
		dest->readNameLength = src->readNameLength;

		/* Allocate memory for the read name */
		dest->readName = malloc(sizeof(int8_t)*(dest->readNameLength+1));
		if(NULL==dest->readName) {
			PrintError(FnName,
					"dest->readName",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		strcpy((char*)dest->readName, (char*)src->readName);
		dest->pairedEnd = src->pairedEnd;
	}
	assert(src->pairedEnd == dest->pairedEnd);

	/* Append the matches */
	RGMatchAppend(&src->matchOne, &dest->matchOne);
	if(1==dest->pairedEnd) {
		RGMatchAppend(&src->matchTwo, &dest->matchTwo);
	}
}

/* TODO */
void RGMatchesFree(RGMatches *m) 
{
	free(m->readName);
	RGMatchFree(&m->matchOne);
	RGMatchFree(&m->matchTwo);
	RGMatchesInitialize(m);
}

/* TODO */
void RGMatchesInitialize(RGMatches *m)
{
	m->pairedEnd = 0;
	m->readNameLength = 0;
	m->readName = NULL;
	RGMatchInitialize(&m->matchOne);
	RGMatchInitialize(&m->matchTwo);
}

/* TODO */
void RGMatchesMirrorPairedEnd(RGMatches *m,
		int32_t pairedEndLength,
		int32_t forceMirroring)
{
	int i;
	int numEntriesOne, numEntriesTwo;

	if(m->pairedEnd == 1) {
		numEntriesOne = m->matchOne.numEntries;
		numEntriesTwo = m->matchTwo.numEntries;
		/* Copy matches from first to second */
		if(forceMirroring == 1 || 
				(m->matchOne.numEntries > 0 && m->matchTwo.numEntries <= 0)) {
			RGMatchReallocate(&m->matchTwo, numEntriesOne + numEntriesTwo);
			for(i=0;i<numEntriesOne;i++) {
				m->matchTwo.contigs[i+numEntriesTwo] = m->matchOne.contigs[i];
				m->matchTwo.strands[i+numEntriesTwo] = m->matchOne.strands[i];
				/* Adjust position */
				m->matchTwo.positions[i+numEntriesTwo] = m->matchOne.positions[i] + m->matchOne.readLength + pairedEndLength;
			}
		}
		/* Copy matches from second to first */
		if(forceMirroring == 1 || 
				(m->matchOne.numEntries <= 0 && m->matchTwo.numEntries > 0)) {
			RGMatchReallocate(&m->matchOne, numEntriesOne + numEntriesTwo);
			for(i=0;i<numEntriesTwo;i++) {
				m->matchOne.contigs[i+numEntriesOne] = m->matchTwo.contigs[i];
				m->matchOne.strands[i+numEntriesOne] = m->matchTwo.strands[i];
				/* Adjust position */
				m->matchOne.positions[i+numEntriesOne] = m->matchTwo.positions[i] - m->matchOne.readLength - pairedEndLength;
			}
		}
		/* Check m */
		/*
		   if(1==RGMATCHES_CHECK) {
		   RGMatchesCheck(m);
		   }
		   */
		RGMatchesRemoveDuplicates(m, INT_MAX);
	}
}

/* TODO */
void RGMatchesCheck(RGMatches *m) 
{
	char *FnName="RGMatchesCheck";
	/* Basic asserts */
	assert(m->pairedEnd == 0 || m->pairedEnd == 1);
	/* Check that the read name length is the same as the length of the read name */
	if(((int)strlen((char*)m->readName)) != m->readNameLength) {
		PrintError(FnName,
				NULL,
				"strlen(m->readName)) != m->readNameLength",
				Exit,
				OutOfRange);
	}
	if(0 == m->pairedEnd) {
		/* Check the first match */
		RGMatchCheck(&m->matchOne);
		/* Make sure the second match is nulled */
		assert(m->matchTwo.readLength == 0);
		assert(m->matchTwo.read == NULL);
		assert(m->matchTwo.maxReached == 0);
		assert(m->matchTwo.numEntries == 0);
		assert(m->matchTwo.contigs == NULL);
		assert(m->matchTwo.positions == NULL);
		assert(m->matchTwo.strands == NULL);
	}
	else {
		/* Check both matches */
		RGMatchCheck(&m->matchOne);
		RGMatchCheck(&m->matchTwo);
	}
}

/* TODO */
void RGMatchesFilterOutOfRange(RGMatches *m,
		int32_t startChr,
		int32_t startPos,
		int32_t endChr,
		int32_t endPos,
		int32_t maxNumMatches)
{
	RGMatchFilterOutOfRange(&m->matchOne,
			startChr,
			startPos,
			endChr,
			endPos,
			maxNumMatches);
	if(m->pairedEnd == 1) {
		RGMatchFilterOutOfRange(&m->matchTwo,
				startChr,
				startPos,
				endChr,
				endPos,
				maxNumMatches);
	}
}
