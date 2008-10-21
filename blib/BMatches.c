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
#include "BMatch.h"
#include "BMatches.h"

#define BMATCHES_CHECK 0

/* TODO */
int32_t BMatchesRead(FILE *fp,
		BMatches *m,
		int32_t pairedEnd,
		int32_t binaryInput)
{
	char *FnName = "BMatchesRead";

	/* Read the matches from the input file */
	if(binaryInput == TextInput) {
		/* Read paired end */
		if(fscanf(fp, "%d", &m->pairedEnd)==EOF) {
			return EOF;
		}
		/* Check paired end */
		if(m->pairedEnd != pairedEnd) {
			PrintError(FnName,
					"pairedEnd",
					"Error.  Paired end did not match",
					Exit,
					OutOfRange);
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
		if(m->pairedEnd != pairedEnd) {
			PrintError(FnName,
					"pairedEnd",
					"Error.  Paired end did not match",
					Exit,
					OutOfRange);
		}
	}
	/* Read read name */
	if(BStringRead(&m->readName, fp, binaryInput)==EOF) {
		PrintError(FnName,
				"m->readName",
				"Could not read name",
				Exit,
				ReadFileError);
	}

	/* Read match one */
	BMatchRead(fp,
			&m->matchOne,
			binaryInput);
	if(1==pairedEnd) {
		/* Read match two if necessary */
		BMatchRead(fp,
				&m->matchTwo,
				binaryInput);
	}

	/* Check m */
	if(1==BMATCHES_CHECK) {
		BMatchesCheck(m);
	}

	return 1;
}

/* TODO */
void BMatchesPrint(FILE *fp,
		BMatches *m,
		int32_t pairedEnd,
		int32_t binaryOutput)
{
	char *FnName = "BMatchesPrint";
	assert(fp!=NULL);

	/* Check m */
	if(1==BMATCHES_CHECK) {
		BMatchesCheck(m);
	}

	/* Print the matches to the output file */
	if(binaryOutput == TextInput) {
		/* Print paired end */
		if(0 > fprintf(fp, "%d %d %s\n",
					m->pairedEnd)) {
			PrintError(FnName,
					NULL,
					"Could not write m->pairedEnd",
					Exit,
					WriteFileError);
		}
	}
	else {
		/* Print paired end, read name length, and read name */
		if(fwrite(&m->pairedEnd, sizeof(int32_t), 1, fp) != 1) {
			PrintError(FnName,
					NULL,
					"Could not write m->pairedEnd",
					Exit,
					WriteFileError);
		}
	}

	/* Print read name */
	if(BStringPrint(&m->readName, fp, binaryOutput)<0) {
		PrintError(FnName,
				"m->readName",
				"Could not write read name",
				Exit,
				WriteFileError);
	}

	/* Print match one */
	BMatchPrint(fp,
			&m->matchOne,
			binaryOutput);
	/* Print match two if necessary */
	if(pairedEnd == 1) {
		BMatchPrint(fp,
				&m->matchTwo,
				binaryOutput);
	}
}

/* TODO */
void BMatchesRemoveDuplicates(BMatches *m,
		int32_t maxNumMatches)
{
	BMatchRemoveDuplicates(&m->matchOne, maxNumMatches);
	if(1==m->pairedEnd) {
		BMatchRemoveDuplicates(&m->matchTwo, maxNumMatches);
		assert(m->matchTwo.numEntries <= maxNumMatches);
	}
	assert(m->matchOne.numEntries <= maxNumMatches);

	/* Check m */
	if(1==BMATCHES_CHECK) {
		BMatchesCheck(m);
	}
}

/* TODO */
/* Merges matches from the same read */
int32_t BMatchesMergeFilesAndOutput(FILE **tempFPs,
		int32_t numFiles,
		FILE *outputFP,
		int32_t pairedEnd,
		int32_t binaryOutput,
		int32_t maxNumMatches)
{
	char *FnName="BMatchesMergeFilesAndOutput";
	int32_t i;
	int32_t counter;
	BMatches matches;
	BMatches tempMatches;
	int32_t numMatches=0;
	int32_t numFinished = 0;

	/* Initialize matches */
	BMatchesInitialize(&matches);
	BMatchesInitialize(&tempMatches);

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
		if(VERBOSE >=0 && counter%BMATCH_MEBE_ROTATE_NUM == 0) {
			fprintf(stderr, "\r[%d]", counter);
		}
		counter++;

		/* Read matches for one read from each file */ 
		for(i=0;i<numFiles;i++) {
			if(BMatchesRead(tempFPs[i],
						&tempMatches,
						pairedEnd,
						binaryOutput)==EOF) {
				numFinished++;
			}
			else {
				if(matches.readName.string != NULL &&
						BStringCompare(&matches.readName, &tempMatches.readName)!=0) {
					PrintError(FnName,
							NULL,
							"Read names do not match",
							Exit,
							OutOfRange);
				}
				/* Append temp matches to matches */
				BMatchesAppend(&tempMatches, &matches);
			}

			/* Free temp matches */
			BMatchesFree(&tempMatches);
		}
		/* We must finish all at the same time */
		assert(numFinished == 0 || numFinished == numFiles);

		if(numFinished == 0) {
			/* Remove duplicates */
			BMatchesRemoveDuplicates(&matches, maxNumMatches);

			/* Print to output file */
			if(matches.matchOne.numEntries > 0 ||
					(1==pairedEnd && matches.matchTwo.numEntries > 0)) {
				numMatches++;
			}
			BMatchesPrint(outputFP,
					&matches,
					pairedEnd,
					binaryOutput);
		}
		/* Free memory */
		BMatchesFree(&matches);
	}

	if(VERBOSE >=0) {
		fprintf(stderr, "\r[%d]... completed.\n", counter-1);
	}

	return numMatches;
}

/* TODO */
/* Merges matches from different reads */
int32_t BMatchesMergeThreadTempFilesIntoOutputTempFile(FILE **threadFPs,
		int32_t numThreads,
		FILE *outputFP,
		int32_t pairedEnd,
		int32_t binaryOutput)
{
	char *FnName = "BMatchesMergeThreadTempFilesIntoOutputTempFile";
	int32_t counter;
	int32_t i;
	BMatches matches;
	int32_t numFinished;
	int *finished=NULL;

	/* Initialize matches */
	BMatchesInitialize(&matches);

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
				if(BMatchesRead(threadFPs[i],
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
					if(VERBOSE >=0 && counter%BMATCH_MEBE_ROTATE_NUM == 0) {
						fprintf(stderr, "\r[%d]", counter);
					}
					counter++;

					BMatchesPrint(outputFP,
							&matches,
							pairedEnd,
							binaryOutput);

				}
				/* Free memory */
				BMatchesFree(&matches);
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
void BMatchesAppend(BMatches *src, BMatches *dest)
{
	char *FnName = "BMatchesAppend";
	/* Check that we are not appending to ourselves */
	assert(src != dest);

	/* Check to see if we need to add in the read name */
	if(dest->readName.length <= 0) {
		BStringCopy(&dest->readName, &src->readName);
		dest->pairedEnd = src->pairedEnd;
	}
	assert(src->pairedEnd == dest->pairedEnd);

	/* Append the matches */
	BMatchAppend(&src->matchOne, &dest->matchOne);
	if(1==dest->pairedEnd) {
		BMatchAppend(&src->matchTwo, &dest->matchTwo);
	}
}

/* TODO */
void BMatchesFree(BMatches *m) 
{
	BStringFree(&m->readName);
	BMatchFree(&m->matchOne);
	BMatchFree(&m->matchTwo);
	BMatchesInitialize(m);
}

/* TODO */
void BMatchesInitialize(BMatches *m)
{
	m->pairedEnd = 0;
	BStringInitialize(&m->readName);
	BMatchInitialize(&m->matchOne);
	BMatchInitialize(&m->matchTwo);
}

/* TODO */
void BMatchesMirrorPairedEnd(BMatches *m,
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
			BMatchReallocate(&m->matchTwo, numEntriesOne + numEntriesTwo);
			for(i=0;i<numEntriesOne;i++) {
				m->matchTwo.contigs[i+numEntriesTwo] = m->matchOne.contigs[i];
				m->matchTwo.strand[i+numEntriesTwo] = m->matchOne.strand[i];
				/* Adjust position */
				m->matchTwo.positions[i+numEntriesTwo] = m->matchOne.positions[i] + m->matchOne.readLength + pairedEndLength;
			}
		}
		/* Copy matches from second to first */
		if(forceMirroring == 1 || 
				(m->matchOne.numEntries <= 0 && m->matchTwo.numEntries > 0)) {
			BMatchReallocate(&m->matchOne, numEntriesOne + numEntriesTwo);
			for(i=0;i<numEntriesTwo;i++) {
				m->matchOne.contigs[i+numEntriesOne] = m->matchTwo.contigs[i];
				m->matchOne.strand[i+numEntriesOne] = m->matchTwo.strand[i];
				/* Adjust position */
				m->matchOne.positions[i+numEntriesOne] = m->matchTwo.positions[i] - m->matchOne.readLength - pairedEndLength;
			}
		}
		/* Check m */
		if(1==BMATCHES_CHECK) {
			BMatchesCheck(m);
		}
		BMatchesRemoveDuplicates(m, INT_MAX);
	}
}

/* TODO */
void BMatchesCheck(BMatches *m) 
{
	char *FnName="BMatchesCheck";
	/* Basic asserts */
	assert(m->pairedEnd == 0 || m->pairedEnd == 1);
	/* Check that the read name length is the same as the length of the read name */
	if(m->readName.length != (int)strlen(m->readName.string)) {
		PrintError(FnName,
				NULL,
				"strlen(m->readName)) != m->readNameLength",
				Exit,
				OutOfRange);
	}
	if(0 == m->pairedEnd) {
		/* Check the first match */
		BMatchCheck(&m->matchOne);
		/* Make sure the second match is nulled */
		assert(m->matchTwo.readLength == 0);
		assert(m->matchTwo.read == NULL);
		assert(m->matchTwo.maxReached == 0);
		assert(m->matchTwo.numEntries == 0);
		assert(m->matchTwo.contigs == NULL);
		assert(m->matchTwo.positions == NULL);
		assert(m->matchTwo.strand == NULL);
	}
	else {
		/* Check both matches */
		BMatchCheck(&m->matchOne);
		BMatchCheck(&m->matchTwo);
	}
}

/* TODO */
void BMatchesFilterOutOfRange(BMatches *m,
		int32_t startChr,
		int32_t startPos,
		int32_t endChr,
		int32_t endPos,
		int32_t maxNumMatches)
{
	BMatchFilterOutOfRange(&m->matchOne,
			startChr,
			startPos,
			endChr,
			endPos,
			maxNumMatches);
	if(m->pairedEnd == 1) {
		BMatchFilterOutOfRange(&m->matchTwo,
				startChr,
				startPos,
				endChr,
				endPos,
				maxNumMatches);
	}
}
