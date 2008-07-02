#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <string.h>
#include "BLibDefinitions.h"
#include "BError.h"
#include "RGMatch.h"

/* TODO */
void RGMatchRemoveDuplicates(RGMatch *s,
		int32_t maxMatches)
{
	int32_t i;
	int32_t prevIndex=0;

	/* Check to see if the max has been reached.  If so free all matches and return */
	if(s->maxReached == 1) {
		RGMatchFree(s);
		s->maxReached=1;
		return;
	}

	if(s->numEntries > 0) {
		/* Quick sort the data structure */
		RGMatchQuickSort(s, 0, s->numEntries-1);

		/* Remove duplicates */
		prevIndex=0;
		for(i=1;i<s->numEntries;i++) {
			if(RGMatchCompareAtIndex(s, prevIndex, s, i)==0) {
				/* ignore */
			}
			else {
				prevIndex++;
				/* Copy to prevIndex (incremented) */
				RGMatchCopyAtIndex(s, i, s, prevIndex);
			}
		}

		/* Reallocate pair */
		/* does not make sense if there are no entries */
		RGMatchReallocate(s, prevIndex+1);

		/* Check to see if we have too many matches */
		if(maxMatches > 0 && s->numEntries > maxMatches) {
			RGMatchFree(s);
			s->maxReached=1;
			return;
		}
		else { 
			s->maxReached = 0;
		}
	}
}

/* TODO */
void RGMatchQuickSort(RGMatch *s, int32_t low, int32_t high)
{
	int32_t i;
	int32_t pivot=-1;
	RGMatch *temp;

	if(low < high) {
		/* Allocate memory for the temp used for swapping */
		temp=malloc(sizeof(RGMatch));
		RGMatchInitialize(temp);
		if(NULL == temp) {
			PrintError("RGMatchQuickSort",
					"temp",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		RGMatchAllocate(temp, 1);

		pivot = (low+high)/2;

		RGMatchCopyAtIndex(s, pivot, temp, 0);
		RGMatchCopyAtIndex(s, high, s, pivot);
		RGMatchCopyAtIndex(temp, 0, s, high);

		pivot = low;

		for(i=low;i<high;i++) {
			if(RGMatchCompareAtIndex(s, i, s, high) <= 0) {
				if(i!=pivot) {
					RGMatchCopyAtIndex(s, i, temp, 0);
					RGMatchCopyAtIndex(s, pivot, s, i);
					RGMatchCopyAtIndex(temp, 0, s, pivot);
				}
				pivot++;
			}
		}
		RGMatchCopyAtIndex(s, pivot, temp, 0);
		RGMatchCopyAtIndex(s, high, s, pivot);
		RGMatchCopyAtIndex(temp, 0, s, high);

		/* Free temp before the recursive call, otherwise we have a worst
		 * case of O(n) space (NOT IN PLACE) 
		 * */
		RGMatchFree(temp);
		free(temp);
		temp=NULL;

		RGMatchQuickSort(s, low, pivot-1);
		RGMatchQuickSort(s, pivot+1, high);
	}
}

/* TODO */
/* Append to to the end of the matches */
int32_t RGMatchRead(FILE *fp,
		char *sequenceName,
		char *sequence,
		char *pairedSequence,
		RGMatch *sequenceMatch,
		RGMatch *pairedSequenceMatch,
		int32_t pairedEnd,
		int32_t binaryInput)
{
	char *FnName = "RGMatchRead";
	int32_t i;
	int32_t tempInt;

	/* Read the matches from the input file */
	if(binaryInput == 0) {
		/* Read sequence name */
		if(fscanf(fp, "%s", sequenceName)==EOF) {
			return EOF;
		}

		/* Read first sequence */
		if(fscanf(fp, "%s", sequence)==EOF) {
			PrintError(FnName,
					"sequence",
					"Could not read in sequence",
					Exit,
					EndOfFile);
		}

		/* Read in if we have reached the maximum number of matches */
		if(fscanf(fp, "%d", &sequenceMatch->maxReached)==EOF) {
			PrintError(FnName,
					"sequenceMatch->maxReached",
					"Could not read in sequenceMatch->maxReached",
					Exit,
					EndOfFile);
		}

		/* Read in the number of matches */
		if(fscanf(fp, "%d", &sequenceMatch->numEntries)==EOF) {
			PrintError(FnName,
					"sequenceMatch->numEntries",
					"Could not read in sequenceMatch->numEntries",
					Exit,
					EndOfFile);
		}
		assert(sequenceMatch->numEntries >= 0);

		/* Allocate memory for the matches */
		RGMatchReallocate(sequenceMatch, sequenceMatch->numEntries);

		/* Read first sequence matches */
		for(i=0;i<sequenceMatch->numEntries;i++) {
			if(fscanf(fp, "%d %d %c", 
						&tempInt,
						&sequenceMatch->positions[i],
						&sequenceMatch->strand[i])==EOF) {
				PrintError(FnName,
						NULL,
						"Could not read in match",
						Exit,
						EndOfFile);
			}
			sequenceMatch->chromosomes[i] = tempInt;
		}

		/* Read Paired end if necessary */
		if(pairedEnd == 1) {
			/* Read paired sequence */
			if(fscanf(fp, "%s", pairedSequence)==EOF) {
				PrintError(FnName,
						"pairedSequence",
						"Could not read in pairedSequence",
						Exit,
						EndOfFile);
			}
			/* Read in if we have reached the maximum number of matches */
			if(fscanf(fp, "%d", &pairedSequenceMatch->maxReached)==EOF) {
				PrintError(FnName,
						"pairedSequenceMatch->maxReached",
						"Could not read in pairedSequenceMatch->maxReached",
						Exit,
						EndOfFile);
			}

			/* Read in the number of matches */
			if(fscanf(fp, "%d", &pairedSequenceMatch->numEntries)==EOF) {
				PrintError(FnName,
						"pairedSequenceMatch->numEntries",
						"Could not read in the number of paired matches",
						Exit,
						EndOfFile);
			}
			assert(pairedSequenceMatch->numEntries >= 0);

			/* Allocate memory for the matches */
			RGMatchReallocate(pairedSequenceMatch, pairedSequenceMatch->numEntries);

			/* Read first pairedSequence matches */
			for(i=0;i<pairedSequenceMatch->numEntries;i++) {
				if(fscanf(fp, "%d %d %c", 
							&tempInt,
							&pairedSequenceMatch->positions[i],
							&pairedSequenceMatch->strand[i])==EOF) {
					PrintError(FnName,
							NULL,
							"Could not read in the paired match",
							Exit,
							EndOfFile);
				}
				pairedSequenceMatch->chromosomes[i] = tempInt;
			}
		}
	}
	else {
		/* Read sequence name */
		tempInt=-1;
		if(fread(&tempInt, sizeof(int32_t), 1, fp)!=1) {
			if(feof(fp) != 0) {
				return EOF;
			}
			else {
				PrintError(FnName,
						"sequenceName lenth",
						"Could not read in sequence name length",
						Exit,
						ReadFileError);
			}
		}
		assert(tempInt>0);
		assert(tempInt <= SEQUENCE_NAME_LENGTH);
		sequenceName[0]='\0';
		if(fread(sequenceName, sizeof(int8_t), tempInt, fp)!=tempInt) {
			PrintError(FnName,
					"sequenceName",
					"Could not read in sequence name",
					Exit,
					ReadFileError);
		}
		sequenceName[tempInt]='\0';

		/* Read first sequence */
		tempInt=0;
		if(fread(&tempInt, sizeof(int32_t), 1, fp)!=1) {
			PrintError(FnName,
					NULL,
					"Could not read in sequence length",
					Exit,
					ReadFileError);
		}
		assert(tempInt>0);
		assert(tempInt <= SEQUENCE_LENGTH);
		sequence[0]='\0';
		if(fread(sequence, sizeof(int8_t), tempInt, fp)!=tempInt) {
			PrintError(FnName,
					"sequence",
					"Could not read in sequence",
					Exit,
					ReadFileError);
		}
		sequence[tempInt]='\0';

		/* Read in if we have reached the maximum number of matches */
		if(fread(&sequenceMatch->maxReached, sizeof(int32_t), 1, fp)!=1) {
			PrintError(FnName,
					"sequenceMatch->maxReached",
					"Could not read in sequenceMatch->maxReached",
					Exit,
					ReadFileError);
		}
		assert(sequenceMatch->maxReached == 0 || sequenceMatch->maxReached == 1);

		/* Read in the number of matches */
		if(fread(&sequenceMatch->numEntries, sizeof(int32_t), 1, fp)!=1) {
			PrintError(FnName,
					"sequenceMatch->numEntries",
					"Could not read in sequenceMatch->numEntries",
					Exit,
					ReadFileError);
		}
		assert(sequenceMatch->numEntries >= 0);

		/* Allocate memory for the matches */
		RGMatchReallocate(sequenceMatch, sequenceMatch->numEntries);

		/* Read first sequence matches */
		if(fread(sequenceMatch->chromosomes, sizeof(uint8_t), sequenceMatch->numEntries, fp)!=sequenceMatch->numEntries) {
			PrintError(FnName,
					"sequenceMatch->chromosomes",
					"Could not read in chromosomes",
					Exit,
					ReadFileError);
		}
		if(fread(sequenceMatch->positions, sizeof(uint32_t), sequenceMatch->numEntries, fp)!=sequenceMatch->numEntries) {
			PrintError(FnName,
					"sequenceMatch->positions",
					"Could not read in positions",
					Exit,
					ReadFileError);
		}
		if(fread(sequenceMatch->strand, sizeof(int8_t), sequenceMatch->numEntries, fp)!=sequenceMatch->numEntries) {
			PrintError(FnName,
					"sequenceMatch->strand",
					"Could not read in strand",
					Exit,
					ReadFileError);
		}

		/* Read Paired end if necessary */
		if(pairedEnd == 1) {

			/* Read first pairedSequence */
			if(fread(&tempInt, sizeof(int32_t), 1, fp)!=1) {
				PrintError(FnName,
						NULL,
						"Could not read in pairedSequence length",
						Exit,
						ReadFileError);
			}
			if(fread(pairedSequence, sizeof(int8_t), tempInt, fp)!=tempInt) {
				PrintError(FnName,
						"pairedSequence",
						"Could not read in pairedSequence",
						Exit,
						ReadFileError);
			}
			pairedSequence[tempInt]='\0';

			/* Read in if we have reached the maximum number of matches */
			if(fread(&pairedSequenceMatch->maxReached, sizeof(int32_t), 1, fp)!=1) {
				PrintError(FnName,
						"pairedSequenceMatch->maxReached",
						"Could not read in pairedSequenceMatch->maxReached",
						Exit,
						ReadFileError);
			}

			/* Read in the number of matches */
			if(fread(&pairedSequenceMatch->numEntries, sizeof(int32_t), 1, fp)!=1) {
				PrintError(FnName,
						"pairedSequenceMatch->numEntries",
						"Could not read in pairedSequenceMatch->numEntries",
						Exit,
						ReadFileError);
			}
			assert(pairedSequenceMatch->numEntries >= 0);

			/* Allocate memory for the matches */
			RGMatchReallocate(pairedSequenceMatch, pairedSequenceMatch->numEntries);

			/* Read first pairedSequence matches */
			if(fread(pairedSequenceMatch->chromosomes, sizeof(uint8_t), pairedSequenceMatch->numEntries, fp)!=pairedSequenceMatch->numEntries) {
				PrintError(FnName,
						NULL,
						"Could not read in chromosome",
						Exit,
						ReadFileError);
			}
			if(fread(pairedSequenceMatch->positions, sizeof(uint32_t), pairedSequenceMatch->numEntries, fp)!=pairedSequenceMatch->numEntries) {
				PrintError(FnName,
						NULL,
						"Could not read in positions",
						Exit,
						ReadFileError);
			}
			if(fread(pairedSequenceMatch->strand, sizeof(int8_t), pairedSequenceMatch->numEntries, fp)!=pairedSequenceMatch->numEntries) {
				PrintError(FnName,
						NULL,
						"Could not read in strand",
						Exit,
						ReadFileError);
			}
		}
	}

	return 1;
}

/* TODO */
void RGMatchPrint(FILE *fp,
		char *sequenceName,
		char *sequence,
		char *pairedSequence,
		RGMatch *sequenceMatch,
		RGMatch *pairedSequenceMatch,
		int32_t pairedEnd,
		int32_t binaryOutput)
{
	char *FnName = "RGMatchPrint";
	int32_t i;
	int32_t tempInt;
	assert(fp!=NULL);
	/* Print the matches to the output file */

	if(binaryOutput == 0) {

		/* Print sequence name */
		fprintf(fp, "%s\n", sequenceName);
		/* Print first sequence */
		fprintf(fp, "%s", sequence);
		/* Print if the maximum number of matches was reached */
		fprintf(fp, "\t%d", sequenceMatch->maxReached);
		/* Print the number of matches */
		fprintf(fp, "\t%d", sequenceMatch->numEntries);
		/* Print first sequence matches */
		for(i=0;i<sequenceMatch->numEntries;i++) {
			fprintf(fp, "\t%d\t%d\t%c", 
					sequenceMatch->chromosomes[i],
					sequenceMatch->positions[i],
					sequenceMatch->strand[i]);
		}
		fprintf(fp, "\n");

		/* Print Paired end if necessary */
		if(pairedEnd == 1) {
			/* Print paired sequence */
			fprintf(fp, "%s", pairedSequence);
			/* Print if the maximum number of matches was reached */
			fprintf(fp, "\t%d", pairedSequenceMatch->maxReached);
			/* Print the number of matches for the paired end */
			fprintf(fp, "\t%d", pairedSequenceMatch->numEntries);
			/* Print first pairedSequence matches */
			for(i=0;i<pairedSequenceMatch->numEntries;i++) {
				fprintf(fp, "\t%d\t%d\t%c", 
						pairedSequenceMatch->chromosomes[i],
						pairedSequenceMatch->positions[i],
						pairedSequenceMatch->strand[i]);
			}
			fprintf(fp, "\n");
		}
	}
	else {
		/* Print sequence name */
		tempInt = strlen(sequenceName);
		if(fwrite(&tempInt, sizeof(int32_t), 1, fp) != 1 ||
				fwrite(sequenceName, sizeof(int8_t), tempInt, fp) != tempInt) {
			PrintError(FnName,
					NULL,
					"Could not write sequence name",
					Exit,
					WriteFileError);
		}

		/* Print first sequence */
		tempInt = strlen(sequence);
		if(fwrite(&tempInt, sizeof(int32_t), 1, fp) != 1 ||
				fwrite(sequence, sizeof(int8_t), tempInt, fp) != tempInt) {
			PrintError(FnName,
					NULL,
					"Could not write sequence",
					Exit,
					WriteFileError);
		}

		/* Print if the maximum number of matches was reached */
		if(fwrite(&sequenceMatch->maxReached, sizeof(int32_t), 1, fp) != 1 ||
				/* Print the number of matches */
				fwrite(&sequenceMatch->numEntries, sizeof(int32_t), 1, fp) != 1 || 
				/* Print first sequence matches */
				fwrite(sequenceMatch->chromosomes, sizeof(uint8_t), sequenceMatch->numEntries, fp) != sequenceMatch->numEntries ||
				fwrite(sequenceMatch->positions, sizeof(uint32_t), sequenceMatch->numEntries, fp) != sequenceMatch->numEntries ||
				fwrite(sequenceMatch->strand, sizeof(int8_t), sequenceMatch->numEntries, fp) != sequenceMatch->numEntries) {
			PrintError(FnName,
					NULL,
					"Could not write RGMatch",
					Exit,
					WriteFileError);
		}

		/* Print Paired end if necessary */
		if(pairedEnd == 1) {
			/* Print first paired sequence */
			tempInt = strlen(pairedSequence);
			if(fwrite(&tempInt, sizeof(int32_t), 1, fp) != 1 ||
					fwrite(pairedSequence, sizeof(int8_t), tempInt, fp) != tempInt) {
				PrintError(FnName,
						NULL,
						"Could not write paired sequence",
						Exit,
						WriteFileError);
			}

			/* Print if the maximum number of matches was reached */
			if(fwrite(&pairedSequenceMatch->maxReached, sizeof(int32_t), 1, fp) != 1 ||
					/* Print the number of matches */
					fwrite(&pairedSequenceMatch->numEntries, sizeof(int32_t), 1, fp) != 1 ||
					/* Print first paired sequence matches */
					fwrite(pairedSequenceMatch->chromosomes, sizeof(uint8_t), pairedSequenceMatch->numEntries, fp) != pairedSequenceMatch->numEntries ||
					fwrite(pairedSequenceMatch->positions, sizeof(uint32_t), pairedSequenceMatch->numEntries, fp) != pairedSequenceMatch->numEntries ||
					fwrite(pairedSequenceMatch->strand, sizeof(int8_t), pairedSequenceMatch->numEntries, fp) != pairedSequenceMatch->numEntries) {
				PrintError(FnName,
						NULL,
						"Could not write paired end RGMatch",
						Exit,
						WriteFileError);
			}
		}
	}
}

/* TODO */
int32_t RGMatchMergeFilesAndOutput(FILE **tempFPs,
		int32_t numFiles,
		FILE *outputFP,
		int32_t pairedEnd,
		int32_t binaryOutput,
		int32_t maxMatches)
{
	char *FnName="RGMatchMergeFilesAndOutput";
	int32_t i;
	int32_t counter;
	RGMatch match;
	RGMatch pairedMatch;
	RGMatch tempMatch;
	RGMatch tempPairedMatch;
	int32_t numFinished = 0;
	char **sequenceNames;
	char **sequences;
	char **pairedSequences;
	int32_t numMatches=0;

	/* Initialize matches */
	RGMatchInitialize(&match);
	RGMatchInitialize(&pairedMatch);
	RGMatchInitialize(&tempMatch);
	RGMatchInitialize(&tempPairedMatch);

	/* Allocate memory for the sequenceNames, sequences and pairedSequences */
	sequenceNames = malloc(sizeof(char*)*numFiles);
	if(NULL == sequenceNames) {
		PrintError(FnName,
				"sequenceNames",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	sequences= malloc(sizeof(char*)*numFiles);
	if(NULL == sequences) {
		PrintError(FnName,
				"sequences",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	pairedSequences = malloc(sizeof(char*)*numFiles);
	if(NULL == pairedSequences) {
		PrintError(FnName,
				"pairedSequences",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	for(i=0;i<numFiles;i++) {
		sequenceNames[i] = malloc(sizeof(char)*SEQUENCE_NAME_LENGTH);
		if(NULL == sequenceNames[i]) {
			PrintError(FnName,
					"sequenceNames[i]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		sequences[i] = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == sequences[i]) {
			PrintError(FnName,
					"sequences[i]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		pairedSequences[i] = malloc(sizeof(char)*SEQUENCE_LENGTH);
		if(NULL == pairedSequences[i]) {
			PrintError(FnName,
					"pairedSequences[i]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}

	/* Seek to the beginning of the files */
	for(i=0;i<numFiles;i++) {
		fseek(tempFPs[i], 0, SEEK_SET);
	}

	/* Read in each sequence/match one at a time */
	counter = 0;
	if(VERBOSE >=0) {
		fprintf(stderr, "\r%d", 0);
	}
	while(numFinished == 0) {
		if(VERBOSE >=0 && counter%RGMATCH_MERGE_ROTATE_NUM == 0) {
			fprintf(stderr, "\r%d", counter);
		}
		counter++;

		/* Read matches for one read from each file */ 
		for(i=0;i<numFiles;i++) {

			if(RGMatchRead(tempFPs[i],
						sequenceNames[i],
						sequences[i],
						pairedSequences[i],
						&tempMatch,
						&tempPairedMatch,
						pairedEnd,
						binaryOutput)==EOF) {
				numFinished++;
			}
			else {
				/* Append temp matches to matches */
				RGMatchAppend(&tempMatch, &match, maxMatches);
				if(pairedEnd == 1) {
					RGMatchAppend(&tempPairedMatch, &pairedMatch, maxMatches);
				}
			}

			RGMatchFree(&tempMatch);
			if(pairedEnd == 1) {
				RGMatchFree(&tempPairedMatch);
			}
		}
		/* We must finish all at the same time */
		assert(numFinished == 0 || numFinished == numFiles);

		if(numFinished == 0) {

			/* Error checking */
			for(i=1;i<numFiles;i++) {
				/* Make sure we are reading the same sequence */
				if(strcmp(sequenceNames[i], sequenceNames[0])!=0) {
					PrintError(FnName,
							NULL,
							"Sequence names are not the same",
							Warn,
							OutOfRange);
					fprintf(stderr, "sequenceNames[%d]:%s\nsequenceNames[%d]:%s\n",
							i,
							sequenceNames[i],
							0,
							sequenceNames[0]);
				}
				assert(strcmp(sequenceNames[i], sequenceNames[0])==0);
			}

			/* Remove duplicates */
			RGMatchRemoveDuplicates(&match, maxMatches);
			if(pairedEnd==1) {
				RGMatchRemoveDuplicates(&pairedMatch, maxMatches);
			}

			/* Print to output file */
			if(match.numEntries > 0) {
				numMatches++;
			}

			RGMatchPrint(outputFP,
					sequenceNames[0],
					sequences[0],
					pairedSequences[0],
					&match,
					&pairedMatch,
					pairedEnd,
					binaryOutput);

		}
		/* Free memory */
		RGMatchFree(&match);
		if(pairedEnd == 1) {
			RGMatchFree(&pairedMatch);
		}
	}
	if(VERBOSE >=0) {
		fprintf(stderr, "\r%d... completed.\n", counter-1);
	}

	/* Free memory */
	for(i=0;i<numFiles;i++) {
		free(sequenceNames[i]);
		free(sequences[i]);
		free(pairedSequences[i]);
	}
	free(sequenceNames);
	free(sequences);
	free(pairedSequences);

	return numMatches;
}

/* TODO */
int32_t RGMatchMergeThreadTempFilesIntoOutputTempFile(FILE **threadFPs,
		int32_t numThreads,
		FILE *outputFP,
		int32_t pairedEnd,
		int32_t binaryOutput)
{
	char *FnName = "RGMatchMergeThreadTempFilesIntoOutputTempFile";
	int32_t counter;
	int32_t i;
	RGMatch match;
	RGMatch pairedMatch;
	int32_t numFinished;
	char *sequenceName;
	char *sequence;
	char *pairedSequence;
	int *finished=NULL;

	/* Initialize matches */
	RGMatchInitialize(&match);
	RGMatchInitialize(&pairedMatch);

	/* Initialize thread file pointers */
	for(i=0;i<numThreads;i++) {
		fseek(threadFPs[i], 0, SEEK_SET);
	}

	/* Allocate memory for the sequenceNames, sequences and pairedSequences */
	sequenceName = malloc(sizeof(char)*SEQUENCE_NAME_LENGTH);
	if(NULL == sequenceName) {
		PrintError(FnName,
				"sequenceName",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	sequence = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL == sequence) {
		PrintError(FnName,
				"sequence",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	pairedSequence = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL == pairedSequence) {
		PrintError(FnName,
				"pairedSequence",
				"Could not allocate memory",
				Exit,
				MallocMemory);
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
	while(numFinished < numThreads) {
		/* For each thread */
		for(i=0;i<numThreads;i++) {
			/* Only try reading from those that are not finished */
			if(0 == finished[i]) {
				if(RGMatchRead(threadFPs[i],
							sequenceName,
							sequence,
							pairedSequence,
							&match,
							&pairedMatch,
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
					counter++;

					RGMatchPrint(outputFP,
							sequenceName,
							sequence,
							pairedSequence,
							&match,
							&pairedMatch,
							pairedEnd,
							binaryOutput);

				}
				/* Free memory */
				RGMatchFree(&match);
				if(pairedEnd == 1) {
					RGMatchFree(&pairedMatch);
				}
			}
		}
	}
	assert(numFinished == numThreads);
	for(i=0;i<numThreads;i++) {
		assert(1==finished[i]);
	}

	/* Free memory */
	free(sequenceName);
	free(sequence);
	free(pairedSequence);
	free(finished);

	return counter;
}

/* TODO */
int32_t RGMatchCompareAtIndex(RGMatch *mOne, int32_t indexOne, RGMatch *mTwo, int32_t indexTwo) 
{
	assert(indexOne >= 0 && indexOne < mOne->numEntries);
	assert(indexTwo >= 0 && indexTwo < mTwo->numEntries);
	if(mOne->chromosomes[indexOne] < mTwo->chromosomes[indexTwo] ||
			(mOne->chromosomes[indexOne] == mTwo->chromosomes[indexTwo] && mOne->positions[indexOne] < mTwo->positions[indexTwo]) ||
			(mOne->chromosomes[indexOne] == mTwo->chromosomes[indexTwo] && mOne->positions[indexOne] == mTwo->positions[indexTwo] && mOne->strand[indexOne] < mTwo->strand[indexTwo])) {
		return -1;
	}
	else if(mOne->chromosomes[indexOne] ==  mTwo->chromosomes[indexTwo] && mOne->positions[indexOne] == mTwo->positions[indexTwo] && mOne->strand[indexOne] == mTwo->strand[indexTwo]) {
		return 0;
	}
	else {
		return 1;
	}
}

/* TODO */
void RGMatchAppend(RGMatch *src, RGMatch *dest, int maxMatches) 
{
	int32_t i, start;

	assert(src != dest);
	start = dest->numEntries;
	RGMatchReallocate(dest, dest->numEntries + src->numEntries);
	assert(dest->numEntries == start + src->numEntries);
	assert(start <= dest->numEntries);
	for(i=start;i<dest->numEntries;i++) {
		RGMatchCopyAtIndex(src, i-start, dest, i);
	}

	if(maxMatches > 0 && dest->numEntries > maxMatches) {
		dest->maxReached = 1;
	}

}

/* TODO */
void RGMatchCopyAtIndex(RGMatch *src, int32_t srcIndex, RGMatch *dest, int32_t destIndex)
{
	assert(srcIndex >= 0 && srcIndex < src->numEntries);
	assert(destIndex >= 0 && destIndex < dest->numEntries);

	if(src != dest || srcIndex != destIndex) {
		dest->positions[destIndex] = src->positions[srcIndex];
		dest->chromosomes[destIndex] = src->chromosomes[srcIndex];
		dest->strand[destIndex] = src->strand[srcIndex];
	}
}

void RGMatchAllocate(RGMatch *m, int32_t numEntries)
{
	assert(m->numEntries==0);
	m->numEntries = numEntries;
	assert(m->positions==NULL);
	m->positions = malloc(sizeof(uint32_t)*numEntries); 
	if(NULL == m->positions) {
		PrintError("RGMatchAllocate",
				"m->positions",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(m->chromosomes==NULL);
	m->chromosomes = malloc(sizeof(uint8_t)*numEntries); 
	if(NULL == m->chromosomes) {
		PrintError("RGMatchAllocate",
				"m->chromosomes",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	assert(m->strand==NULL);
	m->strand = malloc(sizeof(int8_t)*numEntries); 
	if(NULL == m->strand) {
		PrintError("RGMatchAllocate",
				"m->strand",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
}

void RGMatchReallocate(RGMatch *m, int32_t numEntries)
{
	if(numEntries > 0) {
		m->numEntries = numEntries;
		m->positions = realloc(m->positions, sizeof(uint32_t)*numEntries); 
		if(numEntries > 0 && NULL == m->positions) {
			fprintf(stderr, "numEntries:%d\n", numEntries);
			PrintError("RGMatchReaocate",
					"m->positions",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		m->chromosomes = realloc(m->chromosomes, sizeof(uint8_t)*numEntries); 
		if(numEntries > 0 && NULL == m->chromosomes) {
			PrintError("RGMatchReaocate",
					"m->chromosomes",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		m->strand = realloc(m->strand, sizeof(int8_t)*numEntries); 
		if(numEntries > 0 && NULL == m->strand) {
			PrintError("RGMatchReaocate",
					"m->strand",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
	}
	else {
		RGMatchFree(m);
	}
}

void RGMatchFree(RGMatch *m) 
{
	if(m->numEntries>0) {
		free(m->positions);
		free(m->chromosomes);
		free(m->strand);
	}
	RGMatchInitialize(m);
}

void RGMatchInitialize(RGMatch *m)
{
	m->numEntries=0;
	m->maxReached=0;
	m->positions=NULL;
	m->chromosomes=NULL;
	m->strand=NULL;
}

void RGMatchMirrorPairedEnd(RGMatch *src, RGMatch *dest, int length) 
{
	int i;

	assert(dest->numEntries == 0);
	/* Reallocate paired match */
	RGMatchInitialize(dest);
	RGMatchReallocate(dest, src->numEntries);
	/* Copy over */
	for(i=0;i<src->numEntries;i++) {
		dest->chromosomes[i] = src->chromosomes[i];
		dest->strand[i] = src->strand[i];
		/* Adjust position */
		dest->positions[i] = src->positions[i] + length;
	}
}

