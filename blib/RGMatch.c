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

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "In GMatchRemoveDuplicates\n");
	}

	/* Check to see if the max has been reached.  If so free all matches and return */
	if(s->maxReached == 1) {
		RGMatchFree(s);
		s->maxReached=1;
		return;
	}

	if(s->numEntries > 0) {
		/* Quick sort the data structure */
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Quick sorting match\n");
			for(i=0;i<s->numEntries;i++) {
				fprintf(stderr, "%d\t%d\t%c\n",
						s->chromosomes[i],
						s->positions[i],
						s->strand[i]);
			}
		}
		RGMatchQuickSort(s, 0, s->numEntries-1);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Sorted!\n");
			for(i=0;i<s->numEntries;i++) {
				fprintf(stderr, "%d\t%d\t%c\n",
						s->chromosomes[i],
						s->positions[i],
						s->strand[i]);
			}
		}

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
		if(s->numEntries > maxMatches) {
			RGMatchFree(s);
			s->maxReached=1;
			return;
		}

	}
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Exiting RGMatchRemoveDuplicates\n");
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
void RGMatchOutputToFile(FILE *fp,
		char *sequenceName,
		char *sequence,
		char *pairedSequence,
		RGMatch *sequenceMatch,
		RGMatch *pairedSequenceMatch,
		int32_t pairedEnd,
		int32_t binaryOutput)
{
	int32_t i;
	int32_t tempInt;
	assert(fp!=NULL);
	/* Print the matches to the output file */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "In RGMatchOutputToFile.\n");
	}

	if(binaryOutput == 0) {

		/* Print sequence name */
		fprintf(fp, "%s\n", sequenceName);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "sequenceName:%s\n", sequenceName);
		}

		/* Print first sequence */
		fprintf(fp, "%s", sequence);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "sequence:%s", sequence);
		}

		/* Print if the maximum number of matches was reached */
		fprintf(fp, "\t%d", sequenceMatch->maxReached);

		/* Print the number of matches */
		fprintf(fp, "\t%d", sequenceMatch->numEntries);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "\t%d", sequenceMatch->numEntries);
		}

		/* Print first sequence matches */
		for(i=0;i<sequenceMatch->numEntries;i++) {
			fprintf(fp, "\t%d\t%d\t%c", 
					sequenceMatch->chromosomes[i],
					sequenceMatch->positions[i],
					sequenceMatch->strand[i]);
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "\t%d\t%d\t%c", 
						sequenceMatch->chromosomes[i],
						sequenceMatch->positions[i],
						sequenceMatch->strand[i]);
			}
		}
		fprintf(fp, "\n");
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "\n");
		}

		/* Print Paired end if necessary */
		if(pairedEnd == 1) {
			/* Print paired sequence */
			fprintf(fp, "%s", pairedSequence);

			/* Print if the maximum number of matches was reached */
			fprintf(fp, "\t%d", pairedSequenceMatch->maxReached);

			/* Print the number of matches for the paired end */
			fprintf(fp, "\t%d", pairedSequenceMatch->numEntries);

			/* Print first pairedSequence matches */
			fprintf(fp, "\t%d", pairedSequenceMatch->numEntries);
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
		tempInt = strlen(sequenceName)+1;
		assert(tempInt>0);
		fwrite(&tempInt, sizeof(int32_t), 1, fp);
		fwrite(sequenceName, sizeof(int8_t), tempInt, fp);

		/* Print first sequence */
		tempInt = strlen(sequence)+1;
		assert(tempInt>0);
		fwrite(&tempInt, sizeof(int32_t), 1, fp);
		fwrite(sequence, sizeof(int8_t), tempInt, fp);

		/* Print if the maximum number of matches was reached */
		fwrite(&sequenceMatch->maxReached, sizeof(int32_t), 1, fp);

		/* Print the number of matches */
		fwrite(&sequenceMatch->numEntries, sizeof(int32_t), 1, fp);

		/* Print first sequence matches */
		fwrite(sequenceMatch->chromosomes, sizeof(uint8_t), sequenceMatch->numEntries, fp);
		fwrite(sequenceMatch->positions, sizeof(uint32_t), sequenceMatch->numEntries, fp);
		fwrite(sequenceMatch->strand, sizeof(int8_t), sequenceMatch->numEntries, fp);

		/* Print Paired end if necessary */
		if(pairedEnd == 1) {
			/* Print first paired sequence */
			tempInt = strlen(pairedSequence)+1;
			fwrite(&tempInt, sizeof(int32_t), 1, fp);
			fwrite(pairedSequence, sizeof(int8_t), tempInt, fp);

			/* Print if the maximum number of matches was reached */
			fwrite(&pairedSequenceMatch->maxReached, sizeof(int32_t), 1, fp);

			/* Print the number of matches */
			fwrite(&pairedSequenceMatch->numEntries, sizeof(int32_t), 1, fp);

			/* Print first paired sequence matches */
			fwrite(pairedSequenceMatch->chromosomes, sizeof(uint8_t), pairedSequenceMatch->numEntries, fp);
			fwrite(pairedSequenceMatch->positions, sizeof(uint32_t), pairedSequenceMatch->numEntries, fp);
			fwrite(pairedSequenceMatch->strand, sizeof(int8_t), pairedSequenceMatch->numEntries, fp);
		}
	}
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "\nExiting RGMatchOutputToFile.\n");
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
	int32_t continueReading = 1;
	char **sequenceNames;
	char **sequences;
	char **pairedSequences;
	int32_t numMatches=0;

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
	/*
	if(VERBOSE >=0) {
		fprintf(stderr, "%d", 0);
	}
	*/
	while(continueReading == 1) {
		/*
		if(VERBOSE >=0 && counter%RGMATCH_MERGE_ROTATE_NUM == 0) {
			fprintf(stderr, "\r%d", counter);
		}
		*/
		counter++;

		/* Initialize match */
		match.positions=NULL;
		match.chromosomes=NULL;
		match.strand=NULL;
		match.numEntries=0;
		match.maxReached=0;
		pairedMatch.positions=NULL;
		pairedMatch.chromosomes=NULL;
		pairedMatch.strand=NULL;
		pairedMatch.numEntries=0;
		pairedMatch.maxReached=0;

		/* Read matches for one read from each file */ 
		for(i=0;continueReading==1 && i<numFiles;i++) {
			if(RGMatchGetNextFromFile(tempFPs[i],
						sequenceNames[i],
						sequences[i],
						pairedSequences[i],
						&match,
						&pairedMatch,
						pairedEnd,
						binaryOutput)==EOF) {
				continueReading=0;
			}
		}

		if(continueReading==1) {

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
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "Printing to the output file\n");
			}
			if(match.numEntries > 0) {
				numMatches++;
			}

			RGMatchOutputToFile(outputFP,
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
	/*
		if(VERBOSE >=0) {
			fprintf(stderr, "\r%d... completed.\n", counter);
		}
		*/

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
	int32_t i;
	RGMatch match;
	RGMatch pairedMatch;
	int32_t continueReading = 1;
	char *sequenceName;
	char *sequence;
	char *pairedSequence;

	/* Allocate memory for the sequenceNames, sequences and pairedSequences */
	sequenceName = malloc(sizeof(char)*SEQUENCE_NAME_LENGTH);
	if(NULL == sequenceName) {
		PrintError("RGMatchMergeThreadTempFilesIntoOutputTempFile",
				"sequenceName",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	sequence = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL == sequence) {
		PrintError("RGMatchMergeThreadTempFilesIntoOutputTempFile",
				"sequence",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	pairedSequence = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL == pairedSequence) {
		PrintError("RGMatchMergeThreadTempFilesIntoOutputTempFile",
				"pairedSequence",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	continueReading=1;
	while(continueReading == 1) {

		/* For each thread */
		for(i=0;i<numThreads && 1==continueReading;i++) {

			/* Initialize match */
			match.positions=NULL;
			match.chromosomes=NULL;
			match.strand=NULL;
			match.numEntries=0;
			match.maxReached=0;
			pairedMatch.positions=NULL;
			pairedMatch.chromosomes=NULL;
			pairedMatch.strand=NULL;
			pairedMatch.numEntries=0;
			pairedMatch.maxReached=0;

			if(RGMatchGetNextFromFile(threadFPs[i],
						sequenceName,
						sequence,
						pairedSequence,
						&match,
						&pairedMatch,
						pairedEnd,
						binaryOutput)==EOF) {
				continueReading=0;
			}

			if(continueReading==1) {

				RGMatchOutputToFile(outputFP,
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

	/* Free memory */
	free(sequenceName);
	free(sequence);
	free(pairedSequence);

	return 1;
}

/* TODO */
/* Append to to the end of the matches */
int32_t RGMatchGetNextFromFile(FILE *fp,
		char *sequenceName,
		char *sequence,
		char *pairedSequence,
		RGMatch *sequenceMatch,
		RGMatch *pairedSequenceMatch,
		int32_t pairedEnd,
		int32_t binaryInput)
{
	int32_t i;
	int32_t tempInt;
	int32_t tempMaxReached;
	int32_t tempStart;

	/* Read the matches from the input file */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "In RGMatchGetNextFromFile.\n");
	}

	if(binaryInput == 0) {
		/* Read sequence name */
		if(fscanf(fp, "%s", sequenceName)==EOF) {
			return EOF;
		}

		/* Read first sequence */
		if(fscanf(fp, "%s", sequence)==EOF) {
			PrintError("RGMatchGetNextFromFile",
					"sequence",
					"Could not read in sequence",
					Exit,
					EndOfFile);
		}

		/* Read in if we have reached the maximum number of matches */
		if(fscanf(fp, "%d", &tempMaxReached)==EOF) {
			PrintError("RGMatchGetNextFromFile",
					"sequenceMatch->maxReached",
					"Could not read in sequenceMatch->maxReached",
					Exit,
					EndOfFile);
		}
		/* Update max reached */
		if(tempMaxReached==1) {
			sequenceMatch->maxReached = 1;
		}

		/* Read in the number of matches */
		if(fscanf(fp, "%d", &tempInt)==EOF) {
			PrintError("RGMatchGetNextFromFile",
					"sequenceMatch->numEntries",
					"Could not read in sequenceMatch->numEntries",
					Exit,
					EndOfFile);
		}
		assert(tempInt >= 0);
		tempStart = sequenceMatch->numEntries;
		sequenceMatch->numEntries += tempInt;
		assert(sequenceMatch->numEntries >= 0);

		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "sequenceName:%s\nsequence:%s\nnumEntries:%d\n",
					sequenceName,
					sequence,
					sequenceMatch->numEntries);
		}

		/* Allocate memory for the matches */
		RGMatchReallocate(sequenceMatch, sequenceMatch->numEntries);

		/* Read first sequence matches */
		for(i=tempStart;i<sequenceMatch->numEntries;i++) {
			if(fscanf(fp, "%d %d %c", 
						&tempInt,
						&sequenceMatch->positions[i],
						&sequenceMatch->strand[i])==EOF) {
				PrintError("RGMatchGetNextFromFile",
						NULL,
						"Could not read in match",
						Exit,
						EndOfFile);
			}
			sequenceMatch->chromosomes[i] = tempInt;
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "chr%d:%d:%c\t",
						sequenceMatch->chromosomes[i],
						sequenceMatch->positions[i],
						sequenceMatch->strand[i]);
			}
		}
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "\n");
		}

		/* Read Paired end if necessary */
		if(pairedEnd == 1) {
			/* Read paired sequence */
			if(fscanf(fp, "%s", pairedSequence)==EOF) {
				PrintError("RGMatchGetNextFromFile",
						"pairedSequence",
						"Could not read in pairedSequence",
						Exit,
						EndOfFile);
			}
			/* Read in if we have reached the maximum number of matches */
			if(fscanf(fp, "%d", &tempMaxReached)==EOF) {
				PrintError("RGMatchGetNextFromFile",
						"pairedSequenceMatch->maxReached",
						"Could not read in pairedSequenceMatch->maxReached",
						Exit,
						EndOfFile);
			}
			/* Update max reached */
			if(tempMaxReached==1) {
				pairedSequenceMatch->maxReached = 1;
			}

			/* Read in the number of matches */
			if(fscanf(fp, "%d", &tempInt)==EOF) {
				PrintError("RGMatchGetNextFromFile",
						"pairedSequenceMatch->numEntries",
						"Could not read in the number of paired matches",
						Exit,
						EndOfFile);
			}
			assert(tempInt >= 0);
			tempStart = pairedSequenceMatch->numEntries;
			pairedSequenceMatch->numEntries += tempInt;
			assert(pairedSequenceMatch->numEntries >= 0);

			/* Allocate memory for the matches */
			RGMatchReallocate(pairedSequenceMatch, pairedSequenceMatch->numEntries);

			/* Read first pairedSequence matches */
			for(i=tempStart;i<pairedSequenceMatch->numEntries;i++) {
				if(fscanf(fp, "%d %d %c", 
							&tempInt,
							&pairedSequenceMatch->positions[i],
							&pairedSequenceMatch->strand[i])==EOF) {
					PrintError("RGMatchGetNextFromFile",
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
		tempInt=-1;
		if(fread(&tempInt, sizeof(int32_t), 1, fp)!=1) {
			assert(feof(fp)!=0);
			return EOF;
		}
		assert(tempInt>0);
		sequenceName[0]='\0';
		if(fread(sequenceName, sizeof(int8_t), tempInt, fp)==EOF) {
			PrintError("RGMatchGetNextFromFile",
					"sequenceName",
					"Could not read in sequence name",
					Exit,
					EndOfFile);
		}

		/* Read first sequence */
		if(fread(&tempInt, sizeof(int32_t), 1, fp)==EOF) {
			PrintError("RGMatchGetNextFromFile",
					NULL,
					"Could not read in sequence length",
					Exit,
					EndOfFile);
		}
		assert(tempInt>0);
		if(fread(sequence, sizeof(int8_t), tempInt, fp)==EOF) {
			PrintError("RGMatchGetNextFromFile",
					"sequence",
					"Could not read in sequence",
					Exit,
					EndOfFile);
		}

		/* Read in if we have reached the maximum number of matches */
		if(fread(&tempMaxReached, sizeof(int32_t), 1, fp)==EOF) {
			PrintError("RGMatchGetNextFromFile",
					"sequenceMatch->maxReached",
					"Could not read in sequenceMatch->maxReached",
					Exit,
					EndOfFile);
		}
		/* Update max reached */
		if(tempMaxReached==1) {
			sequenceMatch->maxReached = 1;
		}

		/* Read in the number of matches */
		if(fread(&tempInt, sizeof(int32_t), 1, fp)==EOF) {
			PrintError("RGMatchGetNextFromFile",
					"sequenceMatch->numEntries",
					"Could not read in sequenceMatch->numEntries",
					Exit,
					EndOfFile);
		}
		assert(tempInt >= 0);
		tempStart = sequenceMatch->numEntries;
		sequenceMatch->numEntries += tempInt;
		assert(sequenceMatch->numEntries >= 0);

		/* Allocate memory for the matches */
		RGMatchReallocate(sequenceMatch, sequenceMatch->numEntries);

		/* Read first sequence matches */
		for(i=tempStart;i<sequenceMatch->numEntries;i++) {
			if(fread(&tempInt, sizeof(uint8_t), 1, fp)==EOF) {
				PrintError("RGMatchGetNextFromFile",
						NULL,
						"Could not read in chromosome",
						Exit,
						EndOfFile);
			}
			sequenceMatch->chromosomes[i] = tempInt;
		}
		for(i=tempStart;i<sequenceMatch->numEntries;i++) {
			if(fread(&sequenceMatch->positions[i], sizeof(uint32_t), 1, fp)==EOF) {
				PrintError("RGMatchGetNextFromFile",
						NULL,
						"Could not read in position",
						Exit,
						EndOfFile);
			}
		}
		for(i=tempStart;i<sequenceMatch->numEntries;i++) {
			if(fread(&sequenceMatch->strand[i], sizeof(int8_t), 1, fp)==EOF) {
				PrintError("RGMatchGetNextFromFile",
						NULL,
						"Could not read in strand",
						Exit,
						EndOfFile);
			}
		}

		/* Read Paired end if necessary */
		if(pairedEnd == 1) {

			/* Read first pairedSequence */
			if(fread(&tempInt, sizeof(int32_t), 1, fp)==EOF) {
				PrintError("RGMatchGetNextFromFile",
						NULL,
						"Could not read in pairedSequence length",
						Exit,
						EndOfFile);
			}
			if(fread(pairedSequence, sizeof(int8_t), tempInt, fp)==EOF) {
				PrintError("RGMatchGetNextFromFile",
						"pairedSequence",
						"Could not read in pairedSequence",
						Exit,
						EndOfFile);
			}

			/* Read in if we have reached the maximum number of matches */
			if(fread(&tempMaxReached, sizeof(int32_t), 1, fp)==EOF) {
				PrintError("RGMatchGetNextFromFile",
						"pairedSequenceMatch->maxReached",
						"Could not read in pairedSequenceMatch->maxReached",
						Exit,
						EndOfFile);
			}
			/* Update max reached */
			if(tempMaxReached==1) {
				pairedSequenceMatch->maxReached = 1;
			}

			/* Read in the number of matches */
			if(fread(&tempInt, sizeof(int32_t), 1, fp)==EOF) {
				PrintError("RGMatchGetNextFromFile",
						"pairedSequenceMatch->numEntries",
						"Could not read in pairedSequenceMatch->numEntries",
						Exit,
						EndOfFile);
			}
			assert(tempInt >= 0);
			tempStart = pairedSequenceMatch->numEntries;
			pairedSequenceMatch->numEntries += tempInt;
			assert(pairedSequenceMatch->numEntries >= 0);

			/* Allocate memory for the matches */
			RGMatchReallocate(pairedSequenceMatch, pairedSequenceMatch->numEntries);

			/* Read first pairedSequence matches */
			for(i=tempStart;i<pairedSequenceMatch->numEntries;i++) {
				if(fread(&tempInt, sizeof(uint8_t), 1, fp)==EOF) {
					PrintError("RGMatchGetNextFromFile",
							NULL,
							"Could not read in chromosome",
							Exit,
							EndOfFile);
				}
				pairedSequenceMatch->chromosomes[i] = tempInt;
			}
			for(i=tempStart;i<pairedSequenceMatch->numEntries;i++) {
				if(fread(&pairedSequenceMatch->positions[i], sizeof(uint32_t), 1, fp)==EOF) {
					PrintError("RGMatchGetNextFromFile",
							NULL,
							"Could not read in position",
							Exit,
							EndOfFile);
				}
			}
			for(i=tempStart;i<pairedSequenceMatch->numEntries;i++) {
				if(fread(&pairedSequenceMatch->strand[i], sizeof(int8_t), 1, fp)==EOF) {
					PrintError("RGMatchGetNextFromFile",
							NULL,
							"Could not read in strand",
							Exit,
							EndOfFile);
				}
			}
		}
	}

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Exiting RGMatchGetNextFromFile.\n");
	}
	return 1;
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
void RGMatchCopyAtIndex(RGMatch *src, int32_t srcIndex, RGMatch *dest, int32_t destIndex)
{
	if(!(srcIndex >= 0 && srcIndex < src->numEntries)) {
		fprintf(stderr, "Error. srcIndex:%d\tnumEntries:%d\n",
				srcIndex,
				src->numEntries);
	}
	if(!(destIndex >= 0 && destIndex < dest->numEntries)) {
		fprintf(stderr, "Error. destIndex:%d\tnumEntries:%d\n",
				destIndex,
				dest->numEntries);
	}

	assert(srcIndex >= 0 && srcIndex < src->numEntries);
	assert(destIndex >= 0 && destIndex < dest->numEntries);

	dest->positions[destIndex] = src->positions[srcIndex];
	dest->chromosomes[destIndex] = src->chromosomes[srcIndex];
	dest->strand[destIndex] = src->strand[srcIndex];
}

void RGMatchAllocate(RGMatch *m, int32_t numEntries)
{
	m->numEntries = numEntries;
	m->positions = malloc(sizeof(uint32_t)*numEntries); 
	if(NULL == m->positions) {
		PrintError("RGMatchAllocate",
				"m->positions",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	m->chromosomes = malloc(sizeof(uint8_t)*numEntries); 
	if(NULL == m->chromosomes) {
		PrintError("RGMatchAllocate",
				"m->chromosomes",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
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

void RGMatchFree(RGMatch *m) 
{
	if(m->numEntries>0) {
		free(m->positions);
		m->positions=NULL;
		free(m->chromosomes);
		m->chromosomes=NULL;
		free(m->strand);
		m->strand=NULL;
	}
	m->numEntries=0;
	m->maxReached=0;
}
