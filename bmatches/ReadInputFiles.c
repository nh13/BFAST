#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "../blib/BError.h"
#include "../blib/RGMatch.h"
#include "Definitions.h"
#include "ReadInputFiles.h"

/* TODO */
int ReadSequencesToTempFile(FILE *seqFP,
		FILE ***tempSeqFPs, /* One for each thread */ 
		int startReadNum, 
		int endReadNum, 
		int pairedEnd,
		int numThreads,
		int timing)
{
	int i;
	int curSeqFPIndex=0;
	char sequenceName[SEQUENCE_NAME_LENGTH];
	char sequence[SEQUENCE_LENGTH];
	char pairedSequence[SEQUENCE_LENGTH];
	int curReadNum = 1;
	time_t startTime = time(NULL);
	time_t endTime;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "In ReadSequencesToTempFile\n");
	}

	/* open one temporary file, one for each thread */
	for(i=0;i<numThreads;i++) {
		(*tempSeqFPs)[i] = tmpfile();
		assert((*tempSeqFPs)[i]!=NULL);
	}

	/* NOTE: we could implement a counter */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "\n0");
	}
	while((endReadNum<=0 || endReadNum >= curReadNum) && EOF != fscanf(seqFP, "%s", sequenceName)) {
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "\r%d",
					curReadNum);
		}
		/* Get which tmp file to put the read in */
		curSeqFPIndex = (curReadNum-1)%numThreads;

		/* Read sequence label above */
		/* Read sequence */
		if(EOF==fscanf(seqFP, "%s", sequence)) {
			PrintError("ReadSequencesToTempFile",
					"sequence",
					"Could not read sequence",
					Exit,
					EndOfFile);
		}
		SequenceToLower(sequence, strlen(sequence));
		/* Read paired sequence */
		if(pairedEnd==1) {
			if(EOF==fscanf(seqFP, "%s", pairedSequence)) {
				PrintError("ReadSequencesToTempFile",
						"pairedSequence",
						"Could not read sequence",
						Exit,
						EndOfFile);
			}
			SequenceToLower(pairedSequence, strlen(pairedSequence));
		}
		/* Print only if we are within the desired limit */
		if(startReadNum<=0 || curReadNum >= startReadNum) {
			/* Print sequence */
			fprintf((*tempSeqFPs)[curSeqFPIndex], "%s\n", sequenceName);
			/* Print sequence */
			fprintf((*tempSeqFPs)[curSeqFPIndex], "%s\n", sequence);
			/* Print paired sequence */
			if(pairedEnd==1) {
				fprintf((*tempSeqFPs)[curSeqFPIndex], "%s\n", pairedSequence);
			}
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "sequenceName:%s\tsequence:%s",
						sequenceName,
						sequence);
				if(pairedEnd==1) {
					fprintf(stderr, "\tpairedSequence:%s",
							pairedSequence);
				}
				fprintf(stderr, "\n");
			}
		}
		/* Increment sequence number */
		curReadNum++;
	}
	curReadNum--; /* decrement sequence number since we have +1 after the loop */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "\r%d\n",
				curReadNum);
	}

	/* close sequence file */
	fclose(seqFP);

	/* reset pointer to temp files to the beginning of the file */
	for(i=0;i<numThreads;i++) {
		fseek((*tempSeqFPs)[i], 0, SEEK_SET);
	}

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

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "Exiting ReadSequencesToTempFile\n");
	}

	return curReadNum;
}

/* TODO */
/* Read the next sequence from the stream */
int ReadNextSequence(FILE *fp, char **sequenceOne, char **sequenceTwo, char**sequenceName, int pairedEnd)
{
	if(pairedEnd == 1) {
		if(EOF==fscanf(fp, "%s", (*sequenceName)) || EOF==fscanf(fp, "%s", (*sequenceOne)) || EOF==fscanf(fp, "%s", (*sequenceTwo))) {
			return EOF;
		}
		SequenceToLower((*sequenceOne), strlen((*sequenceOne)));
		SequenceToLower((*sequenceTwo), strlen((*sequenceTwo)));
	}
	else {
		if(EOF==fscanf(fp, "%s", (*sequenceName)) || EOF==fscanf(fp, "%s", (*sequenceOne))) {
			return EOF;
		}
		SequenceToLower((*sequenceOne), strlen((*sequenceOne)));
		(*sequenceTwo)=NULL;
	}
	return 1;
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

/* TODO */
/* Go through the temporary output file and output those sequences that have 
 * at least one match to the final output file.  For those sequences that have
 * zero matches, output them to the temporary sequence file *
 * */
void ReadTempSequencesAndOutput(FILE *tempOutputFP,
		FILE *outputFP,
		FILE *tempSeqFP,
		int pairedEnd)
{
	char sequenceName[SEQUENCE_NAME_LENGTH];
	char sequence[SEQUENCE_LENGTH];
	char pairedSequence[SEQUENCE_LENGTH];
	RGMatch sequenceMatch;
	RGMatch pairedSequenceMatch;

	assert(pairedEnd==0);

	/* Initialize match structures */
	sequenceMatch.positions=NULL;
	sequenceMatch.chromosomes=NULL;
	sequenceMatch.strand=NULL;
	sequenceMatch.numEntries=0;
	pairedSequenceMatch.positions=NULL;
	pairedSequenceMatch.chromosomes=NULL;
	pairedSequenceMatch.strand=NULL;
	pairedSequenceMatch.numEntries=0;

	/* Go to the beginning of the temporary output file */
	fseek(tempOutputFP, 0, SEEK_SET);

	while(RGMatchGetNextFromFile(tempOutputFP, 
				sequenceName, 
				sequence, 
				pairedSequence, 
				&sequenceMatch, 
				&pairedSequenceMatch,
				pairedEnd)!=EOF) {
		if(sequenceMatch.numEntries > 0) {
			/* Output to final output file */
			RGMatchOutputToFile(outputFP,
					sequenceName,
					sequence,
					pairedSequence,
					&sequenceMatch,
					&pairedSequenceMatch,
					pairedEnd);
		}
		else {
			/* Put back in the sequence file */
			fprintf(tempSeqFP, "%s\n%s\n",
					sequenceName,
					sequence);
			if(pairedEnd==1) {
				fprintf(tempSeqFP, "%s\n",
						pairedSequence);
			}
		}

		/* Free match memory and reinitialize match structures */
		if(sequenceMatch.numEntries>0) {
			free(sequenceMatch.positions);
			sequenceMatch.positions=NULL;
			free(sequenceMatch.chromosomes);
			sequenceMatch.chromosomes=NULL;
			free(sequenceMatch.strand);
			sequenceMatch.strand=NULL;
		}
		sequenceMatch.numEntries=0;
		if(pairedSequenceMatch.numEntries > 0) {
			pairedSequenceMatch.positions=NULL;
			free(pairedSequenceMatch.positions);
			pairedSequenceMatch.chromosomes=NULL;
			free(pairedSequenceMatch.chromosomes);
			pairedSequenceMatch.strand=NULL;
			free(pairedSequenceMatch.strand);
		}
		pairedSequenceMatch.numEntries=0;
	}
}

/* TODO */
void SequenceToLower(char* sequence, int length)
{
	int i;
	for(i=0;i<length;i++) {
		switch(sequence[i]) {
			case 'A':
				sequence[i] = 'a';
				break;
			case 'C':
				sequence[i] = 'c';
				break;
			case 'G':
				sequence[i] = 'g';
				break;
			case 'T':
				sequence[i] = 't';
				break;
			default:
				/* do nothing */
				break;
		}
	}
}

