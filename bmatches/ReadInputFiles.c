#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "../blib/SRTree.h"
#include "Definitions.h"
#include "ReadInputFiles.h"

/* TODO */
/* Reads the sequences into a SRTree */
void ReadSequences(char *sequencesFileName, SRNode *root, int pairedEnd)
{

	FILE *fp;

	/* NOT FULLY IMPLEMENTED
	   fprintf(stderr, "Error.  ReadSequences not fully implemented.  Terminating!\n");
	   exit(1);
	   */

	/* open file */
	if((fp=fopen(sequencesFileName, "r"))==0) {
		fprintf(stderr, "Error opening %s for reading.  Terminating!\n", sequencesFileName);
		exit(1);
	}

	SRTreeReadFromFile(root, fp);

	/* close file */
	fclose(fp);
}

/* TODO */
void ReadSequencesToTempFile(char *sequenceFileName,
		FILE **tempSeqFP,
		int startReadNum,
		int endReadNum,
		int pairedEnd)
{
	FILE *seqFP=NULL;
	char sequenceName[SEQUENCE_NAME_LENGTH];
	char sequence[SEQUENCE_LENGTH];
	char pairedSequence[SEQUENCE_LENGTH];
	int curReadNum = 1;

	if(VERBOSE > 0) {
	fprintf(stderr, "Reading sequences from %s to a temp file.\n",
		sequenceFileName);
	}

	/* open sequence file */
	if((seqFP=fopen(sequenceFileName, "r"))==0) {
		fprintf(stderr, "Error opening %s for reading.  Terminating!\n", sequenceFileName);
		exit(1);
	}

	/* open a temporary file */
	(*tempSeqFP)=tmpfile();

	/* NOTE: we could implement a counter */
	while((endReadNum<=0 || endReadNum >= curReadNum) && EOF != fscanf(seqFP, "%s", sequenceName)) {

		/* Read sequence label above */
		/* Read sequence */
		if(EOF==fscanf(seqFP, "%s", sequence)) {
			fprintf(stderr, "Error.  Could not read sequence from %s.  Terminating!\n", sequenceFileName);
		}
		SequenceToLower(sequence, strlen(sequence));
		/* Read paired sequence */
		if(pairedEnd==1) {
			if(EOF==fscanf(seqFP, "%s", pairedSequence)) {
				fprintf(stderr, "Error.  Could not read sequence from %s.  Terminating!\n", sequenceFileName);
			}
			SequenceToLower(pairedSequence, strlen(pairedSequence));
		}
		/* Print only if we are within the desired limit */
		if(startReadNum<=0 || curReadNum >= startReadNum) {
			/* Print sequence */
			fprintf((*tempSeqFP), "%s", sequenceName);
			/* Print sequence */
			fprintf((*tempSeqFP), "\t%s", sequence);
			/* Print paired sequence */
			if(pairedEnd==1) {
				fprintf((*tempSeqFP), "\t%s", pairedSequence);
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

	/* close sequence file */
	fclose(seqFP);

	/* reset pointer to temp file to the beginning of the file */
	fseek((*tempSeqFP), 0, SEEK_SET);

	if(VERBOSE > 0) {
		fprintf(stderr, "Read %d reads from %s.\n",
				curReadNum,
				sequenceFileName);
	}
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
/* Reads in a RGTree from file */
void ReadRGTree(char *rgTreeFileName, RGTree *tree)
{
	FILE *fp;

	if(VERBOSE >= 0) {
		fprintf(stderr, "Reading tree from %s.\n",
				rgTreeFileName);
	}

	/* open file */
	if((fp=fopen(rgTreeFileName, "r"))==0) {
		fprintf(stderr, "Error opening %s for reading.  Terminating!\n", rgTreeFileName);
		exit(1);
	}

	/* Read from file */
	RGTreeReadFromFile(tree, fp);

	/* close file */
	fclose(fp);

	if(VERBOSE >= 0) {
		fprintf(stderr, "Reading successful\n");
	}
}

int ReadRGTreeFileNames(char *rgTreeListFileName, char ***rgTreeFileNames, int **offsets, int *numOffsets)
{
	char tempFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fp;
	int numFileNames=0;
	int i;

	if(VERBOSE>0) {
		fprintf(stderr, "Reading in blatter tree file names from %s.\n",
				rgTreeListFileName);
	}

	/* open file */
	if((fp=fopen(rgTreeListFileName, "r"))==0) {
		fprintf(stderr, "Error opening %s for reading.  Terminating!\n", rgTreeListFileName);
		exit(1);
	}

	/* Read offsets */
	if(fscanf(fp, "%d", numOffsets)==EOF) {
		fprintf(stderr, "Error.  Could not read the number offsets in %s.  Terminating!\n", rgTreeListFileName);
		exit(1);
	}
	assert((*numOffsets)>0);
	if(VERBOSE>0) {
		fprintf(stderr, "Will read %d offsets from %s.\n",
				(*numOffsets),
				rgTreeListFileName);
	}
	/* Allocate memory for the offsets */
	(*offsets)=(int*)malloc(sizeof(int)*(*numOffsets));
	for(i=0;i<(*numOffsets);i++) {
		if(fscanf(fp, "%d", &(*offsets)[i])==EOF) {
			fprintf(stderr, "Error.  Could not read %d offsets in %s.  Terminating!\n", 
					(*numOffsets),
					rgTreeListFileName);
			exit(1);
		}
		if(VERBOSE>0) {
			fprintf(stderr, "Offset %d:%d\n", i+1, (*offsets)[i]);
		}
	}

	/* Read in the file names */
	while(fscanf(fp, "%s", tempFileName)!=EOF) {
		numFileNames++;
		(*rgTreeFileNames) = (char**)realloc((*rgTreeFileNames), sizeof(char*)*numFileNames);
		(*rgTreeFileNames)[numFileNames-1] = (char*)malloc(sizeof(char)*MAX_FILENAME_LENGTH);
		strcpy((*rgTreeFileNames)[numFileNames-1], tempFileName);
		if(VERBOSE>0) {
			fprintf(stderr, "blatter tree file name %d:%s\n", 
					numFileNames,
					(*rgTreeFileNames)[numFileNames-1]);
		}
	}

	/* close file */
	fclose(fp);

	if(VERBOSE>0) {
		fprintf(stderr, "Read %d blatter tree file name.\n", numFileNames);
	}

	return numFileNames;
}
