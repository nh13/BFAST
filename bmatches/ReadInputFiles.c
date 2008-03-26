#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
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

	/* open sequence file */
	if((seqFP=fopen(sequenceFileName, "r"))==0) {
		fprintf(stderr, "Error opening %s for reading.  Terminating!\n", sequenceFileName);
		exit(1);
	}

	/* open a temporary file */
	(*tempSeqFP)=tmpfile();

	while((endReadNum<=0 || endReadNum >= curReadNum) && EOF != fscanf(seqFP, "%s", sequenceName)) {

		/* Read sequence label above */
		/* Read sequence */
		if(EOF==fscanf(seqFP, "%s", sequence)) {
			fprintf(stderr, "Error.  Could not read sequence from %s.  Terminating!\n", sequenceFileName);
		}
		/* Read paired sequence */
		if(pairedEnd==1) {
			if(EOF==fscanf(seqFP, "%s", pairedSequence)) {
				fprintf(stderr, "Error.  Could not read sequence from %s.  Terminating!\n", sequenceFileName);
			}
		}
		/* Print only if we are within the desired limit */
		if(startReadNum<=0 || curReadNum >= startReadNum) {
			/* Print sequence */
			fprintf((*tempSeqFP), "%s", sequenceName);
			/* Print sequence */
			fprintf((*tempSeqFP), "%s", sequence);
			/* Print paired sequence */
			if(pairedEnd==1) {
				fprintf((*tempSeqFP), "%s", pairedSequence);
			}
		}
		/* Increment sequence number */
		curReadNum++;
	}

	/* close sequence file */
	fclose(seqFP);

	/* reset pointer to temp file to the beginning of the file */
	fseek((*tempSeqFP), 0, SEEK_SET);

}

/* TODO */
/* Read the next sequence from the stream */
int ReadNextSequence(FILE *fp, char **sequenceOne, char **sequenceTwo, char**sequenceName, int pairedEnd)
{
	if(pairedEnd == 1) {
		if(EOF==fscanf(fp, "%s", (*sequenceName)) || EOF==fscanf(fp, "%s", (*sequenceOne)) || EOF==fscanf(fp, "%s", (*sequenceTwo))) {
			return EOF;
		}
	}
	else {
		if(EOF==fscanf(fp, "%s", (*sequenceName)) || EOF==fscanf(fp, "%s", (*sequenceOne))) {
			return EOF;
		}
		(*sequenceTwo)=NULL;
	}
	return 1;
}

/* TODO */
/* Reads in a RGTree from file */
void ReadRGTree(char *rgTreeFileName, RGTree *tree)
{
	FILE *fp;

	/* open file */
	if((fp=fopen(rgTreeFileName, "r"))==0) {
		fprintf(stderr, "Error opening %s for reading.  Terminating!\n", rgTreeFileName);
		exit(1);
	}

	/* Read from file */
	RGTreeReadFromFile(tree, fp);

	/* close file */
	fclose(fp);
}

int ReadRGTreeFileNames(char *rgTreeListFileName, char ***rgTreeFileNames, int **offsets, int *numOffsets)
{
	char tempFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fp;
	int numFileNames=0;

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
	(*offsets)=(int*)malloc(sizeof(int)*(*numOffsets));
	if(fread((*offsets), sizeof(int), (*numOffsets), fp)==EOF) {
		fprintf(stderr, "Error.  Could not read %d offsets in %s.  Terminating!\n", 
				(*numOffsets),
				rgTreeListFileName);
		exit(1);
	}

	/* Read in the file names */
	while(fscanf(fp, "%s", tempFileName)!=EOF) {
		numFileNames++;
		(*rgTreeFileNames) = (char**)realloc((*rgTreeFileNames), sizeof(char*)*numFileNames);
		(*rgTreeFileNames)[numFileNames-1] = (char*)malloc(sizeof(char)*MAX_FILENAME_LENGTH);
		strcpy((*rgTreeFileNames)[numFileNames-1], tempFileName);
	}

	/* close file */
	fclose(fp);

	return numFileNames;
}
