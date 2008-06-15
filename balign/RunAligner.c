#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include "ReadInputFiles.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "../blib/RGBinary.h"
#include "../blib/RGMatch.h" /* To read in the matches */
#include "../blib/AlignEntry.h" /* For output */
#include "Align.h"
#include "Definitions.h"
#include "RunAligner.h"

/* TODO */
void RunAligner(RGBinary *rgBinary,
		char *matchesFileName,
		char *scoringMatrixFileName,
		int algorithm,
		int offsetLength,
		int maxNumMatches,
		int pairedEnd,
		int binaryInput,
		int numThreads,
		char *outputID,
		char *outputDir,
		char *tmpDir)
{
	int i;
	FILE *outputFP=NULL;
	FILE *matchesFP=NULL;
	FILE *matchFP=NULL;
	char **matchFileNames=NULL;
	int numMatchFileNames=0;
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	char tempFileName[MAX_FILENAME_LENGTH]="\0";

	/* Open matches file */
	if((matchesFP=fopen(matchesFileName, "r"))==0) {
		PrintError("RunAligner",
				matchesFileName,
				"Could not open matchesFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Read in the list of file names*/
	while(fscanf(matchesFP, "%s", tempFileName)!=EOF) {
		numMatchFileNames++;
		/* Reallocate memory */
		matchFileNames = (char**)realloc(matchFileNames, sizeof(char*)*numMatchFileNames);
		if(NULL==matchFileNames) {
			PrintError("RunAligner",
					"matchFileNames",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		/* Allocate memory */
		matchFileNames[numMatchFileNames-1] = (char*)malloc(sizeof(char)*(strlen(tempFileName)+1));
		if(NULL==matchFileNames[numMatchFileNames-1]) {
			PrintError("RunAligner",
					"matchFileNames[numMatchFileNames-1]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Copy over file name */
		strcpy(matchFileNames[numMatchFileNames-1], tempFileName);
	}

	/* Create output file name */
	sprintf(outputFileName, "%sbfast.aligned.file.%s.%d.%s",
			outputDir,
			outputID,
			algorithm,
			BFAST_ALIGN_FILE_EXTENSION);

	/* Open output file */
	if((outputFP=fopen(outputFileName, "w"))==0) {
		PrintError("RunAligner",
				outputFileName,
				"Could not open outputFileName for writing",
				Exit,
				OpenFileError);
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "Will output to %s.\n", outputFileName);
	}

	for(i=0;i<numMatchFileNames;i++) { /* For each match file name */

		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Reading match file #%d from %s.\n",
					i+1,
					matchFileNames[i]);
		}

		/* Open current match file */
		if((matchFP=fopen(matchFileNames[i], "r"))==0) {
			PrintError("RunAligner",
					matchFileNames[i],
					"Could not open matchesFileNames[] for reading",
					Exit,
					OpenFileError);
		}

		/* Run selected algorithm */
		switch(algorithm) {
			case 0:
				RunDynamicProgramming(matchFP,
						rgBinary,
						scoringMatrixFileName,
						offsetLength,
						maxNumMatches,
						pairedEnd,
						binaryInput,
						numThreads,
						tmpDir,
						outputFP);
				break;
			default:
				PrintError("RunAligner",
						NULL,
						"Could not understand algorithm option",
						Exit,
						OutOfRange);
				break;
		}
		/* Close the match file */
		fclose(matchFP);
		matchFP=NULL;
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Close output file */
	fclose(outputFP);

	/* Free memory */
	for(i=0;i<numMatchFileNames;i++) {
		free(matchFileNames[i]);
	}
	free(matchFileNames);
}

/* TODO */
void RunDynamicProgramming(FILE *matchFP,
		RGBinary *rgBinary,
		char *scoringMatrixFileName,
		int offsetLength,
		int maxNumMatches,
		int pairedEnd,
		int binaryInput,
		int numThreads,
		char *tmpDir,
		FILE *outputFP)
{
	/* local variables */
	ScoringMatrix sm;
	char readName[SEQUENCE_NAME_LENGTH]="\0";
	char read[SEQUENCE_LENGTH]="\0";
	char pairedSequence[SEQUENCE_LENGTH]="\0";
	RGMatch readMatch;
	RGMatch pairedReadMatch;
	int i, j;
	int numMatches=0;
	int continueReading=0;
	int numAlignEntries;
	int numAlignments=0;
	AlignEntry aEntry;
	/* Thread specific data */
	ThreadData *data;
	pthread_t *threads=NULL;
	int errCode;
	void *status;

	/* Initialize match */
	RGMatchInitialize(&readMatch);
	RGMatchInitialize(&pairedReadMatch);

	/* Allocate memory for thread arguments */
	data = malloc(sizeof(ThreadData)*numThreads);
	if(NULL==data) {
		PrintError("RunDynamicProgramming",
				"data",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for thread ids */
	threads = malloc(sizeof(pthread_t)*numThreads);
	if(NULL==threads) {
		PrintError("RunDynamicProgramming",
				"threads",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Read in scoring matrix */
	ReadScoringMatrix(scoringMatrixFileName, &sm); 

	/**/
	/* Split the input file into equal temp files for each thread */
	/**/

	/* Open temp files for the threads */
	for(i=0;i<numThreads;i++) {
		data[i].inputFP = OpenTmpFile(tmpDir, &data[i].inputFileName);
		data[i].outputFP = OpenTmpFile(tmpDir, &data[i].outputFileName);
	}

	/* Go through each read in the match file */
	numMatches=0;
	while(EOF!=RGMatchRead(matchFP, 
				readName, 
				read, 
				pairedSequence, 
				&readMatch,
				&pairedReadMatch,
				pairedEnd,
				binaryInput)) {
		/* Get the thread index - do this BEFORE incrementing */
		int threadIndex = numMatches%numThreads;
		/* increment */
		numMatches++;

		/* Print match to temp file */
		RGMatchPrint(data[threadIndex].inputFP,
				readName,
				read,
				pairedSequence,
				&readMatch,
				&pairedReadMatch,
				pairedEnd,
				binaryInput);

		/* Free match */
		RGMatchFree(&readMatch);
		if(pairedEnd==1) {
			RGMatchFree(&pairedReadMatch);
		}
	}

	/* Create thread arguments */
	for(i=0;i<numThreads;i++) {
		fseek(data[i].inputFP, 0, SEEK_SET);
		fseek(data[i].outputFP, 0, SEEK_SET);
		data[i].rgBinary=rgBinary;
		data[i].offsetLength=offsetLength;
		data[i].maxNumMatches=maxNumMatches;
		data[i].pairedEnd=pairedEnd;
		data[i].binaryInput=binaryInput;
		data[i].sm = &sm;
		data[i].threadID = i;
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "Performing alignment...\n");
		fprintf(stderr, "Currently on:\n0");
	}

	/* Create threads */
	for(i=0;i<numThreads;i++) {
		/* Start thread */
		errCode = pthread_create(&threads[i], /* thread struct */
				NULL, /* default thread attributes */
				RunDynamicProgrammingThread, /* start routine */
				&data[i]); /* data to routine */
		if(0!=errCode) {
			PrintError("RunDynamicProgramming",
					"pthread_create: errCode",
					"Could not start thread",
					Exit,
					ThreadError);
		}
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Created threadID:%d\n",
					i);
		}
	}

	/* Wait for the threads to finish */
	for(i=0;i<numThreads;i++) {
		/* Wait for the given thread to return */
		errCode = pthread_join(threads[i],
				&status);
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "Thread returned with errCode:%d\n",
					errCode);
		}
		/* Check the return code of the thread */
		if(0!=errCode) {
			PrintError("FindMatchesInIndexes",
					"pthread_join: errCode",
					"Thread returned an error",
					Exit,
					ThreadError);
		}
		/* Reinitialize file pointer */
		fseek(data[i].outputFP, 0, SEEK_SET);
		assert(NULL!=data[i].outputFP);
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Alignment complete.\n");
	}

	/* Close tmp input files */
	for(i=0;i<numThreads;i++) {
		CloseTmpFile(&data[i].inputFP, &data[i].inputFileName);
	}

	/* Allocate memory for the align entries */
	aEntry.read = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL==aEntry.read) {
		PrintError("RunDynamicProgramming",
				"aEntry.read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	aEntry.reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL==aEntry.reference) {
		PrintError("RunDynamicProgramming",
				"aEntry.reference",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	if(VERBOSE >=0) {
		fprintf(stderr, "Outputting...\n");
	}

	/* Merge all the outputs together */
	numAlignments=0;
	continueReading=1;
	while(continueReading==1) {
		/* Get one align from each thread */
		continueReading=0;
		for(i=0;i<numThreads;i++) {
			/* First get the number of align entries */
			if(EOF != fscanf(data[i].outputFP, "%d", &numAlignEntries)) {
				continueReading=1;
				assert(numAlignEntries >= 0);
				numAlignments++;
				for(j=0;j<numAlignEntries;j++) {
					/* Read in the align entry */
					if(EOF==AlignEntryRead(&aEntry,
								data[i].outputFP)) {
						PrintError("RunDynamicProgramming",
								"AlignEntryRead",
								"Could not read in the align entry",
								Exit,
								EndOfFile);
					}
					/* Print the align entry to file */
					AlignEntryPrint(&aEntry,
							outputFP);
				}
			}
		}
	}

	/* Close tmp output files */
	for(i=0;i<numThreads;i++) {
		CloseTmpFile(&data[i].outputFP, &data[i].outputFileName);
	}

	if(VERBOSE >=0) {
		fprintf(stderr, "Outputted alignments for %d reads.\n", numAlignments);
		fprintf(stderr, "Outputting complete.\n");
	}

	/* Free memory */
	free(data);
	free(threads);
	free(aEntry.read);
	free(aEntry.reference);
	/* Free scores */
	free(sm.key);
	for(i=0;i<ALPHABET_SIZE+1;i++) {
		free(sm.scores[i]);
	}
	free(sm.scores);
}

/* TODO */
void *RunDynamicProgrammingThread(void *arg)
{
	/* Recover arguments */
	ThreadData *data = (ThreadData *)(arg);
	FILE *inputFP=data->inputFP;
	FILE *outputFP=data->outputFP;
	RGBinary *rgBinary=data->rgBinary;
	int offsetLength=data->offsetLength;
	int maxNumMatches=data->maxNumMatches;
	int pairedEnd=data->pairedEnd;
	int binaryInput=data->binaryInput;
	int threadID=data->threadID;
	ScoringMatrix *sm = data->sm;
	/* Local variables */
	AlignEntry *aEntry=NULL;
	int numAlignEntries=0;
	char readName[SEQUENCE_NAME_LENGTH]="\0";
	char read[SEQUENCE_LENGTH]="\0";
	char reverseRead[SEQUENCE_LENGTH]="\0";
	char pairedSequence[SEQUENCE_LENGTH]="\0";
	RGMatch readMatch;
	RGMatch pairedReadMatch;
	int matchLength;
	char *reference=NULL;
	int referenceLength=0;
	int adjustPosition=0;
	int i;
	int position;
	int numMatches=0;
	int numMatchesAligned=0;

	/* Initialize match */
	RGMatchInitialize(&readMatch);
	RGMatchInitialize(&pairedReadMatch);

	/* Allocate memory for the reference */
	referenceLength = 2*offsetLength + SEQUENCE_LENGTH + 1;
	reference = (char*)malloc(sizeof(char)*(referenceLength+1));
	if(NULL==reference) {
		PrintError("RunDynamicProgrammingThread",
				"reference",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	reference[referenceLength] = '\0'; /* Add null terminator */

	/* Go through each read in the match file */
	while(EOF!=RGMatchRead(inputFP, 
				readName, 
				read, 
				pairedSequence, 
				&readMatch,
				&pairedReadMatch,
				pairedEnd,
				binaryInput)) {
		numMatches++;

		numAlignEntries = 0;

		if(VERBOSE >= 0 && numMatches%ALIGN_ROTATE_NUM==0) {
			fprintf(stderr, "\rthread:%d\t%d", threadID, numMatches);
		}

		/* This does not work for paired end */
		if(readMatch.maxReached == 0 && readMatch.numEntries > 0 && (maxNumMatches == 0 || readMatch.numEntries < maxNumMatches)) {

			/* Get the read length */
			matchLength = strlen(read);
			assert(matchLength+2*offsetLength < SEQUENCE_LENGTH);

			/* Allocate memory for the AlignEntries */
			if(readMatch.numEntries > 0) {
				aEntry = (AlignEntry*)malloc(sizeof(AlignEntry)*readMatch.numEntries);
				if(NULL == aEntry) {
					PrintError("RunDynamicProgrammingThread",
							"aEntry",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}
			}

			/* Run the aligner */
			assert(pairedEnd==0);
			if(readMatch.numEntries > 0) {
				numMatchesAligned++;
			}
			for(i=0;i<readMatch.numEntries;i++) { /* For each match */

				/* Get the appropriate reference read */
				RGBinaryGetSequence(rgBinary,
						readMatch.chromosomes[i], 
						readMatch.positions[i],
						FORWARD, /* We have been just reversing the read instead of the reference */
						offsetLength,
						reference,
						matchLength,
						&referenceLength,
						&position);
				assert(referenceLength > 0);

				/* If the direction was reverse, then give the
				 * reverse compliment.
				 * */
				if(readMatch.strand[i] == REVERSE) {
					GetReverseComplimentAnyCase(read, reverseRead, strlen(read));
					/* Get alignment */
					adjustPosition=AlignmentGetScore(reverseRead,
							matchLength,
							reference,
							referenceLength,
							sm,
							&aEntry[i]);
				}
				else {

					/* Get alignment */
					adjustPosition=AlignmentGetScore(read,
							matchLength,
							reference,
							referenceLength,
							sm,
							&aEntry[i]);
				}

				/* Update chromosome, position, strand and sequence name*/
				aEntry[i].chromosome = readMatch.chromosomes[i];
				aEntry[i].position = position+adjustPosition; /* Adjust position */
				aEntry[i].strand = readMatch.strand[i]; 
				strcpy(aEntry[i].readName, readName);

				/* Free memory */
				/* Free AlignEntry */
				if(VERBOSE >= DEBUG) {
					fprintf(stderr, "aligned read:%s\naligned reference:%s\nposition:%d\n",
							aEntry[i].read,
							aEntry[i].reference,
							aEntry[i].position);
				}
				if(VERBOSE >= DEBUG) {
					fprintf(stderr, "Finished entry %d out of %d.\n",
							i+1,
							readMatch.numEntries);
				}
			}
			/* Remove duplicate alignments */
			numAlignEntries=AlignEntryRemoveDuplicates(&aEntry, readMatch.numEntries, AlignEntrySortByAll);
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "Outputting %d aligns.\n",
						numAlignEntries);
			}
			/* Output alignment */
			fprintf(outputFP, "%d\n", numAlignEntries); /* We need this when merging temp files */
			for(i=0;i<numAlignEntries;i++) {
				AlignEntryPrint(&aEntry[i], outputFP);
			}
			/* Free memory */
			for(i=0;i<numAlignEntries;i++) {
				assert(aEntry[i].length>0);
				free(aEntry[i].read);
				free(aEntry[i].reference);
			}
		}
		/* Free match */
		RGMatchFree(&readMatch);
		if(pairedEnd==1) {
			RGMatchFree(&pairedReadMatch);
		}
		if(numAlignEntries > 0) {
			free(aEntry);
			aEntry = NULL;
		}
	}

	/* Free memory */
	free(reference);

	return arg;
}
