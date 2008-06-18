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
	FILE *notAlignedFP=NULL;
	FILE *matchesFP=NULL;
	FILE *matchFP=NULL;
	char **matchFileNames=NULL;
	int numMatchFileNames=0;
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	char notAlignedFileName[MAX_FILENAME_LENGTH]="\0";
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

	/* Close matches file */
	fclose(matchesFP);

	/* Create output file name */
	sprintf(outputFileName, "%sbfast.aligned.file.%s.%d.%s",
			outputDir,
			outputID,
			algorithm,
			BFAST_ALIGN_FILE_EXTENSION);
	/* Create not aligned file name */
	sprintf(notAlignedFileName, "%sbfast.not.aligned.file.%s.%d.%s",
			outputDir,
			outputID,
			algorithm,
			BFAST_NOT_ALIGNED_FILE_EXTENSION);

	/* Open output file */
	if((outputFP=fopen(outputFileName, "w"))==0) {
		PrintError("RunAligner",
				outputFileName,
				"Could not open file for writing",
				Exit,
				OpenFileError);
	}

	/* Open not aligned file */
	if((notAlignedFP=fopen(notAlignedFileName, "w"))==0) {
		PrintError("RunAligner",
				notAlignedFileName,
				"Could not open file for writing",
				Exit,
				OpenFileError);
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "Will output to %s.\n", outputFileName);
		fprintf(stderr, "Not aligned reads outputted to %s.\n", notAlignedFileName);
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
						outputFP,
						notAlignedFP);
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

	/* Close not aligned file */
	fclose(notAlignedFP);

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
		FILE *outputFP,
		FILE *notAlignedFP)
{
	/* local variables */
	ScoringMatrix sm;
	char readName[SEQUENCE_NAME_LENGTH]="\0";
	char read[SEQUENCE_LENGTH]="\0";
	char pairedRead[SEQUENCE_LENGTH]="\0";
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
		data[i].notAlignedFP = OpenTmpFile(tmpDir, &data[i].notAlignedFileName); 
	}

	/* Go through each read in the match file */
	numMatches=0;
	while(EOF!=RGMatchRead(matchFP, 
				readName, 
				read, 
				pairedRead, 
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
				pairedRead,
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

	/* Merge the not aligned tmp files */
	/* Go through each thread */
	for(i=0;i<numThreads;i++) {
		/* Move to the beginning of the file */
		fseek(data[i].notAlignedFP, 0, SEEK_SET);

		continueReading=1;
		while(continueReading==1) {
			/* Initialize the match */
			RGMatchInitialize(&readMatch);
			RGMatchInitialize(&pairedReadMatch);
			/* Read in one match */
			if(RGMatchRead(data[i].notAlignedFP,
					readName,
					read,
					pairedRead,
					&readMatch,
					&pairedReadMatch,
					pairedEnd,
					0) == EOF) {
				continueReading = 0;
			}
			else {
				RGMatchPrint(notAlignedFP,
						readName,
						read,
						pairedRead,
						&readMatch,
						&pairedReadMatch,
						pairedEnd,
						0);

				RGMatchFree(&readMatch);
				RGMatchFree(&pairedReadMatch);
			}
		}
	}

	/* Close tmp output files */
	for(i=0;i<numThreads;i++) {
		CloseTmpFile(&data[i].outputFP, &data[i].outputFileName);
		CloseTmpFile(&data[i].notAlignedFP, &data[i].notAlignedFileName);
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
	FILE *notAlignedFP = data->notAlignedFP;
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
	char pairedRead[SEQUENCE_LENGTH]="\0";
	RGMatch readMatch;
	RGMatch pairedReadMatch;
	int readLength;
	char *reference=NULL;
	int referenceLength=0;
	int adjustPosition=0;
	int i;
	int position;
	int numMatches=0;
	int numMatchesAligned=0;
	int aligned;

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
				pairedRead, 
				&readMatch,
				&pairedReadMatch,
				pairedEnd,
				binaryInput)) {
		numMatches++;

		numAlignEntries = 0;
		aligned = 0;

		if(VERBOSE >= 0 && numMatches%ALIGN_ROTATE_NUM==0) {
			fprintf(stderr, "\rthread:%d\t[%d]", threadID, numMatches);
		}

		/* This does not work for paired end */
		assert(pairedEnd == 0);
		if(readMatch.maxReached == 0 && readMatch.numEntries > 0 && (maxNumMatches == 0 || readMatch.numEntries < maxNumMatches)) {

			/* Get the read length */
			readLength = strlen(read);
			assert(readLength+2*offsetLength < SEQUENCE_LENGTH);

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
						readLength,
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
							readLength,
							reference,
							referenceLength,
							sm,
							&aEntry[i]);
				}
				else {

					/* Get alignment */
					adjustPosition=AlignmentGetScore(read,
							readLength,
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
				/* Adjust based on strand */
				if(aEntry[i].strand == REVERSE) {
					ReverseSequence(aEntry[i].read, aEntry[i].length);
					ReverseSequence(aEntry[i].reference, aEntry[i].length);
				}

			}
			/* Remove duplicate alignments */
			numAlignEntries=AlignEntryRemoveDuplicates(&aEntry, readMatch.numEntries, AlignEntrySortByAll);
			/* Update whether we aligned */
			if(numAlignEntries > 0) {
				aligned = 1;
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
		if(0 == aligned) {
			/* Print out the sequences we could not match */
			RGMatchPrint(notAlignedFP,
					readName,
					read,
					pairedRead,
					&readMatch,
					&pairedReadMatch,
					pairedEnd,
					0);
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
	if(VERBOSE >= 0) {
		fprintf(stderr, "\rthread:%d\t[%d]", threadID, numMatches);
	}

	/* Free memory */
	free(reference);

	return arg;
}
