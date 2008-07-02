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
#include "../blib/AlignEntries.h"
#include "../blib/AlignEntry.h" 
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
		char *tmpDir,
		int *totalAlignTime,
		int *totalFileHandlingTime)
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
	int startTime, endTime;

	/* Start file handling timer */
	startTime = time(NULL);

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
		matchFileNames[numMatchFileNames-1] = malloc(sizeof(char)*(strlen(tempFileName)+1));
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

	/* End file handling timer */
	endTime = time(NULL);
	(*totalFileHandlingTime) += endTime - startTime;

	/* Create output file name */
	sprintf(outputFileName, "%s%s.aligned.file.%s.%d.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			algorithm,
			BFAST_ALIGN_FILE_EXTENSION);
	/* Create not aligned file name */
	sprintf(notAlignedFileName, "%s%s.not.aligned.file.%s.%d.%s",
			outputDir,
			PROGRAM_NAME,
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
		fprintf(stderr, "Will output aligned reads to %s.\n", outputFileName);
		fprintf(stderr, "Will output unaligned reads to %s.\n", notAlignedFileName);
		fprintf(stderr, "%s", BREAK_LINE);
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
						notAlignedFP,
						totalAlignTime,
						totalFileHandlingTime);
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
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
		}
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
		FILE *notAlignedFP,
		int *totalAlignTime,
		int *totalFileHandlingTime)
{
	/* local variables */
	ScoringMatrix sm;
	char readName[SEQUENCE_NAME_LENGTH]="\0";
	char read[SEQUENCE_LENGTH]="\0";
	char pairedRead[SEQUENCE_LENGTH]="\0";
	RGMatch readMatch;
	RGMatch pairedReadMatch;
	int i;
	int continueReading=0;
	int numMatches=0;
	int numAligned=0;
	int numNotAligned=0;
	int startTime, endTime;
	AlignEntries aEntries;
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

	/* Start file handling timer */
	startTime = time(NULL);

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

	/* Go through each read in the match file and partition them for the threads */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Filtering and partitioning matches for threads...\n0");
	}
	i=0;
	while(EOF!=RGMatchRead(matchFP, 
				readName, 
				read, 
				pairedRead, 
				&readMatch,
				&pairedReadMatch,
				pairedEnd,
				binaryInput) 
		 ) {

		if(VERBOSE >= 0 && numMatches%PARTITION_MATCHES_ROTATE_NUM==0) {
			fprintf(stderr, "\r[%d]", numMatches);
		}
		/* increment */
		numMatches++;

		/* Filter those reads we will not be able to align */
		/* line 1 - if both were found to have too many alignments in bmatches */
		/* line 2 - if both have too many alignments in balign */
		if( (readMatch.maxReached == 1 && (pairedEnd == 0 || pairedReadMatch.maxReached == 1)) ||
				(maxNumMatches != 0 && readMatch.numEntries > maxNumMatches && (pairedEnd == 0 || pairedReadMatch.numEntries > maxNumMatches))) {
			numNotAligned++;
			RGMatchPrint(notAlignedFP,
					readName,
					read,
					pairedRead,
					&readMatch,
					&pairedReadMatch,
					pairedEnd,
					0); /* Do not print in binary */
		}
		else {
			/* Print match to temp file */
			RGMatchPrint(data[i].inputFP,
					readName,
					read,
					pairedRead,
					&readMatch,
					&pairedReadMatch,
					pairedEnd,
					binaryInput);
			/* Increment */
			i = (i+1)%numThreads;
		}

		/* Free match */
		RGMatchFree(&readMatch);
		if(pairedEnd==1) {
			RGMatchFree(&pairedReadMatch);
		}
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\r[%d]\n", numMatches);
		fprintf(stderr, "Initially filtered %d out of %d reads.\n",
				numNotAligned,
				numMatches);
	}


	/* End file handling timer */
	endTime = time(NULL);
	(*totalFileHandlingTime) += endTime - startTime;

	/* Create thread arguments */
	for(i=0;i<numThreads;i++) {
		fseek(data[i].inputFP, 0, SEEK_SET);
		fseek(data[i].outputFP, 0, SEEK_SET);
		data[i].rgBinary=rgBinary;
		data[i].offsetLength=offsetLength;
		data[i].pairedEnd=pairedEnd;
		data[i].binaryInput=binaryInput;
		data[i].sm = &sm;
		data[i].threadID = i;
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "Performing alignment...\n");
		fprintf(stderr, "Currently on:\n0");
	}

	/* Start align timer */
	startTime = time(NULL);

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

	/* End align timer */
	endTime = time(NULL);
	(*totalAlignTime) += endTime - startTime;

	/* Start file handling timer */
	startTime = time(NULL);

	/* Close tmp input files */
	for(i=0;i<numThreads;i++) {
		CloseTmpFile(&data[i].inputFP, &data[i].inputFileName);
	}

	/* Merge all the aligned reads from the threads */
	if(VERBOSE >=0) {
		fprintf(stderr, "Merging and outputting aligned reads from threads...\n[0]");
	}
	AlignEntriesInitialize(&aEntries);
	numAligned=0;
	continueReading=1;
	while(continueReading==1) {
		/* Get an align from a thread */
		continueReading=0;
		for(i=0;i<numThreads;i++) {
			/* Read in the align entries */
			if(EOF != AlignEntriesRead(&aEntries, data[i].outputFP)) {
				if(VERBOSE >=0 && numAligned%ALIGN_ROTATE_NUM == 0) {
					fprintf(stderr, "\r[%d]", numAligned);
				}
				continueReading=1;
				/* Update the number that were aligned */
				if(aEntries.numEntriesOne > 0 && (aEntries.pairedEnd == 0 || aEntries.numEntriesTwo > 0)) {
					numAligned++;
				}
				/* Print it out */
				AlignEntriesPrint(&aEntries,
						outputFP);
			}
			AlignEntriesFree(&aEntries);
		}
	}
	if(VERBOSE >=0) {
		fprintf(stderr, "\r[%d]\n", numAligned);
	}

	/* Merge the not aligned tmp files */
	if(VERBOSE >=0) {
		fprintf(stderr, "Merging and outputting unaligned reads from threads and initial filter...\n");
	}
	for(i=0;i<numThreads;i++) {
		fseek(data[i].notAlignedFP, 0, SEEK_SET);

		continueReading=1;
		while(continueReading==1) {
			RGMatchInitialize(&readMatch);
			RGMatchInitialize(&pairedReadMatch);
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
				if(VERBOSE >= 0 && numNotAligned%ALIGN_ROTATE_NUM==0) {
					fprintf(stderr, "\r[%d]", numNotAligned);
				}
				numNotAligned++;
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
	if(VERBOSE >=0) {
		fprintf(stderr, "\r[%d]\n", numNotAligned);
	}

	/* Close tmp output files */
	for(i=0;i<numThreads;i++) {
		CloseTmpFile(&data[i].outputFP, &data[i].outputFileName);
		CloseTmpFile(&data[i].notAlignedFP, &data[i].notAlignedFileName);
	}

	/* Start file handling timer */
	endTime = time(NULL);
	(*totalFileHandlingTime) += endTime - startTime;

	if(VERBOSE >=0) {
		fprintf(stderr, "Outputted alignments for %d reads.\n", numAligned);
		fprintf(stderr, "Outputted %d reads for which there were no alignments.\n", numNotAligned); 
		fprintf(stderr, "Outputting complete.\n");
	}
	assert(numAligned + numNotAligned == numMatches);

	/* Free memory */
	free(data);
	free(threads);
	AlignEntriesFree(&aEntries);
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
	/*
	FILE *notAlignedFP = data->notAlignedFP;
	*/
	RGBinary *rgBinary=data->rgBinary;
	int offsetLength=data->offsetLength;
	int pairedEnd=data->pairedEnd;
	int binaryInput=data->binaryInput;
	int threadID=data->threadID;
	ScoringMatrix *sm = data->sm;
	/* Local variables */
	/*
	   char *FnName = "RunDynamicProgrammingThread";
	   */
	AlignEntries aEntries;
	int numAlignEntries=0;
	char readName[SEQUENCE_NAME_LENGTH]="\0";
	char read[SEQUENCE_LENGTH]="\0";
	char pairedRead[SEQUENCE_LENGTH]="\0";
	RGMatch readMatch;
	RGMatch pairedReadMatch;
	int readLength, pairedReadLength=0;
	int i;
	int numMatches=0;

	/* Initialize */
	RGMatchInitialize(&readMatch);
	RGMatchInitialize(&pairedReadMatch);
	AlignEntriesInitialize(&aEntries);

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

		if(VERBOSE >= 0 && numMatches%ALIGN_ROTATE_NUM==0) {
			fprintf(stderr, "\rthread:%d\t[%d]", threadID, numMatches);
		}

		/* Get the read length */
		readLength = strlen(read);
		assert(readLength+2*offsetLength < SEQUENCE_LENGTH);
		if(pairedEnd == 1) {
			pairedReadLength = strlen(pairedRead);
			assert(pairedReadLength+2*offsetLength < SEQUENCE_LENGTH);
		}

		/* Allocate memory for the AlignEntries */
		AlignEntriesAllocate(&aEntries,
				readMatch.numEntries,
				pairedReadMatch.numEntries,
				pairedEnd);
		/* Copy over read name */
		strcpy(aEntries.readName, readName);

		/* Run the aligner */
		/* First entry */
		assert(readMatch.numEntries == aEntries.numEntriesOne);
		for(i=0;i<readMatch.numEntries;i++) { /* For each match */
			RunDynamicProgrammingThreadHelper(rgBinary,
					readMatch.chromosomes[i],
					readMatch.positions[i],
					readMatch.strand[i],
					read,
					readLength,
					offsetLength,
					sm,
					&aEntries.entriesOne[i]);
		}
		/* Second entry */
		assert(pairedReadMatch.numEntries == aEntries.numEntriesTwo);
		for(i=0;i<pairedReadMatch.numEntries;i++) { /* For each match */
			RunDynamicProgrammingThreadHelper(rgBinary,
					pairedReadMatch.chromosomes[i],
					pairedReadMatch.positions[i],
					pairedReadMatch.strand[i],
					pairedRead,
					pairedReadLength,
					offsetLength,
					sm,
					&aEntries.entriesTwo[i]);
		}

		if(aEntries.numEntriesOne > 0 && aEntries.numEntriesTwo > 0) {
			/* Remove duplicates */
			AlignEntriesRemoveDuplicates(&aEntries,
					AlignEntrySortByAll);
			assert(aEntries.numEntriesOne > 0 && aEntries.numEntriesTwo > 0);
		}

		/* Output alignment */
		AlignEntriesPrint(&aEntries,
				outputFP);

		/* Free alignment entry */
		AlignEntriesFree(&aEntries);
		/* Free match */
		RGMatchFree(&readMatch);
		if(pairedEnd==1) {
			RGMatchFree(&pairedReadMatch);
		}

	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\rthread:%d\t[%d]", threadID, numMatches);
	}


	return arg;
}

/* TODO */
void RunDynamicProgrammingThreadHelper(RGBinary *rgBinary,
		uint8_t chromosome,
		uint32_t position,
		int8_t strand,
		char *read,
		int readLength,
		int offsetLength,
		ScoringMatrix *sm,
		AlignEntry *aEntry)
{
	char *FnName = "RunDynamicProgrammingThreadHelper";
	char reverseRead[SEQUENCE_LENGTH]="\0";
	char *reference=NULL;
	char tmpString[SEQUENCE_LENGTH]="\0";
	int referenceLength=0;
	int referencePosition=0;
	int adjustPosition=0;

	/* Allocate memory for the reference */
	referenceLength = 2*offsetLength + SEQUENCE_LENGTH + 1;
	reference = malloc(sizeof(char)*(referenceLength+1));
	if(NULL==reference) {
		PrintError(FnName,
				"reference",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	reference[referenceLength] = '\0'; /* Add null terminator */

	/* Get the appropriate reference read */
	RGBinaryGetSequence(rgBinary,
			chromosome,
			position,
			FORWARD, /* We have been just reversing the read instead of the reference */
			offsetLength,
			reference,
			readLength,
			&referenceLength,
			&referencePosition);
	assert(referenceLength > 0);

	/* If the direction was reverse, then give the
	 * reverse compliment for the read.
	 * */
	switch(strand) {
		case FORWARD:
			/* Get alignment */
			adjustPosition=AlignmentGetScore(read,
					readLength,
					reference,
					referenceLength,
					sm,
					aEntry);
			break;
		case REVERSE:
			GetReverseComplimentAnyCase(read, reverseRead, readLength);
			/* Get alignment */
			adjustPosition=AlignmentGetScore(reverseRead,
					readLength,
					reference,
					referenceLength,
					sm,
					aEntry);
			/* We must reverse the alignment to match the REVERSE stand */
			GetReverseComplimentAnyCase(aEntry->read, tmpString, aEntry->length);
			strcpy(aEntry->read, tmpString);
			GetReverseComplimentAnyCase(aEntry->reference, tmpString, aEntry->length);
			strcpy(aEntry->reference, tmpString);
			break;
		default:
			PrintError(FnName,
					NULL,
					"Could not understand strand",
					Exit,
					OutOfRange);
			break;
	}
	/* Update adjustPosition based on offsetLength */
	assert(adjustPosition >= 0 && adjustPosition <= referenceLength);

	/* Update chromosome, position, strand and sequence name*/
	aEntry->chromosome = chromosome;
	aEntry->position = referencePosition+adjustPosition; /* Adjust position */
	aEntry->strand = strand;

	/* Check align entry */
	/*
	   AlignEntryCheckReference(&aEntry[i], rgBinary);
	   */
	/* Free memory */
	free(reference);
	reference=NULL;
}
