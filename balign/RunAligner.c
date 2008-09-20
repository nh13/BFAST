#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "../blib/RGBinary.h"
#include "../blib/RGMatch.h" 
#include "../blib/RGMatches.h" 
#include "../blib/AlignEntries.h"
#include "../blib/AlignEntry.h" 
#include "Align.h"
#include "ScoringMatrix.h"
#include "Definitions.h"
#include "RunAligner.h"

/* TODO */
void RunAligner(RGBinary *rgBinary,
		char *matchesFileName,
		char *scoringMatrixFileName,
		int colorSpace,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int offsetLength,
		int maxNumMatches,
		int pairedEnd,
		int binaryInput,
		int numThreads,
		int usePairedEndLength,
		int pairedEndLength,
		char *outputID,
		char *outputDir,
		char *tmpDir,
		int *totalAlignTime,
		int *totalFileHandlingTime)
{
	char *FnName = "RunAligner";
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

	/* Check rg to make sure it is not in color space */
	if(rgBinary->colorSpace == 1) {
		PrintError(FnName,
				"rg->colorSpace",
				"The reference genome must not be in color space",
				Exit,
				OutOfRange);
	}

	/* Adjust start and end based on reference genome */
	AdjustBounds(rgBinary,
			&startChr,
			&startPos,
			&endChr,
			&endPos);

	/* Open matches file */
	if((matchesFP=fopen(matchesFileName, "r"))==0) {
		PrintError(FnName,
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
			PrintError(FnName,
					"matchFileNames",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		/* Allocate memory */
		matchFileNames[numMatchFileNames-1] = malloc(sizeof(char)*(strlen(tempFileName)+1));
		if(NULL==matchFileNames[numMatchFileNames-1]) {
			PrintError(FnName,
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
			colorSpace,
			BFAST_ALIGN_FILE_EXTENSION);
	/* Create not aligned file name */
	sprintf(notAlignedFileName, "%s%s.not.aligned.file.%s.%d.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			colorSpace,
			BFAST_NOT_ALIGNED_FILE_EXTENSION);

	/* Open output file */
	if((outputFP=fopen(outputFileName, "w"))==0) {
		PrintError(FnName,
				outputFileName,
				"Could not open file for writing",
				Exit,
				OpenFileError);
	}

	/* Open not aligned file */
	if((notAlignedFP=fopen(notAlignedFileName, "w"))==0) {
		PrintError(FnName,
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
			PrintError(FnName,
					matchFileNames[i],
					"Could not open matchesFileNames[] for reading",
					Exit,
					OpenFileError);
		}

		RunDynamicProgramming(matchFP,
				rgBinary,
				scoringMatrixFileName,
				colorSpace,
				startChr,
				startPos,
				endChr,
				endPos,
				offsetLength,
				maxNumMatches,
				pairedEnd,
				binaryInput,
				numThreads,
				usePairedEndLength,
				pairedEndLength,
				tmpDir,
				outputFP,
				notAlignedFP,
				totalAlignTime,
				totalFileHandlingTime);

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
		int colorSpace,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int offsetLength,
		int maxNumMatches,
		int pairedEnd,
		int binaryInput,
		int numThreads,
		int usePairedEndLength,
		int pairedEndLength,
		char *tmpDir,
		FILE *outputFP,
		FILE *notAlignedFP,
		int *totalAlignTime,
		int *totalFileHandlingTime)
{
	/* local variables */
	ScoringMatrix sm;
	RGMatches m;
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

	/* Initialize */
	RGMatchesInitialize(&m);
	ScoringMatrixInitialize(&sm);

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
	ScoringMatrixRead(scoringMatrixFileName, &sm, colorSpace); 

	/**/
	/* Split the input file into equal temp files for each thread */
	/* Could put this in a separate function */
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
	while(EOF!=RGMatchesRead(matchFP, 
				&m,
				pairedEnd,
				binaryInput) 
		 ) {

		if(VERBOSE >= 0 && numMatches%PARTITION_MATCHES_ROTATE_NUM==0) {
			fprintf(stderr, "\r[%d]", numMatches);
		}
		/* increment */
		numMatches++;

		/* Filter the matches that are not within the bounds */
		/* Filter one if it has too many entries */
		RGMatchesFilterOutOfRange(&m,
				startChr,
				startPos,
				endChr,
				endPos,
				maxNumMatches);

		/* Filter those reads we will not be able to align */
		/* line 1 - if both were found to have too many alignments in bmatches */
		/* line 3 - if both do not have any possible matches */
		if( (m.matchOne.maxReached == 1 && (pairedEnd == 0 || m.matchTwo.maxReached == 1)) ||
				(m.matchOne.numEntries <= 0 && (pairedEnd == 0 || m.matchTwo.numEntries <= 0))
		  )	{
			/* Print to the not aligned file */
			numNotAligned++;
			RGMatchesPrint(notAlignedFP,
					&m,
					pairedEnd,
					0); /* Do not print in binary */
		}
		else {
			/* Print match to temp file */
			RGMatchesPrint(data[i].inputFP,
					&m,
					pairedEnd,
					binaryInput);
			/* Increment */
			i = (i+1)%numThreads;
		}

		/* Free memory */
		RGMatchesFree(&m);
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
		data[i].colorSpace=colorSpace;
		data[i].offsetLength=offsetLength;
		data[i].pairedEnd=pairedEnd;
		data[i].usePairedEndLength = usePairedEndLength;
		data[i].pairedEndLength = pairedEndLength;
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
			if(EOF != AlignEntriesRead(&aEntries, data[i].outputFP, pairedEnd)) {
				if(VERBOSE >=0 && numAligned%ALIGN_ROTATE_NUM == 0) {
					fprintf(stderr, "\r[%d]", numAligned);
				}
				continueReading=1;
				/* Update the number that were aligned */
				assert(aEntries.numEntriesOne > 0 ||
						(aEntries.pairedEnd == 1 && aEntries.numEntriesTwo > 0));
				numAligned++;
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
			RGMatchesInitialize(&m);
			if(RGMatchesRead(data[i].notAlignedFP,
						&m,
						pairedEnd,
						0) == EOF) {
				continueReading = 0;
			}
			else {
				if(VERBOSE >= 0 && numNotAligned%ALIGN_ROTATE_NUM==0) {
					fprintf(stderr, "\r[%d]", numNotAligned);
				}
				numNotAligned++;
				RGMatchesPrint(notAlignedFP,
						&m,
						pairedEnd,
						0);

				RGMatchesFree(&m);
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
	ScoringMatrixFree(&sm);
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
	int colorSpace=data->colorSpace;
	int offsetLength=data->offsetLength;
	int pairedEnd=data->pairedEnd;
	int usePairedEndLength=data->usePairedEndLength;
	int pairedEndLength=data->pairedEndLength;
	int binaryInput=data->binaryInput;
	int threadID=data->threadID;
	ScoringMatrix *sm = data->sm;
	/* Local variables */
	/*
	   char *FnName = "RunDynamicProgrammingThread";
	   */
	char matchRead[SEQUENCE_LENGTH]="\0";
	int matchReadLength=0;
	AlignEntries aEntries;
	int numAlignEntries=0;
	RGMatches m;
	int i;
	int numMatches=0;

	/* Initialize */
	RGMatchesInitialize(&m);
	AlignEntriesInitialize(&aEntries);

	/* Go through each read in the match file */
	while(EOF!=RGMatchesRead(inputFP, 
				&m,
				pairedEnd,
				binaryInput)) {
		numMatches++;
		numAlignEntries = 0;

		if(VERBOSE >= 0 && numMatches%ALIGN_ROTATE_NUM==0) {
			fprintf(stderr, "\rthread:%d\t[%d]", threadID, numMatches);
		}

		/* Check to see if we should try to align one read with no candidate
		 * locations if the other one has candidate locations.
		 *
		 * This assumes that the the first read is 5'->3' before the second read.
		 * */
		if(pairedEnd == 1 && usePairedEndLength == 1) {
			RGMatchesMirrorPairedEnd(&m, pairedEndLength);
		}

		/* Allocate memory for the AlignEntries */
		AlignEntriesAllocate(&aEntries,
				m.matchOne.numEntries,
				m.matchTwo.numEntries,
				pairedEnd);
		/* Copy over read name */
		strcpy(aEntries.readName, (char*)m.readName);

		/* Run the aligner */
		/* First entry */
		if(m.matchOne.numEntries > 0) {
			strcpy(matchRead, (char*)m.matchOne.read);
			matchReadLength = m.matchOne.readLength;
			/* HERE A1 */
			/*
			   fprintf(stderr, "HERE A1\nmatchRead=%s\nmatchReadLength=%d\n",
			   matchRead,
			   matchReadLength);
			   */
			if(colorSpace == 1) {
				matchReadLength = ConvertReadFromColorSpace(matchRead, matchReadLength);
			}
		}
		assert(m.matchOne.numEntries == aEntries.numEntriesOne);
		for(i=0;i<m.matchOne.numEntries;i++) { /* For each match */
			RunDynamicProgrammingThreadHelper(rgBinary,
					m.matchOne.chromosomes[i],
					m.matchOne.positions[i],
					m.matchOne.strand[i],
					matchRead,
					matchReadLength,
					colorSpace,
					offsetLength,
					sm,
					&aEntries.entriesOne[i]);
		}
		/* Second entry */
		if(m.matchTwo.numEntries > 0) {
			strcpy(matchRead, (char*)m.matchTwo.read);
			matchReadLength = m.matchTwo.readLength;
			if(colorSpace == 1) {
				matchReadLength = ConvertReadFromColorSpace(matchRead, matchReadLength);
			}
		}
		assert(m.matchTwo.numEntries == aEntries.numEntriesTwo);
		for(i=0;i<m.matchTwo.numEntries;i++) { /* For each match */
			RunDynamicProgrammingThreadHelper(rgBinary,
					m.matchTwo.chromosomes[i],
					m.matchTwo.positions[i],
					m.matchTwo.strand[i],
					matchRead,
					matchReadLength,
					colorSpace,
					offsetLength,
					sm,
					&aEntries.entriesTwo[i]);
		}

		if(aEntries.numEntriesOne > 0 ||
				(pairedEnd == 1 && aEntries.numEntriesTwo > 0)) {

			/* Remove duplicates */
			AlignEntriesRemoveDuplicates(&aEntries,
					AlignEntrySortByAll);
		}

		assert(pairedEnd == aEntries.pairedEnd);

		/* Output alignment */
		AlignEntriesPrint(&aEntries,
				outputFP);

		/* Free memory */
		AlignEntriesFree(&aEntries);
		RGMatchesFree(&m);
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
		int colorSpace,
		int offsetLength,
		ScoringMatrix *sm,
		AlignEntry *aEntry)
{
	/*
	   char *FnName = "RunDynamicProgrammingThreadHelper";
	   */
	char *reference=NULL;
	int referenceLength=0;
	int referencePosition=0;
	int adjustPosition=0;

	/* Get the appropriate reference read */
	RGBinaryGetSequence(rgBinary,
			chromosome,
			position,
			FORWARD, /* We have been just reversing the read instead of the reference */
			offsetLength,
			&reference,
			readLength,
			&referenceLength,
			&referencePosition);
	assert(referenceLength > 0);

	/* HERE 41 */
	/*
	   fprintf(stderr, "HERE 41\nread=%s\nreference=%s\n",
	   read,
	   reference);
	   */

	/* Get alignment */
	adjustPosition=Align(read,
			readLength,
			reference,
			referenceLength,
			sm,
			aEntry,
			strand,
			colorSpace);

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
