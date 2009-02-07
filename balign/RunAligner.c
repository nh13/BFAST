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
void RunAligner(RGBinary *rg,
		char *matchFileName,
		char *scoringMatrixFileName,
		int alignmentType,
		int bestOnly,
		int space,
		int scoringType,
		int startContig,
		int startPos,
		int endContig,
		int endPos,
		int offsetLength,
		int maxNumMatches,
		int pairedEnd,
		int binaryInput,
		int numThreads,
		int usePairedEndLength,
		int pairedEndLength,
		int mirroringType,
		int forceMirroring,
		char *outputID,
		char *outputDir,
		char *tmpDir,
		int binaryOutput,
		int *totalAlignTime,
		int *totalFileHandlingTime)
{
	char *FnName = "RunAligner";
	FILE *outputFP=NULL;
	FILE *notAlignedFP=NULL;
	FILE *matchFP=NULL;
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	char notAlignedFileName[MAX_FILENAME_LENGTH]="\0";

	assert(BinaryInput == binaryInput);
	assert(BinaryOutput == binaryOutput);

	/* Check rg to make sure it is in NT Space */
	if(rg->space != NTSpace) {
		PrintError(FnName,
				"rg->space",
				"The reference genome must be in NT space",
				Exit,
				OutOfRange);
	}

	/* Adjust start and end based on reference genome */
	AdjustBounds(rg,
			&startContig,
			&startPos,
			&endContig,
			&endPos);

	/* Create output file name */
	sprintf(outputFileName, "%s%s.aligned.file.%s.%d.%d.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			space,
			pairedEnd,
			BFAST_ALIGNED_FILE_EXTENSION);
	/* Create not aligned file name */
	sprintf(notAlignedFileName, "%s%s.not.aligned.file.%s.%d.%d.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			space,
			pairedEnd,
			BFAST_MATCHES_FILE_EXTENSION);

	/* Open output file */
	if((outputFP=fopen(outputFileName, "wb"))==0) {
		PrintError(FnName,
				outputFileName,
				"Could not open file for writing",
				Exit,
				OpenFileError);
	}

	/* Open not aligned file */
	if((notAlignedFP=fopen(notAlignedFileName, "wb"))==0) {
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

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Reading match file from %s.\n",
				matchFileName);
	}

	/* Open current match file */
	if((matchFP=fopen(matchFileName, "rb"))==0) {
		PrintError(FnName,
				matchFileName,
				"Could not open file for reading",
				Exit,
				OpenFileError);
	}

	RunDynamicProgramming(matchFP,
			rg,
			scoringMatrixFileName,
			alignmentType,
			bestOnly,
			space,
			scoringType,
			startContig,
			startPos,
			endContig,
			endPos,
			offsetLength,
			maxNumMatches,
			pairedEnd,
			binaryInput,
			numThreads,
			usePairedEndLength,
			pairedEndLength,
			mirroringType,
			forceMirroring,
			tmpDir,
			outputFP,
			notAlignedFP,
			binaryOutput,
			totalAlignTime,
			totalFileHandlingTime);

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Close the match file */
	fclose(matchFP);

	/* Close output file */
	fclose(outputFP);

	/* Close not aligned file */
	fclose(notAlignedFP);
}

/* TODO */
void RunDynamicProgramming(FILE *matchFP,
		RGBinary *rg,
		char *scoringMatrixFileName,
		int alignmentType,
		int bestOnly,
		int space,
		int scoringType,
		int startContig,
		int startPos,
		int endContig,
		int endPos,
		int offsetLength,
		int maxNumMatches,
		int pairedEnd,
		int binaryInput,
		int numThreads,
		int usePairedEndLength,
		int pairedEndLength,
		int mirroringType,
		int forceMirroring,
		char *tmpDir,
		FILE *outputFP,
		FILE *notAlignedFP,
		int binaryOutput,
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
	int64_t numLocalAlignments=0;
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
	ScoringMatrixRead(scoringMatrixFileName, &sm, space); 

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
				startContig,
				startPos,
				endContig,
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
					binaryOutput);
		}
		else {
			/* Print match to temp file */
			RGMatchesPrint(data[i].inputFP,
					&m,
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
		data[i].rg=rg;
		data[i].space=space;
		data[i].scoringType=scoringType;
		data[i].offsetLength=offsetLength;
		data[i].pairedEnd=pairedEnd;
		data[i].usePairedEndLength = usePairedEndLength;
		data[i].pairedEndLength = pairedEndLength;
		data[i].mirroringType = mirroringType;
		data[i].forceMirroring = forceMirroring;
		data[i].binaryInput=binaryInput;
		data[i].binaryOutput=binaryOutput;
		data[i].sm = &sm;
		data[i].alignmentType = alignmentType;
		data[i].bestOnly = bestOnly;
		data[i].numLocalAlignments = 0;
		data[i].threadID = i;
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Processing %d reads.\n",
				numMatches - numNotAligned
			   );
		fprintf(stderr, "%s", BREAK_LINE);
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
			PrintError("RunDynamicProgrammingThread",
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
		numLocalAlignments += data[i].numLocalAlignments;
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
			if(EOF != AlignEntriesRead(&aEntries, data[i].outputFP, pairedEnd, space, binaryOutput)) {
				if(VERBOSE >=0 && numAligned%ALIGN_ROTATE_NUM == 0) {
					fprintf(stderr, "\r[%d]", numAligned);
				}
				continueReading=1;
				/* Update the number that were aligned */
				assert(aEntries.numEntriesOne > 0 ||
						(aEntries.pairedEnd == PairedEnd && aEntries.numEntriesTwo > 0));
				numAligned++;
				/* Print it out */
				AlignEntriesPrint(&aEntries,
						outputFP,
						binaryOutput);
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
						binaryInput) == EOF) {
				continueReading = 0;
			}
			else {
				if(VERBOSE >= 0 && numNotAligned%ALIGN_ROTATE_NUM==0) {
					fprintf(stderr, "\r[%d]", numNotAligned);
				}
				numNotAligned++;
				RGMatchesPrint(notAlignedFP,
						&m,
						binaryOutput);

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
		fprintf(stderr, "Performed %lld local alignments.\n", (long long int)numLocalAlignments);
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
	FILE *notAlignedFP = data->notAlignedFP;
	RGBinary *rg=data->rg;
	int space=data->space;
	int scoringType=data->scoringType;
	int offsetLength=data->offsetLength;
	int pairedEnd=data->pairedEnd;
	int usePairedEndLength=data->usePairedEndLength;
	int pairedEndLength=data->pairedEndLength;
	int mirroringType=data->mirroringType;
	int forceMirroring=data->forceMirroring;
	int binaryInput=data->binaryInput;
	int binaryOutput=data->binaryOutput;
	ScoringMatrix *sm = data->sm;
	int alignmentType=data->alignmentType;
	int bestOnly=data->bestOnly;
	int threadID=data->threadID;
	/* Local variables */
	/*
	   char *FnName = "RunDynamicProgrammingThread";
	   */
	AlignEntries aEntries;
	int numAlignEntries=0;
	RGMatches m;
	int numMatches=0;
	int ctrOne=0;
	int ctrTwo=0;

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
		ctrOne=ctrTwo=0;

		if(VERBOSE >= 0 && numMatches%ALIGN_ROTATE_NUM==0) {
			fprintf(stderr, "\rthread:%d\t[%d]", threadID, numMatches);
		}

		/* Update the number of local alignments performed */
		data->numLocalAlignments += AlignRGMatches(&m,
				rg,
				&aEntries,
				space,
				pairedEnd,
				scoringType,
				offsetLength,
				sm,
				alignmentType,
				bestOnly,
				usePairedEndLength,
				pairedEndLength,
				mirroringType,
				forceMirroring);

		if(0 < aEntries.numEntriesOne ||
				(pairedEnd == 1 && 0 < aEntries.numEntriesTwo)) {
			/* Remove duplicates */
			AlignEntriesRemoveDuplicates(&aEntries,
					AlignEntrySortByAll);
		}

		/* Output alignment */
		if(0 < aEntries.numEntriesOne || 
				(PairedEnd == aEntries.pairedEnd && 0 < aEntries.numEntriesTwo)) {
			AlignEntriesPrint(&aEntries,
					outputFP,
					binaryOutput);
		}
		else {
			RGMatchesPrint(notAlignedFP,
					&m,
					binaryOutput);
		}

		/* Free memory */
		AlignEntriesFree(&aEntries);
		RGMatchesFree(&m);
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\rthread:%d\t[%d]", threadID, numMatches);
	}

	return arg;
}
