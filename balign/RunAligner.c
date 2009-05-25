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
#include "../blib/AlignedRead.h"
#include "../blib/AlignedEnd.h" 
#include "../blib/AlignedEntry.h" 
#include "../blib/ScoringMatrix.h"
#include "Align.h"
#include "Definitions.h"
#include "RunAligner.h"

/* TODO */
void RunAligner(RGBinary *rg,
		char *matchFileName,
		char *scoringMatrixFileName,
		int alignmentType,
		int bestOnly,
		int space,
		int startContig,
		int startPos,
		int endContig,
		int endPos,
		int offsetLength,
		int maxNumMatches,
		int avgMismatchQuality,
		int numThreads,
		int usePairedEndLength,
		int pairedEndLength,
		int mirroringType,
		int forceMirroring,
		char *outputID,
		char *outputDir,
		char *tmpDir,
		int *totalAlignedTime,
		int *totalFileHandlingTime)
{
	char *FnName = "RunAligner";
	gzFile outputFP=NULL;
	gzFile notAlignedFP=NULL;
	gzFile matchFP=NULL;
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	char notAlignedFileName[MAX_FILENAME_LENGTH]="\0";

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
	sprintf(outputFileName, "%s%s.aligned.file.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			BFAST_ALIGNED_FILE_EXTENSION);
	/* Create not aligned file name */
	sprintf(notAlignedFileName, "%s%s.not.aligned.file.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			BFAST_MATCHES_FILE_EXTENSION);

	/* Open output file */
	if((outputFP=gzopen(outputFileName, "wb"))==0) {
		PrintError(FnName,
				outputFileName,
				"Could not open file for writing",
				Exit,
				OpenFileError);
	}

	/* Open not aligned file */
	if((notAlignedFP=gzopen(notAlignedFileName, "wb"))==0) {
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
	if((matchFP=gzopen(matchFileName, "rb"))==0) {
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
			startContig,
			startPos,
			endContig,
			endPos,
			offsetLength,
			maxNumMatches,
			avgMismatchQuality,
			numThreads,
			usePairedEndLength,
			pairedEndLength,
			mirroringType,
			forceMirroring,
			tmpDir,
			outputFP,
			notAlignedFP,
			totalAlignedTime,
			totalFileHandlingTime);

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Close the match file */
	gzclose(matchFP);

	/* Close output file */
	gzclose(outputFP);

	/* Close not aligned file */
	gzclose(notAlignedFP);
}

/* TODO */
void RunDynamicProgramming(gzFile matchFP,
		RGBinary *rg,
		char *scoringMatrixFileName,
		int alignmentType,
		int bestOnly,
		int space,
		int startContig,
		int startPos,
		int endContig,
		int endPos,
		int offsetLength,
		int maxNumMatches,
		int avgMismatchQuality,
		int numThreads,
		int usePairedEndLength,
		int pairedEndLength,
		int mirroringType,
		int forceMirroring,
		char *tmpDir,
		gzFile outputFP,
		gzFile notAlignedFP,
		int *totalAlignedTime,
		int *totalFileHandlingTime)
{
	char *FnName="RunDynamicProgramming";
	/* local variables */
	ScoringMatrix sm;
	double mismatchScore;
	RGMatches m;
	int32_t i, j;

	int32_t toAlign, wasAligned;
	int continueReading=0;
	int numMatches=0;
	int numAligned=0;
	int numNotAligned=0;
	int startTime, endTime;
	int64_t numLocalAlignments=0;
	AlignedRead aEntries;
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
		PrintError(FnName,
				"data",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for thread ids */
	threads = malloc(sizeof(pthread_t)*numThreads);
	if(NULL==threads) {
		PrintError(FnName,
				"threads",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Start file handling timer */
	startTime = time(NULL);

	/* Read in scoring matrix */
	ScoringMatrixRead(scoringMatrixFileName, &sm, space); 
	/* Calculate mismatch score */
	/* Assumes all match scores are the same and all substitution scores are the same */
	if(space == NTSpace) {
		mismatchScore = sm.ntMatch - sm.ntMismatch;
	}
	else {
		mismatchScore = sm.colorMatch - sm.colorMismatch;
	}

	/**/
	/* Split the input file into equal temp files for each thread */
	/* Could put this in a separate function */
	/**/

	/* Open temp files for the threads */
	for(i=0;i<numThreads;i++) {
		data[i].inputFP = OpenTmpGZFile(tmpDir, &data[i].inputFileName);
		data[i].outputFP = OpenTmpGZFile(tmpDir, &data[i].outputFileName);
		data[i].notAlignedFP = OpenTmpGZFile(tmpDir, &data[i].notAlignedFileName); 
	}

	/* Go through each read in the match file and partition them for the threads */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Filtering and partitioning matches for threads...\n0");
	}
	i=0;
	while(EOF!=RGMatchesRead(matchFP, 
				&m)) {

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
		/* line 1 - if an end reached its maximum match limit or there are no
		 * candidate alignment locations */
		for(j=toAlign=0;0==toAlign && j<m.numEnds;j++) {
			if(1 != m.ends[j].maxReached && 0 < m.ends[j].numEntries) {
				toAlign = 1;
			}
		}
		if(0 == toAlign) {
			/* Print to the not aligned file */
			numNotAligned++;
			RGMatchesPrint(notAlignedFP,
					&m);
		}
		else {
			/* Print match to temp file */
			RGMatchesPrint(data[i].inputFP,
					&m);
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
		/* Reinitialize file pointer */
		CloseTmpGZFile(&data[i].inputFP, 
				&data[i].outputFileName,
				0);
		if(!(data[i].inputFP=gzopen(data[i].inputFileName, "rb"))) {
			PrintError(FnName,
					data[i].inputFileName,
					"Could not re-open file for reading",
					Exit,
					OpenFileError);
		}
		data[i].rg=rg;
		data[i].space=space;
		data[i].offsetLength=offsetLength;
		data[i].usePairedEndLength = usePairedEndLength;
		data[i].pairedEndLength = pairedEndLength;
		data[i].mirroringType = mirroringType;
		data[i].forceMirroring = forceMirroring;
		data[i].sm = &sm;
		data[i].alignmentType = alignmentType;
		data[i].bestOnly = bestOnly;
		data[i].numLocalAlignments = 0;
		data[i].avgMismatchQuality = avgMismatchQuality;
		data[i].mismatchScore = mismatchScore;
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
			PrintError(FnName,
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
			PrintError(FnName,
					"pthread_join: errCode",
					"Thread returned an error",
					Exit,
					ThreadError);
		}
		/* Reinitialize file pointer */
		CloseTmpGZFile(&data[i].outputFP, 
				&data[i].outputFileName,
				0);
		if(!(data[i].outputFP=gzopen(data[i].outputFileName, "rb"))) {
			PrintError(FnName,
					data[i].outputFileName,
					"Could not re-open file for reading",
					Exit,
					OpenFileError);
		}
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Alignment complete.\n");
	}

	/* End align timer */
	endTime = time(NULL);
	(*totalAlignedTime) += endTime - startTime;

	/* Start file handling timer */
	startTime = time(NULL);

	/* Close tmp input files */
	for(i=0;i<numThreads;i++) {
		CloseTmpGZFile(&data[i].inputFP, &data[i].inputFileName, 1);
		numLocalAlignments += data[i].numLocalAlignments;
	}

	/* Merge all the aligned reads from the threads */
	if(VERBOSE >=0) {
		fprintf(stderr, "Merging and outputting aligned reads from threads...\n[0]");
	}
	AlignedReadInitialize(&aEntries);
	numAligned=0;
	continueReading=1;
	while(continueReading==1) {
		/* Get an align from a thread */
		continueReading=0;
		for(i=0;i<numThreads;i++) {
			/* Read in the align entries */
			if(EOF != AlignedReadRead(&aEntries, data[i].outputFP)) {
				if(VERBOSE >=0 && numAligned%ALIGN_ROTATE_NUM == 0) {
					fprintf(stderr, "\r[%d]", numAligned);
				}
				continueReading=1;
				/* Update the number that were aligned */
				for(j=wasAligned=0;0==wasAligned && j<aEntries.numEnds;j++) {
					if(0 < aEntries.ends[j].numEntries) {
						wasAligned = 1;
					}
				}
				if(1 == wasAligned) {
					numAligned++;
				}
				/* Print it out */
				AlignedReadPrint(&aEntries,
						outputFP);
			}
			AlignedReadFree(&aEntries);
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
		/* Reinitialize file pointer */
		CloseTmpGZFile(&data[i].notAlignedFP, 
				&data[i].notAlignedFileName,
				0);
		if(!(data[i].notAlignedFP=gzopen(data[i].notAlignedFileName, "rb"))) {
			PrintError(FnName,
					data[i].notAlignedFileName,
					"Could not re-open file for reading",
					Exit,
					OpenFileError);
		}

		continueReading=1;
		while(continueReading==1) {
			RGMatchesInitialize(&m);
			if(RGMatchesRead(data[i].notAlignedFP,
						&m) == EOF) {
				continueReading = 0;
			}
			else {
				if(VERBOSE >= 0 && numNotAligned%ALIGN_ROTATE_NUM==0) {
					fprintf(stderr, "\r[%d]", numNotAligned);
				}
				numNotAligned++;
				RGMatchesPrint(notAlignedFP,
						&m);

				RGMatchesFree(&m);
			}
		}
	}
	if(VERBOSE >=0) {
		fprintf(stderr, "\r[%d]\n", numNotAligned);
	}

	/* Close tmp output files */
	for(i=0;i<numThreads;i++) {
		CloseTmpGZFile(&data[i].outputFP, &data[i].outputFileName, 1);
		CloseTmpGZFile(&data[i].notAlignedFP, &data[i].notAlignedFileName, 1);
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
	AlignedReadFree(&aEntries);
}

/* TODO */
void *RunDynamicProgrammingThread(void *arg)
{
	/* Recover arguments */
	ThreadData *data = (ThreadData *)(arg);
	gzFile inputFP=data->inputFP;
	gzFile outputFP=data->outputFP;
	gzFile notAlignedFP = data->notAlignedFP;
	RGBinary *rg=data->rg;
	int space=data->space;
	int offsetLength=data->offsetLength;
	int usePairedEndLength=data->usePairedEndLength;
	int pairedEndLength=data->pairedEndLength;
	int mirroringType=data->mirroringType;
	int forceMirroring=data->forceMirroring;
	ScoringMatrix *sm = data->sm;
	int alignmentType=data->alignmentType;
	int bestOnly=data->bestOnly;
	int threadID=data->threadID;
	int avgMismatchQuality=data->avgMismatchQuality;
	double mismatchScore=data->mismatchScore;
	/* Local variables */
	/*
	   char *FnName = "RunDynamicProgrammingThread";
	   */
	AlignedRead aEntries;
	int32_t i, wasAligned;
	int numAlignedRead=0;
	RGMatches m;
	int numMatches=0;
	int ctrOne=0;
	int ctrTwo=0;

	/* Initialize */
	RGMatchesInitialize(&m);
	AlignedReadInitialize(&aEntries);

	/* Go through each read in the match file */
	while(EOF != RGMatchesRead(inputFP, 
				&m)) {
		numMatches++;
		numAlignedRead = 0;
		ctrOne=ctrTwo=0;

		if(VERBOSE >= 0 && numMatches%ALIGN_ROTATE_NUM==0) {
			fprintf(stderr, "\rthread:%d\t[%d]", threadID, numMatches);
		}

		/* Update the number of local alignments performed */
		data->numLocalAlignments += AlignRGMatches(&m,
				rg,
				&aEntries,
				space,
				offsetLength,
				sm,
				alignmentType,
				bestOnly,
				usePairedEndLength,
				pairedEndLength,
				mirroringType,
				forceMirroring);

		/* Output alignment */
		for(i=wasAligned=0;0==wasAligned && i<aEntries.numEnds;i++) {
			if(0 < aEntries.ends[i].numEntries) {
				wasAligned = 1;
			}
		}
		if(1 == wasAligned) {
			/* Remove duplicates */
			AlignedReadRemoveDuplicates(&aEntries,
					AlignedEntrySortByAll);
			/* Updating mapping quality */
			AlignedReadUpdateMappingQuality(&aEntries, 
					mismatchScore, 
					avgMismatchQuality);
			/* Print */
			AlignedReadPrint(&aEntries,
					outputFP);
		}
		else {
			RGMatchesPrint(notAlignedFP,
					&m);
		}

		/* Free memory */
		AlignedReadFree(&aEntries);
		RGMatchesFree(&m);
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\rthread:%d\t[%d]", threadID, numMatches);
	}

	return arg;
}
