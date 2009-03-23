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
#include "../balign/Align.h"
#include "../balign/Definitions.h"
#include "Repair.h"
#include "RunRepair.h"

/* TODO */
void RunRepair(RGBinary *rg,
		char *unpairedFileName,
		char *scoringMatrixFileName,
		int alignmentType,
		int space,
		int binaryInput,
		int numThreads,
		int minPairedEndDistance,
		int maxPairedEndDistance,
		int mirroringType,
		int strandedness,
		char *outputID,
		char *outputDir,
		char *tmpDir,
		int binaryOutput,
		int *totalAlignedTime,
		int *totalFileHandlingTime)
{
	char *FnName = "RunRepair";
	FILE *outputFP=NULL;
	FILE *notRepairedFP=NULL;
	FILE *unpairedFP=NULL;
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	char notRepairedFileName[MAX_FILENAME_LENGTH]="\0";

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

	/* Create output file name */
	sprintf(outputFileName, "%s%s.repaired.file.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			BFAST_ALIGNED_FILE_EXTENSION);
	/* Create not repaired file name */
	sprintf(notRepairedFileName, "%s%s.not.repaired.file.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			BFAST_ALIGNED_FILE_EXTENSION);

	/* Open output file */
	if((outputFP=fopen(outputFileName, "wb"))==0) {
		PrintError(FnName,
				outputFileName,
				"Could not open file for writing",
				Exit,
				OpenFileError);
	}

	/* Open not repaired file */
	if((notRepairedFP=fopen(notRepairedFileName, "wb"))==0) {
		PrintError(FnName,
				notRepairedFileName,
				"Could not open file for writing",
				Exit,
				OpenFileError);
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "Will output repaired reads to %s.\n", outputFileName);
		fprintf(stderr, "Will output unrepaired reads to %s.\n", notRepairedFileName);
		fprintf(stderr, "%s", BREAK_LINE);
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Reading unpaired file from %s.\n",
				unpairedFileName);
	}

	/* Open current unpaired file */
	if((unpairedFP=fopen(unpairedFileName, "rb"))==0) {
		PrintError(FnName,
				unpairedFileName,
				"Could not open file for reading",
				Exit,
				OpenFileError);
	}

	RunRepairHelper(unpairedFP,
			rg,
			scoringMatrixFileName,
			alignmentType,
			space,
			binaryInput,
			numThreads,
			minPairedEndDistance,
			maxPairedEndDistance,
			mirroringType,
			strandedness,
			tmpDir,
			outputFP,
			notRepairedFP,
			binaryOutput,
			totalAlignedTime,
			totalFileHandlingTime);

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Close the unpaired file */
	fclose(unpairedFP);

	/* Close output file */
	fclose(outputFP);

	/* Close not repaired file */
	fclose(notRepairedFP);
}

/* TODO */
void RunRepairHelper(FILE *unpairedFP,
		RGBinary *rg,
		char *scoringMatrixFileName,
		int alignmentType,
		int space,
		int binaryInput,
		int numThreads,
		int minPairedEndDistance,
		int maxPairedEndDistance,
		int mirroringType,
		int strandedness,
		char *tmpDir,
		FILE *outputFP,
		FILE *notRepairedFP,
		int binaryOutput,
		int *totalAlignedTime,
		int *totalFileHandlingTime)
{
	char *FnName="RunRepairHelper";
	/* local variables */
	ScoringMatrix sm;
	AlignedRead a;
	int32_t i, j;

	int32_t wasAligned;
	int continueReading=0;
	int numUnpaired=0;
	int numRepaired=0;
	int numNotRepaired=0;
	int startTime, endTime;
	int64_t numLocalAlignments=0;
	AlignedRead aEntries;
	/* Thread specific data */
	RepairThreadData *data;
	pthread_t *threads=NULL;
	int errCode;
	void *status;

	/* Initialize */
	AlignedReadInitialize(&a);
	ScoringMatrixInitialize(&sm);

	/* Allocate memory for thread arguments */
	data = malloc(sizeof(RepairThreadData)*numThreads);
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

	/**/
	/* Split the input file into equal temp files for each thread */
	/* Could put this in a separate function */
	/**/

	/* Open temp files for the threads */
	for(i=0;i<numThreads;i++) {
		data[i].inputFP = OpenTmpFile(tmpDir, &data[i].inputFileName);
		data[i].outputFP = OpenTmpFile(tmpDir, &data[i].outputFileName);
		data[i].notRepairedFP = OpenTmpFile(tmpDir, &data[i].notRepairedFileName); 
	}

	/* Go through each read in the unpaired file and partition them for the threads */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Filtering and partitioning unpaired reads for threads...\n0");
	}
	i=0;
	while(EOF!=AlignedReadRead(&a, unpairedFP, binaryInput)) {
		if(VERBOSE >= 0 && numUnpaired%PARTITION_MATCHES_ROTATE_NUM==0) {
			fprintf(stderr, "\r[%d]", numUnpaired);
		}
		/* Print unpaired to temp file */
		AlignedReadPrint(&a, data[i].inputFP, binaryInput);
		/* Increment */
		i = (i+1)%numThreads;
		/* Free memory */
		AlignedReadFree(&a);
		/* increment */
		numUnpaired++;
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\r[%d]\n", numUnpaired);
		fprintf(stderr, "Initially filtered %d out of %d.\n",
				numNotRepaired,
				numUnpaired);
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
		data[i].minPairedEndDistance = minPairedEndDistance;
		data[i].maxPairedEndDistance = maxPairedEndDistance;
		data[i].mirroringType = mirroringType;
		data[i].strandedness = strandedness;
		data[i].binaryInput=binaryInput;
		data[i].binaryOutput=binaryOutput;
		data[i].sm = &sm;
		data[i].alignmentType = alignmentType;
		data[i].numLocalAlignments = 0;
		data[i].threadID = i;
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Processing %d unpaired reads.\n",
				numUnpaired - numNotRepaired
			   );
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Performing repair...\n");
		fprintf(stderr, "Currently on:\n0");
	}

	/* Start align timer */
	startTime = time(NULL);

	/* Create threads */
	for(i=0;i<numThreads;i++) {
		/* Start thread */
		errCode = pthread_create(&threads[i], /* thread struct */
				NULL, /* default thread attributes */
				RunRepairHelperThread, /* start routine */
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
		fseek(data[i].outputFP, 0, SEEK_SET);
		assert(NULL!=data[i].outputFP);
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Repair complete.\n");
	}

	/* End align timer */
	endTime = time(NULL);
	(*totalAlignedTime) += endTime - startTime;

	/* Start file handling timer */
	startTime = time(NULL);

	/* Close tmp input files */
	for(i=0;i<numThreads;i++) {
		CloseTmpFile(&data[i].inputFP, &data[i].inputFileName);
		numLocalAlignments += data[i].numLocalAlignments;
	}

	/* Merge all the aligned reads from the threads */
	if(VERBOSE >=0) {
		fprintf(stderr, "Merging and outputting repaired reads from threads...\n[0]");
	}
	AlignedReadInitialize(&aEntries);
	numRepaired=0;
	continueReading=1;
	while(continueReading==1) {
		/* Get an align from a thread */
		continueReading=0;
		for(i=0;i<numThreads;i++) {
			/* Read in the align entries */
			if(EOF != AlignedReadRead(&aEntries, data[i].outputFP, binaryOutput)) {
				if(VERBOSE >=0 && numRepaired%ALIGN_ROTATE_NUM == 0) {
					fprintf(stderr, "\r[%d]", numRepaired);
				}
				continueReading=1;
				/* Update the number that were aligned */
				for(j=wasAligned=0;0==wasAligned && j<aEntries.numEnds;j++) {
					if(0 < aEntries.ends[j].numEntries) {
						wasAligned = 1;
					}
				}
				if(1 == wasAligned) {
					numRepaired++;
				}
				/* Print it out */
				AlignedReadPrint(&aEntries,
						outputFP,
						binaryOutput);
			}
			AlignedReadFree(&aEntries);
		}
	}
	if(VERBOSE >=0) {
		fprintf(stderr, "\r[%d]\n", numRepaired);
	}

	/* Merge the not repaired tmp files */
	if(VERBOSE >=0) {
		fprintf(stderr, "Merging and outputting unpaired reads from threads and initial filter...\n");
	}
	for(i=0;i<numThreads;i++) {
		fseek(data[i].notRepairedFP, 0, SEEK_SET);

		continueReading=1;
		while(continueReading==1) {
			AlignedReadInitialize(&a);
			if(EOF == AlignedReadRead(&a, data[i].notRepairedFP, binaryInput)) {
				continueReading = 0;
			}
			else {
				if(VERBOSE >= 0 && numNotRepaired%ALIGN_ROTATE_NUM==0) {
					fprintf(stderr, "\r[%d]", numNotRepaired);
				}
				numNotRepaired++;
				AlignedReadPrint(&a, notRepairedFP, binaryOutput);
				AlignedReadFree(&a);
			}
		}
	}
	if(VERBOSE >=0) {
		fprintf(stderr, "\r[%d]\n", numNotRepaired);
	}

	/* Close tmp output files */
	for(i=0;i<numThreads;i++) {
		CloseTmpFile(&data[i].outputFP, &data[i].outputFileName);
		CloseTmpFile(&data[i].notRepairedFP, &data[i].notRepairedFileName);
	}

	/* Start file handling timer */
	endTime = time(NULL);
	(*totalFileHandlingTime) += endTime - startTime;

	if(VERBOSE >=0) {
		fprintf(stderr, "Performed %lld local alignments.\n", (long long int)numLocalAlignments);
		fprintf(stderr, "Repaired %d unpaired reads.\n", numRepaired);
		fprintf(stderr, "Could not repair %d reads.\n", numNotRepaired); 
		fprintf(stderr, "Outputting complete.\n");
	}
	assert(numRepaired + numNotRepaired == numUnpaired);

	/* Free memory */
	free(data);
	free(threads);
	AlignedReadFree(&aEntries);
	/* Free scores */
	ScoringMatrixFree(&sm);
}

/* TODO */
void *RunRepairHelperThread(void *arg)
{
	/* Recover arguments */
	RepairThreadData *data = (RepairThreadData *)(arg);
	FILE *inputFP=data->inputFP;
	FILE *outputFP=data->outputFP;
	FILE *notRepairedFP = data->notRepairedFP;
	RGBinary *rg=data->rg;
	int space=data->space;
	int minPairedEndDistance=data->minPairedEndDistance;
	int maxPairedEndDistance=data->maxPairedEndDistance;
	int mirroringType=data->mirroringType;
	int strandedness=data->strandedness;
	int binaryInput=data->binaryInput;
	int binaryOutput=data->binaryOutput;
	ScoringMatrix *sm = data->sm;
	int alignmentType=data->alignmentType;
	int threadID=data->threadID;
	/* Local variables */
	/*
	   char *FnName = "RunDynamicProgrammingThread";
	   */
	AlignedRead src, dest;
	int32_t i, wasAligned;
	int numRepairedRead=0;
	int numUnpaired=0;
	int ctrOne=0;
	int ctrTwo=0;

	/* Initialize */
	AlignedReadInitialize(&src);
	AlignedReadInitialize(&dest);

	/* Go through each read in the unpaired file */
	while(EOF!=AlignedReadRead(&src, inputFP, binaryInput)) {
		numUnpaired++;
		numRepairedRead = 0;
		ctrOne=ctrTwo=0;

		fprintf(stderr, "\rthread:%d\t[%d]", threadID, numUnpaired);

		/* Update the number of local alignments performed */
		data->numLocalAlignments += Repair(&src,
				rg,
				&dest,
				space,
				sm,
				alignmentType,
				minPairedEndDistance,
				maxPairedEndDistance,
				mirroringType,
				strandedness);

		/* Output alignment */
		wasAligned=1;
		for(i=0;1==wasAligned && i<dest.numEnds;i++) {
			if(0 == dest.ends[i].numEntries) {
				wasAligned = 0;
			}
		}
		if(1 == wasAligned) {
			/* Remove duplicates */
			AlignedReadRemoveDuplicates(&dest,
					AlignedEntrySortByAll);
			/* Print */
			AlignedReadPrint(&dest,
					outputFP,
					binaryOutput);
		}
		else {
			AlignedReadPrint(&src,
					notRepairedFP,
					binaryOutput);
		}

		/* Free memory */
		AlignedReadFree(&src);
		AlignedReadFree(&dest);
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\rthread:%d\t[%d]", threadID, numUnpaired);
	}

	return arg;
}
