#include <stdlib.h>
#include <stdio.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/AlignEntries.h"
#include "Definitions.h"
#include "ReadInputFiles.h"

void ReadInputFilterAndSplitByChr(char *inputFileName,
		int binaryInput,
		int uniqueReads,
		int bestScoreReads,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int pairedEnd,
		int uniqueReadsPaired,
		int bestScoreReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int meanDistancePaired,
		int chrAbPaired,
		int inversionsPaired,
		char *outputID,
		char *outputDir,
		char *tmpDir,
		int timing,
		TmpFP **tmpFPs,
		int *numTmpFPs)
{
	char *FnName="ReadInputFilterAndSplitByChr";
	FILE *fp=NULL;
	int i;
	int64_t counter, numFiltered;
	AlignEntries a;
	char chrAbFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fpChrAb=NULL;
	char inversionsFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fpInversions=NULL;
	cahr *notReportedFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fpNotReported=NULL;

	/* Open the input file */
	if(!(fp=fopen(inputFileName, "rb"))) {
		PrintError(FnName,
				inputFileName,
				"Could not open inputFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Create output file names */
	sprintf(chrAbFileName, "%s%s.inversion.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			BFAST_CHRAB_FILE_EXTENSION);
	sprintf(inversionsFileName, "%s%s.inversion.%s.%s",
			outputDir,
			PROGRAM_NAME,
			outputID,
			BFAST_INVERSIONS_FILE_EXTENSION);
	sprintf(notReportedFileName, "%s%s.not.reported.%s.%s"
			outputDir,
			PROGRAM_NAME,
			outputID,
			BFAST_NOT_REPORTED_FILE_EXTENSION);
	/* Open output files, if necessary */
	if(chrAbPaired == 1) {
		if(!(fpChrAb=fopen(chrAbFileName, "wb"))) {
			PrintError(FnName,
					chrAbFileName,
					"Could not open chrAbFileName for writing",
					Exit,
					OpenFileError);
		}
	}
	if(inversionsPaired == 1) {
		if(!(fpInversions=fopen(inversionsFileName, "wb"))) {
			PrintError(FnName,
					inversionsFileName,
					"Could not open inversionsFileName for writing",
					Exit,
					OpenFileError);
		}
	}
	if(!(fpNotReported=fopen(notReportedFileName, "wb"))) {
		PrintError(FnName,
				notReportedFileName,
				"Could not open notReportedFileName for writing",
				Exit,
				OpenFileError);
	}

	/* Preallocate tmp file pointers */
	(*numTmpFPs) = endChr-startChr+1;
	(*tmpFPs) = malloc(sizeof(TmpFP)*(*numTmpFPs));
	if(NULL==(*tmpFPs)) {
		PrintError(FnName,
				"tmpFPs",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Initialize */
	AlignEntriesInitialize(&a);
	for(i=0;i<(*numTmpFPs);i++) {
		(*tmpFPs)[i].numEntries = 0;
		(*tmpFPs)[i].tmpFileName = NULL;
		(*tmpFPs)[i].tmpFP = NULL;
		/* Open a temporary file */
		(*tmpFPs)[i].tmpFP = OpenTmpFile(tmpDir, &(*tmpFPs)[i].tmpFileName);
	}

	/* Go through each read */
	if(VERBOSE >= 0) {
		fprintf(stderr, "Filtering alignments, currently on:\n0");
	}
	counter = 0;
	numFiltered = 0;
	while(EOF != AlignEntriesRead(&a, fp)) {
		if(VERBOSE >= 0 && counter%ALIGNENTRIES_READ_ROTATE_NUM==0) {
			fprintf(stderr, "\r%lld",
					(long long int)counter);
		}
		/* Filter */
		numFiltered+=Filter(&a,
				algorithmReads,
				minScoreReads,
				startChr,        
				startPos,
				endChr,        
				endPos,
				pairedEnd,        
				algorithmReadsPaired,
				minScoreReadsPaired,
				minDistancePaired,          
				maxDistancePaired,
				meanDistancePaired,
				chrAbPaired,                        
				inversionsPaired,
				tmpFPs,
				(*numTmpFPs),
				fpChrAb,
				fpInversions,
				fpNotReported);

		/* Free memory */
		AlignEntriesFree(&a);
		/* Increment counter */
		counter++;
	}
	if(VERBOSE>=0) {
		fprintf(stderr, "\r%lld",
				(long long int)counter);
		fprintf(stderr, "Filtered %lld out of %lld reads.\n",
				numFiltered,
				counter);
	}

	/* Close output files, if necessary */
	fclose(fpNotReported);
	if(inversionsPaired == 1) {
		fclose(fpInversions);
	}
	if(chrAbPaired == 1) {
		fclose(fpChrAb);
	}

	/* Close the input file */
	fclose(fp);

}
