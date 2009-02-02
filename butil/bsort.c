#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include "../blib/AlignEntry.h"
#include "../blib/AlignEntries.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BLib.h"
#include "../blib/BError.h"
#include "bsort.h"

#define Name "bsort"
#define BSORT_ROTATE_NUM 100000
#define BSORT_MAX_LINE_LENGTH 100

/* Sorts a bfast report file.
 * */

void TmpFileOpen(TmpFile *tmpFile,
		char *tmpDir)
{
	TmpFileInitialize(tmpFile);
	tmpFile->FP = OpenTmpFile(tmpDir, &tmpFile->FileName);
}

void TmpFileClose(TmpFile *tmpFile) 
{
	CloseTmpFile(&tmpFile->FP, &tmpFile->FileName);
	TmpFileInitialize(tmpFile);
}

void TmpFileInitialize(TmpFile *tmpFile)
{
	tmpFile->FP = NULL;
	tmpFile->FileName = NULL;
	tmpFile->startContig = INT_MAX;
	tmpFile->startPos = INT_MAX;
	tmpFile->endContig= 0;
	tmpFile->endPos = 0;
	tmpFile->numEntries = 0;
}

void TmpFileUpdateMetaData(TmpFile *tmpFile,
		AlignEntries *a)
{
	if(SingleEnd == a->pairedEnd ||
			AlignEntryCompareAtIndex(a->entriesOne, 0, a->entriesTwo, 0, AlignEntrySortByContigPos) < 0) {
		TmpFileUpdateMetaDataHelper(tmpFile, &a->entriesOne[0]);
	}
	else {
		TmpFileUpdateMetaDataHelper(tmpFile, &a->entriesTwo[0]);
	}
	tmpFile->numEntries++;
}

void TmpFileUpdateMetaDataHelper(TmpFile *tmpFile,
		AlignEntry *a)
{
	if(a->contig < tmpFile->startContig || 
			(a->contig == tmpFile->startContig && a->position < tmpFile->startPos)) {
		tmpFile->startContig = a->contig;
		tmpFile->startPos = a->position;
	}
	if(tmpFile->endContig < a->contig ||
			(tmpFile->endContig == a->contig && tmpFile->endPos < a->position)) {
		tmpFile->endContig = a->contig;
		tmpFile->endPos = a->position;
	}
}

void MoveAllIntoTmpFile(char *inputFileName, 
		TmpFile *tmpFile,
		char *tmpDir)
{
	char *FnName="MoveAllIntoTmpFile";
	int64_t counter=0;
	FILE *fpIn;
	AlignEntries a;

	/* Open tmp file */
	TmpFileOpen(tmpFile, tmpDir);

	/* Open the input file */
	if(!(fpIn=fopen(inputFileName, "rb"))) {
		PrintError(Name,
				inputFileName,
				"Could not open file for reading",
				Exit,
				OpenFileError);
	}

	/* Move all entries into the tmp file */
	fprintf(stderr, "Moving all entries into a tmp file.  Currently on read:\n0");
	AlignEntriesInitialize(&a);
	while(EOF != AlignEntriesRead(&a, fpIn, PairedEndDoesNotMatter, SpaceDoesNotMatter, BinaryInput)) {
		if(counter%BSORT_ROTATE_NUM==0) {
			fprintf(stderr, "\r%lld",
					(long long int)counter);
		}
		counter++;

		/* Store AlignEntries */
		if(a.numEntriesOne != 1 ||
				(a.pairedEnd == PairedEnd && a.numEntriesTwo != 1)) {
			PrintError(FnName,
					a.readName,
					"Read was not uniquely aligned",
					Exit,
					OutOfRange);
		}
		else {
			assert(a.numEntriesOne == 1 && (SingleEnd == a.pairedEnd || a.numEntriesTwo == 1));
			AlignEntriesPrint(&a, 
					tmpFile->FP,
					BinaryOutput);
			TmpFileUpdateMetaData(tmpFile, 
					&a);
		}
		AlignEntriesFree(&a);
	}
	fprintf(stderr, "\r%lld\n",
			(long long int)counter);
	fprintf(stderr, "Sorting range is contig%lld:%lld to contig%lld:%lld.\n",
			(long long int)tmpFile->startContig,
			(long long int)tmpFile->startPos,
			(long long int)tmpFile->endContig,
			(long long int)tmpFile->endPos
		   );

	/* Close the input file */
	fclose(fpIn);
}

void SplitEntriesAndPrint(FILE *outputFP,
		TmpFile *tmpFile, 
		char *tmpDir,
		int32_t maxNumEntries)
{
	char *FnName="SplitEntriesAndPrint";
	int64_t meanPos;
	int32_t meanContig;
	AlignEntries *entries=NULL;
	AlignEntries **entriesPtr=NULL;
	AlignEntries a;
	int64_t numEntries=0;
	TmpFile belowTmpFile, aboveTmpFile;
	int32_t endContig, endPos;
	int32_t belowMinPos, belowMinContig, belowMaxPos, belowMaxContig;
	int32_t aboveMinPos, aboveMinContig, aboveMaxPos, aboveMaxContig;
	AlignEntry *tmpAlignEntry=NULL;
	int64_t i, numPrinted;

	assert(tmpFile->startContig < tmpFile->endContig ||
			(tmpFile->startContig == tmpFile->endContig && tmpFile->startPos <= tmpFile->endPos));
	assert(0 < tmpFile->numEntries);
	/*
	   if(tmpFile->numEntries <= 0) {
	   return;
	   }
	   */

	/* Move to the beginning of the tmp file */
	fseek(tmpFile->FP, 0, SEEK_SET);

	/* Check if we should print or split */
	if(tmpFile->numEntries <= maxNumEntries) {
		/* Sort and print */
		numPrinted=PrintContigPos(stderr,
				tmpFile->startContig,
				tmpFile->startPos);
		numPrinted+=fprintf(stderr, ": sorting");
		assert(numPrinted<=BSORT_MAX_LINE_LENGTH);
		while(numPrinted<=BSORT_MAX_LINE_LENGTH) {
			fprintf(stderr, " ");
			numPrinted++;
		}
		endContig = tmpFile->endContig;
		endPos = tmpFile->endPos;
		assert(0 < tmpFile->numEntries);

		/* Allocate memory for the entries */
		entriesPtr = malloc(sizeof(AlignEntries*)*tmpFile->numEntries);
		if(NULL == entriesPtr) {
			PrintError(FnName,
					"entriesPtr",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		entries = malloc(sizeof(AlignEntries)*tmpFile->numEntries);
		if(NULL == entries) {
			PrintError(FnName,
					"entries",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Initialize */
		for(i=0;i<tmpFile->numEntries;i++) {
			entriesPtr[i] = &entries[i];
			AlignEntriesInitialize(entriesPtr[i]);
		}

		/* Read in, sort, and print */
		numEntries = 0;
		while(numEntries < tmpFile->numEntries &&
				EOF != AlignEntriesRead(entriesPtr[numEntries],
					tmpFile->FP,
					PairedEndDoesNotMatter,
					SpaceDoesNotMatter,
					BinaryInput)) {
			assert(numEntries < tmpFile->numEntries);
			numEntries++;
		}
		assert(numEntries == tmpFile->numEntries);

		/* Close the file */
		TmpFileClose(tmpFile);
		/* Sort */
		AlignEntriesMergeSortAll(entriesPtr, 
				0,
				numEntries-1);
		/* Print and Free memory */
		for(i=0;i<numEntries;i++) {
			/* Print */
			AlignEntriesPrint(entriesPtr[i],
					outputFP,
					BinaryOutput);
			/* Free memory */
			AlignEntriesFree(entriesPtr[i]);
			entriesPtr[i]=NULL;
		}
		free(entriesPtr);
		entriesPtr=NULL;
		free(entries);
		entries=NULL;
	}
	else if(tmpFile->startContig == tmpFile->endContig && 
			tmpFile->startPos == tmpFile->endPos) {
		PrintError(FnName,
				NULL,
				"Could not split the file any further.  Try increasing your the maximum number of entries.",
				Exit,
				OutOfRange);
	}
	else {
		/* Split and recurse */
		/*
		   fprintf(stderr, "\rSplitting file of size %lld.",
		   (long long int)tmpFile->numEntries);
		   */

		/* Initialize */
		AlignEntriesInitialize(&a);
		TmpFileOpen(&belowTmpFile, tmpDir);
		TmpFileOpen(&aboveTmpFile, tmpDir);
		if(tmpFile->startContig == tmpFile->endContig) {
			meanPos = (tmpFile->startPos + tmpFile->endPos)/2;
			numPrinted=PrintContigPos(stderr,
					tmpFile->startContig,
					tmpFile->startPos);
			numPrinted+=fprintf(stderr, ": splitting file of size %lld by position %lld on contig %lld",
					(long long int)tmpFile->numEntries,
					(long long int)meanPos,
					(long long int)tmpFile->startContig
					);
			belowTmpFile.startContig = tmpFile->startContig;
			belowTmpFile.startPos = tmpFile->startPos;
			belowTmpFile.endContig = tmpFile->endContig;
			belowTmpFile.endPos = meanPos;
			aboveTmpFile.startContig = tmpFile->startContig;
			aboveTmpFile.startPos = meanPos + 1;
			aboveTmpFile.endContig = tmpFile->endContig;
			aboveTmpFile.endPos = tmpFile->endPos;
		}
		else {
			meanContig = (tmpFile->startContig + tmpFile->endContig)/2;
			numPrinted=PrintContigPos(stderr,
					tmpFile->startContig,
					tmpFile->startPos);
			numPrinted+=fprintf(stderr, ": splitting file of size %lld by contig %d",
					(long long int)tmpFile->numEntries,
					meanContig
					);
			belowTmpFile.startContig = tmpFile->startContig;
			belowTmpFile.startPos = tmpFile->startPos;
			belowTmpFile.endContig = meanContig;
			belowTmpFile.endPos = INT_MAX-1;
			aboveTmpFile.startContig = meanContig + 1;
			aboveTmpFile.startPos = 1;
			aboveTmpFile.endContig = tmpFile->endContig;
			aboveTmpFile.endPos = tmpFile->endPos;
		}
		assert(numPrinted<=BSORT_MAX_LINE_LENGTH);
		while(numPrinted<=BSORT_MAX_LINE_LENGTH) {
			fprintf(stderr, " ");
			numPrinted++;
		}

		belowMinPos = INT_MAX;
		belowMinContig = INT_MAX;
		belowMaxPos = 0;
		belowMaxContig = 0;
		aboveMinPos = INT_MAX;
		aboveMinContig = INT_MAX;
		aboveMaxPos = 0;
		aboveMaxContig = 0;

		/* Split */
		while(EOF != AlignEntriesRead(&a,
					tmpFile->FP,
					PairedEndDoesNotMatter,
					SpaceDoesNotMatter,
					BinaryInput)) {

			if(SingleEnd == a.pairedEnd ||
					AlignEntryCompareAtIndex(a.entriesOne, 0, a.entriesTwo, 0, AlignEntrySortByContigPos) < 0) {
				tmpAlignEntry = &a.entriesOne[0];
			}
			else {
				tmpAlignEntry = &a.entriesTwo[0];
			}

			/* Print to the appropriate file */
			if(tmpAlignEntry->contig < belowTmpFile.endContig ||
					(tmpAlignEntry->contig == belowTmpFile.endContig && tmpAlignEntry->position < belowTmpFile.endPos)) {
				/* Print */
				AlignEntriesPrint(&a,
						belowTmpFile.FP, 
						BinaryOutput);
				belowTmpFile.numEntries++;

				/* Update bounds */
				if(tmpAlignEntry->contig < belowMinContig ||
						(tmpAlignEntry->contig == belowMinContig && tmpAlignEntry->position < belowMinPos)) {
					belowMinContig = tmpAlignEntry->contig;
					belowMinPos = tmpAlignEntry->position;
				}
				if(belowMaxContig < tmpAlignEntry->contig || 
						(belowMaxContig == tmpAlignEntry->contig && belowMaxPos < tmpAlignEntry->position)) {
					belowMaxContig = tmpAlignEntry->contig;
					belowMaxPos = tmpAlignEntry->position;
				}
			}
			else {
				/* Print */
				AlignEntriesPrint(&a,
						aboveTmpFile.FP, 
						BinaryOutput);
				aboveTmpFile.numEntries++;

				/* Update bounds */
				if(tmpAlignEntry->contig < aboveMinContig ||
						(tmpAlignEntry->contig == aboveMinContig && tmpAlignEntry->position < aboveMinPos)) {
					aboveMinContig = tmpAlignEntry->contig;
					aboveMinPos = tmpAlignEntry->position;
				}
				if(aboveMaxContig < tmpAlignEntry->contig || 
						(aboveMaxContig == tmpAlignEntry->contig && aboveMaxPos < tmpAlignEntry->position)) {
					aboveMaxContig = tmpAlignEntry->contig;
					aboveMaxPos = tmpAlignEntry->position;
				}
			}
			AlignEntriesFree(&a);
		}

		/* Close tmp file */
		TmpFileClose(tmpFile);

		/* Close tmp files if we are not going to use them */
		if(belowTmpFile.numEntries <= 0) {
			assert(belowTmpFile.startContig < belowMinContig ||
					(belowTmpFile.startContig == belowMinContig && belowTmpFile.startPos <= belowMinPos));
			assert(belowMaxContig < belowTmpFile.endContig ||
					(belowMaxContig == belowTmpFile.endContig && belowMaxPos <= belowTmpFile.endPos));
			belowTmpFile.startContig = belowMinContig;
			belowTmpFile.startPos = belowMinPos;
			belowTmpFile.endContig = belowMaxContig;
			belowTmpFile.endPos = belowMaxPos;
			TmpFileClose(&belowTmpFile);
		}
		if(aboveTmpFile.numEntries <= 0) {
			assert(aboveTmpFile.startContig < belowMinContig ||
					(aboveTmpFile.startContig == belowMinContig && belowTmpFile.startPos <= belowMinPos));
			assert(aboveMaxContig < belowTmpFile.endContig ||
					(aboveMaxContig == belowTmpFile.endContig && belowMaxPos <= belowTmpFile.endPos));
			aboveTmpFile.startContig = aboveMinContig;
			aboveTmpFile.startPos = aboveMinPos;
			aboveTmpFile.endContig = aboveMaxContig;
			aboveTmpFile.endPos = aboveMaxPos;
			TmpFileClose(&aboveTmpFile);
		}

		/* Recurse on the two */
		if(0 < belowTmpFile.numEntries) {
			/* First update based on bounds */
			SplitEntriesAndPrint(outputFP,
					&belowTmpFile,
					tmpDir,
					maxNumEntries);
		}
		if(0 < aboveTmpFile.numEntries) {
			SplitEntriesAndPrint(outputFP,
					&aboveTmpFile,
					tmpDir,
					maxNumEntries);
		}
	}
}

int main(int argc, char *argv[])
{
	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	int32_t maxNumEntries=0;
	char tmpDir[MAX_FILENAME_LENGTH]="\0";
	TmpFile tmpFile;
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *outputFP=NULL;

	if(argc == 4) {
		strcpy(inputFileName, argv[1]);
		maxNumEntries = atoi(argv[2]);
		strcpy(tmpDir, argv[3]);

		/* Move all into a tmp file */
		fprintf(stderr, "%s", BREAK_LINE);
		MoveAllIntoTmpFile(inputFileName, &tmpFile, tmpDir);
		fprintf(stderr, "%s", BREAK_LINE);

		/* Create output file name */ 
		strcpy(outputFileName, inputFileName);
		strcat(outputFileName, ".sorted");

		if(!(outputFP = fopen(outputFileName, "wb"))) {
			PrintError(Name,
					outputFileName,
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}

		/* Split entries and print */
		fprintf(stderr, "Performing initial split, this could take some time.\n");
		SplitEntriesAndPrint(outputFP,
				&tmpFile,
				tmpDir,
				maxNumEntries);
		fprintf(stderr, "\n");

		/* Close files */
		fclose(outputFP);

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully!\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast report file name>\n");
		fprintf(stderr, "\t<maximum number of entries when sorting>\n");
		fprintf(stderr, "\t<tmp directory>\n");
	}
	return 0;
}
