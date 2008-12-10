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
#define BREPORT_ROTATE_NUM 100000
#define MAX_LINE_LENGTH 4028
/* For zero-based, set this to 0, otherwise to 1 */
#define SUBTRACT 0
#define BED "BED"
#define WIG "WIG"
#define DELIMINATORS "\t\n "

/* Converts a bfast .baf file to bfast .bed and .wig files.  The
 * .baf file must be generated from bpostprocess, with each
 * read having a unique alignment.  The files generated files 
 * do not conform with UCSC standards but are nonetheless more 
 * verbose.  There will be two files per contig.
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
	tmpFile->memory = 0;
	tmpFile->numEntries = 0;
}

void TmpFileUpdateMetaData(TmpFile *tmpFile,
		AlignEntries *a)
{
	TmpFileUpdateMetaDataHelper(tmpFile, &a->entriesOne[0]);
	if(PairedEnd == a->pairedEnd) {
		TmpFileUpdateMetaDataHelper(tmpFile, &a->entriesTwo[0]);
	}
	tmpFile->memory += AlignEntriesGetSize(a);
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
	char *FnName="SplitIntoTmpFilesByContig";
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
		if(counter%BREPORT_ROTATE_NUM==0) {
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
			assert(a.numEntriesOne == 1 && a.numEntriesTwo == 1);
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
	fprintf(stderr, "Moved %lld entries into a tmp file\n",
			(long long int)counter);

	/* Close the input file */
	fclose(fpIn);
}

void SplitEntriesAndPrint(FILE *outputFP,
		TmpFile *tmpFile, 
		char *tmpDir,
		int32_t memoryLimit)
{
	char *FnName="SplitEntriesAndPrint";
	int32_t meanPos, meanContig;
	AlignEntries **entries=NULL;
	AlignEntries a;
	int64_t numEntries=0;
	TmpFile belowTmpFile, aboveTmpFile;
	int32_t belowMinPosition, belowMinContig, belowMaxPosition, belowMaxContig;
	int32_t aboveMinPosition, aboveMinContig, aboveMaxPosition, aboveMaxContig;
	int64_t i;

	if(tmpFile->numEntries <= 0) {
		return;
	}

	/* Move to the beginning of the tmp file */
	fseek(tmpFile->FP, 0, SEEK_SET);

	/* Check if we should print or split */
	if(tmpFile->memory <= (memoryLimit/pow(2, 20))) { /* Assumes MEGABYTES */
		/* Sort and print */
		PrintContigPos(stderr,
				tmpFile->startContig,
				tmpFile->startPos);
		assert(tmpFile->numEntries > 0);


		/* Allocate memory for the entries */
		entries = malloc(sizeof(AlignEntries*)*tmpFile->numEntries);
		if(NULL == entries) {
			PrintError(FnName,
					"entries",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		for(i=0;i<tmpFile->numEntries;i++) {
			entries[i] = malloc(sizeof(AlignEntries));
			if(NULL == entries[i]) {
				PrintError(FnName,
						"entries[i]",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
			AlignEntriesInitialize(entries[i]);
		}

		/* Read in, sort, and print */
		numEntries = 0;
		while(EOF != AlignEntriesRead(entries[numEntries],
					tmpFile->FP,
					PairedEndDoesNotMatter,
					SpaceDoesNotMatter,
					BinaryInput)) {
			assert(numEntries < tmpFile->numEntries);
			numEntries++;
		}
		assert(numEntries == tmpFile->numEntries);
		/* Sort */
		AlignEntriesMergeSortAll(entries, 
				0,
				numEntries);
		/* Print and Free memory */
		for(i=0;i<tmpFile->numEntries;i++) {
			/* Print */
			AlignEntriesPrint(entries[i],
					outputFP,
					BinaryOutput);
			/* Free memory */
			AlignEntriesFree(entries[i]);
			free(entries[i]);
			entries[i]=NULL;
		}
		free(entries);
		entries=NULL;
	}
	else if(tmpFile->startContig == tmpFile->endContig && 
			tmpFile->startPos == tmpFile->endPos) {
		PrintError(FnName,
				NULL,
				"Could not split the file any further.  Try increasing your memory limit.",
				Exit,
				OutOfRange);
	}
	else {
		/* Split and recurse */

		/* Initialize */
		AlignEntriesInitialize(&a);
		TmpFileOpen(&belowTmpFile, tmpDir);
		TmpFileOpen(&aboveTmpFile, tmpDir);
		if(tmpFile->startContig == tmpFile->endContig) {
			meanPos = (tmpFile->startPos + tmpFile->endPos)/2;
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
			belowTmpFile.startContig = tmpFile->startContig;
			belowTmpFile.startPos = tmpFile->startPos;
			belowTmpFile.endContig = meanContig;
			belowTmpFile.endPos = INT_MAX;
			aboveTmpFile.startContig = meanContig + 1;
			aboveTmpFile.startPos = 1;
			aboveTmpFile.endContig = tmpFile->endContig;
			aboveTmpFile.endPos = tmpFile->endPos;
		}

		belowMinPosition = INT_MAX;
		belowMinContig = INT_MAX;
		belowMaxPosition = 0;
		belowMaxContig = 0;
		aboveMinPosition = INT_MAX;
		aboveMinContig = INT_MAX;
		aboveMaxPosition = 0;
		aboveMaxContig = 0;

		/* Split */
		while(EOF != AlignEntriesRead(&a,
					tmpFile->FP,
					PairedEndDoesNotMatter,
					SpaceDoesNotMatter,
					BinaryInput)) {

			/* Print to the appropriate file */
			if(a.entriesOne[0].contig < belowTmpFile.endContig ||
					(a.entriesOne[0].contig == belowTmpFile.endContig && a.entriesOne[0].position < belowTmpFile.endPos)) {
				AlignEntriesPrint(&a,
						belowTmpFile.FP, 
						BinaryOutput);
				belowTmpFile.numEntries++;
				if(a.entriesOne[0].contig < belowMinContig ||
						(a.entriesOne[0].contig == belowMinContig && a.entriesOne[0].position < belowMinPosition)) {
					belowMinContig = a.entriesOne[0].contig;
					belowMinPosition = a.entriesOne[0].position;
				}
				if(belowMaxContig < a.entriesOne[0].contig ||
						(belowMaxContig == a.entriesOne[0].contig && belowMaxPosition < a.entriesOne[0].position)) {
					belowMaxContig = a.entriesOne[0].contig;
					aboveMaxPosition = a.entriesOne[0].position;
				}
			}
			else {
				AlignEntriesPrint(&a,
						aboveTmpFile.FP, 
						BinaryOutput);
				aboveTmpFile.numEntries++;
				if(a.entriesOne[0].contig < aboveMinContig ||
						(a.entriesOne[0].contig == aboveMinContig && a.entriesOne[0].position < aboveMinPosition)) {
					aboveMinContig = a.entriesOne[0].contig;
					aboveMinPosition = a.entriesOne[0].position;
				}
				if(aboveMaxContig < a.entriesOne[0].contig ||
						(aboveMaxContig == a.entriesOne[0].contig && aboveMaxPosition < a.entriesOne[0].position)) {
					aboveMaxContig = a.entriesOne[0].contig;
					aboveMaxPosition = a.entriesOne[0].position;
				}
			}
		}

		/* Update ranges */
		belowTmpFile.startContig = belowMinContig;
		belowTmpFile.startPos = belowMinPosition;
		belowTmpFile.endContig = belowMaxContig;
		belowTmpFile.endPos  = belowMaxPosition;
		aboveTmpFile.startContig = aboveMinContig;
		aboveTmpFile.startPos = aboveMinPosition;
		aboveTmpFile.endContig = aboveMaxContig;
		aboveTmpFile.endPos = aboveMaxPosition;

		/* Recurse on the two */
		SplitEntriesAndPrint(outputFP,
				&belowTmpFile,
				tmpDir,
				memoryLimit);
		SplitEntriesAndPrint(outputFP,
				&aboveTmpFile,
				tmpDir,
				memoryLimit);

		/* Close the files */
		TmpFileClose(&belowTmpFile);
		TmpFileClose(&aboveTmpFile);
	}
}

int main(int argc, char *argv[])
{

	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	int32_t memoryLimit=0;
	char tmpDir[MAX_FILENAME_LENGTH]="\0";
	TmpFile tmpFile;
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *outputFP=NULL;

	if(argc == 4) {
		strcpy(inputFileName, argv[1]);
		memoryLimit = atoi(argv[2]);
		strcpy(tmpDir, argv[3]);

		/* Move all into a tmp file */
		fprintf(stderr, "%s", BREAK_LINE);
		MoveAllIntoTmpFile(inputFileName, &tmpFile, tmpDir);
		fprintf(stderr, "%s", BREAK_LINE);

		/* Create output file name 
		 * TODO */
		strcpy(outputFileName, inputFileName);
		strcat(outputFileName, ".sorted");

		if(!(outputFP == fopen(outputFileName, "wb"))) {
			PrintError(Name,
					outputFileName,
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}

		/* Split entries and print */
		SplitEntriesAndPrint(outputFP,
				&tmpFile,
				tmpDir,
				memoryLimit);

		/* Close files */
		fclose(outputFP);
		TmpFileClose(&tmpFile);

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully!\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast report file name>\n");
		fprintf(stderr, "\t<memory limit in MB>\n");
		fprintf(stderr, "\t<tmp directory>\n");
	}
	return 0;
}
