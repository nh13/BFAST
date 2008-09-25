#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "../blib/AlignEntry.h"
#include "../blib/AlignEntries.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BLib.h"
#include "../blib/BError.h"
#include "breport.h"

#define Name "breport"
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
		char *tmpDir,
		int contig) 
{
	TmpFileInitialize(tmpFile);
	tmpFile->FP = OpenTmpFile(tmpDir, &tmpFile->FileName);
	tmpFile->contig = contig;
}

void TmpFileClose(TmpFile *tmpFile) 
{
	CloseTmpFile(&tmpFile->FP, &tmpFile->FileName);
	TmpFileInitialize(tmpFile);
}

void TmpFileInitialize(TmpFile *tmpFile)
{
	tmpFile->FP = NULL;
	tmpFile->contig = 0;
	tmpFile->minPos = INT_MAX;
	tmpFile->maxPos = 0;
	tmpFile->numEntries = 0;
}

void PrintEntriesToBedAndWig(AlignEntries *a,
		int contig,
		int64_t startPos,
		int64_t endPos,
		FILE *bedFP,
		FILE *wigFP)
{
	char *FnName="PrintEntriesToBedAndWig";
	char reference[SEQUENCE_LENGTH]="\0";
	char read[SEQUENCE_LENGTH]="\0";
	int64_t curPos;
	char referenceBase='n';
	char indel[SEQUENCE_LENGTH]="\0";
	char type='N';
	int fCounts[5] = {0,0,0,0,0};
	int rCounts[5] = {0,0,0,0,0};
	int total, totalF, totalR;
	int32_t start, cur, i, j, tempJ;
	int numEntries = a->numEntriesOne;

	if(numEntries <= 0) {
		return;
	}

	/* Initialize */
	for(curPos=startPos, start=0;
			start<numEntries &&
			curPos <= endPos;
			curPos++) { /* For each position */
		/* Initialize */
		fCounts[0] = fCounts[1] = fCounts[2] = fCounts[3] = fCounts[4] = 0;
		rCounts[0] = rCounts[1] = rCounts[2] = rCounts[3] = rCounts[4] = 0;
		total = totalF = totalR = 0;
		referenceBase = 'n';

		if(curPos < a->entriesOne[start].position) {
			curPos = a->entriesOne[start].position;
		}
		assert(a->entriesOne[start].position <= curPos &&
				curPos <= a->entriesOne[start].position + a->entriesOne[start].referenceLength -1);

		/* Go through every entry at that overlaps this position */
		for(cur = start;
				cur < numEntries &&
				a->entriesOne[cur].position <= curPos;
				cur++) {
			/* Only use if it is within bounds */ 
			if(a->entriesOne[cur].position <= curPos &&
					curPos <= a->entriesOne[cur].position + a->entriesOne[cur].referenceLength - 1) {
				assert(a->entriesOne[cur].contig == contig);
				/* Copy over reference and read.  Adjust if they are on the - strand */
				switch(a->entriesOne[cur].strand) {
					case FORWARD:
						strcpy(reference, a->entriesOne[cur].reference);
						strcpy(read, a->entriesOne[cur].read);
						break;
					case REVERSE:
						GetReverseComplimentAnyCase(a->entriesOne[cur].reference, reference, a->entriesOne[cur].length);
						GetReverseComplimentAnyCase(a->entriesOne[cur].read, read, a->entriesOne[cur].length);
						break;
					default:
						PrintError(FnName,
								"a->entriesOne[cur].strand",
								"Could not understand strand",
								Exit,
								OutOfRange);
				}
				/* Add to counts */
				/* Move to correct position */
				i=0;
				j=0;
				while(a->entriesOne[cur].position + i < curPos) {
					assert(j < a->entriesOne[cur].length);
					/* Only move our position if there is not gap */
					if(reference[j] != GAP) {
						i++;
					}
					j++;
				}
				tempJ = j; /* Save this for bed */
				assert(a->entriesOne[cur].position + i == curPos);
				/************/
				/* WIG FILE */
				/************/
				/* Skip over gaps in the reference */
				while(j < a->entriesOne[cur].length &&
						reference[j] == GAP) {
					j++;
				}
				/* We should not end with a gap so this should be true */
				assert(j < a->entriesOne[cur].length);
				assert(reference[j] != GAP);
				/* Update based on the current base */
				if(referenceBase == 'n') {
					referenceBase = reference[j];
					assert(referenceBase != GAP);
				}
				/* Update counts */
				switch(a->entriesOne[cur].strand) {
					case FORWARD:
						totalF++;
						total++;
						switch(read[j]) {
							case 'a':
							case 'A':
								fCounts[0]++;
								break;
							case 'c':
							case 'C':
								fCounts[1]++;
								break;
							case 'g':
							case 'G':
								fCounts[2]++;
								break;
							case 't':
							case 'T':
								fCounts[3]++;
								break;
							case GAP:
								fCounts[4]++;
								break;
							default:
								PrintError(FnName,
										"read[j]",
										"Could not understand base forward",
										Exit,
										OutOfRange);
						}
						break;
					case REVERSE:
						totalR++;
						total++;
						switch(read[j]) {
							case 'a':
							case 'A':
								rCounts[0]++;
								break;
							case 'c':
							case 'C':
								rCounts[1]++;
								break;
							case 'g':
							case 'G':
								rCounts[2]++;
								break;
							case 't':
							case 'T':
								rCounts[3]++;
								break;
							case GAP:
								rCounts[4]++;
								break;
							default:
								PrintError(FnName,
										"read[j]",
										"Could not understand base reverse",
										Exit,
										OutOfRange);
						}
						break;
					default:
						PrintError(FnName,
								"a->entriesOne[cur].strand",
								"Could not understand strand",
								Exit,
								OutOfRange);
						break;
				}
				/************/
				/* WIG FILE */
				/************/
				j = tempJ;
				/* Don't output if we are on the last letter */
				/* Don't output if we are in a gap already in either the reference or the read */
				if(j < a->entriesOne[cur].length - 1 &&
						reference[j] != GAP &&
						read[j] != GAP) { 
					/* The indel starts after the current base */
					j=j+1;
					i=0;
					type='N';
					/* Check gap in reference */
					while(j<a->entriesOne[cur].length && 
							reference[j] == GAP) {
						/* Get the insertion in the read */
						indel[i] = read[j];
						type = 'I';
						i++;
						j++;
					}
					/* Check gap in read */
					while(j<a->entriesOne[cur].length && 
							read[j] == GAP) {
						/* Get the deletion from the read */
						indel[i] = reference[j];
						type = 'D';
						i++;
						j++;
					}
					/* Print out */
					indel[i]='\0';
					switch(type) {
						case 'N':
							/* Do nothing */
							break;
						case 'I':
						case 'D':
							assert(i>0);
							if(0>fprintf(bedFP, "contig%d %lld %lld %c %c %s\n",
										a->entriesOne[cur].contig,
										(long long int)(curPos-SUBTRACT),
										(long long int)(curPos+i-1-SUBTRACT),
										a->entriesOne[cur].strand,
										type,
										indel)) {
								PrintError(FnName,
										"bedFP",
										"Could not write to file",
										Exit,
										WriteFileError);
							}
							break;
						default:
							PrintError(FnName,
									"type",
									"Could not understand type",
									Exit,
									OutOfRange);
					}
				}
			}
		}
		/* Print out counts */
		if(total > 0) {
			i=0;
			switch(referenceBase) {
				case 'a':
				case 'A':
					i=0;
					break;
				case 'c':
				case 'C':
					i=1;
					break;
				case 'g':
				case 'G':
					i=2;
					break;
				case 't':
				case 'T':
					i=3;
					break;
				case GAP:
					i=4;
					break;
				default:
					PrintError(FnName,
							"referenceBase",
							"Could not understand base",
							Exit,
							OutOfRange);
			}
			assert(fCounts[i] + rCounts[i] <= total);
			if(0>fprintf(wigFP, "contig%d %lld %lld %d %d %d %d %d %d %d %d %d %d %d %3.2lf %3.2lf %3.2lf\n",
						contig,
						(long long int)(curPos-SUBTRACT),
						(long long int)(curPos+1-SUBTRACT),
						total,
						fCounts[0],
						fCounts[1],
						fCounts[2],
						fCounts[3],
						fCounts[4],
						rCounts[0],
						rCounts[1],
						rCounts[2],
						rCounts[3],
						rCounts[4],
						(100.0*fCounts[i])/((double)totalF),
						(100.0*rCounts[i])/((double)totalR),
						(100.0*(fCounts[i] + rCounts[i]))/((double)total))) {
				PrintError(FnName,
						"wigFP",
						"Could not write to file",
						Exit,
						WriteFileError);
			}
		}
		/* Update next start index */
		while(start < numEntries && 
				a->entriesOne[start].position + a->entriesOne[start].referenceLength - 1 < curPos + 1) {
			start++;
		}
	}
}

int SplitIntoTmpFilesByContig(char *inputFileName,
		TmpFile **tmpFiles,
		char *tmpDir,
		int startContig,
		int endContig)
{
	char *FnName="SplitIntoTmpFilesByContig";
	int i;
	int64_t counter=0;
	FILE *fpIn;
	int numFiles = endContig - startContig + 1;
	AlignEntries a;

	/* Open the input file */
	if(!(fpIn=fopen(inputFileName, "r"))) {
		PrintError(Name,
				inputFileName,
				"Could not open file for reading",
				Exit,
				OpenFileError);
	}

	/* Create tmp files */
	(*tmpFiles) = malloc(sizeof(TmpFile)*numFiles);
	if(NULL == (*tmpFiles)) {
		PrintError(FnName,
				"tmpFiles",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Open tmp files */
	for(i=0;i<numFiles;i++) {
		TmpFileOpen(&(*tmpFiles)[i], tmpDir, i+1);
	}

	/* Split into temporary files */
	fprintf(stderr, "Splitting by contig.  Currently on read:\n0");
	AlignEntriesInitialize(&a);
	while(EOF != AlignEntriesRead(&a, fpIn, PairedEndDoesNotMatter, SpaceDoesNotMatter, BinaryInput)) {
		if(counter%BREPORT_ROTATE_NUM==0) {
			fprintf(stderr, "\r%lld",
					(long long int)counter);
		}
		counter++;

		/* Print to the appropriate file both end.  We store
		 * "AlignEntry"s not "AlignEntries" because we split 
		 * the paired end */
		/* Print read one */
		if(a.numEntriesOne > 1) {
			PrintError(FnName,
					a.readName,
					"Read one was not uniquely aligned",
					Exit,
					OutOfRange);
		}
		else if(a.numEntriesOne == 1) {
			assert(a.entriesOne[0].contig > 0 && a.entriesOne[0].contig <= numFiles);
			AlignEntryPrint(&a.entriesOne[0], 
					(*tmpFiles)[a.entriesOne[0].contig-1].FP,
					NTSpace, /* Dont print color space information */
					BinaryOutput);
			/* Update meta-data */
			(*tmpFiles)[a.entriesOne[0].contig-1].numEntries++;
			if(a.entriesOne[0].position < (*tmpFiles)[a.entriesOne[0].contig-1].minPos) {
				(*tmpFiles)[a.entriesOne[0].contig-1].minPos = a.entriesOne[0].position;
			}
			if(a.entriesOne[0].position + a.entriesOne[0].referenceLength - 1 > (*tmpFiles)[a.entriesOne[0].contig-1].maxPos) {
				(*tmpFiles)[a.entriesOne[0].contig-1].maxPos = a.entriesOne[0].position + a.entriesOne[0].referenceLength - 1;
			}
		}
		else {
			/* Ignore */
		}
		/* Print read two */
		if(a.pairedEnd == PairedEnd) {
			if(a.numEntriesTwo > 1) {
				PrintError(FnName,
						a.readName,
						"Read two was not uniquely aligned",
						Exit,
						OutOfRange);
			}
			else if(a.numEntriesTwo == 1) {
				assert(a.entriesTwo[0].contig > 0 && a.entriesTwo[0].contig <= numFiles);
				AlignEntryPrint(&a.entriesTwo[0], 
						(*tmpFiles)[a.entriesTwo[0].contig-1].FP,
						NTSpace, /* Dont print color space information */
						BinaryOutput);
				/* Update meta-data */
				(*tmpFiles)[a.entriesTwo[0].contig-1].numEntries++;
				if(a.entriesTwo[0].position < (*tmpFiles)[a.entriesTwo[0].contig-1].minPos) {
					(*tmpFiles)[a.entriesTwo[0].contig-1].minPos = a.entriesTwo[0].position;
				}
				if(a.entriesTwo[0].position + a.entriesTwo[0].referenceLength - 1 > (*tmpFiles)[a.entriesTwo[0].contig-1].maxPos) {
					(*tmpFiles)[a.entriesTwo[0].contig-1].maxPos = a.entriesTwo[0].position + a.entriesTwo[0].referenceLength - 1;
				}
			}
		else {
			/* Ignore */
		}
		}
		AlignEntriesFree(&a);
	}
	fprintf(stderr, "\r%lld\n",
			(long long int)counter);

	/* Close the input file */
	fclose(fpIn);

	return numFiles;
}

void SplitEntriesAndPrint(FILE *bedFP,
		FILE *wigFP,
		TmpFile *tmpFile, /* Should all be from the same contig */
		char *tmpDir,
		int maxNumEntries)
{
	/*
	char *FnName="SplitEntriesAndPrint";
	*/
	int meanPos;
	AlignEntry a;;
	AlignEntries entries;
	int numEntries=0;
	TmpFile belowTmpFile, aboveTmpFile;

	if(tmpFile->numEntries <= 0) {
		return;
	}

	/* Move to the beginning of the tmp file */
	fseek(tmpFile->FP, 0, SEEK_SET);

	/* Check if we should print or split */
	if(tmpFile->numEntries <= maxNumEntries) {
		fprintf(stderr, "\r[%10d-%-10d]", 
				tmpFile->minPos,
				tmpFile->maxPos);
		assert(tmpFile->numEntries > 0);
		/* Initialize */
		AlignEntriesInitialize(&entries);
		/* Allocate memory for the entries */
		AlignEntriesAllocate(&entries, 
				"DUMMY",
				tmpFile->numEntries,
				0,
				SingleEnd,
				SpaceDoesNotMatter);

			/* Read in, sort, and print */
			numEntries = 0;
		while(EOF != AlignEntryRead(&entries.entriesOne[numEntries], 
					tmpFile->FP, 
					NTSpace,
					BinaryInput)) {
			numEntries++;
		}
		assert(numEntries == tmpFile->numEntries);
		/* Sort */
		AlignEntriesMergeSort(&entries, 
				AlignEntrySortByContigPos,
				0);
		/* Print Out */
		PrintEntriesToBedAndWig(&entries, 
				tmpFile->contig, 
				tmpFile->minPos, 
				tmpFile->maxPos, 
				bedFP, 
				wigFP);
		/* Free memory */
		AlignEntriesFree(&entries);
	}
	else {
		/* Split and recurse */

		/* Initialize */
		AlignEntryInitialize(&a);
		TmpFileOpen(&belowTmpFile, tmpDir, tmpFile->contig);
		TmpFileOpen(&aboveTmpFile, tmpDir, tmpFile->contig);

		/* Where will we split */
		meanPos = (tmpFile->maxPos + tmpFile->minPos)/2;

		/* Update meta-data */
		belowTmpFile.minPos = tmpFile->minPos; 
		belowTmpFile.maxPos = meanPos - 1; 
		aboveTmpFile.minPos = meanPos; 
		aboveTmpFile.maxPos = tmpFile->maxPos; 

		/* Split */
		while(EOF != AlignEntryRead(&a,
					tmpFile->FP,
					SpaceDoesNotMatter,
					BinaryInput)) {
			/* Print to the appropriate file */
			assert(a.contig == tmpFile->contig);

			/* Print */
			if(a.position < meanPos) {
				AlignEntryPrint(&a, 
						belowTmpFile.FP, 
						NTSpace,
						BinaryOutput);
				belowTmpFile.numEntries++;
			}
			if(meanPos <= a.position ||
					meanPos <= a.position + a.referenceLength - 1 ) {
				AlignEntryPrint(&a,
						aboveTmpFile.FP, 
						NTSpace,
						BinaryOutput);
				aboveTmpFile.numEntries++;
			}
		}

		/* Recurse on the two */
		SplitEntriesAndPrint(bedFP,
				wigFP,
				&belowTmpFile,
				tmpDir,
				maxNumEntries);
		SplitEntriesAndPrint(bedFP,
				wigFP,
				&aboveTmpFile,
				tmpDir,
				maxNumEntries);

		/* Close the files */
		TmpFileClose(&belowTmpFile);
		TmpFileClose(&aboveTmpFile);
	}
}

int main(int argc, char *argv[])
{

	char inputFileName[MAX_FILENAME_LENGTH]="\0";
	char tmpDir[MAX_FILENAME_LENGTH]="\0";
	TmpFile *tmpFiles=NULL;
	int numTmpFiles=0;
	int i;
	int maxNumEntries;
	FILE *bedFP=NULL;
	char bedFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *wigFP=NULL;
	char wigFileName[MAX_FILENAME_LENGTH]="\0";
	/* Hard coded */
	int startContig = 1;
	int endContig = 22;

	if(argc == 4) {
		strcpy(inputFileName, argv[1]);
		maxNumEntries = atoi(argv[2]);
		strcpy(tmpDir, argv[3]);

		/* Split AlignEntries by contig */
		fprintf(stderr, "%s", BREAK_LINE);
		numTmpFiles=SplitIntoTmpFilesByContig(inputFileName,
				&tmpFiles,
				tmpDir,
				startContig,
				endContig);
		fprintf(stderr, "%s", BREAK_LINE);

		/* Output bed and wig files for each contig */
		for(i=0;i<numTmpFiles;i++) {
			assert(tmpFiles[i].contig == i+1);
			if(tmpFiles[i].numEntries > 0) {
				fprintf(stderr, "%s", BREAK_LINE);
				/* Create output file names */
				sprintf(bedFileName, "%s.contig%d.bed",
						inputFileName,
						tmpFiles[i].contig);
				sprintf(wigFileName, "%s.contig%d.wig",
						inputFileName,
						tmpFiles[i].contig);
				/* Open output files */
				if(!(bedFP = fopen(bedFileName, "wb"))) {
					PrintError(Name,
							"bedFP",
							"Could not open file for writing",
							Exit,
							OpenFileError);
				}
				if(!(wigFP = fopen(wigFileName, "wb"))) {
					PrintError(Name,
							"wigFP",
							"Could not open file for writing",
							Exit,
							OpenFileError);
				}
				/* Print */
				fprintf(stderr, "On contig%d:%d-%d.  Currently on position:\n0", 
						tmpFiles[i].contig,
						tmpFiles[i].minPos,
						tmpFiles[i].maxPos
					   );
				SplitEntriesAndPrint(bedFP,
						wigFP,
						&tmpFiles[i],
						tmpDir,
						maxNumEntries);
				fprintf(stderr, "\n");
				/* Close files */
				fclose(bedFP);
				fclose(wigFP);
			}
			/* Close */
			TmpFileClose(&tmpFiles[i]);
		}

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully!\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast aligned file name>\n");
		fprintf(stderr, "\t<maximum number of entries when sorting>\n");
		fprintf(stderr, "\t<tmp file directory>\n");
	}
	return 0;
}
