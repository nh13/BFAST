#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "../blib/AlignedEntry.h"
#include "../blib/AlignedEnd.h"
#include "../blib/AlignedRead.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BLib.h"
#include "../blib/BError.h"
#include "breport.h"

#define Name "breport"
#define BREPORT_ROTATE_NUM 100000
#define MAX_LINE_LENGTH 4028
/* For zero-based, set this to 1, otherwise to 0 */
#define SUBTRACT 1
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

void PrintEntriesToBedAndWig(AlignedEnd *a,
		RGBinary *rg,
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
	int fCounts[6] = {0,0,0,0,0,0};
	int rCounts[6] = {0,0,0,0,0,0};
	int total, totalF, totalR;
	int32_t start, cur, i, j, tempJ;
	int numEntries = a->numEntries;

	if(numEntries <= 0) {
		return;
	}

	/* Initialize */
	for(curPos=startPos, start=0;
			start<numEntries &&
			curPos <= endPos;
			curPos++) { /* For each position */
		/* Initialize */
		fCounts[0] = fCounts[1] = fCounts[2] = fCounts[3] = fCounts[4] = fCounts[5] = 0;
		rCounts[0] = rCounts[1] = rCounts[2] = rCounts[3] = rCounts[4] = fCounts[5] = 0;
		total = totalF = totalR = 0;

		if(curPos < a->entries[start].position) {
			curPos = a->entries[start].position;
		}
		assert(a->entries[start].position <= curPos &&
				curPos <= a->entries[start].position + a->entries[start].referenceLength -1);

		/* Go through every entry at that overlaps this position */
		for(cur = start;
				cur < numEntries &&
				a->entries[cur].position <= curPos;
				cur++) {
			/* Only use if it is within bounds */ 
			if(a->entries[cur].position <= curPos &&
					curPos <= a->entries[cur].position + a->entries[cur].referenceLength - 1) {
				assert(a->entries[cur].contig == contig);
				/* Copy over reference and read.  Adjust if they are on the - strand */
				switch(a->entries[cur].strand) {
					case FORWARD:
						strcpy(reference, a->entries[cur].reference);
						strcpy(read, a->entries[cur].read);
						break;
					case REVERSE:
						GetReverseComplimentAnyCase(a->entries[cur].reference, reference, a->entries[cur].length);
						GetReverseComplimentAnyCase(a->entries[cur].read, read, a->entries[cur].length);
						break;
					default:
						PrintError(FnName,
								"a->entries[cur].strand",
								"Could not understand strand",
								Exit,
								OutOfRange);
				}
				/* Add to counts */
				/* Move to correct position */
				i=0;
				j=0;
				while(a->entries[cur].position + i < curPos) {
					assert(j < a->entries[cur].length);
					/* Only move our position if there is not gap */
					if(reference[j] != GAP) {
						i++;
					}
					j++;
				}
				tempJ = j; /* Save this for bed */
				assert(a->entries[cur].position + i == curPos);
				/************/
				/* WIG FILE */
				/************/
				/* Skip over gaps in the reference */
				while(j < a->entries[cur].length &&
						reference[j] == GAP) {
					j++;
				}
				/* We should not end with a gap so this should be true */
				assert(j < a->entries[cur].length);
				assert(reference[j] != GAP);
				/* Update counts */
				switch(a->entries[cur].strand) {
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
							case 'N':
							case 'n':
								fCounts[5]++;
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
								rCounts[3]++;
								break;
							case 'c':
							case 'C':
								rCounts[2]++;
								break;
							case 'g':
							case 'G':
								rCounts[1]++;
								break;
							case 't':
							case 'T':
								rCounts[0]++;
								break;
							case 'N':
							case 'n':
								rCounts[5]++;
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
								"a->entries[cur].strand",
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
				if(j < a->entries[cur].length - 1 &&
						reference[j] != GAP &&
						read[j] != GAP) { 
					/* The indel starts after the current base */
					j=j+1;
					i=0;
					type='N';
					/* Check gap in reference */
					while(j<a->entries[cur].length && 
							reference[j] == GAP) {
						/* Get the insertion in the read */
						indel[i] = read[j];
						type = 'I';
						i++;
						j++;
					}
					/* Check gap in read */
					while(j<a->entries[cur].length && 
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
										a->entries[cur].contig,
										(long long int)(curPos-SUBTRACT+2),
										(long long int)(curPos+i-SUBTRACT+1),
										a->entries[cur].strand,
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
		/* Get reference base */
		referenceBase = (char)RGBinaryGetBase(rg,
				contig,
				curPos);

		/* Print out counts */
		if(total > 0) {
			i=0;
			j=0;
			switch(referenceBase) {
				case 'a':
				case 'A':
					i=0;
					j=3;
					break;
				case 'c':
				case 'C':
					i=1;
					j=2;
					break;
				case 'g':
				case 'G':
					i=2;
					j=1;
					break;
				case 't':
				case 'T':
					i=3;
					j=0;
					break;
				case 'n':
				case 'N':
					i=4;
					j=4;
					break;
				case GAP:
					i=5;
					j=5;
					break;
				default:
					fprintf(stderr, "\nreferenceBase=%c.\n", referenceBase);
					PrintError(FnName,
							"referenceBase",
							"Could not understand base",
							Exit,
							OutOfRange);
			}
			assert(fCounts[i] + rCounts[j] <= total);
			if(0>fprintf(wigFP, "contig%d %lld %lld %d f[%c %d %d %d %d %d %d %3.2lf] r[%c %d %d %d %d %d %d %3.2lf] t[%3.2lf]\n",
						contig,
						(long long int)(curPos-SUBTRACT),
						(long long int)(curPos-SUBTRACT+1),
						total,
						referenceBase,
						fCounts[0],
						fCounts[1],
						fCounts[2],
						fCounts[3],
						fCounts[4],
						fCounts[5],
						(100.0*fCounts[i])/((double)totalF),
						GetReverseComplimentAnyCaseBase(referenceBase),
						rCounts[0],
						rCounts[1],
						rCounts[2],
						rCounts[3],
						rCounts[4],
						rCounts[5],
						(100.0*rCounts[j])/((double)totalR),
						(100.0*(fCounts[i] + rCounts[j]))/((double)total))) {
							PrintError(FnName,
									"wigFP",
									"Could not write to file",
									Exit,
									WriteFileError);
						}
		}
		/* Update next start index */
		while(start < numEntries && 
				a->entries[start].position + a->entries[start].referenceLength - 1 < curPos + 1) {
			start++;
		}
	}
}

void SplitIntoTmpFilesByContig(char *inputFileName,
		TmpFile **tmpFiles,
		int *numFiles,
		char *tmpDir,
		int startContig,
		int endContig,
		int number,
		int total)
{
	char *FnName="SplitIntoTmpFilesByContig";
	int32_t i;
	int64_t counter=0;
	FILE *fpIn;
	AlignedRead a;

	/* Create tmp files */
	if((*numFiles) <= 0) {
		(*numFiles) = endContig - startContig + 1;

		(*tmpFiles) = malloc(sizeof(TmpFile)*(*numFiles));
		if(NULL == (*tmpFiles)) {
			PrintError(FnName,
					"tmpFiles",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Open tmp files */
		for(i=0;i<(*numFiles);i++) {
			TmpFileOpen(&(*tmpFiles)[i], tmpDir, i+1);
		}
	}

	/* Open the input file */
	fprintf(stderr, "(%d/%d) Splitting %s by contig.\n",
			number,
			total,
			inputFileName);
	if(!(fpIn=fopen(inputFileName, "r"))) {
		PrintError(Name,
				inputFileName,
				"Could not open file for reading",
				Exit,
				OpenFileError);
	}

	/* Split into temporary files */
	fprintf(stderr, "Currently on read:\n0");
	AlignedReadInitialize(&a);
	while(EOF != AlignedReadRead(&a, fpIn, BinaryInput)) {
		if(counter%BREPORT_ROTATE_NUM==0) {
			fprintf(stderr, "\r%lld",
					(long long int)counter);
		}
		counter++;

		/* Print to the appropriate file both end.  We store
		 * "AlignedEntry"s not "AlignedRead" because we split 
		 * the paired end */
		/* Print each end */
		for(i=0;i<a.numEnds;i++) {
			if(a.ends[i].numEntries > 1) {
				PrintError(FnName,
						a.readName,
						"Read i was not uniquely aligned",
						Exit,
						OutOfRange);
			}
			else if(a.ends[i].numEntries == 1) {
				assert(a.ends[i].entries[0].contig > 0 && a.ends[i].entries[0].contig <= (*numFiles));
				AlignedEntryPrint(&a.ends[i].entries[0], 
						(*tmpFiles)[a.ends[i].entries[0].contig-1].FP,
						NTSpace, /* Dont print color space information */
						BinaryOutput);
				/* Update meta-data */
				(*tmpFiles)[a.ends[i].entries[0].contig-1].numEntries++;
				if(a.ends[i].entries[0].position < (*tmpFiles)[a.ends[i].entries[0].contig-1].minPos) {
					(*tmpFiles)[a.ends[i].entries[0].contig-1].minPos = a.ends[i].entries[0].position;
				}
				if(a.ends[i].entries[0].position + a.ends[i].entries[0].referenceLength - 1 > (*tmpFiles)[a.ends[i].entries[0].contig-1].maxPos) {
					(*tmpFiles)[a.ends[i].entries[0].contig-1].maxPos = a.ends[i].entries[0].position + a.ends[i].entries[0].referenceLength - 1;
				}
			}
			else {
				/* Ignore */
			}
		}
		AlignedReadFree(&a);
	}
	fprintf(stderr, "\r%lld\n",
			(long long int)counter);

	/* Close the input file */
	fclose(fpIn);
}

void SplitEntriesAndPrint(RGBinary *rg,
		FILE *bedFP,
		FILE *wigFP,
		TmpFile *tmpFile, /* Should all be from the same contig */
		char *tmpDir,
		int maxNumEntries)
{
	/*
	   char *FnName="SplitEntriesAndPrint";
	   */
	int meanPos;
	AlignedEntry a;
	AlignedEnd entries;
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
		AlignedEndInitialize(&entries);
		/* Allocate memory for the entries */
		AlignedEndAllocate(&entries, 
				tmpFile->numEntries);

		/* Read in, sort, and print */
		numEntries = 0;
		while(EOF != AlignedEntryRead(&entries.entries[numEntries], 
					tmpFile->FP, 
					NTSpace,
					BinaryInput)) {
			numEntries++;
		}
		assert(numEntries == tmpFile->numEntries);
		/* Sort */
		AlignedEndMergeSort(&entries, 
				AlignedEntrySortByContigPos,
				0);
		/* Print Out */
		PrintEntriesToBedAndWig(&entries, 
				rg,
				tmpFile->contig, 
				tmpFile->minPos, 
				tmpFile->maxPos, 
				bedFP, 
				wigFP);
		/* Free memory */
		AlignedEndFree(&entries);
	}
	else {
		/* Split and recurse */

		/* Initialize */
		AlignedEntryInitialize(&a);
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
		while(EOF != AlignedEntryRead(&a,
					tmpFile->FP,
					SpaceDoesNotMatter,
					BinaryInput)) {
			/* Print to the appropriate file */
			assert(a.contig == tmpFile->contig);

			/* Print */
			if(a.position < meanPos) {
				AlignedEntryPrint(&a, 
						belowTmpFile.FP, 
						NTSpace,
						BinaryOutput);
				belowTmpFile.numEntries++;
			}
			if(meanPos <= a.position ||
					meanPos <= a.position + a.referenceLength - 1 ) {
				AlignedEntryPrint(&a,
						aboveTmpFile.FP, 
						NTSpace,
						BinaryOutput);
				aboveTmpFile.numEntries++;
			}
			AlignedEntryFree(&a);
		}

		/* Recurse on the two */
		SplitEntriesAndPrint(rg,
				bedFP,
				wigFP,
				&belowTmpFile,
				tmpDir,
				maxNumEntries);
		SplitEntriesAndPrint(rg,
				bedFP,
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
	RGBinary rg;
	char rgFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *bedFP=NULL;
	char bedFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *wigFP=NULL;
	char wigFileName[MAX_FILENAME_LENGTH]="\0";
	char outputID[MAX_FILENAME_LENGTH]="\0";
	char outputDir[MAX_FILENAME_LENGTH]="\0";

	if(7 <= argc) {
		strcpy(rgFileName, argv[1]);
		maxNumEntries = atoi(argv[2]);
		strcpy(outputID, argv[3]);
		strcpy(outputDir, argv[4]);
		strcpy(tmpDir, argv[5]);

		/* Read in rg file */
		RGBinaryReadBinary(&rg, rgFileName);

		/* Split AlignedRead by contig */
		fprintf(stderr, "%s", BREAK_LINE);
		for(i=6;i<argc;i++) {
			strcpy(inputFileName, argv[i]);
			SplitIntoTmpFilesByContig(inputFileName,
					&tmpFiles,
					&numTmpFiles,
					tmpDir,
					1,
					rg.numContigs,
					i-5,
					argc-6
					);
		}
		fprintf(stderr, "%s", BREAK_LINE);

		/* Output bed and wig files for each contig */
		for(i=0;i<numTmpFiles;i++) {
			assert(tmpFiles[i].contig == i+1);
			if(tmpFiles[i].numEntries > 0) {
				fprintf(stderr, "%s", BREAK_LINE);
				/* Create output file names */
				sprintf(bedFileName, "%sbfast.%s.contig%d.bed",
						outputDir,
						outputID,
						tmpFiles[i].contig);
				sprintf(wigFileName, "%sbfast.%s.contig%d.wig",
						outputDir,
						outputID,
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
				SplitEntriesAndPrint(&rg,
						bedFP,
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
		fprintf(stderr, "\t<bfast reference genome file>\n");
		fprintf(stderr, "\t<maximum number of entries when sorting>\n");
		fprintf(stderr, "\t<output ID>\n");
		fprintf(stderr, "\t<output directory>\n");
		fprintf(stderr, "\t<tmp file directory>\n");
		fprintf(stderr, "\t<bfast report file names>\n");
	}
	return 0;
}
