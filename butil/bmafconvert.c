#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/BLib.h"
#include "../blib/BError.h"
#include "bmafconvert.h"

#define Name "bmafconvert"
#define BMFCONVERT_ROTATE_NUM 100000
#define MAX_LINE_LENGTH 4028
/* For zero-based, set this to 0, otherwise to 1 */
#define SUBTRACT 0
#define BED "BED"
#define WIG "WIG"
#define DELIMINATORS "\t\n "

/* Converts a bfast .maf file to bfast .bed and .wig files.  The
 * generated files do not conform with USCS standards but are 
 * nonetheless more verbose.
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

void MAFPrint(FILE *fp,
		MAF *m)
{
	char *FnName="MAFPrint";
	if(0>fprintf(fp, "a score=%lf paired-end=0 read=1\n",
				m->score) ||
			0>fprintf(fp, "s\tcontig%u\t%u\t%d\t%c\t%d\t%s\n",
				m->contig,
				m->position,
				m->alignmentLength,
				m->strand,
				m->referenceLength,
				m->reference)) {
		PrintError(FnName,
				NULL,
				"Could not write MAF",
				Exit,
				WriteFileError);
	}
	if(0>fprintf(fp, "s\t%s\t%d\t%d\t%c\t%d\t%s\n\n",
				"NAME",
				0,
				m->alignmentLength,
				m->strand,
				m->readLength,
				m->read)) {
		PrintError(FnName,
				NULL,
				"Could not write MAF",
				Exit,
				WriteFileError);
	}
	assert(strlen(m->read) == m->alignmentLength);
	assert(strlen(m->reference) == m->alignmentLength);
	assert(m->alignmentLength > 0);
	assert(m->referenceLength > 0);
	assert(m->readLength > 0);
}

/* Only works for bfast's maf files */
int MAFRead(FILE *fp,
		MAF *m)
{
	char *FnName="MAFRead";
	int tempIntOne, tempIntTwo;
	char line[MAX_LINE_LENGTH];
	char name[MAX_LINE_LENGTH];
	char tempC;
	fpos_t fpos;
	int newReport = 0;
	int aCounter=0;
	int sCounter=0;

	fgetpos(fp, &fpos); /* store current position in case we need to move back */ 
	/* Read in each line */
	while(NULL != fgets(line, MAX_LINE_LENGTH, fp) &&
			newReport == 0) {
		switch(line[0]) {
			case '#':
				break;
			case 'a':
				if(aCounter == 0) {
					if(0>sscanf(line, "a score=%lf",
								&m->score)) {
						PrintError(FnName,
								"a line",
								"Could not get score from a line",
								Exit,
								ReadFileError);
					}
				}
				aCounter++;
				break;
			case 's':
				if(sCounter==0) {
					if(6>sscanf(line, "s %u %u %d %c %d %s",
								&m->contig,
								&m->position,
								&m->alignmentLength,
								&m->strand,
								&m->referenceLength,
								m->reference)) {
						fprintf(stderr, "line=%s\n", line);
						PrintError(FnName,
								"s line",
								"Could not read in first s line",
								Exit,
								ReadFileError);
					}
				}
				else if(sCounter==1) {
					if(6>sscanf(line, "s %s %d %d %c %d %s",
								name,
								&tempIntOne,
								&tempIntTwo,
								&tempC,
								&m->readLength,
								m->read)) {
						PrintError(FnName,
								"s line",
								"Could not read in second s line",
								Exit,
								ReadFileError);
					}
				}
				sCounter++;
				break;
			case '\n':
				newReport = 1;
				break;
			default:
				fprintf(stderr, "[%c,%d]\n",
						line[0],
						(int)line[0]);
				PrintError(FnName,
						"line[0]",
						"Could not understand beginning character in the line",
						Exit,
						OutOfRange);
				break;
		}
		fgetpos(fp, &fpos); /* store current position in case we need to move back */ 
	}

	if(newReport == 0) {
		return EOF;
	}

	assert(m->alignmentLength > 0);
	assert(m->referenceLength > 0);
	assert(m->readLength > 0);
	return 1;
}

void MAFMergeSort(MAF *entries,
		int low,
		int high,
		int sortOrder)
{
	char *FnName="MAFMergeSort";
	int i, ctr;
	int mid = (low + high)/2;
	int startLower = low;
	int endLower = mid;
	int startUpper = mid+1;
	int endUpper = high;
	MAF *tmpEntries=NULL;

	if(low >= high) {
		return;
	}

	/* Partition the list */
	MAFMergeSort(entries,
			low,
			mid,
			sortOrder);
	MAFMergeSort(entries,
			mid+1,
			high,
			sortOrder);

	/* Allocate memory for tmp */
	tmpEntries = malloc(sizeof(MAF)*(high - low + 1));
	if(NULL == tmpEntries) {
		PrintError(FnName,
				"tmpEntries",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Merge the two lists */
	ctr=0;
	while( (startLower <= endLower) && (startUpper <= endUpper)) {
		if(MAFCompare(&entries[startLower], &entries[startUpper], sortOrder) <= 0) {
			MAFCopy(&entries[startLower], &tmpEntries[ctr]);
			startLower++;
		}
		else {
			MAFCopy(&entries[startUpper], &tmpEntries[ctr]);
			startUpper++;
		}
		ctr++;
	}
	while(startLower <= endLower) {
		MAFCopy(&entries[startLower], &tmpEntries[ctr]);
		startLower++;
		ctr++;
	}
	while(startUpper <= endUpper) {
		MAFCopy(&entries[startUpper], &tmpEntries[ctr]);
		startUpper++;
		ctr++;
	}
	/* Copy back */
	for(i=low, ctr=0;
			i<=high;
			i++, ctr++) {
		MAFCopy(&tmpEntries[ctr], &entries[i]);
	}

	/* Free memory */
	free(tmpEntries);
	tmpEntries = NULL;
}

int MAFCompare(MAF *m1, MAF *m2, int sortOrder) 
{
	char *FnName="MAFCompare";
	int cmp[2]={0,0};
	int top=1;
	int i;
	assert(m1 != NULL);
	assert(m2 != NULL);

	switch(sortOrder) {
		case MAFSortByContigPos:
			cmp[0] = (m1->contig <= m2->contig)?( (m1->contig == m2->contig)?0:-1):1;
			cmp[1] = (m1->position <= m2->position)?( (m1->position==m2->position)?0:-1):1;
			/* We want longer reference lengths first */
			top=2;
			break;
		default:
			PrintError(FnName,
					"sortOrder",
					"Could not understand sort order",
					Exit,
					OutOfRange);
	}
	for(i=0;i<top;i++) {
		if(cmp[i] < 0) {
			return -1;
		}
		else if(cmp[i] > 0) {
			return 1;
		}
	}
	return 0;
}

void MAFCopy(MAF *src, MAF *dest)
{
	dest->score = src->score;
	dest->contig = src->contig;
	dest->position = src->position;
	dest->strand = src->strand;
	dest->alignmentLength = src->alignmentLength;
	dest->referenceLength = src->referenceLength;
	dest->readLength = src->readLength;
	strcpy(dest->read, src->read);
	strcpy(dest->reference, src->reference);
}

void MAFInitialize(MAF *m) 
{
	m->score = 0.0;
	m->contig = 0;
	m->position = 0;
	m->strand = 0;
	m->alignmentLength = 0;
	m->referenceLength = 0;
	m->readLength = 0;
	strcpy(m->read, "\0");
	strcpy(m->reference, "\0");
}

void MAFPrintToBedAndWig(MAF *m,
		int numEntries,
		int contig,
		int64_t startPos,
		int64_t endPos,
		FILE *bedFP,
		FILE *wigFP)
{
	char *FnName="MAFPrintToWig";
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

		/* HERE */
		/*
		   fprintf(stderr, "Before\tcurPos=%lld\tm[%d].position=%d\tm[start].position + m[start].referenceLength - 1=%d\n",
		   (long long int)curPos,
		   start,
		   m[start].position,
		   m[start].position + m[start].referenceLength - 1);
		   */

		if(curPos < m[start].position) {
			curPos = m[start].position;
		}
		assert(m[start].position <= curPos &&
				curPos <= m[start].position + m[start].referenceLength -1);

		/* Go through every entry at that overlaps this position */
		for(cur = start;
				cur < numEntries &&
				m[cur].position <= curPos;
				cur++) {
			/* Only use if it is within bounds */ 
			if(m[cur].position <= curPos &&
					curPos <= m[cur].position + m[cur].referenceLength - 1) {
				assert(m[cur].contig == contig);
				/* Copy over reference and read.  Adjust if they are on the - strand */
				switch(m[cur].strand) {
					case FORWARD:
						strcpy(reference, m[cur].reference);
						strcpy(read, m[cur].read);
						break;
					case REVERSE:
						GetReverseComplimentAnyCase(m[cur].reference, reference, m[cur].alignmentLength);
						GetReverseComplimentAnyCase(m[cur].read, read, m[cur].alignmentLength);
						break;
					default:
						PrintError(FnName,
								"m[cur].strand",
								"Could not understand strand",
								Exit,
								OutOfRange);
				}
				/* Add to counts */
				/* Move to correct position */
				i=0;
				j=0;
				while(m[cur].position + i < curPos) {
					assert(j < m[cur].alignmentLength);
					/* Only move our position if there is not gap */
					if(reference[j] != GAP) {
						i++;
					}
					j++;
				}
				tempJ = j; /* Save this for bed */
				/* HERE */
				/*
				   fprintf(stderr, "HERE curPos=%lld m[%d].[start,end]=[%d,%d] i=%d j=%d\n",
				   (long long int)curPos,
				   cur,
				   m[cur].position,
				   m[cur].position + m[cur].referenceLength - 1,
				   i,
				   j);
				   */
				assert(m[cur].position + i == curPos);
				/************/
				/* WIG FILE */
				/************/
				/* Skip over gaps in the reference */
				while(j < m[cur].alignmentLength &&
						reference[j] == GAP) {
					j++;
				}
				/* We should not end with a gap so this should be true */
				assert(j < m[cur].alignmentLength);
				assert(reference[j] != GAP);
				/* Update based on the current base */
				if(referenceBase == 'n') {
					referenceBase = reference[j];
					assert(referenceBase != GAP);
				}
				/* Update counts */
				switch(m[cur].strand) {
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
								"m[cur].strand",
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
				if(j < m[cur].alignmentLength - 1 &&
						reference[j] != GAP &&
						read[j] != GAP) { 
					/* The indel starts after the current base */
					j=j+1;
					i=0;
					type='N';
					/* Check gap in reference */
					while(j<m[cur].alignmentLength && 
							reference[j] == GAP) {
						/* Get the insertion in the read */
						indel[i] = read[j];
						type = 'I';
						i++;
						j++;
					}
					/* Check gap in read */
					while(j<m[cur].alignmentLength && 
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
										m[cur].contig,
										(long long int)(curPos-SUBTRACT),
										(long long int)(curPos+i-1-SUBTRACT),
										m[cur].strand,
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
				m[start].position + m[start].referenceLength - 1 < curPos + 1) {
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
	MAF m;

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
	MAFInitialize(&m);
	while(EOF != MAFRead(fpIn, &m)) {
		if(counter%BMFCONVERT_ROTATE_NUM==0) {
			fprintf(stderr, "\r%lld",
					(long long int)counter);
		}
		counter++;
		/* Print to the appropriate file */
		if(!(m.contig > 0 && m.contig <= numFiles)) {
			MAFPrint(stderr, &m);
		}
		assert(m.contig > 0 && m.contig <= numFiles);
		MAFPrint((*tmpFiles)[m.contig-1].FP, &m);
		/* Update meta-data */
		(*tmpFiles)[m.contig-1].numEntries++;
		if(m.position < (*tmpFiles)[m.contig-1].minPos) {
			(*tmpFiles)[m.contig-1].minPos = m.position; 
		}
		if(m.position + m.referenceLength - 1 > (*tmpFiles)[m.contig-1].maxPos) {
			(*tmpFiles)[m.contig-1].maxPos = m.position + m.referenceLength - 1; 
		}
		MAFInitialize(&m);
	}
	fprintf(stderr, "\r%lld\n",
			(long long int)counter);

	/* Close the input file */
	fclose(fpIn);

	return numFiles;
}

void SplitMAFAndPrint(FILE *bedFP,
		FILE *wigFP,
		TmpFile *tmpFile, /* Should all be from the same contig */
		char *tmpDir,
		int maxNumEntries)
{
	char *FnName="SplitMAFAndPrint";
	int meanPos;
	MAF m;
	MAF *entries=NULL;
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
		/* Read in, sort, and print */
		numEntries = 0;
		while(EOF != MAFRead(tmpFile->FP, &m)) {
			numEntries++;
			/* Reallocate memory */
			entries = realloc(entries, sizeof(MAF)*numEntries);
			if(entries == NULL) {
				PrintError(FnName,
						"entries",
						"Could not allocate memory",
						Exit,
						MallocMemory);
			}
			/* Copy MAF */
			MAFCopy(&m, &entries[numEntries-1]);
		}
		assert(numEntries == tmpFile->numEntries);
		/* Sort */
		MAFMergeSort(entries, 0, numEntries-1, MAFSortByContigPos);
		/* Print Out */
		MAFPrintToBedAndWig(entries, numEntries, tmpFile->contig, tmpFile->minPos, tmpFile->maxPos, bedFP, wigFP);
		/* Free memory */
		free(entries);
		entries = NULL;
	}
	else {
		/* Split and recurse */

		/* Initialize */
		TmpFileOpen(&belowTmpFile, tmpDir, tmpFile->contig);
		TmpFileOpen(&aboveTmpFile, tmpDir, tmpFile->contig);

		/* Where will we split */
		meanPos = (tmpFile->maxPos + tmpFile->minPos)/2;

		/* Split */
		MAFInitialize(&m);
		while(EOF != MAFRead(tmpFile->FP, &m)) {
			/* Print to the appropriate file */
			assert(m.contig == tmpFile->contig);

			/* Update meta-data */
			belowTmpFile.minPos = tmpFile->minPos; 
			belowTmpFile.maxPos = meanPos - 1; 
			aboveTmpFile.minPos = meanPos; 
			aboveTmpFile.maxPos = tmpFile->maxPos; 

			/* Print */
			if(m.position < meanPos) {
				MAFPrint(belowTmpFile.FP, &m);
				belowTmpFile.numEntries++;
			}
			if(meanPos <= m.position ||
					meanPos <= m.position + m.referenceLength - 1 ) {
				MAFPrint(aboveTmpFile.FP, &m);
				aboveTmpFile.numEntries++;
			}
		}

		/* Recurse on the two */
		SplitMAFAndPrint(bedFP,
				wigFP,
				&belowTmpFile,
				tmpDir,
				maxNumEntries);
		SplitMAFAndPrint(bedFP,
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

		/* Split MAF by contig */
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
				SplitMAFAndPrint(bedFP,
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
		fprintf(stderr, "%s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast maf file name>\n");
		fprintf(stderr, "\t<maximum number of entries>\n");
		fprintf(stderr, "\t<tmp file directory>\n");
	}
	return 0;
}
