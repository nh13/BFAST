#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "../blib/AlignEntry.h"
#include "Definitions.h"
#include "FilterAlignments.h"
#include "PrintOutputFiles.h"

/* TODO */
int PrintAlignEntriesToTempFilesByChr(FILE* fp,
		RGFiles *rgFiles,
		int uniqueMatches,
		int bestScore,
		int minScore,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		char *tmpDir)
{
	int i, numEntries, chrIndex, numRead, numFiltered; 
	char *FnName = "PrintAlignEntriesToTempFilesByChr";
	AlignEntry *entries=NULL;

	/* Initialize */
	rgFiles->startChr = startChr;
	rgFiles->endChr =  endChr;
	rgFiles->chrFiles = NULL;
	rgFiles->chrFileNames = NULL;

	/* Allocate memory */
	rgFiles->chrFiles = malloc(sizeof(FILE*)*(endChr - startChr + 1));
	if(NULL == rgFiles->chrFiles) {
		PrintError(FnName,
				"rgFiles->chrFiles",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	rgFiles->chrFileNames = malloc(sizeof(char*)*(endChr - startChr + 1));
	if(NULL == rgFiles->chrFileNames) {
		PrintError(FnName,
				"rgFiles->chrFileNames",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Open tmp files */
	for(i=0;i<rgFiles->endChr - rgFiles->startChr + 1;i++) {
		rgFiles->chrFiles[i] = OpenTmpFile(tmpDir, &rgFiles->chrFileNames[i]);
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "Binning into temp files by chromosome and filtering reads.  Currently on [read, filtered]:\n0");
	}

	/* Read in each read and output to the correct tmp file */
	for(numEntries = AlignEntryGetOneRead(&entries, fp), numRead=0, numFiltered=0;
			numEntries > 0;
			numEntries = AlignEntryGetOneRead(&entries, fp)) {
		numRead++;

		/* Apply filtering */
		numEntries = FilterEntries(&entries,
				numEntries,
				uniqueMatches,
				bestScore,
				minScore,
				startChr,
				startPos,
				endChr,
				endPos);
		assert(numEntries == 1 || numEntries == 0);

		/* Output if there is only one entry left after filtering */
		if(numEntries==1) {

			/* Get the correct file to which to print for the start of the read */
			assert(entries[0].chromosome >= startChr && entries[0].chromosome <= endChr);
			chrIndex = entries[0].chromosome - startChr;

			/* Print to the appropriate file */
			AlignEntryPrint(&entries[0], rgFiles->chrFiles[chrIndex]);

			/* Free entries */
			AlignEntryFree(&entries[0]);
			free(entries);
			entries=NULL;
		}
		else {
			numFiltered++;
			/* We should have zero entries */
			assert(numEntries==0);
		}
		if(VERBOSE >= 0 && numRead%ROTATE_BINNING == 0) {
			fprintf(stderr, "\r[%d,%d]", numRead, numFiltered);
		}
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\r[%d,%d]\nBinning complete.\n", numRead, numFiltered);
	}

	return numRead-numFiltered;
}

/* TODO */
/* Does not filter */
int PrintAlignEntriesToTempFilesWithinChr(FILE *fp,
		int uniqueMatches,
		int bestScore,
		int minScore,
		int curChr,
		int startPos,
		int endPos,
		int regionLength,
		char *tmpDir,
		ChrFiles *chrFiles)
{
	int i, numEntries, regionIndex, startIndex, numRead;
	char *FnName = "PrintAlignEntriesToTempFilesWithinChr";
	AlignEntry *entries=NULL;

	/* Initialize */
	chrFiles->files = NULL;
	chrFiles->fileNames = NULL;
	chrFiles->numFiles = 0;

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Binning into temp files for chromosome %d.  Currently on:\n0", 
				curChr);
	}

	/* Read in each read and output to the correct tmp file */
	numRead = 0;
	for(numEntries = AlignEntryGetOneRead(&entries, fp);
			numEntries > 0;
			numEntries = AlignEntryGetOneRead(&entries, fp)) {
		numRead++;

		if(VERBOSE >= 0 && numRead%ROTATE_BINNING == 0) {
			fprintf(stderr, "\r[%d]", numRead);
		}

		/* Output if there is only one entry left after filtering */
		assert(numEntries==1);

		/* Get the correct file to which to print for the start of the read */
		assert(entries[0].chromosome == curChr);
		regionIndex = (entries[0].position)/regionLength;

		/* Reallocate files if the region file does not exist*/
		if(chrFiles->numFiles < regionIndex + 1) {
			/* Reallocate memory */
			startIndex = chrFiles->numFiles;
			chrFiles->numFiles = regionIndex + 1;
			chrFiles->files = realloc(chrFiles->files, sizeof(FILE*)*(chrFiles->numFiles));
			if(NULL==chrFiles->files) {
				PrintError(FnName,
						"chrFiles->files",
						"Could not reallocate memory",
						Exit,
						ReallocMemory);
			}
			chrFiles->fileNames = realloc(chrFiles->fileNames, sizeof(char*)*(chrFiles->numFiles));
			if(NULL==chrFiles->fileNames) {
				PrintError(FnName,
						"chrFiles->fileNames",
						"Could not reallocate memory",
						Exit,
						ReallocMemory);
			}
			/* Open new temp files */
			for(i=startIndex;i<chrFiles->numFiles;i++) {
				chrFiles->files[i] = OpenTmpFile(tmpDir, &chrFiles->fileNames[i]);
				if(chrFiles->files[i] == NULL) {
					fprintf(stderr, "i:%d\n", i);
					PrintError(FnName,
							"chrFiles->files[i]",
							"Could not open file",
							Exit,
							OpenFileError);
				}
			}
		}
		/* Output the entry to the tmp file*/
		assert(regionIndex < chrFiles->numFiles);
		assert(chrFiles->files[regionIndex] != NULL);
		AlignEntryPrint(&entries[0], chrFiles->files[regionIndex]);

		/* If the read overlaps regions, then we output it to both regions */
		if(regionIndex +1 == (entries[0].position + entries[0].referenceLength-1)/regionLength) {
			/* Get the correct file to which to print */
			regionIndex++;

			/* Reallocate files if the region file does not exist*/
			if(chrFiles->numFiles < regionIndex + 1) {
				/* Reallocate memory */
				startIndex = chrFiles->numFiles;
				chrFiles->numFiles = regionIndex + 1;
				chrFiles->files = realloc(chrFiles->files, sizeof(FILE*)*(chrFiles->numFiles));
				if(NULL==chrFiles->files) {
					PrintError(FnName,
							"chrFiles->files",
							"Could not reallocate memory",
							Exit,
							ReallocMemory);
				}
				chrFiles->fileNames = realloc(chrFiles->fileNames, sizeof(char*)*(chrFiles->numFiles));
				if(NULL==chrFiles->fileNames) {
					PrintError(FnName,
							"chrFiles->fileNames",
							"Could not reallocate memory",
							Exit,
							ReallocMemory);
				}
				/* Open new temp files */
				for(i=startIndex;i<chrFiles->numFiles;i++) {
					chrFiles->files[i] = OpenTmpFile(tmpDir, &chrFiles->fileNames[i]);
					if(chrFiles->files[i] == NULL) {
						fprintf(stderr, "i:%d\n", i);
						PrintError(FnName,
								"chrFiles->files[i]",
								"Could not open file",
								Exit,
								OpenFileError);
					}
				}
			}
			/* Output the entry to the tmp file*/
			assert(regionIndex < chrFiles->numFiles);
			assert(regionIndex < chrFiles->numFiles);
			assert(chrFiles->files[regionIndex] != NULL);
			AlignEntryPrint(&entries[0], chrFiles->files[regionIndex]);
		}
		else {
			assert(regionIndex == (entries[0].position + entries[0].referenceLength-1)/regionLength);
		}
		/* Free entries */
		AlignEntryFree(&entries[0]);
		free(entries);
		entries=NULL;
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\r[%d]\nBinning complete.\n", numRead);
	}

	return numRead;
}

/* TODO */
void PrintAlignEntries(ChrFiles *chrFile,
		int regionLength,
		int outputFormat,
		char *outputDir,
		char *outputID,
		int curChr)
{
	char *FnName = "PrintAlignEntries";
	int i, j, numEntries;
	double curPercent;
	char **outputFileNames=NULL;
	AlignEntry *entries=NULL;
	FILE **outputFPs=NULL;
	int numOutputFPs=-1;
	int numReads=0;

	/* Get the number of output files */
	switch(outputFormat) {
		case BAlignFile:
			PrintError(FnName,
					"outputFormat",
					"Output format not supported",
					Exit,
					OutOfRange);
			break;
		case WigFile:
			/* Create only one output file */
			numOutputFPs=1;
			break;
		case BedFile:
			numOutputFPs=3; /* One for mismatches, one for insertions, one for deletions */
			break;
		case BedAndWigFile:
			numOutputFPs=4;
			break;
		default:
			PrintError(FnName,
					"outputFormat",
					"Output format not supported",
					Exit,
					OutOfRange);
			break;
	}

	/* Allocate memory */
	outputFileNames = malloc(sizeof(char*)*numOutputFPs);
	if(NULL==outputFileNames) {
		PrintError(FnName,
				"outputFileNames",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	outputFPs = malloc(sizeof(FILE*)*numOutputFPs);
	if(NULL==outputFPs) {
		PrintError(FnName,
				"outputFPs",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	for(i=0;i<numOutputFPs;i++) {
		outputFileNames[i] = malloc(sizeof(char)*MAX_FILENAME_LENGTH);
		if(NULL==outputFileNames[i]) {
			PrintError(FnName,
					"outputFileNames[i]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
	}

	/* Create the output file names */
	switch(outputFormat) {
		case BAlignFile:
			PrintError(FnName,
					"outputFormat",
					"Output format not supported",
					Exit,
					OutOfRange);
			break;
		case WigFile:
			/* Create only one output file */
			assert(numOutputFPs==1);
			/* Create output file name */
			sprintf(outputFileNames[0], "%s%s.%s.%d.%s",
					outputDir,
					PROGRAM_NAME,
					outputID,
					curChr,
					"wig");
			break;
		case BedFile:
			assert(numOutputFPs==3);
			sprintf(outputFileNames[0], "%s%s.%s.%s.%d.%s",
					outputDir,
					PROGRAM_NAME,
					outputID,
					"mismatches",
					curChr,
					"bed");
			sprintf(outputFileNames[1], "%s%s.%s.%s.%d.%s",
					outputDir,
					PROGRAM_NAME,
					outputID,
					"insertions",
					curChr,
					"bed");
			sprintf(outputFileNames[2], "%s%s.%s.%s.%d.%s",
					outputDir,
					PROGRAM_NAME,
					outputID,
					"deletions",
					curChr,
					"bed");
			break;
		case BedAndWigFile:
			/* Create only one output file */
			assert(numOutputFPs==4);
			/* Create output file name */
			sprintf(outputFileNames[0], "%s%s.%s.%d.%s",
					outputDir,
					PROGRAM_NAME,
					outputID,
					curChr,
					"wig");
			sprintf(outputFileNames[1], "%s%s.%s.%s.%d.%s",
					outputDir,
					PROGRAM_NAME,
					outputID,
					"mismatches",
					curChr,
					"bed");
			sprintf(outputFileNames[2], "%s%s.%s.%s.%d.%s",
					outputDir,
					PROGRAM_NAME,
					outputID,
					"insertions",
					curChr,
					"bed");
			sprintf(outputFileNames[3], "%s%s.%s.%s.%d.%s",
					outputDir,
					PROGRAM_NAME,
					outputID,
					"deletions",
					curChr,
					"bed");
			break;
		default:
			PrintError(FnName,
					"outputFormat",
					"Output format not supported",
					Exit,
					OutOfRange);
			break;
	}

	for(i=0;i<numOutputFPs;i++) {
		if(0==(outputFPs[i]=fopen(outputFileNames[i], "w"))) {
			PrintError(FnName,
					outputFileNames[i],
					"Could not open output file for writing",
					Exit,
					OpenFileError);
		}
	}

	/* Go through each region */
	for(i=0;i<chrFile->numFiles;i++) {
		/* Get all the entries from the tmp file */
		/*
		   if(VERBOSE >= 0) {
		   fprintf(stderr, "\rGetting entries for [chr,startPos-endPos]=[%d,%d-%d]...",
		   curChr,
		   i*regionLength+1,
		   (i+1)*regionLength);
		   }
		   */

		numEntries=AlignEntryGetAll(&entries, chrFile->files[i]);

		numReads += numEntries;

		/* Close the file */
		CloseTmpFile(&chrFile->files[i], &chrFile->fileNames[i]);

		/* Output the entries */
		if(numEntries > 0) {
			/* Sort the entries */
			/*
			   if(VERBOSE >= 0) {
			   fprintf(stderr, "\rSorting %d entries for [chr,startPos-endPos]=[%d,%d-%d]...\n%3.2lf",
			   numEntries,
			   curChr,
			   i*regionLength+1,
			   (i+1)*regionLength,
			   0.0);
			   }
			   */
			curPercent = 0.0;
			if(VERBOSE >= 0) {
				fprintf(stderr, "\rSorting and outputting %d entries for [chr,startPos-endPos]=[%d,%d-%d]:    ",
						numEntries,
						curChr,
						i*regionLength+1,
						(i+1)*regionLength);
				fprintf(stderr, "%3.2lf percent complete", 0.0);
			}
			AlignEntryMergeSort(&entries, 
					0, 
					numEntries-1, 
					AlignEntrySortByChrPos,
					1,
					&curPercent,
					numEntries);
			/*
			   AlignEntryQuickSort(&entries, 
			   0, 
			   numEntries-1, 
			   AlignEntrySortByChrPos,
			   1,
			   &curPercent,
			   numEntries);
			   */
			if(VERBOSE >= 0) {
				PrintPercentCompleteShort(100.0);
			}

			switch(outputFormat) {
				case BAlignFile:
					PrintError(FnName,
							"outputFormat",
							"Output format not supported",
							Exit,
							OutOfRange);
					break;
				case WigFile:
					/* Create only one output file */
					assert(numOutputFPs==1);
					/* Print the entries to file */
					PrintSortedAlignEntriesToWig(entries,
							numEntries,
							curChr,
							outputFPs[0]);
					break;
				case BedFile:
					assert(numOutputFPs==3);
					PrintSortedAlignEntriesToBed(entries,
							numEntries,
							curChr,
							outputFPs,
							0,
							numOutputFPs);
					break;
				case BedAndWigFile:
					assert(numOutputFPs==4);
					/* Print the entries to file */
					PrintSortedAlignEntriesToWig(entries,
							numEntries,
							curChr,
							outputFPs[0]);
					PrintSortedAlignEntriesToBed(entries,
							numEntries,
							curChr,
							outputFPs,
							1,
							numOutputFPs);
					break;
				default:
					PrintError(FnName,
							"outputFormat",
							"Output format not supported",
							Exit,
							OutOfRange);
					break;
			}

			/* Free the entries */
			for(j=0;j<numEntries;j++) {
				AlignEntryFree(&entries[j]);
			}
			free(entries);
			entries=NULL;
			/*
			   if(VERBOSE >= 0) {
			   fprintf(stderr, "\rOutputted %d entries for [chr,startPos-endPos]=[%d,%d-%d].\n",
			   numEntries,
			   curChr,
			   i*regionLength+1,
			   (i+1)*regionLength);
			   }
			   */
		}
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\rSorted and outputted %d entries for [chr,startPos-endPos]=[%d,%d-%d]:    %3.2lf percent complete",
				numReads,
				curChr,
				1,
				(chrFile->numFiles)*regionLength,
				100.00);
		/* Just for good measure */
		fprintf(stderr, "                  \n");
		fprintf(stderr, "Outputted %d entries.\n", numReads);
	}

	/* Close the output files */
	for(i=0;i<numOutputFPs;i++) {
		fclose(outputFPs[i]);
	}

	/* Free memory */
	for(i=0;i<numOutputFPs;i++) {
		free(outputFileNames[i]);
	}
	free(outputFileNames);
	free(outputFPs);
}

/* TODO */
void PrintSortedAlignEntriesToWig(AlignEntry *entries,
		int numEntries,
		int curChr,
		FILE *outputFP)
{
	char *FnName = "PrintSortedAlignEntriesToWig";
	char tmpRead[SEQUENCE_LENGTH]="\0";
	char tmpReference[SEQUENCE_LENGTH]="\0";
	int i, j;
	int curPos, coverage;
	int startIndex, curIndex;
	char referenceBase = '\0';
	assert(numEntries > 0);

	/* Go through each position */
	curPos = entries[0].position;
	startIndex = 0;
	while(startIndex < numEntries) {

		/* Get the coverage for that position */
		for(curIndex=startIndex, coverage=0;
				/* Do not go past the number of entries */
				curIndex < numEntries && 
				/* The start position of the entry must be within bounds */
				entries[curIndex].position <= curPos; 
				curIndex++) {

			/* Check to see that the end position is within bounds */
			if(curPos <= entries[curIndex].position + entries[curIndex].referenceLength - 1) {

				/* Alignment and coverage depend on strandedness */
				switch(entries[curIndex].strand) {
					case FORWARD:
						strcpy(tmpRead, entries[curIndex].read);
						strcpy(tmpReference, entries[curIndex].reference);
						break;
					case REVERSE:
						/*
						   assert(strlen(entries[curIndex].reference) == entries[curIndex].length);
						   assert(strlen(entries[curIndex].read) == entries[curIndex].length);
						   assert(strlen(entries[curIndex].reference) == entries[curIndex].length);
						   */
						GetReverseComplimentAnyCase(entries[curIndex].read, tmpRead, entries[curIndex].length);
						GetReverseComplimentAnyCase(entries[curIndex].reference, tmpReference, entries[curIndex].length);
						break;
					default:
						fprintf(stderr, "[%c]\t[%d]\n",
								(int)entries[curIndex].strand,
								entries[curIndex].strand);
						PrintError(FnName,
								"entries[curIndex].strand",
								"Could not understand strand",
								Exit,
								OutOfRange);
						break;
				}
				/* Convert all to upper case */
				for(i=0;i<entries[curIndex].length;i++) {
					tmpRead[i] = ToUpper(tmpRead[i]);
					tmpReference[i] = ToUpper(tmpReference[i]);
				}
				/* Now everthing should be relative to the + strand */

				/* Go to the base in the reference for this entry */
				i = entries[curIndex].position;
				j = 0;
				while(i<curPos && j<entries[curIndex].length) {
					/* Update position if it is not a gap */
					if(tmpReference[j] != GAP) {
						i++;
					}
					j++;
				}
				while(j<entries[curIndex].length && tmpReference[j] == GAP) {
					j++;
				}
				assert(j<entries[curIndex].length); /* by definition that we went into the for loop */
				assert(i==curPos);
				/* Check that the reference sequence matches */
				if(coverage==0) {
					referenceBase = tmpReference[j];
					assert(referenceBase == 'A' ||
							referenceBase == 'C' || 
							referenceBase == 'G' ||
							referenceBase == 'T');
				}
				else {
					if(!(referenceBase == tmpReference[j])) {
						fprintf(stderr, "%c\t%c\n",
								referenceBase,
								tmpReference[j]);
						fprintf(stderr, "%s\n%s\n",
								tmpReference,
								tmpRead);
						AlignEntryPrint(&entries[curIndex], stderr);
					}
					assert(referenceBase == tmpReference[j]);
				}
				coverage++;
			}
		}

		/* Print the coverage if it is > 0 */
		if(coverage > 0) {
			fprintf(outputFP, "chr%d\t%d\t%d\t%d\n",
					curChr,
					curPos-1,
					curPos,
					coverage);
		}

		/* update position */
		curPos++;
		/* Update the start index */
		while(startIndex < numEntries &&
				entries[startIndex].position + entries[startIndex].referenceLength - 1 < curPos) {
			startIndex++;
		}
	}
}


/* TODO */
void PrintSortedAlignEntriesToBed(AlignEntry *entries,
		int numEntries,
		int curChr,
		FILE **outputFPs,
		int startOutputFPs,
		int numOutputFPs)
{
	char *FnName="PrintSortedAlignEntriesToBed";
	int i, j;
	int curPos, coverage;
	int startIndex, curIndex;
	int numF[5]={0,0,0,0,0}; /* A, C, G, T, Deletion */
	int numR[5]={0,0,0,0,0}; /* A, C, G, T, Deletion */
	char tmpRead[SEQUENCE_LENGTH]="\0";
	char tmpReference[SEQUENCE_LENGTH]="\0";
	char referenceBase = '\0';
	assert(numEntries > 0);
	assert(numOutputFPs - startOutputFPs == 3);

	/* Go through each position */
	curPos = entries[0].position;
	startIndex = 0;
	while(startIndex < numEntries) {
		/* Initialize */
		numF[0]=numF[1]=numF[2]=numF[3]=numF[4]=0;
		numR[0]=numR[1]=numR[2]=numR[3]=numR[4]=0;
		referenceBase = '\0';
		coverage = 0;
		/* Go through each entry within bounds */
		for(curIndex=startIndex, coverage=0;
				curIndex < numEntries &&
				entries[curIndex].position <= curPos;
				curIndex++) {
			/* Check to see that the end position is within bounds */
			if(curPos <= entries[curIndex].position + entries[curIndex].referenceLength - 1) {
				/* Alignment and coverage depend on strandedness */
				switch(entries[curIndex].strand) {
					case FORWARD:
						strcpy(tmpRead, entries[curIndex].read);
						strcpy(tmpReference, entries[curIndex].reference);
						break;
					case REVERSE:
						assert(strlen(entries[curIndex].read) == entries[curIndex].length);
						assert(strlen(entries[curIndex].reference) == entries[curIndex].length);
						GetReverseComplimentAnyCase(entries[curIndex].read, tmpRead, entries[curIndex].length);
						GetReverseComplimentAnyCase(entries[curIndex].reference, tmpReference, entries[curIndex].length);
						break;
					default:
						PrintError(FnName,
								"entries[curIndex].strand",
								"Could not understand strand",
								Exit,
								OutOfRange);
						break;
				}
				/* Convert all to upper case */
				for(i=0;i<entries[curIndex].length;i++) {
					tmpRead[i] = ToUpper(tmpRead[i]);
					tmpReference[i] = ToUpper(tmpReference[i]);
				}
				/* Now everthing should be relative to the + strand */

				/* Go the base in the reference */
				i = entries[curIndex].position;
				j = 0;
				while(i<curPos && j<entries[curIndex].length) {
					if(tmpReference[j] == GAP) {
						/* Skip over */
						j++;
					}
					else {
						i++;
						j++;
					}
				}
				while(j<entries[curIndex].length && tmpReference[j] == GAP) {
					j++;
				}
				assert(j<entries[curIndex].length); /* by definition that we are in the for loop */
				assert(i==curPos);
				/* Check that the reference sequence matches */
				if(coverage==0) {
					referenceBase = tmpReference[j];
					assert(referenceBase == 'A' ||
							referenceBase == 'C' || 
							referenceBase == 'G' ||
							referenceBase == 'T');
				}
				else {
					assert(referenceBase == tmpReference[j]);
				}
				/* Get mismatches and deletions */
				switch(tmpRead[j]) {
					case 'A':
						(entries[curIndex].strand==FORWARD)?(numF[0]++):(numR[0]++);
						break;
					case 'C':
						(entries[curIndex].strand==FORWARD)?(numF[1]++):(numR[1]++);
						break;
					case 'G':
						(entries[curIndex].strand==FORWARD)?(numF[2]++):(numR[2]++);
						break;
					case 'T':
						(entries[curIndex].strand==FORWARD)?(numF[3]++):(numR[3]++);
						break;
					case GAP:
						(entries[curIndex].strand==FORWARD)?(numF[4]++):(numR[4]++);
						break;
					default:
						PrintError(FnName,
								entries[curIndex].read,
								"Illegal base in read",
								Exit,
								OutOfRange);
						break;
				}
				/* Get insertions */
				if(tmpReference[j] !=GAP && /* No gap at the current position */
						j<entries[curIndex].referenceLength-1 && /* There are more bases in the alignment */
						tmpReference[j+1] == GAP) { /* The next base is a gap */
					/* We started a gap */
					fprintf(outputFPs[startOutputFPs+1], "chr%d\t%d\t%d\t",
							curChr,
							curPos-1,
							curPos);
					/* Print the insertion */
					for(i=j+1;
							i<entries[curIndex].referenceLength &&
							tmpReference[i] == GAP;
							i++) {
						fprintf(outputFPs[startOutputFPs+1], "%c",
								entries[curIndex].read[i]);
					}
					fprintf(outputFPs[startOutputFPs+1], "\n");
				}
				coverage++;
			}
		}

		/* Print the coverage if it is > 0 */
		if(coverage > 0) {
			int i = -1;
			switch(referenceBase) {
				case 'A':
					i = 0;
					break;
				case 'C':
					i = 1;
					break;
				case 'G':
					i = 2;
					break;
				case 'T':
					i = 3;
					break;
				default:
					fprintf(stderr, "\n[%c]\n[%d]\n",
							referenceBase,
							(int)referenceBase);
					PrintError(FnName,
							"referenceBase ",
							"Could not understand reference base",
							Exit,
							OutOfRange);
			}
			/* Print the mismatches */
			if(numF[i] + numR[i] != numF[0]+numF[1]+numF[2]+numF[3]+numR[0]+numR[1]+numR[2]+numR[3]) {
				fprintf(outputFPs[startOutputFPs], "chr%d\t%d\t%d\t%c(%d:%d:%d:%d:%3.2lf[F:%d:%d:%d:%d:%3.2lf][R:%d:%d:%d:%d:%3.2lf])\n",
						curChr,
						curPos-1,
						curPos,
						referenceBase,
						numF[0]+numR[0],
						numF[1]+numR[1],
						numF[2]+numR[2],
						numF[3]+numR[3],
						(100.0*(numF[i]+numR[i]))/(numF[0]+numF[1]+numF[2]+numF[3]+numR[0]+numR[1]+numR[2]+numR[3]),
						numF[0],
						numF[1],
						numF[2],
						numF[3],
						(100.0*(numF[i]))/(numF[0]+numF[1]+numF[2]+numF[3]),
						numR[0],
						numR[1],
						numR[2],
						numR[3],
						(100.0*(numR[i]))/(numR[0]+numR[1]+numR[2]+numR[3]));
			}
			/* Print deletions */
			if(numF[4] > 0 || numF[4] > 0) {
				fprintf(outputFPs[startOutputFPs+2], "chr%d\t%d\t%d\t%c(%d:%d:%d:%d:%d:%3.2lf[F:%d:%d:%d:%d:%d:%3.2lf][R:%d:%d:%d:%d:%d:%3.2lf])\n",
						curChr,
						curPos-1,
						curPos,
						referenceBase,
						numF[0]+numR[0],
						numF[1]+numR[1],
						numF[2]+numR[2],
						numF[3]+numR[3],
						numF[4]+numF[4],
						(100.0*(numF[4]+numR[4]))/(numF[0]+numF[1]+numF[2]+numF[3]+numF[4]+numR[0]+numR[1]+numR[2]+numR[3]+numR[4]),
						numF[0],
						numF[1],
						numF[2],
						numF[3],
						numF[4],
						(100.0*(numF[4]))/(numF[0]+numF[1]+numF[2]+numF[3]+numR[4]),
						numR[0],
						numR[1],
						numR[2],
						numR[3],
						numR[4],
						(100.0*(numR[4]))/(numR[0]+numR[1]+numR[2]+numR[3]+numR[4]));
			}
		}

		/* update position */
		curPos++;
		/* Update the start index */
		while(startIndex < numEntries &&
				entries[startIndex].position + entries[startIndex].referenceLength - 1 < curPos) {
			startIndex++;
		}
	}
}
















