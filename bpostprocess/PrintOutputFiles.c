#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "../blib/AlignEntry.h"
#include "Definitions.h"
#include "FilterAlignments.h"
#include "PrintOutputFiles.h"

/* TODO */
void PrintAlignEntriesToTempFiles(FILE *fp,
		int uniqueMatches,
		int bestScore,
		int minScore,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		ChrFiles *chrFiles)
{
	int i, numEntries, chrIndex, regionIndex, startIndex, curRead;
	char *FnName = "PrintAlignEntriesToTempFiles";
	AlignEntry *entries=NULL;

	/* Initialize */
	chrFiles->files = NULL;
	chrFiles->numFiles = 0;

	if(VERBOSE >= 0) {
		fprintf(stderr, "Binning into temp files.  Currently on:\n0");
	}

	/* Read in each read and output to the correct tmp file */
	curRead = 0;
	for(numEntries = AlignEntryGetOneRead(&entries, fp);
			numEntries > 0;
			numEntries = AlignEntryGetOneRead(&entries, fp)) {
		curRead++;

		if(VERBOSE >= 0 && curRead%ROTATE_BINNING == 0) {
			fprintf(stderr, "\r%d", curRead);
		}

		/* Apply filtering */
			fprintf(stderr, "HERE 1\t%d\n", curRead);
		numEntries = FilterEntries(&entries,
				numEntries,
				uniqueMatches,
				bestScore,
				minScore,
				startChr,
				startPos,
				endChr,
				endPos);
			fprintf(stderr, "HERE 2\t%d\n", curRead);

		/* Output if there is only one entry left after filtering */
		if(numEntries==1) {
			/* Get the correct file to which to print for the start of the read */
			chrIndex = entries[0].chromosome - startChr;
			assert(chrIndex==0);
			regionIndex = (entries[0].position)/REGION_LENGTH;

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
				/* Open new temp files */
				for(i=startIndex;i<chrFiles->numFiles;i++) {
					chrFiles->files[i] = tmpfile();
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
			if(regionIndex +1 == (entries[0].position + entries[0].referenceLength-1)/REGION_LENGTH) {
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
					/* Open new temp files */
					for(i=startIndex;i<chrFiles->numFiles;i++) {
						chrFiles->files[i] = tmpfile();
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
				assert(regionIndex == (entries[0].position + entries[0].referenceLength-1)/REGION_LENGTH);
			}
			/* Free entries */
			free(entries[0].read);
			free(entries[0].reference);
			free(entries);
			entries=NULL;
		}
		else {
			/* We should have zero entries */
			assert(numEntries==0);
		}
	}
	if(VERBOSE >= 0) {
		fprintf(stderr, "\r%d\nBinning complete.\n", curRead);
	}
}

/* TODO */
void PrintAlignEntries(ChrFiles *chrFile,
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
			sprintf(outputFileNames[0], "%sblatter.%s.%d.%s",
					outputDir,
					outputID,
					curChr,
					"wig");
			break;
		case BedFile:
			assert(numOutputFPs==3);
			sprintf(outputFileNames[0], "%sblatter.%s.%s.%d.%s",
					outputDir,
					outputID,
					"mismatches",
					curChr,
					"bed");
			sprintf(outputFileNames[1], "%sblatter.%s.%s.%d.%s",
					outputDir,
					outputID,
					"insertions",
					curChr,
					"bed");
			sprintf(outputFileNames[2], "%sblatter.%s.%s.%d.%s",
					outputDir,
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
		if(VERBOSE >= 0) {
			fprintf(stderr, "\rGetting entries for [chr,startPos-endPos[%d,%d-%d]...",
					curChr,
					i*REGION_LENGTH+1,
					(i+1)*REGION_LENGTH);
		}
		numEntries=AlignEntryGetAll(&entries, chrFile->files[i]);

		/* Close the file */
		fclose(chrFile->files[i]);
		chrFile->files[i] = NULL;

		/* Output the entries */
		if(numEntries > 0) {
			/* Sort the entries */
			if(VERBOSE >= 0) {
				fprintf(stderr, "\rSorting entries for [chr,startPos-endPos[%d,%d-%d]...\n%3.2lf",
						curChr,
						i*REGION_LENGTH+1,
						(i+1)*REGION_LENGTH,
						0.0);
			}
			curPercent = 0.0;
			AlignEntryQuickSort(&entries, 
					0, 
					numEntries-1, 
					AlignEntrySortByChrPos,
					1,
					&curPercent,
					numEntries);

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
				free(entries[j].read);
				free(entries[j].reference);
			}
			if(VERBOSE >= 0) {
				fprintf(stderr, "\rOutputted entries for [chr,startPos-endPos[%d,%d-%d].\n",
						curChr,
						i*REGION_LENGTH+1,
						(i+1)*REGION_LENGTH);
			}
		}
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
	int i, j;
	int curPos, coverage;
	int startIndex, curIndex;
	char reference = '\0';
	assert(numEntries > 0);

	/* Go through each position */
	curPos = entries[0].position;
	startIndex = 0;
	while(startIndex < numEntries) {
		/* Get the coverage */
		for(curIndex=startIndex, coverage=0;
				curIndex < numEntries &&
				entries[curIndex].position <= curPos &&
				curPos <= entries[curIndex].position + entries[curIndex].referenceLength - 1;
				curIndex++, coverage++) {
			/* Go the base in the reference */
			i = entries[curIndex].position;
			j = 0;
			while(i<curPos && j<entries[curIndex].length) {
				/* Update position if it is not a gap */
				if(entries[curIndex].reference[j] != GAP) {
					i++;
				}
				j++;
			}
			while(j<entries[curIndex].length && entries[curIndex].reference[j] == GAP) {
				j++;
			}
			assert(j<entries[curIndex].length); /* by definition that we are in the for loop */
			assert(i==curPos);
			/* Check that the reference sequence matches */
			if(coverage==0) {
				reference = ToUpper(entries[curIndex].reference[j]);
			}
			else {
				assert(reference == ToUpper(entries[curIndex].reference[j]));
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
		int numOutputFPs)
{
	char *FnName="PrintSortedAlignEntriesToBed";
	int i, j;
	int curPos, coverage;
	int startIndex, curIndex;
	int numF[5]={0,0,0,0,0}; /* A, C, G, T, Deletion */
	int numR[5]={0,0,0,0,0}; /* A, C, G, T, Deletion */
	char reference = '\0';
	assert(numEntries > 0);

	/* Go through each position */
	curPos = entries[0].position;
	startIndex = 0;
	while(startIndex < numEntries) {
		/* Initialize */
		numF[0]=numF[1]=numF[2]=numF[3]=numF[4]=0;
		numR[0]=numR[1]=numR[2]=numR[3]=numR[4]=0;
		reference = '\0';
		coverage = 0;
		/* Go through each entry within bounds */
		for(curIndex=startIndex, coverage=0;
				curIndex < numEntries &&
				entries[curIndex].position <= curPos &&
				curPos <= entries[curIndex].position + entries[curIndex].referenceLength - 1;
				curIndex++, coverage++) {
			/* Go the base in the reference */
			i = entries[curIndex].position;
			j = 0;
			while(i<curPos && j<entries[curIndex].length) {
				if(entries[curIndex].reference[j] == GAP) {
					/* Skip over */
					j++;
				}
				else {
					i++;
					j++;
				}
			}
			while(j<entries[curIndex].length && entries[curIndex].reference[j] == GAP) {
				j++;
			}
			assert(j<entries[curIndex].length); /* by definition that we are in the for loop */
			assert(i==curPos);
			/* Check that the reference sequence matches */
			if(coverage==0) {
				reference = ToUpper(entries[curIndex].reference[j]);
			}
			else {
				assert(reference == ToUpper(entries[curIndex].reference[j]));
			}
			/* Get mismatches and deletions */
			switch(entries[curIndex].read[j]) {
				case 'A':
				case 'a':
					(entries[curIndex].strand==FORWARD)?(numF[0]++):(numR[0]++);
					break;
				case 'C':
				case 'c':
					(entries[curIndex].strand==FORWARD)?(numF[1]++):(numR[1]++);
					break;
				case 'G':
				case 'g':
					(entries[curIndex].strand==FORWARD)?(numF[2]++):(numR[2]++);
					break;
				case 'T':
				case 't':
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
			if(entries[curIndex].reference[j] !=GAP && /* No gap at the current position */
					j<entries[curIndex].referenceLength-1 && /* There are more bases in the alignment */
					entries[curIndex].reference[j+1] == GAP) { /* The next base is a gap */
				/* We started a gap */
				fprintf(outputFPs[1], "chr%d\t%d\t%d\t",
						curChr,
						curPos-1,
						curPos);
				/* Print the insertion */
				for(i=j+1;
						i<entries[curIndex].referenceLength &&
						entries[curIndex].reference[i] == GAP;
						i++) {
					fprintf(outputFPs[1], "%c",
							entries[curIndex].read[i]);
				}
				fprintf(outputFPs[1], "\n");
			}
		}

		/* Print the coverage if it is > 0 */
		if(coverage > 0) {
			int i = -1;
			switch(reference) {
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
					PrintError(FnName,
							"reference",
							"Could not understand reference",
							Exit,
							OutOfRange);
			}
			/* Print the mismatches */
			if(numF[i] + numR[i] != numF[0]+numF[1]+numF[2]+numF[3]+numR[0]+numR[1]+numR[2]+numR[3]) {
				fprintf(outputFPs[0], "chr%d\t%d\t%d\t%c(%d:%d:%d:%d:%3.2lf[F:%d:%d:%d:%d:%3.2lf][R:%d:%d:%d:%d:%3.2lf])\n",
						curChr,
						curPos-1,
						curPos,
						reference,
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
				fprintf(outputFPs[2], "chr%d\t%d\t%d\t%c(%d:%d:%d:%d:%d:%3.2lf[F:%d:%d:%d:%d:%d:%3.2lf][R:%d:%d:%d:%d:%d:%3.2lf])\n",
						curChr,
						curPos-1,
						curPos,
						reference,
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
















