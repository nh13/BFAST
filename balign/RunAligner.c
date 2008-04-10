#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include "ReadInputFiles.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/RGMatch.h" /* To read in the matches */
#include "../blib/RGSeqPair.h" /* To get the reverse compliment */
#include "../blib/AlignEntry.h" /* For output */
#include "Align.h"
#include "Definitions.h"
#include "RunAligner.h"

/* TODO */
void RunAligner(RGBinary *rgBinary,
		char *matchesFileName,
		char *scoringMatrixFileName,
		int algorithm,
		int offsetLength,
		int maxNumMatches,
		int pairedEnd,
		int numThreads,
		char *outputID,
		char *outputDir)
{
	int i;
	FILE *outputFP=NULL;
	FILE *matchesFP=NULL;
	FILE *matchFP=NULL;
	char **matchFileNames=NULL;
	int numMatchFileNames=0;
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	char tempFileName[MAX_FILENAME_LENGTH]="\0";

	/* Open matches file */
	if((matchesFP=fopen(matchesFileName, "r"))==0) {
		PrintError("RunAligner",
				matchesFileName,
				"Could not open matchesFileName for reading",
				Exit,
				OpenFileError);
	}

	/* Read in the list of file names*/
	while(fscanf(matchesFP, "%s", tempFileName)!=EOF) {
		numMatchFileNames++;
		/* Reallocate memory */
		matchFileNames = (char**)realloc(matchFileNames, sizeof(char*)*numMatchFileNames);
		if(NULL==matchFileNames) {
			PrintError("RunAligner",
					"matchFileNames",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		/* Allocate memory */
		matchFileNames[numMatchFileNames-1] = (char*)malloc(sizeof(char)*(strlen(tempFileName+1)));
		if(NULL==matchFileNames[numMatchFileNames-1]) {
			PrintError("RunAligner",
					"matchFileNames[numMatchFileNames-1]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Copy over file name */
		strcpy(matchFileNames[numMatchFileNames-1], tempFileName);
	}

	/* Create output file name */
	sprintf(outputFileName, "%sblatter.aligned.file.%s.%d.%s",
			outputDir,
			outputID,
			algorithm,
			BLATTER_ALIGN_FILE_EXTENSION);

	/* Open output file */
	if((outputFP=fopen(outputFileName, "w"))==0) {
		PrintError("RunAligner",
				outputFileName,
				"Could not open outputFileName for writing",
				Exit,
				OpenFileError);
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "Will output to %s.\n", outputFileName);
	}

	for(i=0;i<numMatchFileNames;i++) { /* For each match file name */

		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Reading match file #%d from %s.\n",
					i+1,
					matchFileNames[i]);
		}

		/* Open current match file */
		if((matchFP=fopen(matchFileNames[i], "r"))==0) {
			PrintError("RunAligner",
					matchFileNames[i],
					"Could not open matchesFileNames[] for reading",
					Exit,
					OpenFileError);
		}

		/* Run selected algorithm */
		switch(algorithm) {
			case 0:
				RunDynamicProgramming(matchFP,
						rgBinary,
						scoringMatrixFileName,
						offsetLength,
						maxNumMatches,
						pairedEnd,
						numThreads,
						outputFP);
				break;
			default:
				PrintError("RunAligner",
						NULL,
						"Could not understand algorithm option",
						Exit,
						OutOfRange);
				break;
		}
		/* Close the match file */
		fclose(matchFP);
		matchFP=NULL;
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "%s", BREAK_LINE);
	}

	/* Close output file */
	fclose(outputFP);

	/* Free memory */
	for(i=0;i<numMatchFileNames;i++) {
		free(matchFileNames[i]);
	}
	free(matchFileNames);
}

/* TODO */
void RunDynamicProgramming(FILE *matchFP,
		RGBinary *rgBinary,
		char *scoringMatrixFileName,
		int offsetLength,
		int maxNumMatches,
		int pairedEnd,
		int numThreads,
		FILE *outputFP)
{
	/* local variables */
	ScoringMatrix sm;
	char readName[SEQUENCE_NAME_LENGTH]="\0";
	char read[SEQUENCE_LENGTH]="\0";
	char pairedSequence[SEQUENCE_LENGTH]="\0";
	RGMatch readMatch;
	RGMatch pairedSequenceMatch;
	int i, j;
	int numMatches=0;
	int continueReading=0;
	int numAlignEntries;
	int numAlignments=0;
	AlignEntry aEntry;
	/* Thread specific data */
	ThreadData *data;
	pthread_t *threads=NULL;
	int errCode;
	void *status;

	/* Allocate memory for thread arguments */
	data = malloc(sizeof(ThreadData)*numThreads);
	if(NULL==data) {
		PrintError("RunDynamicProgramming",
				"data",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	/* Allocate memory for thread ids */
	threads = malloc(sizeof(pthread_t)*numThreads);
	if(NULL==threads) {
		PrintError("RunDynamicProgramming",
				"threads",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Read in scoring matrix */
	ReadScoringMatrix(scoringMatrixFileName, &sm); 

	/**/
	/* Split the input file into equal temp files for each thread */
	/**/

	/* Open temp files for the threads */
	for(i=0;i<numThreads;i++) {
		data[i].inputFP = tmpfile();
		data[i].outputFP = tmpfile();
	}

	/* Initialize match */
	readMatch.positions=NULL;
	readMatch.chromosomes=NULL;
	readMatch.strand=NULL;
	readMatch.numEntries=0;
	pairedSequenceMatch.positions=NULL;
	pairedSequenceMatch.chromosomes=NULL;
	pairedSequenceMatch.strand=NULL;
	pairedSequenceMatch.numEntries=0;

	/* Go through each read in the match file */
	numMatches=0;
	while(EOF!=RGMatchGetNextFromFile(matchFP, 
				readName, 
				read, 
				pairedSequence, 
				&readMatch,
				&pairedSequenceMatch,
				pairedEnd)) {
		/* Get the thread index - do this BEFORE incrementing */
		int threadIndex = numMatches%numThreads;
		/* increment */
		numMatches++;

		/* Print match to temp file */
		RGMatchOutputToFile(data[threadIndex].inputFP,
				readName,
				read,
				pairedSequence,
				&readMatch,
				&pairedSequenceMatch,
				pairedEnd);

		/* Free match */
		if(readMatch.numEntries > 0) {
			/* Free AlignEntry */
			free(readMatch.positions);
			free(readMatch.chromosomes);
			free(readMatch.strand);
		}
		if(pairedSequenceMatch.numEntries > 0) {
			free(pairedSequenceMatch.positions);
			free(pairedSequenceMatch.chromosomes);
			free(pairedSequenceMatch.strand);
		}
		/* Initialize match */
		readMatch.positions=NULL;
		readMatch.chromosomes=NULL;
		readMatch.strand=NULL;
		readMatch.numEntries=0;
		pairedSequenceMatch.positions=NULL;
		pairedSequenceMatch.chromosomes=NULL;
		pairedSequenceMatch.strand=NULL;
		pairedSequenceMatch.numEntries=0;
	}

	/* Create thread arguments */
	for(i=0;i<numThreads;i++) {
		fseek(data[i].inputFP, 0, SEEK_SET);
		fseek(data[i].outputFP, 0, SEEK_SET);
		data[i].rgBinary=rgBinary;
		data[i].offsetLength=offsetLength;
		data[i].maxNumMatches=maxNumMatches;
		data[i].pairedEnd=pairedEnd;
		data[i].sm = &sm;
	}

	if(VERBOSE >= 0) {
		fprintf(stderr, "Performing alignment...\n");
	}

	/* Create threads */
	for(i=0;i<numThreads;i++) {
		/* Start thread */
		errCode = pthread_create(&threads[i], /* thread struct */
				NULL, /* default thread attributes */
				RunDynamicProgrammingThread, /* start routine */
				&data[i]); /* data to routine */
		if(0!=errCode) {
			PrintError("RunDynamicProgramming",
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
			PrintError("FindMatchesInIndexes",
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
		fprintf(stderr, "Alignment complete.\n");
	}

	/* Allocate memory for the align entries */
	aEntry.read = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL==aEntry.read) {
		PrintError("RunDynamicProgramming",
				"aEntry.read",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	aEntry.reference = malloc(sizeof(char)*SEQUENCE_LENGTH);
	if(NULL==aEntry.reference) {
		PrintError("RunDynamicProgramming",
				"aEntry.reference",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	
	if(VERBOSE >=0) {
		fprintf(stderr, "Outputting...\n");
	}

	/* Merge all the outputs together */
	numAlignments=0;
	continueReading=1;
	while(continueReading==1) {
		/* Get one align from each thread */
		continueReading=0;
		for(i=0;i<numThreads;i++) {
			/* First get the number of align entries */
			if(EOF != fscanf(data[i].outputFP, "%d", &numAlignEntries)) {
				continueReading=1;
				assert(numAlignEntries >= 0);
				numAlignments++;
				for(j=0;j<numAlignEntries;j++) {
					/* Read in the align entry */
					if(EOF==AlignEntryRead(&aEntry,
								data[i].outputFP)) {
						PrintError("RunDynamicProgramming",
								"AlignEntryRead",
								"Could not read in the align entry",
								Exit,
								EndOfFile);
					}
					/* Print the align entry to file */
					AlignEntryPrint(&aEntry,
							outputFP);
				}
			}
		}
	}
	
	if(VERBOSE >=0) {
		fprintf(stderr, "Outputted alignments for %d reads.\n", numAlignments);
		fprintf(stderr, "Outputting complete.\n");
	}

	/* Free memory */
	free(data);
	free(threads);
	free(aEntry.read);
	free(aEntry.reference);
	/* Free scores */
	free(sm.key);
	for(i=0;i<ALPHABET_SIZE+1;i++) {
		free(sm.scores[i]);
	}
	free(sm.scores);
}

/* TODO */
void *RunDynamicProgrammingThread(void *arg)
{
	/* Recover arguments */
	ThreadData *data = (ThreadData *)(arg);
	FILE *inputFP=data->inputFP;
	FILE *outputFP=data->outputFP;
	RGBinary *rgBinary=data->rgBinary;
	int offsetLength=data->offsetLength;
	int maxNumMatches=data->maxNumMatches;
	int pairedEnd=data->pairedEnd;
	ScoringMatrix *sm = data->sm;
	/* Local variables */
	AlignEntry *aEntry=NULL;
	int numAlignEntries=0;
	char readName[SEQUENCE_NAME_LENGTH]="\0";
	char read[SEQUENCE_LENGTH]="\0";
	char pairedSequence[SEQUENCE_LENGTH]="\0";
	RGMatch readMatch;
	RGMatch pairedSequenceMatch;
	int matchLength;
	char *reference=NULL;
	int referenceLength=0;
	int adjustPosition=0;
	int i;
	int position;
	int numMatches=0;
	int numMatchesAligned=0;

	/* Allocate memory for the reference */
	referenceLength = 2*offsetLength + SEQUENCE_LENGTH + 1;
	reference = (char*)malloc(sizeof(char)*(referenceLength+1));
	if(NULL==reference) {
		PrintError("RunDynamicProgrammingThread",
				"reference",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	reference[referenceLength] = '\0'; /* Add null terminator */

	/* Initialize match */
	readMatch.positions=NULL;
	readMatch.chromosomes=NULL;
	readMatch.strand=NULL;
	readMatch.numEntries=0;
	pairedSequenceMatch.positions=NULL;
	pairedSequenceMatch.chromosomes=NULL;
	pairedSequenceMatch.strand=NULL;
	pairedSequenceMatch.numEntries=0;

	/* Go through each read in the match file */
	while(EOF!=RGMatchGetNextFromFile(inputFP, 
				readName, 
				read, 
				pairedSequence, 
				&readMatch,
				&pairedSequenceMatch,
				pairedEnd)) {
		numMatches++;

		numAlignEntries = 0;

		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "\nreadMatch.numEntries=%d.\n",
					readMatch.numEntries);
		}
		if(readMatch.numEntries > 0 && (maxNumMatches == 0 || readMatch.numEntries < maxNumMatches)) {

			/* Get the read length */
			matchLength = strlen(read);
			assert(matchLength+2*offsetLength < SEQUENCE_LENGTH);

			/* Allocate memory for the AlignEntries */
			if(readMatch.numEntries > 0) {
				aEntry = (AlignEntry*)malloc(sizeof(AlignEntry)*readMatch.numEntries);
				if(NULL == aEntry) {
					PrintError("RunDynamicProgrammingThread",
							"aEntry",
							"Could not allocate memory",
							Exit,
							MallocMemory);
				}
			}

			/* Run the aligner */
			assert(pairedEnd==0);
			if(readMatch.numEntries > 0) {
				numMatchesAligned++;
			}
			for(i=0;i<readMatch.numEntries;i++) { /* For each match */
				if(VERBOSE >= DEBUG) {
					fprintf(stderr, "On entry %d out of %d.  chr%d:%d %c\n",
							i+1,
							readMatch.numEntries,
							readMatch.chromosomes[i],
							readMatch.positions[i],
							readMatch.strand[i]);
				}

				/* Get the appropriate reference read */
				GetSequenceFromReferenceGenome(rgBinary, 
						readMatch.chromosomes[i], 
						readMatch.positions[i],
						readMatch.strand[i],
						offsetLength,
						reference,
						matchLength,
						&referenceLength,
						&position);
				assert(referenceLength > 0);

				/* Get alignment */
				adjustPosition=AlignmentGetScore(read,
						matchLength,
						reference,
						referenceLength,
						sm,
						&aEntry[i]);

				/* Update chromosome, position, strand and sequence name*/
				aEntry[i].chromosome = readMatch.chromosomes[i];
				aEntry[i].position = position; /* remember to use the update position */
				aEntry[i].position = (readMatch.strand[i] == FORWARD)?(aEntry[i].position+adjustPosition):(aEntry[i].position-adjustPosition); /* Adjust position */
				aEntry[i].strand = readMatch.strand[i]; 
				strcpy(aEntry[i].readName, readName);


				/* Free memory */
				/* Free AlignEntry */
				if(VERBOSE >= DEBUG) {
					fprintf(stderr, "aligned read:%s\naligned reference:%s\nposition:%d\n",
							aEntry[i].read,
							aEntry[i].reference,
							aEntry[i].position);
				}
				if(VERBOSE >= DEBUG) {
					fprintf(stderr, "Finished entry %d out of %d.\n",
							i+1,
							readMatch.numEntries);
				}
			}
			/* Remove duplicate alignments */
			numAlignEntries=AlignEntryRemoveDuplicates(&aEntry, readMatch.numEntries);
			if(VERBOSE >= DEBUG) {
				fprintf(stderr, "Outputting %d aligns.\n",
						numAlignEntries);
			}
			/* Output alignment */
			fprintf(outputFP, "%d\n", numAlignEntries); /* We need this when merging temp files */
			for(i=0;i<numAlignEntries;i++) {
				AlignEntryPrint(&aEntry[i], outputFP);
			}
			/* Free memory */
			for(i=0;i<numAlignEntries;i++) {
				assert(aEntry[i].length>0);
				free(aEntry[i].read);
				free(aEntry[i].reference);
			}
		}
		/* Free match */
		if(readMatch.numEntries > 0) {
			/* Free AlignEntry */
			free(readMatch.positions);
			free(readMatch.chromosomes);
			free(readMatch.strand);
		}
		if(pairedSequenceMatch.numEntries > 0) {
			free(pairedSequenceMatch.positions);
			free(pairedSequenceMatch.chromosomes);
			free(pairedSequenceMatch.strand);
		}
		if(numAlignEntries > 0) {
			free(aEntry);
			aEntry = NULL;
		}
		/* Initialize match */
		readMatch.positions=NULL;
		readMatch.chromosomes=NULL;
		readMatch.strand=NULL;
		readMatch.numEntries=0;
		pairedSequenceMatch.positions=NULL;
		pairedSequenceMatch.chromosomes=NULL;
		pairedSequenceMatch.strand=NULL;
		pairedSequenceMatch.numEntries=0;
	}

	/* Free memory */
	free(reference);

	return arg;
}

/* TODO */
void GetSequenceFromReferenceGenome(RGBinary *rgBinary,
		int chromosome,
		int position,
		char strand,
		int offsetLength,
		char *reference,
		int matchLength,
		int *returnReferenceLength,
		int *returnPosition)
{
	int i;
	int chrIndex;
	int posIndex;
	int byteIndex;
	int curByteIndex;
	unsigned char curChar;
	char c;
	unsigned char repeat;
	char *reverseCompliment;
	int numCharsPerByte;
	int startPos, endPos;
	int referenceLength = 2*offsetLength + matchLength;
	/* We assume that we can hold 2 [acgt] (nts) in each byte */
	assert(ALPHABET_SIZE==4);
	numCharsPerByte=ALPHABET_SIZE/2;

	/* Get the start and end positions */
	startPos=-1;
	endPos=-1;
	if(FORWARD == strand) {
		startPos = position - offsetLength;
		endPos = position + matchLength - 1 + offsetLength;
	}
	else if(REVERSE == strand) {
		startPos = position - matchLength + 1 - offsetLength;
		endPos = position + offsetLength;
	}
	else {
		PrintError("GetSequenceFromReferenceGenome",
				NULL,
				"Could not recognize strand",
				Exit,
				OutOfRange);
	}

	/* Get chromosome index in rgBinary */
	chrIndex = chromosome - rgBinary->startChr;

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "In GetSequenceFromReferenceGenome: user chromosome [%d] position [%d] strand [%c].\n",
				chromosome,
				position,
				strand);
		fprintf(stderr, "In GetSequenceFromReferenceGenome: chromosome [%d] with range [%d,%d].\n",
				chromosome,
				rgBinary->startChr,
				rgBinary->endChr);
	}

	/* Check chromosome bounds */
	if(chromosome < rgBinary->startChr || chromosome > rgBinary->endChr) {
		PrintError("GetSequenceFromReferenceGenome",
				NULL,
				"Chromosome is out of range",
				Exit,
				OutOfRange);
	}

	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "In GetSequenceFromReferenceGenome: start position [%d] range [%d,%d] and adjusted range [%d,%d]\n",
				position,
				rgBinary->chromosomes[chrIndex].startPos,
				rgBinary->chromosomes[chrIndex].endPos,
				startPos,
				endPos);
	}

	/* Check position bounds */
	if(startPos < rgBinary->chromosomes[chrIndex].startPos || endPos > rgBinary->chromosomes[chrIndex].endPos) {
		/* Adjust position bounds if possible */
		if(startPos < rgBinary->chromosomes[chrIndex].startPos) {
			/* Check that we have enough sequence */
			startPos = rgBinary->chromosomes[chrIndex].startPos;
			assert(startPos + matchLength - 1 + 2*offsetLength <= rgBinary->chromosomes[chrIndex].endPos);
			endPos = startPos + matchLength - 1 + 2*offsetLength;
		}
		else if(endPos > rgBinary->chromosomes[chrIndex].endPos) {
			/* Check that we have enough sequence */
			endPos = rgBinary->chromosomes[chrIndex].endPos;
			assert(endPos - matchLength + 1 - 2*offsetLength >= rgBinary->chromosomes[chrIndex].startPos);
			startPos = endPos - matchLength + 1 - 2*offsetLength;
		}
	}

	if(position < rgBinary->chromosomes[chrIndex].startPos || position > rgBinary->chromosomes[chrIndex].endPos) {
		PrintError("GetSequenceFromReferenceGenome",
				NULL,
				"Position out of range",
				Exit,
				OutOfRange);
	}

	/* Get the reference sequence */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "startPos:%d\tendPos:%d\n",
				startPos,
				endPos);
	}
	assert(startPos <= endPos);
	assert(startPos >= 1);
	for(i=startPos;i<=endPos;i++) {
		/* Get position index in rgBinary */
		posIndex = i - rgBinary->chromosomes[chrIndex].startPos;
		byteIndex = posIndex%numCharsPerByte; /* Which bits in the byte */
		curByteIndex = (posIndex - byteIndex)/numCharsPerByte; /* Get which byte */

		/* Get the sequence */
		curChar = rgBinary->chromosomes[chrIndex].sequence[curByteIndex];
		repeat = 0;
		c = 'E';
		switch(byteIndex) {
			case 0:
				/* left-most 2-bits */
				repeat = curChar & 0xC0; /* zero out the irrelevant bits */
				switch(repeat) {
					case 0x00:
						repeat = 0; 
						break;
					case 0x40:
						repeat = 1;
						break;
					case 0x80:
						repeat = 2;
						break;
					default:
						PrintError("GetSequenceFromReferenceGenome",
								NULL,
								"Could not understand case 0 repeat",
								Exit,
								OutOfRange);
						break;
				}
				/* third and fourth bits from the left */
				curChar = curChar & 0x30; /* zero out the irrelevant bits */
				switch(curChar) {
					case 0x00:
						c = 'a';
						break;
					case 0x10:
						c = 'c';
						break;
					case 0x20:
						c = 'g';
						break;
					case 0x30:
						c = 't';
						break;
					default:
						PrintError("GetSequenceFromReferenceGenome",
								NULL,
								"Could not understand case 0 base",
								Exit,
								OutOfRange);
						break;
				}
				break;
			case 1:
				/* third and fourth bits from the right */
				repeat = curChar & 0x0C; /* zero out the irrelevant bits */
				switch(repeat) {
					case 0x00:
						repeat = 0;
						break;
					case 0x04:
						repeat = 1;
						break;
					case 0x08:
						repeat = 2;
						break;
					default:
						PrintError("GetSequenceFromReferenceGenome",
								NULL,
								"Could not understand case 1 repeat",
								Exit,
								OutOfRange);
						break;
				}
				/* right-most 2-bits */
				curChar = curChar & 0x03; /* zero out the irrelevant bits */
				switch(curChar) {
					case 0x00:
						c = 'a';
						break;
					case 0x01:
						c = 'c';
						break;
					case 0x02:
						c = 'g';
						break;
					case 0x03:
						c = 't';
						break;
					default:
						PrintError("GetSequenceFromReferenceGenome",
								NULL,
								"Could not understand case 1 base",
								Exit,
								OutOfRange);
						break;
				}
				break;
			default:
				PrintError("GetSequenceFromReferenceGenome",
						"byteIndex",
						"Could not understand byteIndex",
						Exit,
						OutOfRange);
		}
		/* Update based on repeat */
		switch(repeat) {
			case 0:
				/* ignore, not a repeat */
				break;
			case 1:
				/* repeat, convert char to upper */
				c=ToUpper(c);
				break;
			case 2:
				/* N character */
				c='N';
				break;
			default:
				fprintf(stderr, "Error.  In GetSequenceFromReferenceGenome, could not understand repeat indexed [%d].  Terminating!\n",
						repeat);
				break;
		}
		/* Error check */
		switch(c) {
			case 'a':
			case 'c':
			case 'g':
			case 't':
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case 'N':
				break;
			default:
				PrintError("GetSequenceFromReferenceGenome",
						NULL,
						"Could not understand base",
						Exit,
						OutOfRange);
		}
		/* Update reference */
		reference[i-startPos] = c;
	}
	reference[i-startPos] = '\0';
	/* Get the reverse compliment if necessary */
	if(strand == FORWARD) {
		/* ignore */
	}
	else if(strand == REVERSE) {
		/* Get the reverse compliment */
		reverseCompliment = (char*)malloc(sizeof(char)*(referenceLength+1));
		if(NULL == reverseCompliment) {
			PrintError("GetSequenceFromReferenceGenome",
					"reverseCompliment",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		GetReverseComplimentAnyCase(reference, reverseCompliment, referenceLength);
		strcpy(reference, reverseCompliment);
		free(reverseCompliment);
	}
	else {
		PrintError("GetSequenceFromReferenceGenome",
				"strand",
				"Could not understand strand",
				Exit,
				OutOfRange);
	}
	/* Update start pos and reference length */
	(*returnReferenceLength) = referenceLength;
	(*returnPosition) = startPos;
	if(VERBOSE>=DEBUG) {
		fprintf(stderr, "Exiting GetSequenceFromReferenceGenome:[%s] length [%d] referenceLength [%d] and startPos [%d].\n",
				reference,
				(int)strlen(reference),
				referenceLength,
				startPos);
	}
}

void GetReverseComplimentAnyCase(char *s, 
		char *r,
		int length)
{
	int i;
	/* Get reverse compliment sequence */
	for(i=length-1;i>=0;i--) {
		switch(s[length-1-i]) {
			case 'a':
				r[i] = 't';
				break;
			case 'c':
				r[i] = 'g';
				break;
			case 'g':
				r[i] = 'c';
				break;
			case 't':
				r[i] = 'a';
				break;
			case 'A':
				r[i] = 'T';
				break;
			case 'C':
				r[i] = 'G';
				break;
			case 'G':
				r[i] = 'C';
				break;
			case 'T':
				r[i] = 'A';
				break;
			default:
				PrintError("GetReverseComplimentAnyCase",
						NULL,
						"Could not understand sequence base",
						Exit,
						OutOfRange);
				break;
		}
	}
	r[length]='\0';
}
