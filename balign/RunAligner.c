#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "ReadInputFiles.h"
#include "../blib/BLibDefinitions.h"
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
		int pairedEnd,
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
		fprintf(stderr, "Error.  Could not open %s for reading.  Terminating!\n", matchesFileName);
		exit(1);
	}

	/* Read in the list of file names*/
	while(fscanf(matchesFP, "%s", tempFileName)!=EOF) {
		numMatchFileNames++;
		/* Reallocate memory */
		matchFileNames = (char**)realloc(matchFileNames, sizeof(char*)*numMatchFileNames);
		/* Allocate memory */
		matchFileNames[numMatchFileNames-1] = (char*)malloc(sizeof(char)*(strlen(tempFileName+1)));
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
		fprintf(stderr, "Error.  Could not open %s for writing.  Terminating!\n",
				outputFileName);
		exit(1);
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
			fprintf(stderr, "Error.  Could not open %s for reading.  Terminating!\n",
					matchFileNames[i]);
			exit(1);
		}

		/* Run selected algorithm */
		switch(algorithm) {
			case 0:
				RunDynamicProgramming(matchFP,
						rgBinary,
						scoringMatrixFileName,
						algorithm,
						offsetLength,
						pairedEnd,
						outputFP);
				break;
			default:
				fprintf(stderr, "Error.  Could not understand algorithm option %d.  Terminating!\n",
						algorithm);
				exit(1);
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
		int algorithm,
		int offsetLength,
		int pairedEnd,
		FILE *outputFP)
{
	AlignEntry *aEntry=NULL;
	int numAlignEntries=0;
	ScoringMatrix sm;
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
	reference[referenceLength] = '\0'; /* Add null terminator */

	/* Read in scoring matrix */
	ReadScoringMatrix(scoringMatrixFileName, &sm); 

	/* Initialize match */
	readMatch.positions=NULL;
	readMatch.chromosomes=NULL;
	readMatch.strand=NULL;
	readMatch.numEntries=0;
	pairedSequenceMatch.positions=NULL;
	pairedSequenceMatch.chromosomes=NULL;
	pairedSequenceMatch.strand=NULL;
	pairedSequenceMatch.numEntries=0;

	if(VERBOSE >= 0) {
		fprintf(stderr, "Currently on match:\n0");
	}
	/* Go through each read in the match file */
	while(EOF!=RGMatchGetNextFromFile(matchFP, 
				readName, 
				read, 
				pairedSequence, 
				&readMatch,
				&pairedSequenceMatch,
				pairedEnd)) {
		numMatches++;

		if(VERBOSE >= 0 && numMatches%ALIGN_ROTATE_NUM==0) {
			fprintf(stderr, "\r%d", numMatches);
		}
		if(VERBOSE >= DEBUG) {
			fprintf(stderr, "\nreadMatch.numEntries=%d.\n",
					readMatch.numEntries);
		}

		/* Get the read length */
		matchLength = strlen(read);
		assert(matchLength+2*offsetLength < SEQUENCE_LENGTH);

		/* Allocate memory for the AlignEntries */
		if(readMatch.numEntries > 0) {
			aEntry = (AlignEntry*)malloc(sizeof(AlignEntry)*readMatch.numEntries);
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
					&sm,
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
		for(i=0;i<numAlignEntries;i++) {
			AlignEntryPrint(&aEntry[i], outputFP);
		}
		/* Free memory */
		for(i=0;i<numAlignEntries;i++) {
			assert(aEntry[i].length>0);
			free(aEntry[i].read);
			free(aEntry[i].reference);
		}
		/* Free match */
		if(readMatch.numEntries > 0) {
			/* Free AlignEntry */
			free(aEntry);
			aEntry = NULL;
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
	if(VERBOSE >= 0) {
		fprintf(stderr, "\n");
		fprintf(stderr, "There were %d empty matches.\n", numMatches - numMatchesAligned);
		fprintf(stderr, "Aligned %d matches.\n", numMatchesAligned);
	}

	/* Free memory */
	free(reference);
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
	if(FORWARD == strand) {
		startPos = position - offsetLength;
		endPos = position + matchLength - 1 + offsetLength;
	}
	else if(REVERSE == strand) {
		startPos = position - matchLength + 1 - offsetLength;
		endPos = position + offsetLength;
	}
	else {
		fprintf(stderr, "Error.  Could not recognize strand [%c].  Terminating!\n",
				strand);
		exit(1);
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
		fprintf(stderr, "Error.  In GetSequenceFromReferenceGenome: chromosome [%d] out of range [%d,%d]\n",
				chromosome,
				rgBinary->startChr,
				rgBinary->endChr);
		exit(1);
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
		fprintf(stderr, "Error.  In GetSequenceFromReferenceGenome: end position [%d] out of range [%d,%d] for chromosome [%d].  Terminating!\n",
				endPos,
				rgBinary->chromosomes[chrIndex].startPos,
				rgBinary->chromosomes[chrIndex].endPos,
				chromosome);
		exit(1);
	}

	/* Get the reference sequence */
	if(VERBOSE >= DEBUG) {
		fprintf(stderr, "startPos:%d\tendPos:%d\n",
				startPos,
				endPos);
	}
	assert(startPos <= endPos);
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
						fprintf(stderr, "Error.  In GetSequenceFromReferenceGenome, could not understand repeat [%c].  Terminating!\n",
								repeat);
						exit(1);
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
						fprintf(stderr, "Error.  Could not get third and fourth bits from the left of %x.  Terminating!\n",
								curChar);
						exit(1);
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
						fprintf(stderr, "Error.  In GetSequenceFromReferenceGenome, could not understand repeat [%c].  Terminating!\n",
								repeat);
						exit(1);
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
						fprintf(stderr, "Error.  Could not get right-most 2-bits of %x.  Terminating!\n",
								curChar);
						exit(1);
						break;
				}
				break;
			default:
				fprintf(stderr, "Error.  In GetSequenceFromReferenceGenome, could not understand byteIndex [%d].  Terminating!\n",
						byteIndex);
				exit(1);
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
				fprintf(stderr, "Error.  In GetSequenceFromReferenceGenome, could not understand final char [%c].  Terminating!\n",
						c);
				exit(1);
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
		GetReverseComplimentAnyCase(reference, reverseCompliment, referenceLength);
		strcpy(reference, reverseCompliment);
		free(reverseCompliment);
	}
	else {
		fprintf(stderr, "Error.  Could not understand strandedness [%c].  Terminating!\n", strand);
		exit(1);
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
				fprintf(stderr, "Error.  In GetReverseComplimentAnyCase, could not understand %c.  Terminating!\n", s[length-1-i]);
				exit(1);
				break;
		}
	}
	r[length]='\0';
}
