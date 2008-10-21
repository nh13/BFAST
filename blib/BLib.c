#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "BLibDefinitions.h"
#include "BIndex.h"
#include "BError.h"
#include "BLib.h"

char DNA[5] = "ACGTN";
char COLORS[5] = "01234";

/* TODO */
int32_t GetFastaHeaderLine(FILE *fp,
		BString *header);
{
	char *FnName="GetFastaHeaderLine";
	int i, length;

	/* Read in the line */
	if(BStringGetLine(header, fp, MAX_CONTIG_NAME_LENGTH) == EOF) {
		return EOF;
	}

	/* Check that the first character is a ">" */
	if(header->string[0] != '>') {
		PrintError(FnName,
				"header",
				"Header of a contig must start with a '>'",
				Exit,
				OutOfRange);
	}

	/* Shift over to remove the ">" and remove whitespaces */
	BStringShiftLeftAndReplaceWhitespace(header, 1, '_');

	return 1;
}

/* TODO */
int8_t ToLower(int8_t a) 
{
	switch(a) {
		case 'A':
			return 'a';
			break;
		case 'C':
			return 'c';
			break;
		case 'G':
			return 'g';
			break;
		case 'T':
			return 't';
			break;
		case 'N':
			return 'n';
			break;
		default:
			return a;
	}
}

/* TODO */
int8_t ToUpper(int8_t a) 
{
	switch(a) {
		case 'a':
			return 'A';
			break;
		case 'c':
			return 'C';
			break;
		case 'g':
			return 'G';
			break;
		case 't':
			return 'T';
			break;
		case 'n':
			return 'N';
			break;
		default:
			return a;
	}
}

int8_t GetReverseComplimentAnyCaseBase(int8_t a) 
{
	char *FnName = "GetReverseComplimentAnyCaseBase";
	switch(a) {
		case 'a':
			return 't';
			break;
		case 'c':
			return 'g';
			break;
		case 'g':
			return 'c';
			break;
		case 't':
			return 'a';
			break;
		case 'n':
			return 'n';
			break;
		case 'A':
			return 'T';
			break;
		case 'C':
			return 'G';
			break;
		case 'G':
			return 'C';
			break;
		case 'T':
			return 'A';
			break;
		case 'N':
			return 'N';
			break;
		case GAP:
			return GAP;
			break;
		default:
			fprintf(stderr, "\n[%c]\t[%d]\n",
					a,
					(int)a);
			PrintError(FnName,
					NULL,
					"Could not understand sequence base",
					Exit,
					OutOfRange);
			break;
	}
	PrintError(FnName,
			NULL,
			"Control should not reach here",
			Exit,
			OutOfRange);
	return '0';
}

/* TODO */
int32_t ValidateBasePair(int8_t c) {
	switch(c) {
		case 'a':
		case 'c':
		case 'g':
		case 't':
		case 'A':
		case 'C':
		case 'G':
		case 'T':
		case 'n':
		case 'N':
			return 1;
			break;
		default:
			return 0;
			break;
	}
}

int32_t IsAPowerOfTwo(unsigned int a) {
	int i;

	for(i=0;i<8*sizeof(unsigned int);i++) {
		/*
		   fprintf(stderr, "i:%d\ta:%d\tshifted:%d\tres:%d\n",
		   i,
		   a,
		   a>>i,
		   (a >> i)%2);
		   */
		if( (a >> i) == 2) {
			return 1;
		}
		else if( (a >> i)%2 != 0) {
			return 0;
		}
	}
	return 1;
}

/* TODO */
uint32_t Log2(uint32_t num) 
{
	char *FnName = "Log2";
	int i;

	if(IsAPowerOfTwo(num)==0) {
		PrintError(FnName,
				"num",
				"Num is not a power of 2",
				Exit,
				OutOfRange);
	}
	/* Not the most efficient but we are not going to use this often */
	for(i=0;num>1;i++,num/=2) {
	}
	return i;
}

int8_t TransformFromIUPAC(int8_t a) 
{
	switch(a) {
		case 'U':
			return 'T';
			break;
		case 'u':
			return 't';
			break;
		case 'R':
		case 'Y':
		case 'M':
		case 'K':
		case 'W':
		case 'S':
		case 'B':
		case 'D':
		case 'H':
		case 'V':
			return 'N';
			break;
		case 'r':
		case 'y':
		case 'm':
		case 'k':
		case 'w':
		case 's':
		case 'b':
		case 'd':
		case 'h':
		case 'v':
			return 'n';
			break;
		default:
			return a;
			break;
	}
}

void CheckBIndexes(BString *mainFileNames,
		int numMainFileNames,
		BString *secondaryFileNames,
		int numSecondaryFileNames,
		int binaryInput,
		int32_t *startContig,
		int32_t *startPos,
		int32_t *endContig,
		int32_t *endPos,
		int32_t space)
{
	int i;
	int32_t mainStartContig, mainStartPos, mainEndContig, mainEndPos;
	int32_t secondaryStartContig, secondaryStartPos, secondaryEndContig, secondaryEndPos;
	int32_t mainColorSpace=space;
	int32_t secondaryColorSpace=space;
	int32_t mainContigType=0;
	int32_t secondaryContigType=0;
	mainStartContig = mainStartPos = mainEndContig = mainEndPos = 0;
	secondaryStartContig = secondaryStartPos = secondaryEndContig = secondaryEndPos = 0;

	BIndex tempIndex;
	FILE *fp;

	/* Read in main indexes */
	for(i=0;i<numMainFileNames;i++) {
		/* Open file */
		if((fp=fopen(mainFileNames[i], "r"))==0) {
			PrintError("CheckBIndexes",
					mainFileNames[i],
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}

		/* Get the header */
		BIndexReadHeader(fp, &tempIndex, binaryInput); 

		assert(tempIndex.startContig < tempIndex.endContig ||
				(tempIndex.startContig == tempIndex.endContig && tempIndex.startPos <= tempIndex.endPos));

		if(i==0) {
			mainContigType = tempIndex.contigType;
			mainStartContig = tempIndex.startContig;
			mainStartPos = tempIndex.startPos;
			mainEndContig = tempIndex.endContig;
			mainEndPos = tempIndex.endPos;
			mainColorSpace = tempIndex.space;
		}
		else {
			/* Update bounds if necessary */
			assert(mainContigType == tempIndex.contigType);
			if(tempIndex.startContig < mainStartContig ||
					(tempIndex.startContig == mainStartContig && tempIndex.startPos < mainStartPos)) {
				mainStartContig = tempIndex.startContig;
				mainStartPos = tempIndex.startPos;
			}
			if(tempIndex.endContig > mainEndContig ||
					(tempIndex.endContig == mainEndContig && tempIndex.endPos > mainEndPos)) {
				mainEndContig = tempIndex.endContig;
				mainEndPos = tempIndex.endPos;
			}
			assert(mainColorSpace == tempIndex.space);
		}

		/* Free masks */
		free(tempIndex.mask);
		tempIndex.mask=NULL;

		/* Close file */
		fclose(fp);
	}
	/* Read in secondary indexes */
	for(i=0;i<numSecondaryFileNames;i++) {
		/* Open file */
		if((fp=fopen(secondaryFileNames[i], "r"))==0) {
			PrintError("CheckBIndexes",
					"secondaryFileNames[i]",
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}

		/* Get the header */
		BIndexReadHeader(fp, &tempIndex, binaryInput); 

		assert(tempIndex.startContig < tempIndex.endContig ||
				(tempIndex.startContig == tempIndex.endContig && tempIndex.startPos <= tempIndex.endPos));

		if(i==0) {
			secondaryContigType = tempIndex.contigType;
			secondaryStartContig = tempIndex.startContig;
			secondaryStartPos = tempIndex.startPos;
			secondaryEndContig = tempIndex.endContig;
			secondaryEndPos = tempIndex.endPos;
			secondaryColorSpace = tempIndex.space;
		}
		else {
			/* Update bounds if necessary */
			assert(secondaryContigType == tempIndex.contigType);
			if(tempIndex.startContig < secondaryStartContig ||
					(tempIndex.startContig == secondaryStartContig && tempIndex.startPos < secondaryStartPos)) {
				secondaryStartContig = tempIndex.startContig;
				secondaryStartPos = tempIndex.startPos;
			}
			if(tempIndex.endContig > secondaryEndContig ||
					(tempIndex.endContig == secondaryEndContig && tempIndex.endPos > secondaryEndPos)) {
				secondaryEndContig = tempIndex.endContig;
				secondaryEndPos = tempIndex.endPos;
			}
			assert(secondaryColorSpace == tempIndex.space);
		}

		/* Free masks */
		free(tempIndex.mask);
		tempIndex.mask=NULL;

		/* Close file */
		fclose(fp);
	}

	/* Check the bounds between main and secondary indexes */
	assert(numSecondaryFileNames == 0 ||
			mainContigType == secondaryContigType);
	if(mainStartContig != secondaryStartContig ||
			mainStartPos != secondaryStartPos ||
			mainEndContig != secondaryEndContig ||
			mainEndPos != secondaryEndPos ||
			mainColorSpace != secondaryColorSpace) {
		PrintError("CheckBIndexes",
				NULL,
				"The ranges between main and secondary indexes differ",
				Exit,
				OutOfRange);
	}

	(*startContig) = mainStartContig;
	(*startPos) = mainStartPos;
	(*endContig) = mainEndContig;
	(*endPos) = mainEndPos;

	assert(mainColorSpace == space);
	assert(secondaryColorSpace == space);
}

/* TODO */
FILE *OpenTmpFile(BString *tmpDir,
		BString *tmpFileName) 
{
	char *FnName = "OpenTmpFile";
	FILE *fp;


	/* Create the templated file name */
	/* Allocate memory */
	BStringAllocate(tmpFileName, MAX_FILENAME_LENGTH);
	/* Copy over tmp directory */
	strcpy(tmpFileName->string, tmpDir->string);
	/* Copy over the tmp name */
	strcat(tmpFileName->string, BFAST_TMP_TEMPLATE);
	/* Reallocate memory */
	BStringReallocate(tmpFileName, strlen(tmpFileName->string));

	/* Create a new tmp file name */
	if(NULL == mktemp(tmpFileName->string)) {
		PrintError(FnName,
				tmpFileName->string,
				"Could not create a tmp file name",
				Exit,
				IllegalFileName);
	}

	/* Open a new file */
	if(!(fp = fopen(tmpFileName->string, "wb+"))) {
		PrintError(FnName,
				tmpFileName->string,
				"Could not open temporary file",
				Exit,
				OpenFileError);
	}

	return fp;
}

/* TODO */
void CloseTmpFile(FILE **fp,
		BString *tmpFileName)
{
	char *FnName="CloseTmpFile";

	/* Close the file */
	assert((*fp)!=NULL);
	fclose((*fp));
	(*fp)=NULL;

	/* Remove the file */
	assert(tmpFileName->string!=NULL);
	if(0!=remove(tmpFileName->string)) {
		PrintError(FnName,
				tmpFileName->string,
				"Could not delete temporary file",
				Exit,
				DeleteFileError);
	}

	/* Free file name */
	BStringFree(tmpFileName);
}

/* TODO */
void PrintPercentCompleteShort(double percent)
{
	/* back space the " percent complete" */
	fprintf(stderr, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	/* back space the "%3.2lf" */
	if(percent < 10.0) {
		/* Only really %1.2lf */
		fprintf(stderr, "\b\b\b\b");
	}
	else if(percent < 100.0) {
		/* Only really %2.2lf */
		fprintf(stderr, "\b\b\b\b\b");
	}
	else {
		fprintf(stderr, "\b\b\b\b\b\b");
	}
	fprintf(stderr, "%3.2lf percent complete", percent);
}

/* TODO */
void PrintPercentCompleteLong(double percent)
{
	/* back space the " percent complete" */
	fprintf(stderr, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	/* back space the "%3.3lf" */
	if(percent < 10.0) {
		/* Only really %1.3lf */
		fprintf(stderr, "\b\b\b\b\b");
	}
	else if(percent < 100.0) {
		/* Only really %2.3lf */
		fprintf(stderr, "\b\b\b\b\b\b");
	}
	else {
		fprintf(stderr, "\b\b\b\b\b\b\b");
	}
	fprintf(stderr, "%3.3lf percent complete", percent);
}

void PrintContigPos(FILE *fp,
		int32_t contig,
		int32_t position)
{
	int i;
	int contigLog10 = (int)floor(log10(contig));
	int positionLog10 = (int)floor(log10(position));

	contigLog10 = (contig < 1)?0:((int)floor(log10(contig)));
	positionLog10 = (position < 1)?0:((int)floor(log10(position)));

	assert(contigLog10 <= MAX_CONTIG_LOG_10);
	assert(positionLog10 <= MAX_POSITION_LOG_10);

	fprintf(fp, "\r[");
	for(i=0;i<(MAX_CONTIG_LOG_10-contigLog10);i++) {
		fprintf(fp, "-");
	}
	fprintf(fp, "%d,",
			contig);
	for(i=0;i<(MAX_POSITION_LOG_10-positionLog10);i++) {
		fprintf(fp, "-");
	}
	fprintf(fp, "%d",
			position);
	fprintf(fp, "]");
}

/* TODO */
int32_t UpdateRead(BString *read)
{
	int i;

	if(read->length <= 0) {
		/* Ignore zero length reads */
		return 0;
	}

	/* Update the read if possible to lower case, 
	 * if we encounter a base we do not recognize, 
	 * return 0 
	 * */
	for(i=0;i<read->length;i++) {
		switch(read->string[i]) {
			/* Do nothing for a, c, g, and t */
			case 'a':
			case 'c':
			case 'g':
			case 't':
			case 'n':
			case '0':
			case '1':
			case '2':
			case '3':
			case '4':
				break;
			case 'A':
				read->string[i] = 'a';
				break;
			case 'C':
				read->string[i] = 'c';
				break;
			case 'G':
				read->string[i] = 'g';
				break;
			case 'T':
				read->string[i] = 't';
				break;
			case 'N':
				read->string[i] = 'n';
				break;
			default:
				return 0;
				break;
		}
	}
	return 1;
}

/* TODO */
/* Debugging function */
int32_t CheckReadAgainstIndex(BIndex *index,
		BString *read)
{
	char *FnName = "CheckReadAgainstIndex";
	int i;
	for(i=0;i<index->width;i++) {
		assert(i<read->length);
		switch(CheckReadBase(read->string[i])) {
			case 0:
				return 0;
				break;
			case 1:
				break;
			default:
				PrintError(FnName,
						NULL,
						"Could not understand return value of CheckReadBase",
						Exit,
						OutOfRange);
				break;
		}
	}
	return 1;
}

/* TODO */
/* Debugging function */
int32_t CheckReadBase(int8_t base) 
{
	/* Do not include "n"s */
	switch(base) {
		case 'a':
		case 'c':
		case 'g':
		case 't':
			return 1;
			break;
		default:
			return 0;
			break;
	}
}

/* TODO */
/* Two bases */
/* If either of the two bases is an "N" or an "n", then
 * we return the color code 4 */
uint8_t ConvertBaseToColorSpace(uint8_t A, 
		uint8_t B)
{
	char *FnName = "ConvertBaseToColorSpace";
	int start=0;
	int by=0;
	int result=0;

	switch(A) {
		case 'A':
		case 'a':
			start = 0;
			by = 1;
			break;
		case 'C':
		case 'c':
			start = 1;
			by = -1;
			break;
		case 'G':
		case 'g':
			start = 2;
			by = 1;
			break;
		case 'T':
		case 't':
			start = 3;
			by = -1;
			break;
		case 'N':
		case 'n':
			return 4;
			break;
		default:
			fprintf(stderr, "\nA=(%c,%d)\n",
					A,
					(int)A);
			PrintError(FnName,
					"A",
					"Could not understand base",
					Exit,
					OutOfRange);
			break;
	}

	switch(B) {
		case 'A':
		case 'a':
			result = start;
			break;
		case 'C':
		case 'c':
			result = start + by;
			break;
		case 'G':
		case 'g':
			result = start + 2*by;
			break;
		case 'T':
		case 't':
			result = start + 3*by;
			break;
		case 'N':
		case 'n':
			return 4;
			break;
		default:
			fprintf(stderr, "\nB=(%c,%d)\n",
					B,
					(int)B);
			PrintError(FnName,
					"B",
					"Could not understand base",
					Exit,
					OutOfRange);
			break;
	}

	if(result < 0) {
		result =  ALPHABET_SIZE - ( (-1*result)% ALPHABET_SIZE);
	}
	else {
		result = (result%ALPHABET_SIZE);
	}

	return result;
}

/* TODO */
/* color must be an integer, and a base a character */
uint8_t ConvertBaseAndColor(uint8_t base, uint8_t color) 
{
	/* sneaky */
	return (uint8_t)DNA[(int)ConvertBaseToColorSpace(base, DNA[(int)color])];
}

/* TODO */
/* Include the first letter adaptor */
/* Does not reallocate memory */
void ConvertReadFromColorSpace(BString *read) 
{
	int i, index;

	/* Convert character numbers to 8-bit ints */
	for(i=0;i<read->length;i++) {
		switch(read->string[i]) {
			case '0':
				read->string[i] = 0;
				break;
			case '1':
				read->string[i] = 1;
				break;
			case '2':
				read->string[i] = 2;
				break;
			case '3':
				read->string[i] = 3;
				break;
			default:
				/* Ignore */
				break;
		}
	}

	for(i=0;i<read->length-1;i++) {
		if(0==i) {
			index = 0;
		}
		else {
			index = i-1; 
		}
		read->string[i] = ConvertBaseAndColor(read->string[index], read->string[i+1]); 
	}
	read->string[readLength-1] = '\0';
	BStringReallocate(read, read->length-1);
}

/* TODO */
/* NT read to color space */
void ConvertReadToColorSpace(BString *read) 
{
	char *FnName="ConvertReadToColorSpace";
	int i;
	BString tempRead;

	assert((*readLength) < SEQUENCE_LENGTH);

	/* Initialize */
	BStringInitialize(&tempRead);
	BStringAllocate(&tempread, read->length+1);
	tempRead->string[0] =  COLOR_SPACE_START_NT;
	tempRead->string[1] = ConvertBaseToColorSpace(tempRead->string[0], read->string[0]);

	/* Convert to colors represented as integers */
	for(i=1;i<(*readLength);i++) {
		tempRead->string[i+1] = ConvertBaseToColorSpace(read->string[i-1], read->string[i]);
	}

	/* Convert integers to characters */
	for(i=1;i<(*readLength)+1;i++) {
		assert(0<=tempRead->string[i] && tempRead->string[i] <= 4);
		tempRead->string[i] = COLORS[(int)(tempRead->string[i])];
	}
	tempRead->string[(*readLength)+1]='\0';

	/* Reallocate read to make sure */
	BStrinCopy(read, &tempRead);
}

/* TODO */
/* Takes in a NT read, converts to color space,
 * and then converts back to NT space using the
 * start NT */
void NormalizeRead(BString *read,
		char startNT)
{
	int i;
	char prevOldBase, prevNewBase;
	uint8_t tempColor;

	prevOldBase = startNT;
	prevNewBase = COLOR_SPACE_START_NT;
	for(i=0;i<read->length;i++) {
		/* Convert to color space using the old previous NT and current old NT */
		tempColor = ConvertBaseToColorSpace(prevOldBase, read->string[i]);
		prevOldBase = read->string[i];
		/* Convert to NT space but using the new previous NT and current color */
		read->string[i] = ConvertBaseAndColor(prevNewBase, tempColor);;
		prevNewBase = read->string[i];
	}
}

/* TODO */
void ConvertColorsToStorage(BString *colors)
{
	int i;
	for(i=0;i<colors->length;i++) {
		colors->string[i] = ConvertColorToStorage(colors->string[i]);
	}
}

/* TODO */
int8_t ConvertColorToStorage(int8_t c)
{
	switch(c) {
		case 0:
		case '0':
			c = 'A';
			break;
		case 1:
		case '1':
			c = 'C';
			break;
		case 2:
		case '2':
			c = 'G';
			break;
		case 3:
		case '3':
			c = 'T';
			break;
		default:
			c = 'N';
			break;
	}
	return c;
}

/* TODO */
void AdjustBounds(RGBinary *rg,
		int32_t *startContig,
		int32_t *startPos,
		int32_t *endContig,
		int32_t *endPos
		)
{
	char *FnName = "AdjustBounds";

	/* Adjust start and end based on reference genome */
	/* Adjust start */
	if((*startContig) <= 0) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Warning: startContig was less than zero.\n");
			fprintf(stderr, "Defaulting to contig=%d and position=%d.\n",
					1,
					1);
			fprintf(stderr, "%s", BREAK_LINE);
		}
		(*startContig) = 1;
		(*startPos) = 1;
	}
	else if((*startPos) <= 0) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Warning: startPos was less than zero.\n");
			fprintf(stderr, "Defaulting to position=%d.\n",
					1);
		}
		(*startPos) = 1;
	}

	/* Adjust end */
	if((*endContig) > rg->numContigs) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Warning: endContig was greater than the number of contigs in the reference genome.\n");
			fprintf(stderr, "Defaulting to reference genome's end contig=%d and position=%d.\n",
					rg->numContigs,
					rg->contigs[rg->numContigs-1].sequenceLength);
			fprintf(stderr, "%s", BREAK_LINE);
		}
		(*endContig) = rg->numContigs;
		(*endPos) = rg->contigs[rg->numContigs-1].sequenceLength;
	}
	else if((*endContig) <= rg->numContigs && 
			(*endPos) > rg->contigs[(*endContig)-1].sequenceLength) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Warning: endPos was greater than reference genome's contig %d end position.\n",
					(*endContig));
			fprintf(stderr, "Defaulting to reference genome's contig %d end position: %d.\n",
					(*endContig),
					rg->contigs[(*endContig)-1].sequenceLength);
			fprintf(stderr, "%s", BREAK_LINE);
		}
		(*endPos) = rg->contigs[(*endContig)-1].sequenceLength;
	}

	/* Check that the start and end bounds are ok */
	if((*startContig) > (*endContig)) {
		PrintError(FnName,
				NULL,
				"The start contig is greater than the end contig",
				Exit,
				OutOfRange);
	}
	else if((*startContig) == (*endContig) &&
			(*startPos) > (*endPos)) {
		PrintError(FnName,
				NULL,
				"The start position is greater than the end position on the same contig",
				Exit,
				OutOfRange);
	}
}

/* TODO */
int32_t WillGenerateValidKey(BIndex *index,
		BString *read) 
{
	int i;

	for(i=0;i<index->width;i++) {
		if(i >= read->length ||
				(1 == index->mask[i] && 1==RGBinaryIsBaseN(read->string[i]))) {
			return 0;
		}
	}
	return 1;
}

/* TODO */
int ValidateFileName(BString *Name)
{
	/* 
	 *        Checking that strings are good: FileName = [a-zA-Z_0-9][a-zA-Z0-9-.]+
	 *               FileName can start with only [a-zA-Z_0-9]
	 *                      */

	char *ptr=Name->string;
	int counter=0;
	/*   fprintf(stderr, "Validating FileName %s with length %d\n", ptr, strlen(Name));  */

	assert(ptr!=0);

	while(*ptr) {
		if((isalnum(*ptr) || (*ptr=='_') || (*ptr=='+') ||
					((*ptr=='.') /* && (counter>0)*/) || /* FileNames can't start  with . or - */
					((*ptr=='/')) || /* Make sure that we can navigate through folders */
					((*ptr=='-') && (counter>0)))) {
			ptr++;
			counter++;
		}
		else return 0;
	}
	return 1;
}
