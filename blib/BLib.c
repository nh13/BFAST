#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "BLibDefinitions.h"
#include "RGIndex.h"
#include "BError.h"
#include "BLib.h"

char DNA[5] = "ACGTN";

/* TODO */
char ToLower(char a) 
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
char ToUpper(char a)
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
		default:
			return a;
	}
}

/* TODO */
void ReverseRead(char *s,
		char *r,
		int length)
{       
	int i;
	/* Get reverse */
	for(i=length-1;i>=0;i--) {
		r[i] = s[length-1-i];
	}
	r[length]='\0';
}

/* TODO */
void GetReverseComplimentAnyCase(char *s,
		char *r,
		int length)
{       
	int i;
	/* Get reverse compliment sequence */
	for(i=length-1;i>=0;i--) {
		r[i] = GetReverseComplimentAnyCaseBase(s[length-1-i]);
	}
	r[length]='\0';
}

char GetReverseComplimentAnyCaseBase(char a) {
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
			PrintError("GetReverseComplimentAnyCaseBase",
					NULL,
					"Could not understand sequence base",
					Exit,
					OutOfRange);
			break;
	}
	PrintError("GetReverseComplimentAnyCaseBase",
			NULL,
			"Control should not reach here",
			Exit,
			OutOfRange);
	return '0';
}

/* TODO */
int ValidateBasePair(char c) {
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

int IsAPowerOfTwo(unsigned int a) {
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

char TransformFromIUPAC(char a) 
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

void CheckRGIndexes(char **mainFileNames,
		int numMainFileNames,
		char **secondaryFileNames,
		int numSecondaryFileNames,
		int binaryInput,
		int32_t *startChr,
		int32_t *startPos,
		int32_t *endChr,
		int32_t *endPos,
		int32_t colorSpace)
{
	int i;
	int32_t mainStartChr, mainStartPos, mainEndChr, mainEndPos;
	int32_t secondaryStartChr, secondaryStartPos, secondaryEndChr, secondaryEndPos;
	int32_t mainColorSpace, secondaryColorSpace;
	mainStartChr = mainStartPos = mainEndChr = mainEndPos = 0;
	secondaryStartChr = secondaryStartPos = secondaryEndChr = secondaryEndPos = 0;

	RGIndex tempIndex;
	FILE *fp;

	/* Read in main indexes */
	for(i=0;i<numMainFileNames;i++) {
		/* Open file */
		if((fp=fopen(mainFileNames[i], "r"))==0) {
			PrintError("CheckRGIndexes",
					mainFileNames[i],
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}

		/* Get the header */
		RGIndexReadHeader(fp, &tempIndex, binaryInput); 

		assert(tempIndex.startChr < tempIndex.endChr ||
				(tempIndex.startChr == tempIndex.endChr && tempIndex.startPos <= tempIndex.endPos));

		if(i==0) {
			mainStartChr = tempIndex.startChr;
			mainStartPos = tempIndex.startPos;
			mainEndChr = tempIndex.endChr;
			mainEndPos = tempIndex.endPos;
			mainColorSpace = tempIndex.colorSpace;
		}
		else {
			/* Update bounds if necessary */
			if(tempIndex.startChr < mainStartChr ||
					(tempIndex.startChr == mainStartChr && tempIndex.startPos < mainStartPos)) {
				mainStartChr = tempIndex.startChr;
				mainStartPos = tempIndex.startPos;
			}
			if(tempIndex.endChr > mainEndChr ||
					(tempIndex.endChr == mainEndChr && tempIndex.endPos > mainEndPos)) {
				mainEndChr = tempIndex.endChr;
				mainEndPos = tempIndex.endPos;
			}
			assert(mainColorSpace == tempIndex.colorSpace);
		}

		/* Close file */
		fclose(fp);
	}
	/* Read in secondary indexes */
	for(i=0;i<numSecondaryFileNames;i++) {
		/* Open file */
		if((fp=fopen(secondaryFileNames[i], "r"))==0) {
			PrintError("CheckRGIndexes",
					"secondaryFileNames[i]",
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}

		/* Get the header */
		RGIndexReadHeader(fp, &tempIndex, binaryInput); 

		assert(tempIndex.startChr < tempIndex.endChr ||
				(tempIndex.startChr == tempIndex.endChr && tempIndex.startPos <= tempIndex.endPos));

		if(i==0) {
			secondaryStartChr = tempIndex.startChr;
			secondaryStartPos = tempIndex.startPos;
			secondaryEndChr = tempIndex.endChr;
			secondaryEndPos = tempIndex.endPos;
			secondaryColorSpace = tempIndex.colorSpace;
		}
		else {
			/* Update bounds if necessary */
			if(tempIndex.startChr < secondaryStartChr ||
					(tempIndex.startChr == secondaryStartChr && tempIndex.startPos < secondaryStartPos)) {
				secondaryStartChr = tempIndex.startChr;
				secondaryStartPos = tempIndex.startPos;
			}
			if(tempIndex.endChr > secondaryEndChr ||
					(tempIndex.endChr == secondaryEndChr && tempIndex.endPos > secondaryEndPos)) {
				secondaryEndChr = tempIndex.endChr;
				secondaryEndPos = tempIndex.endPos;
			}
			assert(secondaryColorSpace == tempIndex.colorSpace);
		}

		/* Close file */
		fclose(fp);
	}

	/* Check the bounds between main and secondary indexes */
	if(mainStartChr != secondaryStartChr ||
			mainStartPos != secondaryStartPos ||
			mainEndChr != secondaryEndChr ||
			mainEndPos != secondaryEndPos ||
			mainColorSpace != secondaryColorSpace) {
		PrintError("CheckRGIndexes",
				NULL,
				"The ranges between main and secondary indexes differ",
				Exit,
				OutOfRange);
	}

	(*startChr) = mainStartChr;
	(*startPos) = mainStartPos;
	(*endChr) = mainEndChr;
	(*endPos) = mainEndPos;

	assert(mainColorSpace == colorSpace);
	assert(secondaryColorSpace == colorSpace);
}

/* TODO */
FILE *OpenTmpFile(char *tmpDir,
		char **tmpFileName)
{
	char *FnName = "OpenTmpFile";
	FILE *fp;

	/* Allocate memory */
	(*tmpFileName) = malloc(sizeof(char)*MAX_FILENAME_LENGTH);
	if(NULL == (*tmpFileName)) {
		PrintError(FnName,
				"tmpFileName",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	/* Create the templated */
	/* Copy over tmp directory */
	strcpy((*tmpFileName), tmpDir);
	/* Copy over the tmp name */
	strcat((*tmpFileName), BFAST_TMP_TEMPLATE);

	/* Create a new tmp file name */
	if(NULL == mktemp((*tmpFileName))) {
		PrintError(FnName,
				(*tmpFileName),
				"Could not create a tmp file name",
				Exit,
				IllegalFileName);
	}

	/* Open a new file */
	if(!(fp = fopen((*tmpFileName), "wb+"))) {
		PrintError(FnName,
				(*tmpFileName),
				"Could not open temporary file",
				Exit,
				OpenFileError);
	}

	return fp;
}

/* TODO */
void CloseTmpFile(FILE **fp,
		char **tmpFileName)
{
	char *FnName="CloseTmpFile";

	/* Close the file */
	assert((*fp)!=NULL);
	fclose((*fp));
	(*fp)=NULL;

	/* Remove the file */
	assert((*tmpFileName)!=NULL);
	if(0!=remove((*tmpFileName))) {
		PrintError(FnName,
				(*tmpFileName),
				"Could not delete temporary file",
				Exit,
				DeleteFileError);
	}

	/* Free file name */
	free((*tmpFileName));
	(*tmpFileName) = NULL;
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

/* TODO */
int UpdateRead(char *read, int readLength) 
{
	int i;

	if(readLength <= 0) {
		/* Ignore zero length reads */
		return 0;
	}

	/* Update the read if possible to lower case, 
	 * if we encounter a base we do not recognize, 
	 * return 0 
	 * */
	for(i=0;i<readLength;i++) {
		switch(read[i]) {
			/* Do nothing for a, c, g, and t */
			case 'a':
			case 'c':
			case 'g':
			case 't':
			case '0':
			case '1':
			case '2':
			case '3':
				break;
			case 'A':
				read[i] = 'a';
				break;
			case 'C':
				read[i] = 'c';
				break;
			case 'G':
				read[i] = 'g';
				break;
			case 'T':
				read[i] = 't';
				break;
			default:
				return 0;
				break;
		}
	}
	return 1;
}

/* TODO */
int CheckReadAgainstIndex(RGIndex *index,
		char *read,
		int readLength)
{
	char *FnName = "CheckReadAgainstIndex";
	int curPos, curTile, curTilePos;
	for(curPos = 0, curTile = 0;curTile < index->numTiles;curTile++) {
		for(curTilePos=0;curTilePos < index->tileLengths[curTile];curTilePos++) {
			switch(CheckReadBase(read[curPos])) {
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
			/* Update position */
			curPos++;
		}
		if(curTile < index->numTiles-1) {
			curPos += index->gaps[curTile];
		}
	}
	return 1;
}

/* TODO */
int CheckReadBase(char base) 
{
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
int ConvertReadFromColorSpace(char *read,
		int readLength)
{
	int i, index;

	/* HERE A2 */
	/*
	   fprintf(stderr, "HERE A2\nread=%s\n",
	   read);
	   fprintf(stderr, "readLength=%d\n",
	   (int)strlen(read));
	   fprintf(stderr, "act length=%d\n",
	   readLength);
	   */

	/* Convert character numbers to 8-bit ints */
	for(i=0;i<readLength;i++) {
		switch(read[i]) {
			case '0':
				read[i] = 0;
				break;
			case '1':
				read[i] = 1;
				break;
			case '2':
				read[i] = 2;
				break;
			case '3':
				read[i] = 3;
				break;
			default:
				/* Ignore */
				break;
		}
	}

	for(i=0;i<readLength-1;i++) { 
		if(0==i) {
			index = 0;
		}
		else {
			index = i-1; 
		}
		read[i] = ConvertBaseAndColor(read[index], read[i+1]); 
	}
	read[readLength-1] = '\0';
	readLength--;

	return readLength;
}

/* TODO */
/* Must reallocate memory */
void ConvertReadToColorSpace(char **read,
		int *readLength)
{
	char *FnName="ConvertReadToColorSpace";
	int i;
	char tempRead[SEQUENCE_LENGTH]="\0";

	/* Initialize */
	tempRead[0] =  COLOR_SPACE_START_NT;

	/* Go through read and convert to colors */
	for(i=0;i<(*readLength);i++) {
		tempRead[i+1] = ConvertBaseToColorSpace(tempRead[i], (*read)[i]);
		switch(tempRead[i+1]) {
			case 0:
				tempRead[i+1] = '0';
				break;
			case 1:
				tempRead[i+1] = '1';
				break;
			case 2:
				tempRead[i+1] = '2';
				break;
			case 3:
				tempRead[i+1] = '3';
				break;
			default:
				PrintError(FnName,
						"result",
						"Could not understand result",
						Exit,
						OutOfRange);
		}
	}
	tempRead[(*readLength)+1]='\0';

	/* Reallocate read to make sure */
	(*read) = realloc((*read), sizeof(char)*SEQUENCE_LENGTH);
	if((*read)==NULL) {
		PrintError(FnName,
				"(*read)",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	strcpy((*read), tempRead);
	(*readLength)++;
}

/* TODO */
/* Takes in a NT read, converts to color space,
 * and then converts back to NT space using the
 * start NT */
void NormalizeRead(char **read,
		int *readLength,
		char startNT)
{
	int i;
	char prevOldBase, prevNewBase;
	uint8_t tempColor;

	/* HERE 23 */
	/*
	   fprintf(stderr, "Before=%s\n",
	   (*read));
	   */

	prevOldBase = startNT;
	prevNewBase = COLOR_SPACE_START_NT;
	for(i=0;i<(*readLength);i++) {
		/* Convert to color space using the old previous NT and current old NT */
		tempColor = ConvertBaseToColorSpace(prevOldBase, (*read)[i]);
		prevOldBase = (*read)[i];
		/* Convert to NT space but using the new previous NT and current color */
		(*read)[i] = ConvertBaseAndColor(prevNewBase, tempColor);;
		prevNewBase = (*read)[i];
	}

	/* HERE 24 */
	/*
	   fprintf(stderr, "\n");
	   fprintf(stderr, "After=%s\n", 
	   (*read));
	   */
}

/* TODO */
void ConvertColorsToStorage(char *colors, int length)
{
	int i;
	for(i=0;i<length;i++) {
		colors[i] = ConvertColorToStorage(colors[i]);
	}
}

/* TODO */
char ConvertColorToStorage(char c)
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
		int32_t *startChr,
		int32_t *startPos,
		int32_t *endChr,
		int32_t *endPos
		)
{
	char *FnName = "AdjustBounds";

	/* Adjust start and end based on reference genome */
	/* Adjust start */
	if((*startChr) < rg->startChr) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Warning: startChr was less than reference genome's start chromosome.\n");
			fprintf(stderr, "Defaulting to reference genome's start chromosome and position: chr%d:%d.\n",
					rg->startChr,
					rg->startPos);
		}
		(*startChr) = rg->startChr;
		(*startPos) = rg->startPos;
	}
	else if((*startChr) == rg->startChr &&
			(*startPos) < rg->startPos) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Warning: startPos was less than reference genome's start position.\n");
			fprintf(stderr, "Defaulting to reference genome's start position: %d.\n",
					rg->startPos);
		}
		(*startPos) = rg->startPos;
	}
	/* Adjust end */
	if((*endChr) > rg->endChr) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Warning: endChr was greater than reference genome's end chromosome.\n");
			fprintf(stderr, "Defaulting to reference genome's end chromosome and position: chr%d:%d.\n",
					rg->endChr,
					rg->endPos);
		}
		(*endChr) = rg->endChr;
		(*endPos) = rg->endPos;
	}
	else if((*endChr) == rg->endChr &&
			(*endPos) > rg->endPos) {
		if(VERBOSE >= 0) {
			fprintf(stderr, "%s", BREAK_LINE);
			fprintf(stderr, "Warning: endPos was greater than reference genome's end position.\n");
			fprintf(stderr, "Defaulting to reference genome's end position: %d.\n",
					rg->endPos);
		}
		(*endPos) = rg->endPos;
	}
	/* Check that the start and end bounds are ok */
	if((*startChr) > (*endChr)) {
		PrintError(FnName,
				NULL,
				"The start chromosome is greater than the end chromosome",
				Exit,
				OutOfRange);
	}
	else if((*startChr) == (*endChr) &&
			(*startPos) > (*endPos)) {
		PrintError(FnName,
				NULL,
				"The start position is greater than the end position on the same chromosome",
				Exit,
				OutOfRange);
	}
}

