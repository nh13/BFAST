#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "BLib.h"
#include "BLibDefinitions.h"
#include "BError.h"
#include "RGMatches.h"
#include "AlignedEntry.h"
#include "ScoringMatrix.h"
#include "Align.h"
#include "AlignColorSpace.h"

// Do we need to enforce NEGATIVE_INFINITY as a lower bound when adding?

/* TODO */
int32_t AlignColorSpaceUngapped(char *read,
		char *mask,
		int readLength,
		char *reference,
		int referenceLength,
		int unconstrained,
		ScoringMatrix *sm,
		AlignedEntry *a,
		int offset,
		int32_t position,
		char strand)
{
	/* read goes on the rows, reference on the columns */
	char *FnName = "AlignColorSpaceUngapped";
	int i, j, k, l;

	int offsetAligned=-1;
	int32_t prevScore[ALPHABET_SIZE+1];
	int prevNT[ALPHABET_SIZE+1][SEQUENCE_LENGTH];
	int32_t maxScore = NEGATIVE_INFINITY;
	int maxNT[SEQUENCE_LENGTH];
	char DNA[ALPHABET_SIZE+1] = "ACGTN";
	char readAligned[SEQUENCE_LENGTH]="\0";
	char referenceAligned[SEQUENCE_LENGTH]="\0";
	char colorErrorAligned[SEQUENCE_LENGTH]="\0";
	int32_t alphabetSize = ALPHABET_SIZE+1;

	assert(readLength <= referenceLength);
	assert(2*offset + readLength <= referenceLength); 

	alphabetSize = AlignColorSpaceGetAlphabetSize(read, readLength, reference, referenceLength);

	for(i=offset;i<referenceLength-readLength-offset+1;i++) { /* Starting position */
		/* Initialize */
		for(j=0;j<alphabetSize;j++) {
			if(DNA[j] == COLOR_SPACE_START_NT) { 
				prevScore[j] = 0.0;
			}
			else {
			}
		}
		for(j=0;j<readLength;j++) { /* Position in the alignment */
			char curColor;
			int32_t curReadBaseInt = -1;

			/* Get the current color for the read */
			if(0 == ConvertBaseToColorSpace((j==0)?COLOR_SPACE_START_NT:read[j-1], read[j], &curColor)) {
				PrintError(FnName, "curColor", "Could not convert base to color space", Exit, OutOfRange);
			}
			if(Constrained == unconstrained) {
				curReadBaseInt = BaseToInt(read[j]);
			}
			int32_t nextScore[ALPHABET_SIZE+1];
			char nextNT[ALPHABET_SIZE+1];
			for(k=0;k<alphabetSize;k++) { /* To NT */

				/* Get the best score to this NT */
				int32_t bestScore = NEGATIVE_INFINITY;
				int bestNT=-1;
				char bestColor = 'X';

				if(Constrained == unconstrained && '1' == mask[j]) { // If we are to use the constraint and it exists
					AlignColorSpaceUngappedGetBest(sm, 
							prevScore[curReadBaseInt],
							curColor,
							reference[i+j],
							k,
							curReadBaseInt, // Use the current base (as an integer)
							alphabetSize,
							&bestScore,
							&bestNT,
							&bestColor);
				}
				else { // Ignore constraint, go through all possible transitions
					for(l=0;l<alphabetSize;l++) { /* From NT */
						AlignColorSpaceUngappedGetBest(sm, 
								prevScore[l],
								curColor,
								reference[i+j],
								k,
								l,
								alphabetSize,
								&bestScore,
								&bestNT,
								&bestColor);
					}
				}
				nextScore[k] = bestScore;
				nextNT[k] = bestNT;
			}

			for(k=0;k<alphabetSize;k++) { /* To NT */
				prevScore[k] = nextScore[k];
				prevNT[k][j] = nextNT[k];
			}
		}
		/* Check if the score is better than the max */
		k=0;
		for(j=0;j<alphabetSize;j++) { /* To NT */
			if(prevScore[k] < prevScore[j]) {
				k=j;
			}
		}
		if(maxScore < prevScore[k]) {
			maxScore = prevScore[k];
			/* TO GET COLORS WE NEED TO BACKTRACK */
			for(j=readLength-1;0<=j;j--) {
				maxNT[j] = k;
				k=prevNT[k][j];
			}
			offsetAligned = i;
		}
	}

	if(NEGATIVE_INFINITY < maxScore) {
		for(i=0;i<readLength;i++) {
			char c[2];
			readAligned[i] = DNA[maxNT[i]];
			referenceAligned[i] = reference[i+offsetAligned];
			ConvertBaseToColorSpace((i==0)?COLOR_SPACE_START_NT:read[i-1],
					read[i],
					&c[0]);
			ConvertBaseToColorSpace((i==0)?COLOR_SPACE_START_NT:readAligned[i-1],
					readAligned[i],
					&c[1]);
			colorErrorAligned[i] = ConvertIntColorToCharColor((c[0] == c[1])?GAP:c[0]); /* Keep original color */
		}
		readAligned[readLength]=referenceAligned[readLength]=colorErrorAligned[readLength]='\0';

		/* Copy over */
		AlignedEntryUpdateAlignment(a,
				(FORWARD==strand) ? (position + offsetAligned) : (position + referenceLength - readLength - offsetAligned),
				maxScore,
				readLength,
				readLength,
				readAligned,
				referenceAligned,
				colorErrorAligned);
		return 1;
	}
	else {
		/* For a case where this actually occurs, think reads at the beginning 
		 * or end of a contig, with adaptor sequence!
		 * */
		return 0;
	}
}

void AlignColorSpaceUngappedGetBest(
		ScoringMatrix *sm,
		int32_t curScore, // previous score (prevScore[l])
		char curColor, // observed color
		char refBase, // reference base (reference[i+j])
		int32_t k, // To NT
		int32_t l, // From NT
		int32_t alphabetSize,
		int32_t *bestScore,
		int32_t *bestNT,
		char *bestColor)
{
	char *FnName="AlignColorSpaceUngappedGetBest";

	char convertedColor='X';
	/* Get color */
	if(0 == ConvertBaseToColorSpace(DNA[l], DNA[k], &convertedColor)) {
		PrintError(FnName, "convertedColor", "Could not convert base to color space", Exit, OutOfRange);
	}
	/* Add score for color error, if any */
	curScore += ScoringMatrixGetColorScore(curColor,
			convertedColor,
			sm);
	/* Add score for NT */
	curScore += ScoringMatrixGetNTScore(refBase, DNA[k], sm);

	if(curScore < NEGATIVE_INFINITY/2) {
		curScore = NEGATIVE_INFINITY;
	}

	if((*bestScore)< curScore) {
		(*bestScore) = curScore;
		(*bestNT) = l;
		(*bestColor) = convertedColor;
	}
}

void AlignColorSpaceGappedBounded(char *read,
		int readLength,
		char *reference,
		int referenceLength,
		ScoringMatrix *sm,
		AlignedEntry *a,
		AlignMatrix *matrix,
		int32_t position,
		char strand,
		int32_t maxH,
		int32_t maxV)
{
	char *FnName = "AlignColorSpaceGappedBounded";
	int i, j;
	int alphabetSize=ALPHABET_SIZE;

	assert(0 < readLength);
	assert(0 < referenceLength);

	alphabetSize = AlignColorSpaceGetAlphabetSize(read, readLength, reference, referenceLength);

	AlignColorSpaceInitializeAtStart(read, matrix, sm, readLength, referenceLength, alphabetSize);

	/* Fill in the matrix according to the recursive rules */
	for(i=0;i<readLength;i++) { /* read/rows */
		/* Get the current color */
		char curColor;
		if(0 == ConvertBaseToColorSpace((i==0)?COLOR_SPACE_START_NT:read[i-1], read[i], &curColor)) {
			PrintError(FnName, "curColor", "Could not convert base to color space", Exit, OutOfRange);
		}

		for(j=GETMAX(0, i - maxV);
				j <= GETMIN(referenceLength-1, referenceLength - (readLength - maxH) + i);
				j++) { /* reference/columns */
			assert(i-maxV <= j && j <= referenceLength - (readLength - maxH) + i);
			AlignColorSpaceFillInCell(read, readLength, reference, referenceLength, sm, matrix, i, j, curColor, maxH, maxV, alphabetSize);
		}
	}

	AlignColorSpaceRecoverAlignmentFromMatrix(a, matrix, read, readLength, reference, referenceLength, readLength - maxV, position, strand, alphabetSize, 0);
}

void AlignColorSpaceGappedConstrained(char *read,
		char *mask,
		int readLength,
		char *reference,
		int referenceLength,
		ScoringMatrix *sm,
		AlignedEntry *a,
		AlignMatrix *matrix,
		int32_t referenceOffset,
		int32_t readOffset,
		int32_t position,
		char strand)
{
	char *FnName = "AlignColorSpaceGappedConstrained";
	int i, j, k;
	int alphabetSize=ALPHABET_SIZE;
	int32_t endRowStepOne, endColStepOne, endRowStepTwo, endColStepTwo;

	assert(0 < readLength);
	assert(0 < referenceLength);

	alphabetSize = AlignColorSpaceGetAlphabetSize(read, readLength, reference, referenceLength);

	/* Get where to transition */
	endRowStepOne = endColStepOne = endRowStepTwo = endColStepTwo = -1;
	i=0;
	while(i<readLength) {
		if('1' == mask[i]) {
			endRowStepOne=i;
			endColStepOne=referenceOffset+i;
			break;
		}
		i++;
	}
	i=readLength;
	while(0<=i) {
		if('1' == mask[i]) {
			endRowStepTwo=i;
			endColStepTwo=referenceOffset+i;
			break;
		}
		i--;
	}

	/* Adjust based off of the read offset */
	if(0 < readOffset) {
		endRowStepOne += readOffset;
		endRowStepTwo += readOffset;
	}

	assert(0 <= endRowStepOne && 0 <= endColStepOne);
	assert(0 <= endRowStepTwo && 0 <= endColStepTwo);


	// HERE
	/*
	fprintf(stderr, "%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
			FnName,
			reference,
			read,
			mask,
			reference + referenceOffset,
			read + readOffset,
			mask + readOffset);
	fprintf(stderr, "[%d,%d]\n",
			referenceOffset,
			readOffset);
	fprintf(stderr, "[%d,%d]\n[%d,%d]\n",
			endRowStepOne, endColStepOne,
			endRowStepTwo, endColStepTwo);
			*/

	/* Step 1 - upper left */
	AlignColorSpaceInitializeAtStart(read, matrix, sm, endRowStepOne, endColStepTwo, alphabetSize);
	for(i=0;i<endRowStepOne;i++) { /* read/rows */
		/* Get the current color for the read */
		char curColor;
		if(0 == ConvertBaseToColorSpace((i==0)?COLOR_SPACE_START_NT:read[i-1], read[i], &curColor)) {
			PrintError(FnName, "curColor - step 1", "Could not convert base to color space", Exit, OutOfRange);
		}
		for(j=0;j<endColStepOne;j++) { /* reference/columns */
			AlignColorSpaceFillInCell(read, readLength, reference, referenceLength, sm, matrix, i, j, curColor, INT_MAX, INT_MAX, alphabetSize);
		}
	}

	/* Step 2 - align along the mask */

	// Must consider ins, del, and match on first "color"

	/* Step 2 - align along the mask */
	// HERE XY
	// TODO
	for(i=endRowStepOne,j=endColStepOne;
			i<endRowStepTwo && j<endColStepTwo;
			i++,j++) {
		/*
		fprintf(stderr, "%s: filling in [%d,%d]\n",
				FnName,
				i+1, j+1);
				*/
		char curColor, curReferenceColor;
		if(0 == ConvertBaseToColorSpace((i==0)?COLOR_SPACE_START_NT:read[i-1], read[i], &curColor)) {
			PrintError(FnName, "curColor - step 2", "Could not convert base to color space", Exit, OutOfRange);
		}
		/* Get the current color for the reference */
		assert(0 < j); // this could be a problem since we need more reference
		if(0 == ConvertBaseToColorSpace(reference[j-1], reference[j], &curReferenceColor)) {
			PrintError(FnName, "curReferenceColor", "Could not convert base to color space", Exit, OutOfRange);
		}
		/*
		fprintf(stderr, "[%c,%c]\n",
				"01234"[(int)curColor],
				"01234"[(int)curReferenceColor]);
				*/
		/* Check the colors match */
		if('1' == mask[i] && curColor != curReferenceColor) {
			char c;
			fprintf(stderr, "%s\n", reference + referenceOffset);
			for(i=referenceOffset;i<referenceLength;i++) {
				assert(0 < i);
				if(0 == ConvertBaseToColorSpace(reference[i-1], reference[i], &c)) {
					PrintError(FnName, "c", "Could not convert base to color space", Exit, OutOfRange);
				}
				fprintf(stderr, "%c", "01234"[(int)c]);
			}
			fprintf(stderr, "\n%s\n", mask + readOffset);
			for(i=readOffset;i<readLength;i++) {
				if(0 == ConvertBaseToColorSpace((i==0)?COLOR_SPACE_START_NT:read[i-1], read[i], &c)) {
					PrintError(FnName, "c", "Could not convert base to color space", Exit, OutOfRange);
				}
				fprintf(stderr, "%c", "01234"[(int)c]);
			}
			fprintf(stderr, "\n");
			fprintf(stderr, "%s\n", read + readOffset);
			fprintf(stderr, "%d:%c\n", position, strand);
			PrintError(FnName, NULL, "read and reference did not match", Exit, OutOfRange);
		}
		for(k=0;k<alphabetSize;k++) { /* To NT */
			char fromNT;
			int32_t fromNTInt;

			/* Get the from base */
			if(0 == ConvertBaseAndColor(DNA[k], curColor, &fromNT)) { 
				fprintf(stderr, "\n[%c,%c]\n", DNA[k], "01234"[(int)curColor]);
				PrintError(FnName, "fromNT", "Could not convert base and color to base", Exit, OutOfRange);
			}
			fromNTInt=BaseToInt(fromNT);

			// Add color score
			matrix->cells[i+1][j+1].s.score[k] = ScoringMatrixGetColorScore(curColor, curColor, sm);
			/* Add score for NT */
			//Assertion above handles this ... assert(0 < j); // this could be a problem
			matrix->cells[i+1][j+1].s.score[k] += ScoringMatrixGetNTScore(reference[j], DNA[k], sm);
			matrix->cells[i+1][j+1].s.from[k] = fromNTInt + 1 + (ALPHABET_SIZE + 1); 
			matrix->cells[i+1][j+1].s.length[k] = matrix->cells[i][j].s.length[fromNTInt] + 1;
			matrix->cells[i+1][j+1].s.colorError[k] = GAP;
			
			// Consider from an indel on the first extension
			if(i == endRowStepOne && j == endColStepOne) {
				int32_t curScore = ScoringMatrixGetColorScore(curColor, curColor, sm);
				curScore += ScoringMatrixGetNTScore(reference[j], DNA[k], sm);

				/* From Horizontal - Deletion */
				if(matrix->cells[i+1][j+1].s.score[k] < curScore + matrix->cells[i][j].h.score[fromNTInt]) { 
					matrix->cells[i+1][j+1].s.score[k] = curScore + matrix->cells[i][j].h.score[fromNTInt];
					matrix->cells[i+1][j+1].s.from[k] = fromNTInt + 1;
					matrix->cells[i+1][j+1].s.length[k] = matrix->cells[i][j].h.length[fromNTInt] + 1;
					matrix->cells[i+1][j+1].s.colorError[k] = GAP;
				}

				/* From Vertical - Insertion */
				if(matrix->cells[i+1][j+1].s.score[k] < curScore + matrix->cells[i][j].v.score[fromNTInt]) { 
					matrix->cells[i+1][j+1].s.score[k] = curScore + matrix->cells[i][j].v.score[fromNTInt];
					matrix->cells[i+1][j+1].s.from[k] = fromNTInt + 1 + 2*(ALPHABET_SIZE + 1);
					matrix->cells[i+1][j+1].s.length[k] = matrix->cells[i][j].v.length[fromNTInt] + 1;
					matrix->cells[i+1][j+1].s.colorError[k] = GAP;
				}
			}
		}
	}
	for(k=0;k<alphabetSize;k++) {
		assert(1 + ALPHABET_SIZE < matrix->cells[endRowStepTwo][endColStepTwo].s.from[k] &&
				matrix->cells[endRowStepTwo][endColStepTwo].s.from[k] <= 2*(ALPHABET_SIZE + 1));
	}

	/* Step 3 - lower right */
	AlignColorSpaceInitializeToExtend(read, matrix, sm, readLength, referenceLength, endRowStepTwo, endColStepTwo, alphabetSize);
	// Note: we ignore any cells on row==endRowStepTwo or col==endRowStepTwo
	// since we assumed they were filled in by the previous re-initialization
	for(i=endRowStepTwo;i<readLength;i++) { /* read/rows */
		/* Get the current color for the read */
		char curColor;
		if(0 == ConvertBaseToColorSpace((i==0)?COLOR_SPACE_START_NT:read[i-1], read[i], &curColor)) {
			fprintf(stderr, "\n");
			fprintf(stderr, "readLength=%d\n", readLength);
			fprintf(stderr, "i=%d\n", i);
			fprintf(stderr, "%s\n", read);
			PrintError(FnName, "curColor - step 3", "Could not convert base to color space", Exit, OutOfRange);
		}
		for(j=endColStepTwo;j<referenceLength;j++) { /* reference/columns */
			AlignColorSpaceFillInCell(read, readLength, reference, referenceLength, sm, matrix, i, j, curColor, INT_MAX, INT_MAX, alphabetSize);
		}
	}

	/* Step 4 - recover alignment */
	AlignColorSpaceRecoverAlignmentFromMatrix(a, matrix, read, readLength, reference, referenceLength, endColStepTwo, position, strand, alphabetSize, 0);
}

void AlignColorSpaceRecoverAlignmentFromMatrix(AlignedEntry *a,
		AlignMatrix *matrix,
		char *read,
		int readLength,
		char *reference,
		int referenceLength,
		int toExclude,
		int32_t position,
		char strand,
		int alphabetSize,
		int debug)
{
	char *FnName="AlignColorSpaceRecoverAlignmentFromMatrix";
	int curRow, curCol, startRow, startCol, startCell;
	char curReadBase;
	int nextRow, nextCol, nextCell, nextFrom;
	char nextReadBase;
	int curFrom=-1;
	double maxScore;
	int i, j;
	int offset;
	char readAligned[SEQUENCE_LENGTH]="\0";
	char referenceAligned[SEQUENCE_LENGTH]="\0";
	char colorErrorAligned[SEQUENCE_LENGTH]="\0";
	int32_t referenceLengthAligned, length;

	curReadBase = nextReadBase = 'X';
	nextRow = nextCol = nextCell = -1;

	/* Get the best alignment.  We can find the best score in the last row and then
	 * trace back.  We choose the best score from the last row since we want to 
	 * align the read completely and only locally to the reference. */
	startRow=-1;
	startCol=-1;
	startCell=-1;
	maxScore = NEGATIVE_INFINITY-1;
	//fprintf(stderr, "Checking\n");
	for(i=toExclude+1;i<referenceLength+1;i++) {
		for(j=0;j<alphabetSize;j++) {
			/*
			fprintf(stderr, "Checking [%d,%d,%d]\n",
					readLength, i, j);
					*/
			/* Don't end with a Deletion in the read */

			/* End with a Match/Mismatch */
			if(maxScore < matrix->cells[readLength][i].s.score[j]) {
				maxScore = matrix->cells[readLength][i].s.score[j];
				startRow = readLength;
				startCol = i;
				startCell = j + 1 + (ALPHABET_SIZE + 1);
				/*
				fprintf(stderr, "Ending with a mismatch\n");
				fprintf(stderr, "length=%d\n",
						matrix->cells[readLength][i].s.length[j]);
						*/
			}

			/* End with an Insertion */
			if(maxScore < matrix->cells[readLength][i].v.score[j]) {
				maxScore = matrix->cells[readLength][i].v.score[j];
				startRow = readLength;
				startCol = i;
				startCell = j + 1 + 2*(ALPHABET_SIZE + 1);
				/*
				fprintf(stderr, "Ending with an insertion\n");
				*/
			}
		}
	}
	assert(startRow >= 0 && startCol >= 0 && startCell >= 0);

	/* Initialize variables for the loop */
	curRow=startRow;
	curCol=startCol;
	curFrom=startCell;

	// HERE
	/*
	fprintf(stderr, "curFrom=%d\tcurRow=%d\tcurCol=%d\n",
			curFrom, curRow, curCol);
			*/
	assert(0 < curFrom);
	assert(curFrom <= 3*(ALPHABET_SIZE + 1));

	referenceLengthAligned=0;
	/* Init */
	if(curFrom <= (ALPHABET_SIZE + 1)) {
		PrintError(FnName, "curFrom", "Cannot end with a deletion", Exit, OutOfRange);
		length = matrix->cells[curRow][curCol].h.length[(curFrom - 1) % (ALPHABET_SIZE + 1)];
	}
	else if(2*(ALPHABET_SIZE + 1) < curFrom) {
		length = matrix->cells[curRow][curCol].v.length[(curFrom - 1) % (ALPHABET_SIZE + 1)];
	}
	else {
		length = matrix->cells[curRow][curCol].s.length[(curFrom - 1) % (ALPHABET_SIZE + 1)];
	}
	// HERE
	/*
	fprintf(stderr, "\nlength=%d\nreadLength=%d\n", length, readLength);
	if(length < readLength) {
		fprintf(stderr, "\nlength=%d\nreadLength=%d\n", length, readLength);
	}
	*/
	assert(readLength <= length);
	i=length-1;

	/* Now trace back the alignment using the "from" member in the matrix */
	while(0 <= i) {
		/*
		fprintf(stderr, "i=%d\tcurFrom=%d\tcurRow=%d\tcurCol=%d\n",
				i,
				curFrom, curRow, curCol);
				*/
		assert(0 <= curRow && 0 <= curCol);
		/* Where did the current cell come from */
		/* Get if there was a color error */
		if(curFrom <= (ALPHABET_SIZE + 1)) {
			/*
			   fprintf(stderr, "\ni=%d\ncurFrom=%d\nh.length=%d\n%s",
			   i,
			   curFrom,
			   matrix->cells[curRow][curCol].h.length[(curFrom - 1) % (ALPHABET_SIZE + 1)],
			   BREAK_LINE);
			   assert(i + 1 == matrix->cells[curRow][curCol].h.length[(curFrom - 1) % (ALPHABET_SIZE + 1)]);
			   */
			nextFrom = matrix->cells[curRow][curCol].h.from[(curFrom - 1) % (ALPHABET_SIZE + 1)];
			colorErrorAligned[i] = matrix->cells[curRow][curCol].h.colorError[(curFrom - 1) % (ALPHABET_SIZE + 1)];
		}
		else if(2*(ALPHABET_SIZE + 1) < curFrom) {
			/*
			   fprintf(stderr, "\ni=%d\ncurFrom=%d\nv.length=%d\n%s",
			   i,
			   curFrom,
			   matrix->cells[curRow][curCol].v.length[(curFrom - 1) % (ALPHABET_SIZE + 1)],
			   BREAK_LINE);
			   assert(i + 1 == matrix->cells[curRow][curCol].v.length[(curFrom - 1) % (ALPHABET_SIZE + 1)]);
			   */
			nextFrom = matrix->cells[curRow][curCol].v.from[(curFrom - 1) % (ALPHABET_SIZE + 1)];
			colorErrorAligned[i] = matrix->cells[curRow][curCol].v.colorError[(curFrom - 1) % (ALPHABET_SIZE + 1)];
		}
		else {
			/*
			   fprintf(stderr, "\ni=%d\ncurFrom=%d\ns.length=%d\n%s",
			   i,
			   curFrom,
			   matrix->cells[curRow][curCol].s.length[(curFrom - 1) % (ALPHABET_SIZE + 1)],
			   BREAK_LINE);
			   assert(i + 1 == matrix->cells[curRow][curCol].s.length[(curFrom - 1) % (ALPHABET_SIZE + 1)]);
			   */
			nextFrom = matrix->cells[curRow][curCol].s.from[(curFrom - 1) % (ALPHABET_SIZE + 1)];
			colorErrorAligned[i] = matrix->cells[curRow][curCol].s.colorError[(curFrom - 1) % (ALPHABET_SIZE + 1)];
		}
		colorErrorAligned[i] = ConvertIntColorToCharColor(colorErrorAligned[i]);

		switch(curFrom) {
			case MatchA:
			case InsertionA:
				readAligned[i] = 'A';
				break;
			case MatchC:
			case InsertionC:
				readAligned[i] = 'C';
				break;
			case MatchG:
			case InsertionG:
				readAligned[i] = 'G';
				break;
			case MatchT:
			case InsertionT:
				readAligned[i] = 'T';
				break;
			case MatchN:
			case InsertionN:
				readAligned[i] = 'N';
				break;
			case DeletionA:
			case DeletionC:
			case DeletionG:
			case DeletionT:
			case DeletionN:
				readAligned[i] = GAP;
				break;
			default:
				fprintf(stderr, "curFrom=%d\n",
						curFrom);
				PrintError(FnName, "curFrom", "Could not understand curFrom", Exit, OutOfRange);
		}

		switch(curFrom) {
			case InsertionA:
			case InsertionC:
			case InsertionG:
			case InsertionT:
			case InsertionN:
				referenceAligned[i] = GAP;
				break;
			default:
				referenceAligned[i] = reference[curCol-1];
				referenceLengthAligned++;
				break;
		}

		assert(readAligned[i] != GAP || readAligned[i] != referenceAligned[i]);

		/* Update next row/col */
		if(curFrom <= (ALPHABET_SIZE + 1)) {
			nextRow = curRow;
			nextCol = curCol-1;
		}
		else if(2*(ALPHABET_SIZE +1) < curFrom) {
			nextRow = curRow-1;
			nextCol = curCol;
		}
		else {
			nextRow = curRow-1;
			nextCol = curCol-1;
		}

		/* Update for next loop iteration */
		curFrom = nextFrom;
		curRow = nextRow;
		curCol = nextCol;
		i--;
	} /* End loop */
	assert(-1==i);
	assert(referenceLengthAligned <= length);

	offset = curCol;
	readAligned[length]='\0';
	referenceAligned[length]='\0';
	colorErrorAligned[length]='\0';

	/* Copy over */
	AlignedEntryUpdateAlignment(a,
			(FORWARD==strand) ? (position + offset) : (position + referenceLength - referenceLengthAligned - offset),
			maxScore,
			referenceLengthAligned,
			length,
			readAligned,
			referenceAligned,
			colorErrorAligned);

}

void AlignColorSpaceInitializeAtStart(char *read,
		AlignMatrix *matrix,
		ScoringMatrix *sm, 
		int32_t endRow, 
		int32_t endCol,
		int32_t alphabetSize)
{
	//char *FnName="AlignColorSpaceInitializeAtStart";
	int32_t i, j, k;

	// HERE
	/*
	fprintf(stderr, "%s: filling in [%d-%d,%d-%d]\n",
			FnName,
			0, endRow,
			0, endCol);
			*/

	/* Normal initialization */
	/* Allow the alignment to start anywhere in the reference */
	for(j=0;j<endCol+1;j++) {
		for(k=0;k<alphabetSize;k++) {
			matrix->cells[0][j].h.score[k] = NEGATIVE_INFINITY;
			matrix->cells[0][j].h.from[k] = StartCS;
			matrix->cells[0][j].h.length[k] = 0;
			matrix->cells[0][j].h.colorError[k] = GAP;

			/* Assumes both DNA and COLOR_SPACE_START_NT are upper case */
			if(DNA[k] == COLOR_SPACE_START_NT) { 
				/* Starting adaptor NT */
				matrix->cells[0][j].s.score[k] = 0;
			}
			else {
				matrix->cells[0][j].s.score[k] = NEGATIVE_INFINITY;
			}
			matrix->cells[0][j].s.from[k] = StartCS;
			matrix->cells[0][j].s.length[k] = 0;
			matrix->cells[0][j].s.colorError[k] = GAP;

			matrix->cells[0][j].v.score[k] = NEGATIVE_INFINITY;
			matrix->cells[0][j].v.from[k] = StartCS;
			matrix->cells[0][j].v.length[k] = 0;
			matrix->cells[0][j].v.colorError[k] = GAP;
		}
	}
	/* Row i (i>0) column 0 should be negative infinity since we want to
	 * align the full read */
	for(i=1;i<endRow+1;i++) {
		for(k=0;k<alphabetSize;k++) {
			matrix->cells[i][0].h.score[k] = NEGATIVE_INFINITY;
			matrix->cells[i][0].h.from[k] = StartCS;
			matrix->cells[i][0].h.length[k] = 0;
			matrix->cells[i][0].h.colorError[k] = GAP;

			matrix->cells[i][0].s.score[k] = NEGATIVE_INFINITY;
			matrix->cells[i][0].s.from[k] = StartCS;
			matrix->cells[i][0].s.length[k] = 0;
			matrix->cells[i][0].s.colorError[k] = GAP;

			// Allow an insertion
			if(DNA[k] == read[i-1]) { // Must be consistent with the read (no color errors please)
				if(i == 1) { // Allow for an insertion start
					matrix->cells[i][0].v.score[k] = matrix->cells[i-1][0].s.score[COLOR_SPACE_START_NT_INT] + sm->gapOpenPenalty;
					matrix->cells[i][0].v.from[k] = COLOR_SPACE_START_NT_INT + 1 + (ALPHABET_SIZE + 1); /* see the enum */
					matrix->cells[i][0].v.length[k] = matrix->cells[i-1][0].s.length[COLOR_SPACE_START_NT_INT] + 1;
					matrix->cells[i][0].v.colorError[k] = GAP;
				}
				else { // Allow for an insertion extension
					int32_t fromNT = BaseToInt(read[i-2]); // previous NT
					matrix->cells[i][0].v.score[k] = matrix->cells[i-1][0].s.score[fromNT] + sm->gapOpenPenalty;
					matrix->cells[i][0].v.from[k] = fromNT + 1 + 2*(ALPHABET_SIZE + 1); /* see the enum */
					matrix->cells[i][0].v.length[k] = matrix->cells[i-1][0].s.length[fromNT] + 1;
					matrix->cells[i][0].v.colorError[k] = GAP;

				}
			}
			else {
				matrix->cells[i][0].v.score[k] = NEGATIVE_INFINITY;
				matrix->cells[i][0].v.from[k] = StartCS;
				matrix->cells[i][0].v.length[k] = 0;
				matrix->cells[i][0].v.colorError[k] = GAP;
			}
		}
	}
}

void AlignColorSpaceInitializeToExtend(char *read,
		AlignMatrix *matrix,
		ScoringMatrix *sm, 
		int32_t readLength,
		int32_t referenceLength,
		int32_t startRow, 
		int32_t startCol,
		int32_t alphabetSize)
{
	char *FnName="AlignColorSpaceInitializeToExtend";
	int32_t i, j, k, endRow, endCol;

	assert(0 < startRow && 0 < startCol);

	endRow = readLength;
	endCol = referenceLength;

	// HERE
	/*
	fprintf(stderr, "%s: filling in [%d-%d,%d-%d]\n",
			FnName,
			startRow, endRow,
			startCol, endCol);
			*/

	/* Special initialization */

	/* Initialize the corner cell */
	// Check that the match has been filled in 
	for(k=0;k<alphabetSize;k++) {
		assert(1 + ALPHABET_SIZE < matrix->cells[startRow][startCol].s.from[k] &&
				matrix->cells[startRow][startCol].s.from[k] <= 2*(ALPHABET_SIZE + 1));
		// HERE
		/*
		fprintf(stderr, "startRow=%d <= matrix->cells[startRow][startCol].s.length[k]=%d\n",
				startRow,
				matrix->cells[startRow][startCol].s.length[k]);
		fprintf(stderr, "startRow=%d\tstartCol=%d\tk=%d\n",
				startRow,
				startCol,
				k);
		*/
		assert(startRow <= matrix->cells[startRow][startCol].s.length[k]);
		// Do not allow a deletion or insertion
		matrix->cells[startRow][startCol].h.score[k] = matrix->cells[startRow][startCol].v.score[k] = NEGATIVE_INFINITY;
		matrix->cells[startRow][startCol].h.from[k] = matrix->cells[startRow][startCol].v.from[k] = StartNT;
		matrix->cells[startRow][startCol].h.length[k] = matrix->cells[startRow][startCol].v.length[k] = 0;
		matrix->cells[startRow][startCol].h.colorError[k] = matrix->cells[startRow][startCol].v.colorError[k] = GAP;
	}

	// TODO
	for(j=startCol+1;j<endCol+1;j++) { // Columns
		for(k=0;k<alphabetSize;k++) { // To NT
			if(j == startCol + 1) { // Allow for a deletion start
				matrix->cells[startRow][j].h.score[k] = matrix->cells[startRow][j-1].s.score[k] + sm->gapOpenPenalty;
				matrix->cells[startRow][j].h.length[k] = matrix->cells[startRow][j-1].s.length[k] + 1;
				matrix->cells[startRow][j].h.from[k] = k + 1 + (ALPHABET_SIZE + 1); /* see the enum */ 
				matrix->cells[startRow][j].h.colorError[k] = GAP;
			}
			else { // Allow for a deletion extension
				matrix->cells[startRow][j].h.score[k] = matrix->cells[startRow][j-1].h.score[k] + sm->gapExtensionPenalty;
				matrix->cells[startRow][j].h.length[k] = matrix->cells[startRow][j-1].h.length[k] + 1;
				matrix->cells[startRow][j].h.from[k] = k + 1;
				matrix->cells[startRow][j].h.colorError[k] = GAP;
			}

			// Do not allow for a match or an insertion
			matrix->cells[startRow][j].s.score[k] = matrix->cells[startRow][j].v.score[k] = NEGATIVE_INFINITY;
			matrix->cells[startRow][j].s.from[k] = matrix->cells[startRow][j].v.from[k] = StartNT;
			matrix->cells[startRow][j].s.length[k] = matrix->cells[startRow][j].v.length[k] = 0;
			matrix->cells[startRow][j].s.colorError[k] = matrix->cells[startRow][j].v.colorError[k] = GAP;
		}
	}

	/* Align the full read */
	for(i=startRow+1;i<endRow+1;i++) {
		char curColor, base;
		int32_t fromNT;
		/* Get the current color for the read */
		assert(1 < i); // Otherwise we should use the COLOR_SPACE_START_NT
		if(0 == ConvertBaseToColorSpace(read[i-2], read[i-1], &curColor)) {
			PrintError(FnName, "curColor", "Could not convert base to color space", Exit, OutOfRange);
		}
		for(k=0;k<alphabetSize;k++) {
			// Do not allow for a match or a deletion
			matrix->cells[i][startCol].h.score[k] = matrix->cells[i][startCol].s.score[k] = NEGATIVE_INFINITY;
			matrix->cells[i][startCol].h.from[k] = matrix->cells[i][startCol].s.from[k] = StartNT;
			matrix->cells[i][startCol].h.length[k] = matrix->cells[i][startCol].s.length[k] = 0;
			matrix->cells[i][startCol].h.colorError[k] = matrix->cells[i][startCol].s.colorError[k] = GAP;

			/* Get from base for extending an insertion */
			if(0 == ConvertBaseAndColor(DNA[k], curColor, &base)) {
				PrintError(FnName, NULL, "Could not convert base and color", Exit, OutOfRange);
			}
			fromNT=BaseToInt(base);

			if(i == startRow + 1) { // Allow for an insertion start
				matrix->cells[i][startCol].v.score[k] = matrix->cells[i-1][startCol].s.score[fromNT] + sm->gapOpenPenalty;
				matrix->cells[i][startCol].v.length[k] = matrix->cells[i-1][startCol].s.length[fromNT] + 1;
				matrix->cells[i][startCol].v.from[k] = fromNT + 1 + (ALPHABET_SIZE + 1);
				matrix->cells[i][startCol].v.colorError[k] = GAP;
			}
			else { // Allow for an insertion extension
				matrix->cells[i][startCol].v.score[k] = matrix->cells[i-1][startCol].v.score[fromNT] + sm->gapExtensionPenalty;
				matrix->cells[i][startCol].v.length[k] = matrix->cells[i-1][startCol].v.length[fromNT] + 1;
				matrix->cells[i][startCol].v.from[k] = fromNT + 1 + 2*(ALPHABET_SIZE + 1);
				matrix->cells[i][startCol].v.colorError[k] = GAP;
			}
		}
	}
}

void AlignColorSpaceFillInCell(char *read,
		int32_t readLength,
		char *reference,
		int32_t referenceLength,
		ScoringMatrix *sm,
		AlignMatrix *matrix,
		int32_t row,
		int32_t col,
		char curColor,
		int32_t maxH,
		int32_t maxV,
		int32_t alphabetSize)
{
	char *FnName = "AlignColorSpaceFillInCell";
	int32_t k, l;

	// HERE
	/*
	fprintf(stderr, "%s: filling in [%d,%d]\n",
			FnName,
			row+1,
			col+1);
			*/

	/* Deletion */
	if(maxV <= row - col) { // Out of bounds, do not consider
		for(k=0;k<alphabetSize;k++) { /* To NT */
			/* Update */
			matrix->cells[row+1][col+1].h.score[k] = NEGATIVE_INFINITY-1;
			matrix->cells[row+1][col+1].h.from[k] = NoFromCS;
			matrix->cells[row+1][col+1].h.colorError[k] = GAP;
			matrix->cells[row+1][col+1].h.length[k] = INT_MIN;
		}
	}
	else {
		for(k=0;k<alphabetSize;k++) { /* To NT */
			int32_t maxScore = NEGATIVE_INFINITY-1;
			int maxFrom = -1;
			char maxColorError = GAP;
			int maxLength = 0;

			int32_t curScore=NEGATIVE_INFINITY;
			int curLength=-1;

			/* Deletion starts or extends from the same base */

			/* New deletion */
			curLength = matrix->cells[row+1][col].s.length[k] + 1;
			/* Deletion - previous column */
			/* Ignore color error since one color will span the entire
			 * deletion.  We will consider the color at the end of the deletion.
			 * */
			curScore = matrix->cells[row+1][col].s.score[k] + sm->gapOpenPenalty;
			/* Make sure we aren't below infinity */
			if(curScore < NEGATIVE_INFINITY/2) {
				curScore = NEGATIVE_INFINITY;
				assert(curScore < 0);
			}
			if(curScore > maxScore) {
				maxScore = curScore;
				maxFrom = k + 1 + (ALPHABET_SIZE + 1); /* see the enum */ 
				maxColorError = GAP;
				maxLength = curLength;
			}

			/* Extend current deletion */
			curLength = matrix->cells[row+1][col].h.length[k] + 1;
			/* Deletion - previous column */
			curScore = matrix->cells[row+1][col].h.score[k] + sm->gapExtensionPenalty;
			/* Ignore color error since one color will span the entire
			 * deletion.  We will consider the color at the end of the deletion.
			 * */
			/* Make sure we aren't below infinity */
			if(curScore < NEGATIVE_INFINITY/2) {
				curScore = NEGATIVE_INFINITY;
			}
			if(curScore > maxScore) {
				maxScore = curScore;
				maxFrom = k + 1; /* see the enum */ 
				maxColorError = GAP;
				maxLength = curLength;
			}
			/* Update */
			matrix->cells[row+1][col+1].h.score[k] = maxScore;
			matrix->cells[row+1][col+1].h.from[k] = maxFrom;
			matrix->cells[row+1][col+1].h.colorError[k] = maxColorError;
			matrix->cells[row+1][col+1].h.length[k] = maxLength;
		}
	}

	/* Match/Mismatch */
	for(k=0;k<alphabetSize;k++) { /* To NT */
		int32_t maxScore = NEGATIVE_INFINITY-1;
		int maxFrom = -1;
		char maxColorError = GAP;
		int maxLength = 0;

		for(l=0;l<alphabetSize;l++) { /* From NT */
			int32_t curScore=NEGATIVE_INFINITY;
			int curLength=-1;
			char convertedColor='X';
			int32_t scoreNT, scoreColor;

			/* Get color */
			if(0 == ConvertBaseToColorSpace(DNA[l], DNA[k], &convertedColor)) {
				fprintf(stderr, "DNA[l=%d]=%c\tDNA[k=%d]=%c\n", l, DNA[l], k, DNA[k]);
				PrintError(FnName, "convertedColor", "Could not convert base to color space", Exit, OutOfRange);
			}
			/* Get NT and Color scores */
			scoreNT = ScoringMatrixGetNTScore(reference[col], DNA[k], sm);
			scoreColor = ScoringMatrixGetColorScore(curColor,
					convertedColor,
					sm);

			/* From Horizontal - Deletion */
			curLength = matrix->cells[row][col].h.length[l] + 1;
			/* Add previous with current NT */
			curScore = matrix->cells[row][col].h.score[l] + scoreNT;
			/* Add score for color error, if any */
			curScore += scoreColor;
			/* Make sure we aren't below infinity */
			if(curScore < NEGATIVE_INFINITY/2) {
				curScore = NEGATIVE_INFINITY;
			}
			if(curScore > maxScore) {
				maxScore = curScore;
				maxFrom = l + 1; /* see the enum */ 
				maxColorError = (curColor == convertedColor)?GAP:curColor; /* Keep original color */
				maxLength = curLength;
			}

			/* From Vertical - Insertion */
			curLength = matrix->cells[row][col].v.length[l] + 1;
			/* Add previous with current NT */
			curScore = matrix->cells[row][col].v.score[l] + scoreNT;
			/* Add score for color error, if any */
			curScore += scoreColor;
			/* Make sure we aren't below infinity */
			if(curScore < NEGATIVE_INFINITY/2) {
				curScore = NEGATIVE_INFINITY;
			}
			if(curScore > maxScore) {
				maxScore = curScore;
				maxFrom = l + 1 + 2*(ALPHABET_SIZE + 1); /* see the enum */ 
				maxColorError = (curColor == convertedColor)?GAP:curColor; /* Keep original color */
				maxLength = curLength;
			}

			/* From Diagonal - Match/Mismatch */
			curLength = matrix->cells[row][col].s.length[l] + 1;
			/* Add previous with current NT */
			curScore = matrix->cells[row][col].s.score[l] + scoreNT;
			/* Add score for color error, if any */
			curScore += scoreColor;
			/* Make sure we aren't below infinity */
			if(curScore < NEGATIVE_INFINITY/2) {
				curScore = NEGATIVE_INFINITY;
			}
			if(curScore > maxScore) {
				maxScore = curScore;
				maxFrom = l + 1 + (ALPHABET_SIZE + 1); /* see the enum */ 
				maxColorError = (curColor == convertedColor)?GAP:curColor; /* Keep original color */
				maxLength = curLength;
			}
		}
		/* Update */
		matrix->cells[row+1][col+1].s.score[k] = maxScore;
		matrix->cells[row+1][col+1].s.from[k] = maxFrom;
		matrix->cells[row+1][col+1].s.colorError[k] = maxColorError;
		matrix->cells[row+1][col+1].s.length[k] = maxLength;
	}

	/* Insertion */
	if(maxH <= col - referenceLength + readLength + row) {
		/* We are on the boundary, do not consider an insertion */
		for(k=0;k<alphabetSize;k++) { /* To NT */
			/* Update */
			matrix->cells[row+1][col+1].v.score[k] = NEGATIVE_INFINITY-1;
			matrix->cells[row+1][col+1].v.from[k] = NoFromCS;
			matrix->cells[row+1][col+1].v.colorError[k] = GAP;
			matrix->cells[row+1][col+1].v.length[k] = INT_MIN;
		}
	}
	else {
		for(k=0;k<alphabetSize;k++) { /* To NT */
			int32_t maxScore = NEGATIVE_INFINITY-1;
			int maxFrom = -1;
			char maxColorError = GAP;
			int maxLength = 0;

			int32_t curScore=NEGATIVE_INFINITY;
			int curLength=-1;
			char B;
			int fromNT=-1;

			/* Get from base for extending an insertion */
			if(0 == ConvertBaseAndColor(DNA[k], curColor, &B)) {
				PrintError(FnName, NULL, "Could not convert base and color", Exit, OutOfRange);
			}
			fromNT=BaseToInt(B);

			/* New insertion */
			curScore=NEGATIVE_INFINITY;
			curLength=-1;
			/* Get NT and Color scores */
			curLength = matrix->cells[row][col+1].s.length[fromNT] + 1;
			curScore = matrix->cells[row][col+1].s.score[fromNT] + sm->gapOpenPenalty;
			/*
			   curScore += ScoringMatrixGetColorScore(curColor,
			   convertedColor,
			   sm);
			   */
			/* Make sure we aren't below infinity */
			if(curScore < NEGATIVE_INFINITY/2) {
				curScore = NEGATIVE_INFINITY;
			}
			if(curScore > maxScore) {
				maxScore = curScore;
				maxFrom = fromNT + 1 + (ALPHABET_SIZE + 1); /* see the enum */ 
				maxColorError = GAP;
				maxLength = curLength;
			}

			/* Extend current insertion */
			curLength = matrix->cells[row][col+1].v.length[fromNT] + 1;
			/* Insertion - previous row */
			curScore = matrix->cells[row][col+1].v.score[fromNT] + sm->gapExtensionPenalty;
			curScore += ScoringMatrixGetColorScore(curColor,
					curColor,
					sm);
			/* Make sure we aren't below infinity */
			if(curScore < NEGATIVE_INFINITY/2) {
				curScore = NEGATIVE_INFINITY;
			}
			if(curScore > maxScore) {
				maxScore = curScore;
				maxFrom = fromNT + 1 + 2*(ALPHABET_SIZE + 1); /* see the enum */ 
				maxColorError = GAP;
				maxLength = curLength;
			}

			/* Update */
			matrix->cells[row+1][col+1].v.score[k] = maxScore;
			matrix->cells[row+1][col+1].v.from[k] = maxFrom;
			matrix->cells[row+1][col+1].v.colorError[k] = maxColorError;
			matrix->cells[row+1][col+1].v.length[k] = maxLength;
		}
	}
}

int32_t AlignColorSpaceGetAlphabetSize(char *read,
		int32_t readLength,
		char *reference,
		int32_t referenceLength) 
{
	int32_t i;
	int32_t alphabetSize=ALPHABET_SIZE;
	/* Check if there are any Ns or 0s */
	for(i=0;i<readLength;i++) {
		if(1 == RGBinaryIsBaseN(read[i])) {
			alphabetSize++; // include room for the N
			return alphabetSize;
		}
	}
	for(i=0;i<referenceLength && ALPHABET_SIZE==alphabetSize;i++) {
		if(1 == RGBinaryIsBaseN(reference[i])) {
			alphabetSize++; // include room for the N
			return alphabetSize;
		}
	}
	return alphabetSize;
}
