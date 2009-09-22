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

/* TODO */
void AlignColorSpaceUngapped(char *read,
		char *mask,
		int readLength,
		char *reference,
		int referenceLength,
		int unconstrained,
		ScoringMatrix *sm,
		AlignedEntry *a,
		int offset,
		uint32_t position,
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

	/* Check if there are any Ns or 0s */
	alphabetSize=ALPHABET_SIZE;
	for(i=0;i<readLength;i++) {
		if(1 == RGBinaryIsBaseN(read[i])) {
			alphabetSize++; // include room for the N
			break;
		}
	}
	for(i=offset;i<referenceLength-offset+1 && ALPHABET_SIZE==alphabetSize;i++) {
		if(1 == RGBinaryIsBaseN(reference[i])) {
			alphabetSize++; // include room for the N
			break;
		}
	}

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
			char curReadBase = read[j];
			char prevReadBase = (j==0)?COLOR_SPACE_START_NT:read[j-1];

			/* Get the current color for the read */
			if(0 == ConvertBaseToColorSpace(prevReadBase, curReadBase, &curColor)) {
				fprintf(stderr, "read=%s\n", read);
				fprintf(stderr, "prevReadBase=%c\tcurReadBase=%c\n",
						prevReadBase,
						curReadBase);
				PrintError(FnName, "curColor", "Could not convert base to color space", Exit, OutOfRange);
			}
			int32_t nextScore[ALPHABET_SIZE+1];
			char nextNT[ALPHABET_SIZE+1];
			for(k=0;k<alphabetSize;k++) { /* To NT */

				/* Get the best score to this NT */
				int32_t bestScore = NEGATIVE_INFINITY;
				int bestNT=-1;
				char bestColor = 'X';
				int allowedNT[ALPHABET_SIZE+1];

				if(Constrained == unconstrained) {
					for(l=0;l<alphabetSize;l++) { /* From NT */
						allowedNT[l]=0;
					}
					char base;
					assert(1 == ConvertBaseAndColor(DNA[k], curColor, &base));
					switch(ToLower(base)) {
						case 'a':
							allowedNT[0] = 1; break;
						case 'c':
							allowedNT[1] = 1; break;
						case 'g':
							allowedNT[2] = 1; break;
						case 't':
							allowedNT[3] = 1; break;
						case 'n':
							assert(ALPHABET_SIZE < alphabetSize);
							allowedNT[4] = 1; break;
						default:
							PrintError(FnName, "base", "Could not understand base", Exit, OutOfRange);
							break;
					}
				}
				else {
					for(l=0;l<alphabetSize;l++) { /* From NT */
						allowedNT[l]=1;
					}
				}

				for(l=0;l<alphabetSize;l++) { /* From NT */
					if(1 == allowedNT[l]) {
						char convertedColor='X';
						int32_t curScore = prevScore[l];
						/* Get color */
						if(0 == ConvertBaseToColorSpace(DNA[l], DNA[k], &convertedColor)) {
							fprintf(stderr, "DNA[l=%d]=%c\tDNA[k=%d]=%c\n",
									l,
									DNA[l],
									k,
									DNA[k]);
							PrintError(FnName, "convertedColor", "Could not convert base to color space", Exit, OutOfRange);
						}
						/* Add score for color error, if any */
						curScore += ScoringMatrixGetColorScore(curColor,
								convertedColor,
								sm);
						/* Add score for NT */
						curScore += ScoringMatrixGetNTScore(reference[i+j], DNA[k], sm);

						if(curScore < NEGATIVE_INFINITY/2) {
							curScore = NEGATIVE_INFINITY;
						}

						if(bestScore < curScore) {
							bestScore = curScore;
							bestNT = l;
							bestColor = convertedColor;
						}
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
}

/* TODO */
void AlignColorSpaceFull(char *read,
		int readLength,
		char *reference,
		int referenceLength,
		ScoringMatrix *sm,
		AlignedEntry *a,
		AlignMatrixCS **matrix,
		uint32_t position,
		char strand)
{
	/* read goes on the rows, reference on the columns */
	char *FnName = "AlignColorSpaceFull";
	int i, j, k, l;
	int alphabetSize=ALPHABET_SIZE;

	/* Check if there are any Ns or 0s */
	alphabetSize=ALPHABET_SIZE;
	for(i=0;i<readLength;i++) {
		if(1 == RGBinaryIsBaseN(read[i])) {
			alphabetSize++; // include room for the N
			break;
		}
	}
	for(i=0;i<referenceLength && ALPHABET_SIZE==alphabetSize;i++) {
		if(1 == RGBinaryIsBaseN(reference[i])) {
			alphabetSize++; // include room for the N
			break;
		}
	}

	/* Initialize "the matrix" */
	/* Row i (i>0) column 0 should be negative infinity since we want to
	 * align the full read */
	for(i=1;i<readLength+1;i++) {
		for(k=0;k<alphabetSize;k++) {
			matrix[i][0].h.score[k] = NEGATIVE_INFINITY;
			matrix[i][0].h.from[k] = StartCS;
			matrix[i][0].h.length[k] = 0;
			matrix[i][0].h.colorError[k] = GAP;

			matrix[i][0].s.score[k] = NEGATIVE_INFINITY;
			matrix[i][0].s.from[k] = StartCS;
			matrix[i][0].s.length[k] = 0;
			matrix[i][0].s.colorError[k] = GAP;

			matrix[i][0].v.score[k] = NEGATIVE_INFINITY;
			matrix[i][0].v.from[k] = StartCS;
			matrix[i][0].v.length[k] = 0;
			matrix[i][0].v.colorError[k] = GAP;
		}
	}
	/* Row 0 column j should be zero since we want to find the best
	 * local alignment within the reference */
	for(j=0;j<referenceLength+1;j++) {
		for(k=0;k<alphabetSize;k++) {
			matrix[0][j].h.score[k] = NEGATIVE_INFINITY;
			matrix[0][j].h.from[k] = StartCS;
			matrix[0][j].h.length[k] = 0;
			matrix[0][j].h.colorError[k] = GAP;

			/* Assumes both DNA and COLOR_SPACE_START_NT are upper case */
			if(DNA[k] == COLOR_SPACE_START_NT) { 
				/* Starting adaptor NT */
				matrix[0][j].s.score[k] = 0;
			}
			else {
				matrix[0][j].s.score[k] = NEGATIVE_INFINITY;
			}
			matrix[0][j].s.from[k] = StartCS;
			matrix[0][j].s.length[k] = 0;
			matrix[0][j].s.colorError[k] = GAP;

			matrix[0][j].v.score[k] = NEGATIVE_INFINITY;
			matrix[0][j].v.from[k] = StartCS;
			matrix[0][j].v.length[k] = 0;
			matrix[0][j].v.colorError[k] = GAP;
		}
	}

	/* Fill in the matrix according to the recursive rules */
	for(i=0;i<readLength;i++) { /* read/rows */
		/* Get the current color */
		char curColor;
		char curReadBase, prevReadBase;
		/* In color space, the first color is determined by the adapter NT */
		curReadBase = read[i];
		prevReadBase = (i==0)?COLOR_SPACE_START_NT:read[i-1];

		/* Get the current color for the read */
		if(0 == ConvertBaseToColorSpace(prevReadBase, curReadBase, &curColor)) {
			fprintf(stderr, "prevReadBase=%c\tcurReadBase=%c\n",
					prevReadBase,
					curReadBase);
			PrintError(FnName, "curColor", "Could not convert base to color space", Exit, OutOfRange);
		}

		for(j=0;j<referenceLength;j++) { /* reference/columns */

			/* Deletion */
			for(k=0;k<alphabetSize;k++) { /* To NT */
				int32_t maxScore = NEGATIVE_INFINITY-1;
				int maxFrom = -1;
				char maxColorError = GAP;
				int maxLength = 0;

				int32_t curScore=NEGATIVE_INFINITY;
				int curLength=-1;

				/* Deletion starts or extends from the same base */

				/* New deletion */
				curLength = matrix[i+1][j].s.length[k] + 1;
				/* Deletion - previous column */
				/* Ignore color error since one color will span the entire
				 * deletion.  We will consider the color at the end of the deletion.
				 * */
				curScore = matrix[i+1][j].s.score[k] + sm->gapOpenPenalty;
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
				curLength = matrix[i+1][j].h.length[k] + 1;
				/* Deletion - previous column */
				curScore = matrix[i+1][j].h.score[k] + sm->gapExtensionPenalty;
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
				matrix[i+1][j+1].h.score[k] = maxScore;
				matrix[i+1][j+1].h.from[k] = maxFrom;
				matrix[i+1][j+1].h.colorError[k] = maxColorError;
				matrix[i+1][j+1].h.length[k] = maxLength;
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
						fprintf(stderr, "DNA[l=%d]=%c\tDNA[k=%d]=%c\n",
								l,
								DNA[l],
								k,
								DNA[k]);
						PrintError(FnName, "convertedColor", "Could not convert base to color space", Exit, OutOfRange);
					}
					/* Get NT and Color scores */
					scoreNT = ScoringMatrixGetNTScore(reference[j], DNA[k], sm);
					scoreColor = ScoringMatrixGetColorScore(curColor,
							convertedColor,
							sm);

					/* From Horizontal - Deletion */
					curLength = matrix[i][j].h.length[l] + 1;
					/* Add previous with current NT */
					curScore = matrix[i][j].h.score[l] + scoreNT;
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
					curLength = matrix[i][j].v.length[l] + 1;
					/* Add previous with current NT */
					curScore = matrix[i][j].v.score[l] + scoreNT;
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
					curLength = matrix[i][j].s.length[l] + 1;
					/* Add previous with current NT */
					curScore = matrix[i][j].s.score[l] + scoreNT;
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
				matrix[i+1][j+1].s.score[k] = maxScore;
				matrix[i+1][j+1].s.from[k] = maxFrom;
				matrix[i+1][j+1].s.colorError[k] = maxColorError;
				matrix[i+1][j+1].s.length[k] = maxLength;
			}

			/* Insertion */
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
				switch(B) {
					case 'a':
					case 'A':
						fromNT=0;
						break;
					case 'c':
					case 'C':
						fromNT=1;
						break;
					case 'g':
					case 'G':
						fromNT=2;
						break;
					case 't':
					case 'T':
						fromNT=3;
						break;
					default:
						fromNT=4;
						break;
				}

				/* New insertion */
				curScore=NEGATIVE_INFINITY;
				curLength=-1;
				/* Get NT and Color scores */
				curLength = matrix[i][j+1].s.length[fromNT] + 1;
				curScore = matrix[i][j+1].s.score[fromNT] + sm->gapOpenPenalty;
				curScore += ScoringMatrixGetColorScore(curColor,
						curColor,
						sm);
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
				curLength = matrix[i][j+1].v.length[fromNT] + 1;
				/* Insertion - previous row */
				curScore = matrix[i][j+1].v.score[fromNT] + sm->gapExtensionPenalty;
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
				matrix[i+1][j+1].v.score[k] = maxScore;
				matrix[i+1][j+1].v.from[k] = maxFrom;
				matrix[i+1][j+1].v.colorError[k] = maxColorError;
				matrix[i+1][j+1].v.length[k] = maxLength;
			}
		}
	}

	FillAlignedEntryFromMatrixColorSpace(a,
			matrix,
			read,
			readLength,
			reference,
			referenceLength,
			0,
			position,
			strand,
			alphabetSize,
			0);
}

/* TODO */
void AlignColorSpaceFullWithBound(char *read,
		int readLength,
		char *reference,
		int referenceLength,
		ScoringMatrix *sm,
		AlignedEntry *a,
		int32_t maxH,
		int32_t maxV,
		AlignMatrixCS **matrix,
		uint32_t position,
		char strand)
{
	char *FnName = "AlignColorSpaceFullWithBound";
	int i, j, k, l;
	int alphabetSize=ALPHABET_SIZE;

	assert(0 < readLength);
	assert(0 < referenceLength);

	/* Check if there are any Ns or 0s */
	alphabetSize=ALPHABET_SIZE;
	for(i=0;i<readLength;i++) {
		if(1 == RGBinaryIsBaseN(read[i])) {
			alphabetSize++; // include room for the N
			break;
		}
	}
	for(i=0;i<referenceLength && ALPHABET_SIZE==alphabetSize;i++) {
		if(1 == RGBinaryIsBaseN(reference[i])) {
			alphabetSize++; // include room for the N
			break;
		}
	}

	/* Initialize "the matrix" */
	/* Row i (i>0) column 0 should be negative infinity since we want to
	 * align the full read */
	for(i=1;i<readLength+1;i++) {
		for(k=0;k<alphabetSize;k++) {
			matrix[i][0].h.score[k] = NEGATIVE_INFINITY;
			matrix[i][0].h.from[k] = StartCS;
			matrix[i][0].h.length[k] = 0;
			matrix[i][0].h.colorError[k] = GAP;

			matrix[i][0].s.score[k] = NEGATIVE_INFINITY;
			matrix[i][0].s.from[k] = StartCS;
			matrix[i][0].s.length[k] = 0;
			matrix[i][0].s.colorError[k] = GAP;

			matrix[i][0].v.score[k] = NEGATIVE_INFINITY;
			matrix[i][0].v.from[k] = StartCS;
			matrix[i][0].v.length[k] = 0;
			matrix[i][0].v.colorError[k] = GAP;
		}
	}
	/* Row 0 column j should be zero since we want to find the best
	 * local alignment within the reference */
	for(j=0;j<referenceLength+1;j++) {
		for(k=0;k<alphabetSize;k++) {
			matrix[0][j].h.score[k] = NEGATIVE_INFINITY;
			matrix[0][j].h.from[k] = StartCS;
			matrix[0][j].h.length[k] = 0;
			matrix[0][j].h.colorError[k] = GAP;

			/* Assumes both DNA and COLOR_SPACE_START_NT are upper case */
			if(DNA[k] == COLOR_SPACE_START_NT) { 
				/* Starting adaptor NT */
				matrix[0][j].s.score[k] = 0;
			}
			else {
				matrix[0][j].s.score[k] = NEGATIVE_INFINITY;
			}
			matrix[0][j].s.from[k] = StartCS;
			matrix[0][j].s.length[k] = 0;
			matrix[0][j].s.colorError[k] = GAP;

			matrix[0][j].v.score[k] = NEGATIVE_INFINITY;
			matrix[0][j].v.from[k] = StartCS;
			matrix[0][j].v.length[k] = 0;
			matrix[0][j].v.colorError[k] = GAP;
		}
	}

	/* Fill in the matrix according to the recursive rules */
	for(i=0;i<readLength;i++) { /* read/rows */
		/* Get the current color */
		char curColor;
		char curReadBase, prevReadBase;
		/* In color space, the first color is determined by the adapter NT */
		curReadBase = read[i];
		prevReadBase = (i==0)?COLOR_SPACE_START_NT:read[i-1];

		/* Get the current color for the read */
		if(0 == ConvertBaseToColorSpace(prevReadBase, curReadBase, &curColor)) {
			fprintf(stderr, "prevReadBase=%c\tcurReadBase=%c\n",
					prevReadBase,
					curReadBase);
			PrintError(FnName, "curColor", "Could not convert base to color space", Exit, OutOfRange);
		}

		for(j=GETMAX(0, i - maxV);
				j <= GETMIN(referenceLength-1, referenceLength - (readLength - maxH) + i);
				j++) { /* reference/columns */
			assert(i-maxV <= j && j <= referenceLength - (readLength - maxH) + i);

			/* Deletion */
			if(maxV == i - j) {
				/* We are on the boundary, do not consider a deletion */
				for(k=0;k<alphabetSize;k++) { /* To NT */
					/* Update */
					matrix[i+1][j+1].h.score[k] = NEGATIVE_INFINITY-1;
					matrix[i+1][j+1].h.from[k] = NoFromCS;
					matrix[i+1][j+1].h.colorError[k] = GAP;
					matrix[i+1][j+1].h.length[k] = INT_MIN;
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
					curLength = matrix[i+1][j].s.length[k] + 1;
					/* Deletion - previous column */
					/* Ignore color error since one color will span the entire
					 * deletion.  We will consider the color at the end of the deletion.
					 * */
					curScore = matrix[i+1][j].s.score[k] + sm->gapOpenPenalty;
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
					curLength = matrix[i+1][j].h.length[k] + 1;
					/* Deletion - previous column */
					curScore = matrix[i+1][j].h.score[k] + sm->gapExtensionPenalty;
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
					matrix[i+1][j+1].h.score[k] = maxScore;
					matrix[i+1][j+1].h.from[k] = maxFrom;
					matrix[i+1][j+1].h.colorError[k] = maxColorError;
					matrix[i+1][j+1].h.length[k] = maxLength;
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
						fprintf(stderr, "DNA[l=%d]=%c\tDNA[k=%d]=%c\n",
								l,
								DNA[l],
								k,
								DNA[k]);
						PrintError(FnName, "convertedColor", "Could not convert base to color space", Exit, OutOfRange);
					}
					/* Get NT and Color scores */
					scoreNT = ScoringMatrixGetNTScore(reference[j], DNA[k], sm);
					scoreColor = ScoringMatrixGetColorScore(curColor,
							convertedColor,
							sm);

					/* From Horizontal - Deletion */
					curLength = matrix[i][j].h.length[l] + 1;
					/* Add previous with current NT */
					curScore = matrix[i][j].h.score[l] + scoreNT;
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
					curLength = matrix[i][j].v.length[l] + 1;
					/* Add previous with current NT */
					curScore = matrix[i][j].v.score[l] + scoreNT;
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
					curLength = matrix[i][j].s.length[l] + 1;
					/* Add previous with current NT */
					curScore = matrix[i][j].s.score[l] + scoreNT;
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
				matrix[i+1][j+1].s.score[k] = maxScore;
				matrix[i+1][j+1].s.from[k] = maxFrom;
				matrix[i+1][j+1].s.colorError[k] = maxColorError;
				matrix[i+1][j+1].s.length[k] = maxLength;
			}

			/* Insertion */
			if(j == referenceLength - (readLength - maxH) + i) {
				/* We are on the boundary, do not consider an insertion */
				for(k=0;k<alphabetSize;k++) { /* To NT */
					/* Update */
					matrix[i+1][j+1].v.score[k] = NEGATIVE_INFINITY-1;
					matrix[i+1][j+1].v.from[k] = NoFromCS;
					matrix[i+1][j+1].v.colorError[k] = GAP;
					matrix[i+1][j+1].v.length[k] = INT_MIN;
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
					switch(B) {
						case 'a':
						case 'A':
							fromNT=0;
							break;
						case 'c':
						case 'C':
							fromNT=1;
							break;
						case 'g':
						case 'G':
							fromNT=2;
							break;
						case 't':
						case 'T':
							fromNT=3;
							break;
						default:
							fromNT=4;
							break;
					}

					/* New insertion */
					curScore=NEGATIVE_INFINITY;
					curLength=-1;
					/* Get NT and Color scores */
					curLength = matrix[i][j+1].s.length[fromNT] + 1;
					curScore = matrix[i][j+1].s.score[fromNT] + sm->gapOpenPenalty;
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
					curLength = matrix[i][j+1].v.length[fromNT] + 1;
					/* Insertion - previous row */
					curScore = matrix[i][j+1].v.score[fromNT] + sm->gapExtensionPenalty;
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
					matrix[i+1][j+1].v.score[k] = maxScore;
					matrix[i+1][j+1].v.from[k] = maxFrom;
					matrix[i+1][j+1].v.colorError[k] = maxColorError;
					matrix[i+1][j+1].v.length[k] = maxLength;
				}
			}
		}
	}

	FillAlignedEntryFromMatrixColorSpace(a,
			matrix,
			read,
			readLength,
			reference,
			referenceLength,
			readLength - maxV,
			position,
			strand,
			alphabetSize,
			0);

	/* Debug code */
	/*
	   AlignedEntry tmp;
	   AlignedEntryInitialize(&tmp);
	   AlignColorSpaceFull(read,
	   readLength,
	   reference,
	   referenceLength,
	   sm,
	   &tmp,
	   strand,
	   position);
	   if(a->score < tmp.score ||
	   tmp.score < a->score ||
	   !(a->length == tmp.length) ||
	   !(a->referenceLength == tmp.referenceLength)) {
	   fprintf(stderr, "\nreferenceLength=%d\n", referenceLength);
	   fprintf(stderr, "\nstrand=%c\n", strand);
	   AlignedEntryPrint(a,
	   stderr,
	   ColorSpace,
	   TextOutput);
	   AlignedEntryPrint(&tmp,
	   stderr,
	   ColorSpace,
	   TextOutput);
	   PrintError(FnName, NULL, "Alignments did not match", Exit, OutOfRange);
	   }
	   AlignedEntryFree(&tmp);
	   */
}

void FillAlignedEntryFromMatrixColorSpace(AlignedEntry *a,
		AlignMatrixCS **matrix,
		char *read,
		int readLength,
		char *reference,
		int referenceLength,
		int toExclude,
		uint32_t position,
		char strand,
		int alphabetSize,
		int debug)
{
	char *FnName="FillAlignedEntryFromMatrixColorSpace";
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
	for(i=toExclude+1;i<referenceLength+1;i++) {
		for(j=0;j<alphabetSize;j++) {
			/* Don't end with a Deletion in the read */

			/* End with a Match/Mismatch */
			if(maxScore < matrix[readLength][i].s.score[j]) {
				maxScore = matrix[readLength][i].s.score[j];
				startRow = readLength;
				startCol = i;
				startCell = j + 1 + (ALPHABET_SIZE + 1);
			}

			/* End with an Insertion */
			if(maxScore < matrix[readLength][i].v.score[j]) {
				maxScore = matrix[readLength][i].v.score[j];
				startRow = readLength;
				startCol = i;
				startCell = j + 1 + 2*(ALPHABET_SIZE + 1);
			}
		}
	}
	assert(startRow >= 0 && startCol >= 0 && startCell >= 0);

	/* Initialize variables for the loop */
	curRow=startRow;
	curCol=startCol;
	curFrom=startCell;

	referenceLengthAligned=0;
	/* Init */
	if(curFrom <= (ALPHABET_SIZE + 1)) {
		PrintError(FnName, "curFrom", "Cannot end with a deletion", Exit, OutOfRange);
		length = matrix[curRow][curCol].h.length[(curFrom - 1) % (ALPHABET_SIZE + 1)];
	}
	else if(2*(ALPHABET_SIZE + 1) < curFrom) {
		length = matrix[curRow][curCol].v.length[(curFrom - 1) % (ALPHABET_SIZE + 1)];
	}
	else {
		length = matrix[curRow][curCol].s.length[(curFrom - 1) % (ALPHABET_SIZE + 1)];
	}
	if(length < readLength) {
		fprintf(stderr, "\nlength=%d\nreadLength=%d\n", length, readLength);
	}
	assert(readLength <= length);
	i=length-1;

	/* Now trace back the alignment using the "from" member in the matrix */
	while(curRow > 0 && curCol > 0) {
		assert(i>=0);
		/* Where did the current cell come from */
		/* Get if there was a color error */
		if(curFrom <= (ALPHABET_SIZE + 1)) {
			/*
			   fprintf(stderr, "\ni=%d\ncurFrom=%d\nh.length=%d\n%s",
			   i,
			   curFrom,
			   matrix[curRow][curCol].h.length[(curFrom - 1) % (ALPHABET_SIZE + 1)],
			   BREAK_LINE);
			   assert(i + 1 == matrix[curRow][curCol].h.length[(curFrom - 1) % (ALPHABET_SIZE + 1)]);
			   */
			nextFrom = matrix[curRow][curCol].h.from[(curFrom - 1) % (ALPHABET_SIZE + 1)];
			colorErrorAligned[i] = matrix[curRow][curCol].h.colorError[(curFrom - 1) % (ALPHABET_SIZE + 1)];
		}
		else if(2*(ALPHABET_SIZE + 1) < curFrom) {
			/*
			   fprintf(stderr, "\ni=%d\ncurFrom=%d\nv.length=%d\n%s",
			   i,
			   curFrom,
			   matrix[curRow][curCol].v.length[(curFrom - 1) % (ALPHABET_SIZE + 1)],
			   BREAK_LINE);
			   assert(i + 1 == matrix[curRow][curCol].v.length[(curFrom - 1) % (ALPHABET_SIZE + 1)]);
			   */
			nextFrom = matrix[curRow][curCol].v.from[(curFrom - 1) % (ALPHABET_SIZE + 1)];
			colorErrorAligned[i] = matrix[curRow][curCol].v.colorError[(curFrom - 1) % (ALPHABET_SIZE + 1)];
		}
		else {
			/*
			   fprintf(stderr, "\ni=%d\ncurFrom=%d\ns.length=%d\n%s",
			   i,
			   curFrom,
			   matrix[curRow][curCol].s.length[(curFrom - 1) % (ALPHABET_SIZE + 1)],
			   BREAK_LINE);
			   assert(i + 1 == matrix[curRow][curCol].s.length[(curFrom - 1) % (ALPHABET_SIZE + 1)]);
			   */
			nextFrom = matrix[curRow][curCol].s.from[(curFrom - 1) % (ALPHABET_SIZE + 1)];
			colorErrorAligned[i] = matrix[curRow][curCol].s.colorError[(curFrom - 1) % (ALPHABET_SIZE + 1)];
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
		assert(curFrom > 0);
		curRow = nextRow;
		curCol = nextCol;
		i--;
	} /* End loop */
	if(-1!=i) {
		fprintf(stderr, "i=%d\n", i);
	}
	assert(-1==i);
	assert(length >= referenceLengthAligned);

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
