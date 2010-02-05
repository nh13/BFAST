#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>

#include "config.h"
#include "AlignedRead.h"
#include "AlignedEnd.h"
#include "AlignedEntry.h"
#include "BLib.h"
#include "BLibDefinitions.h"
#include "BError.h"
#include "AlignedReadConvert.h"

/* TODO */
void AlignedReadConvertPrintHeader(FILE *fp,
		RGBinary *rg,
		int32_t outputFormat,
		char *readGroup
		) 
{
	char *FnName = "AlignedReadConvertPrintHeader";
	int32_t i;

	switch(outputFormat) {
		case BAF:
			/* Do nothing */
			break;
		case MAF:
			if(0>fprintf(fp, "##maf version=%s scoring=%s\n",
						PACKAGE_VERSION,
						PROGRAM_NAME)) {
				PrintError(FnName, "header", "Could not write to file", Exit, WriteFileError);
			}
			break;
		case GFF:
			if(0>fprintf(fp, "##gff version=%s scoring=%s\n",
						PACKAGE_VERSION,
						PROGRAM_NAME)) {
				PrintError(FnName, "header", "Could not write to file", Exit, WriteFileError);
			}
			break;
		case SAM:
			/* Header */
			if(0>fprintf(fp, "@HD\tVN:%s\tSO:unsorted\tGO:none\n",
						BFAST_SAM_VERSION)) {
				PrintError(FnName, "header", "Could not write to file", Exit, WriteFileError);
			}
			/* Sequence dictionary */
			for(i=0;i<rg->numContigs;i++) {
				if(0>fprintf(fp, "@SQ\tSN:%s\tLN:%d\n",
							rg->contigs[i].contigName,
							rg->contigs[i].sequenceLength)) {
					PrintError(FnName, "header", "Could not write to file", Exit, WriteFileError);
				}
			}
			/* Print read group */
			if(NULL != readGroup && 0>fprintf(fp, "%s\n", readGroup)) {
				PrintError(FnName, "header", "Could not write to file", Exit, WriteFileError);
			}
			/* Program */
			if(0>fprintf(fp, "@PG\tID:%s\tVN:%s\n",
						PACKAGE_NAME,
						PACKAGE_VERSION)) {
				PrintError(FnName, "header", "Could not write to file", Exit, WriteFileError);
			}
			break;
		default:
			PrintError(FnName, "outputFormat", "Could not understand outputFormat", Exit, OutOfRange);
			break;
	}
}

/* TODO */
void AlignedReadConvertPrintOutputFormat(AlignedRead *a, 
		RGBinary *rg,
		FILE *fp,
		gzFile fpGZ,
		char *outputID,
		char *readGroupString,
		int32_t postprocessAlgorithm,
		int32_t *numOriginalEntries,
		int32_t outputFormat,
		int32_t binaryOutput)
{
	char *FnName = "AlignedReadConvertPrintOutputFormat";
	switch(outputFormat) {
		case BAF:
			if(BinaryOutput == binaryOutput) {
				AlignedReadPrint(a, fpGZ);
			}
			else {
				AlignedReadPrintText(a, fp);
			}
			break;
		case MAF:
			AlignedReadConvertPrintMAF(a, rg, fp);
			break;
		case GFF:
			AlignedReadConvertPrintGFF(a, fp);
			break;
		case SAM:
			AlignedReadConvertPrintSAM(a, postprocessAlgorithm, numOriginalEntries, outputID, readGroupString, fp);
			break;
		default:
			PrintError(FnName, "outputFormat", "Could not understand outputFormat", Exit, OutOfRange);
			break;
	}
}

/* TODO */
void AlignedReadConvertPrintMAF(AlignedRead *a,
		RGBinary *rg,
		FILE *fp)
{
	int32_t i, j;

	/* Get Data */
	for(i=0;i<a->numEnds;i++) {
		for(j=0;j<a->ends[i].numEntries;j++) {
			AlignedReadConvertPrintAlignedEntryToMAF(&a->ends[i].entries[j], 
					rg, 
					a->readName, 
					a->ends[i].qual,
					i+1,
					a->space, 
					j+1,
					fp); 
		}
	}
}

/* TODO */
void AlignedReadConvertPrintAlignedEntryToMAF(AlignedEntry *a,
		RGBinary *rg,
		char *readName,
		char *qual,
		int32_t whichEnd,
		int32_t space,
		int32_t alignmentNum,
		FILE *fp)
{
	char *FnName="AlignedReadConvertPrintAlignedEntryToMAF";
	int32_t i;
	int32_t originalReferenceLength=0;
	int32_t originalReadLength=0; 

	/* Recover original lengths */
	for(i=0;i<a->length;i++) {
		if(a->reference[i] != GAP) {
			originalReferenceLength++;
		}
		if(a->read[i] != GAP) {
			originalReadLength++;
		}
	}
	assert(originalReferenceLength == a->referenceLength);

	/* Print the score */
	if(space == ColorSpace) {
		if(0>fprintf(fp, "a score=%lf which-end=%d alignment-num=%d color-errors=%s contig-index=%d\n",
					a->score,
					whichEnd,
					alignmentNum,
					a->colorError,
					a->contig)) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}
	else {
		assert(space == NTSpace);
		if(0>fprintf(fp, "a score=%lf which-end=%d alignment-num=%d contig-index=%d\n",
					a->score,
					whichEnd,
					alignmentNum,
					a->contig)) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}

	/* Make sure the contig reported is within bounds of the reference genome */
	assert(1 <= a->contig && a->contig <= rg->numContigs);

	/* Print the reference */
	if(0>fprintf(fp, "s %s %u %d %c %d %s\n",
				a->contigName,
				a->position-1, /* zero based */
				originalReferenceLength,
				a->strand,
				rg->contigs[a->contig-1].sequenceLength, /* original contig length */
				a->reference)) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
	/* Print the read */
	if(0>fprintf(fp, "s %s %u %d %c %d %s\n\n", /* Include a blank line */
				readName,
				0,
				originalReadLength, /* We align the full read */
				a->strand,
				originalReadLength, /* original read length */
				a->read)) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
	/* Print the qualities */
	if(0>fprintf(fp, "q %s ",
				a->contigName)) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
	for(i=0;i<strlen(qual);i++) {
		if(0>fprintf(fp, "%1d",
					QUAL_TO_MAF_QUAL(CHAR2QUAL(qual[i])))) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}
	if(0>fprintf(fp, "\n")) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
}

/* TODO */
void AlignedReadConvertPrintGFF(AlignedRead *a,
		FILE *fp)
{
	int32_t i, j;

	/* Get Data */
	for(i=0;i<a->numEnds;i++) {
		for(j=0;j<a->ends[i].numEntries;j++) {
			/* Get Data */
			AlignedReadConvertPrintAlignedEntryToGFF(&a->ends[i].entries[j], 
					a->readName, 
					a->ends[i].qual,
					i+1,
					a->space, 
					j+1,
					fp); 
		}
	}
}

/* TODO */
void AlignedReadConvertPrintAlignedEntryToGFF(AlignedEntry *a,
		char *readName,
		char *qual,
		int32_t whichEnd,
		int32_t space,
		int32_t alignmentNum,
		FILE *fp)
{
	char *FnName="AlignedReadConvertPrintAlignedEntryToGFF";
	int32_t i;
	int32_t originalReferenceLength=0;
	int32_t originalReadLength=0; 
	int32_t initialized=0;
	char string[SEQUENCE_LENGTH]="\0";
	char tempString[SEQUENCE_LENGTH]="\0";
	char color;
	char prevBase;

	/* Recover original lengths */
	for(i=0;i<a->length;i++) {
		if(a->reference[i] != GAP) {
			originalReferenceLength++;
		}
		if(a->read[i] != GAP) {
			originalReadLength++;
		}
	}
	assert(originalReferenceLength == a->referenceLength);

	/* Write fields */
	if(0>fprintf(fp, "%s\t%s\t%s\t%d\t%d\t%lf\t%c\t.\twhich_end=%d;",
				readName,
				PROGRAM_NAME,
				"read",
				a->position, /* one based */
				a->position + a->referenceLength-1,
				a->score,
				a->strand,
				alignmentNum)) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}

	/* Write attributes */ 
	if(NTSpace == space) {
		/* b attribute - the base-space representation of the read */
		/* r attribute - the base-space representation of the reference */
		if(0>fprintf(fp, "b=%s;r=%s",
					a->read,
					a->reference)) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
		/* Print the qualities */
		if(0>fprintf(fp, ";q=")) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
		for(i=0;i<strlen(qual);i++) {
			if(0>fprintf(fp, "%d",
						QUAL_TO_MAF_QUAL(CHAR2QUAL(qual[i])))) {
				PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
			}
			if(i<strlen(qual)-1) {
				if(0>fprintf(fp, ",")) {
					PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
				}
			}
		}
	}
	else {
		/* Write ABI SOLiD attributes */

		/* b attribute - the corrected base-space representation of the read. */
		assert(a->length < SEQUENCE_LENGTH);
		for(i=0;i<a->length;i++) {
			if(a->read[i] == GAP ||
					ToLower(a->read[i]) == ToLower(a->reference[i])) {
				string[i] = ToUpper(a->read[i]);
			}
			else {
				string[i] = ToLower(a->read[i]);
			}
		}
		string[a->length]='\0';
		if(0>fprintf(fp, "b=%s", string)) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}

		/* c - ignore, since we may not have unique alignments */

		/* g - the color-sace string for this read, written from 5' to 3'.  Also,
		 * convert the adaptor and color to base space.  This is really stupid, since you
		 * should keep the same format as the ABI csfasta files, but typical ABI: don't be
		 * consistent.
		 * */
		prevBase = string[0] = a->read[0];
		for(i=1;i<a->length;i++) {
			if(GAP == a->read[i]) {
				string[i] = GAP;
			}
			else {
				if(1!=ConvertBaseToColorSpace(prevBase, a->read[i], &color)) {
					PrintError(FnName, NULL, "Could not convert bases to color space", Exit, OutOfRange);
				}
				string[i] = ConvertIntColorToCharColor(color);
				prevBase = a->read[i];
			}
		}
		string[a->length]='\0';
		if(0>fprintf(fp, ";g=%s", string)) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}

		/* i - ignore, since we only use one reference */

		/* p - ignore, since we do not define the mappability */

		/* q - qualities */
		if(0>fprintf(fp, ";q=")) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
		for(i=0;i<strlen(qual);i++) {
			if(0>fprintf(fp, "%d",
						QUAL_TO_MAF_QUAL(CHAR2QUAL(qual[i])))) {
				PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
			}
			if(i<strlen(qual)-1) {
				if(0>fprintf(fp, ",")) {
					PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
				}
			}
		}

		/* r - a comma separated list of {position}_{ref_color} for all of the color
		 * calls in the read sequence that differ from the reference sequence. The position
		 * is 1-based */
		initialized=0;
		string[0]='\0';
		for(i=0;i<a->length;i++) {
			if(GAP != a->colorError[i]) {
				if(0 != initialized) {
					strcat(string, ",");
				}
				if(0>sprintf(tempString, "%d_%c", 
							i+1,
							a->colorError[i])) {
					PrintError(FnName, "tempString", "Could not write to string", Exit, OutOfRange);
				}
				strcat(string, tempString);
				initialized=1;
			}
		}
		if(1==initialized) {
			if(0>fprintf(fp, ";r=%s", string)) {
				PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
			}
		}

		/* s - ignore, since we don't use color rules to align our data */

		/* u - ignore, since we don't define mappability */
	}
	/* Print new-line */
	if(0>fprintf(fp, "\n")) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
}

/* TODO */
void AlignedReadConvertPrintSAM(AlignedRead *a,
		int32_t postprocessAlgorithm,
		int32_t *numOriginalEntries,
		char *outputID,
		char *readGroupString,
		FILE *fp)
{
	char *FnName="AlignedReadConvertPrintSAM";
	int32_t i, j;
	/* Assumes that one end is mapped */

	/* SAM can't deal with generalized multi-end reads */
	if(2 < a->numEnds) {
		PrintError(FnName, NULL, "Outputting reads with greater than two ends to SAM format not supported. Skipping...", Warn, OutOfRange);
		return;
	}

	/* Get Data */
	for(i=0;i<a->numEnds;i++) {
		if(0 == a->ends[i].numEntries) { /* Unmapped read */
			AlignedReadConvertPrintAlignedEntryToSAM(a,
					i,
					-1,
					postprocessAlgorithm,
					numOriginalEntries,
					outputID,
					readGroupString,
					fp);
		}
		else {
			for(j=0;j<a->ends[i].numEntries;j++) {
				AlignedReadConvertPrintAlignedEntryToSAM(a,
						i,
						j,
						postprocessAlgorithm,
						numOriginalEntries,
						outputID,
						readGroupString,
						fp);
			}
		}
	}
}

/* TODO */
void AlignedReadConvertPrintAlignedEntryToSAM(AlignedRead *a,
		int32_t endIndex,
		int32_t entriesIndex,
		int32_t postprocessAlgorithm,
		int32_t *numOriginalEntries,
		char *outputID,
		char *readGroupString,
		FILE *fp) 
{
	char *FnName="AlignedReadConvertPrintAlignedEntryToSAM";
	int32_t i, j;
	uint64_t flag;
	int32_t mateEndIndex, mateEntriesIndex, mapq;
	int32_t numEdits=0;
	char read[SEQUENCE_LENGTH]="\0";
	char readRC[SEQUENCE_LENGTH]="\0";
	char qual[SEQUENCE_LENGTH]="\0";
	char qualRC[SEQUENCE_LENGTH]="\0";
	char colorError[SEQUENCE_LENGTH]="\0";
	char MD[SEQUENCE_LENGTH]="*";

	/* Get mate end and mate index if they exist */
	mateEndIndex=mateEntriesIndex=-1;
	for(i=0;mateEndIndex < 0 && i < a->numEnds;i++) { /* Try other ends */
		if(endIndex != i && 0 <a->ends[i].numEntries) {
			mateEndIndex=i;
			mateEntriesIndex=0;
		}
	}

	/* QNAME */
	assert(strlen(outputID) + strlen(a->readName) < BFAST_SAM_MAX_QNAME); /* One less for separator */
	if(0 < strlen(outputID)) {
		if(0>fprintf(fp, "%s%s%s",
					outputID,
					BFAST_SAM_MAX_QNAME_SEPARATOR,
					a->readName)) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}
	else {
		if(0>fprintf(fp, "%s",
					a->readName)) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}
	/* FLAG */
	flag = 0;
	if(2 == a->numEnds) {
		flag |= 0x0001; /* Paired end */
		flag |= 0x0002; /* Always a proper pair */
		if(mateEndIndex < 0) {
			/* Other end is unmapped */
			flag |= 0x0008;
		}
		else {
			/* Other end is mapped */
			flag |= (REVERSE == a->ends[mateEndIndex].entries[mateEntriesIndex].strand)?0x0020:0x0000; /* Strand of the mate */
		}
		flag |= (0 == endIndex)?0x0040:0x0080; /* Which end */
	}
	if(entriesIndex < 0) { /* Unmapped */
		flag |= 0x0004;
	}
	if(0 < entriesIndex ) {
		flag |= 0x0100; /* This read is not primary */
	}
	if(0 <= entriesIndex) { /* Mapped */
		flag |= (REVERSE==a->ends[endIndex].entries[entriesIndex].strand)?0x0010:0x0000;
	}
	if(0>fprintf(fp, "\t%llu",
				(unsigned long long int)flag)) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
	/* RNAME and POS */
	if(entriesIndex < 0) { /* Current is unmapped */
		/* Use mate */
		if(0 <= mateEndIndex) {
			if(0>fprintf(fp, "\t%s\t%d",
						a->ends[mateEndIndex].entries[mateEntriesIndex].contigName,
						a->ends[mateEndIndex].entries[mateEntriesIndex].position)) {
				PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
			}
		}
		else {
			/* Make absent */ 
			if(0>fprintf(fp, "\t*\t0")) {
				PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
			}
		}
	}
	else {
		if(0>fprintf(fp, "\t%s\t%d",
					a->ends[endIndex].entries[entriesIndex].contigName,
					a->ends[endIndex].entries[entriesIndex].position)) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}
	/* MAPQ */
	if(entriesIndex < 0) {
		mapq = 0;
	}
	else {
		mapq = (int32_t)a->ends[endIndex].entries[entriesIndex].mappingQuality;
	}
	if(mapq < 0) mapq = 0;
	if(mapq > MAXIMUM_MAPPING_QUALITY) mapq = MAXIMUM_MAPPING_QUALITY;
	if(0>fprintf(fp, "\t%d", mapq)) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
	/* CIGAR - no alignment */
	if(entriesIndex < 0) { /* Unmapped */
		if(0>fprintf(fp, "\t*")) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}
	else {
		AlignedReadConvertPrintAlignedEntryToCIGAR(&a->ends[endIndex].entries[entriesIndex], a->space, colorError, MD, &numEdits, fp);
	}
	/* MRNM and MPOS */
	if(2 == a->numEnds) {
		if(0 <= mateEndIndex) {
			if(0 <= entriesIndex &&
					0 == strcmp(a->ends[mateEndIndex].entries[mateEntriesIndex].contigName, a->ends[endIndex].entries[entriesIndex].contigName)) {
			if(0>fprintf(fp, "\t=\t%d",
						a->ends[mateEndIndex].entries[mateEntriesIndex].position)) {
				PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
			}
			}
			else {
			if(0>fprintf(fp, "\t%s\t%d",
						a->ends[mateEndIndex].entries[mateEntriesIndex].contigName,
						a->ends[mateEndIndex].entries[mateEntriesIndex].position)) {
				PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
			}
			}
		}
		else {
			/* Use contig current */ 
			if(entriesIndex < 0) { /* Current is unmapped */
				/* Make absent */ 
				if(0>fprintf(fp, "\t*\t0")) {
					PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
				}
			}
			else { /* Current is mapped */
				if(0>fprintf(fp, "\t=\t%d",
							a->ends[endIndex].entries[entriesIndex].position)) {
					PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
				}
			}
		}
	}
	else {
		if(0>fprintf(fp, "\t*\t0")) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}
	/* ISIZE */
	if(entriesIndex < 0 || /* Unmapped */
			mateEndIndex < 0 || /* Mate is unmapped */
			a->ends[endIndex].entries[entriesIndex].contig != a->ends[mateEndIndex].entries[mateEntriesIndex].contig) {
		if(0>fprintf(fp, "\t0")) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}
	else {
		if(0>fprintf(fp, "\t%d",
					a->ends[mateEndIndex].entries[mateEntriesIndex].position -
					a->ends[endIndex].entries[entriesIndex].position)) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}
	/* SEQ and QUAL */
	if(NTSpace == a->space) {
		if(0 <= entriesIndex && /* Was mapped */
				REVERSE == a->ends[endIndex].entries[entriesIndex].strand) {
			/* Reverse compliment */
			GetReverseComplimentAnyCase(a->ends[endIndex].read,
					read,
					strlen(a->ends[endIndex].read));
			ReverseRead(a->ends[endIndex].qual,
					qual,
					strlen(a->ends[endIndex].qual));
		}
		else {
			strcpy(read, a->ends[endIndex].read);
			strcpy(qual, a->ends[endIndex].qual);
		}
		assert(strlen(qual) == strlen(read));
	}
	else {
		/* Convert read to NT space */
		if(entriesIndex < 0) { /* Unmapped */
			/* Just decode original color space read */
			strcpy(read, a->ends[endIndex].read);
			assert(0 < ConvertReadFromColorSpace(read, strlen(read)));
			/* Convert quals to NT Space */
			for(i=0;i<strlen(a->ends[endIndex].qual);i++) {
				if(0 == i) {
					qual[i] = CHAR2QUAL(a->ends[endIndex].qual[i]);
				}
				else {
					/* How do we determine this? This does not make sense but for now
					 * SAM requires it. For now we will take the average */
					if(0 == CHAR2QUAL(a->ends[endIndex].qual[i-1]) ||
							0 == CHAR2QUAL(a->ends[endIndex].qual[i])) {
						qual[i] = 0; // Default to 0 even though we may be able to recover?
					}
					else {
						qual[i] = (int8_t)(-10*(AddLog10(CHAR2QUAL(a->ends[endIndex].qual[i-1])/-10.0, 
									CHAR2QUAL(a->ends[endIndex].qual[i])/-10.0) - log10(2.0)) + 0.5);
						qual[i] = QUAL2CHAR(qual[i]);
					}
				}
				if(qual[i] <= 0) {
					qual[i] = QUAL2CHAR(0);
				}
				else if(qual[i] > 63) {
					qual[i] = QUAL2CHAR(63);
				}
				else {
					qual[i] = QUAL2CHAR(qual[i]);
				}
			}
			qual[i]='\0';
		}
		else { /* Mapped */
			/* Remove gaps from the read (deletions) */
			for(i=j=0;i<strlen(a->ends[endIndex].entries[entriesIndex].read);i++) {
				if(GAP != a->ends[endIndex].entries[entriesIndex].read[i]) {
					read[j] = a->ends[endIndex].entries[entriesIndex].read[i];
					j++;
				}
			}
			read[j]='\0';
			/* Convert quals to NT Space - use MAQ 0.7.1 conversion */
			for(i=j=0;i<a->ends[endIndex].entries[entriesIndex].length;i++) {
				if(GAP != a->ends[endIndex].entries[entriesIndex].read[i]) { /* Not a deletion */
					if(a->ends[endIndex].entries[entriesIndex].length - 1 == i) { /* At the end of the alignment */
						assert(j==a->ends[endIndex].qualLength-1);
						qual[j] = CHAR2QUAL(a->ends[endIndex].qual[j]);
					}
					else if(GAP == a->ends[endIndex].entries[entriesIndex].colorError[i] &&
							GAP == a->ends[endIndex].entries[entriesIndex].colorError[i+1]) {
						qual[j] = CHAR2QUAL(a->ends[endIndex].qual[j]) + 
							CHAR2QUAL(a->ends[endIndex].qual[j+1]) + 10;
					}
					else if(GAP == a->ends[endIndex].entries[entriesIndex].colorError[i]) {
						qual[j] = CHAR2QUAL(a->ends[endIndex].qual[j]) - 
							CHAR2QUAL(a->ends[endIndex].qual[j+1]);
					}
					else if(GAP == a->ends[endIndex].entries[entriesIndex].colorError[i+1]) {
						qual[j] = CHAR2QUAL(a->ends[endIndex].qual[j+1]) - 
							CHAR2QUAL(a->ends[endIndex].qual[j]);
					}
					else {
						qual[j] = 0;
					}
					/* Round */
					if(qual[j] <= 0) {
						qual[j] = QUAL2CHAR(1);
					}
					else if(qual[j] > 63) {
						qual[j] = QUAL2CHAR(63);
					}
					else {
						qual[j] = QUAL2CHAR(qual[j]);
					}
					j++;
				}
			}
			qual[j]='\0';
			if(REVERSE == a->ends[endIndex].entries[entriesIndex].strand) {
				/* Reverse compliment */
				GetReverseComplimentAnyCase(read, /* src */
						readRC, /* dest */
						strlen(read));
				strcpy(read, readRC);
				ReverseRead(qual, /* src */
						qualRC, /* dest */
						strlen(qual));
				strcpy(qual, qualRC);
			}
		}
		assert(strlen(qual) == strlen(read));
	}
	if(0>fprintf(fp, "\t%s\t%s",
				read,
				qual)) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
	/* RG - optional field */
	/* LB - optional field */
	/* PU - optional field */
	if(NULL != readGroupString && 0>fprintf(fp, "%s", readGroupString)) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
	/* PG - optional field */
	if(0>fprintf(fp, "\tPG:Z:%s", PACKAGE_NAME)) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
	/* AS - optional field */
	if(entriesIndex < 0) { /* Unmapped */
		if(0>fprintf(fp, "\tAS:i:%d", INT_MIN)) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}
	else {
		if(0>fprintf(fp, "\tAS:i:%d", (int32_t)a->ends[endIndex].entries[entriesIndex].score)) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}
	/* MQ - optional field */
	if(2 == a->numEnds && 0 <= mateEndIndex) {
		if(0>fprintf(fp, "\tMQ:i:%d",
					a->ends[mateEndIndex].entries[mateEntriesIndex].mappingQuality)) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}
	/* NM - optional field */
	if(0>fprintf(fp, "\tNM:i:%d", numEdits)) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
	/* NH - optional field */
	if(0>fprintf(fp, "\tNH:i:%d",
				(NULL == numOriginalEntries) ? ((entriesIndex < 0) ? 1:a->ends[endIndex].numEntries) : numOriginalEntries[endIndex])) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
	/* IH - optional field */
	if(0>fprintf(fp, "\tIH:i:%d",
				(entriesIndex < 0)?1:a->ends[endIndex].numEntries)) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
	/* HI - optional field */
	if(0>fprintf(fp, "\tHI:i:%d",
				(entriesIndex < 0)?1:(entriesIndex+1))) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
	/* MD - optional field */
	if(0>fprintf(fp, "\tMD:Z:%s", MD)) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
	/* CS, CQ and CM - optional fields */
	if(ColorSpace == a->space) {
		int32_t numCM=0;
		if(0 <=entriesIndex) {
			for(i=0;i<a->ends[endIndex].entries[entriesIndex].length;i++) {
				if(GAP != a->ends[endIndex].entries[entriesIndex].colorError[i]) numCM++;
			}
		}
		if(0>fprintf(fp, "\tCS:Z:%s\tCQ:Z:%s\tCM:i:%d",
					a->ends[endIndex].read,
					a->ends[endIndex].qual,
					numCM)) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}
	/* CC - optional field */
	/* CP - optional field */
	if(entriesIndex < 0 || /* Unmapped */
			entriesIndex == a->ends[endIndex].numEntries-1) { /* Last hit */
		/* Leave empty */
	}
	else {
		if(0>fprintf(fp, "\tCC:Z:%s\tCP:i:%d",
					a->ends[endIndex].entries[entriesIndex+1].contigName,
					a->ends[endIndex].entries[entriesIndex+1].position)) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}
	/* BFAST specific fields */
	if(0 <= postprocessAlgorithm && 0>fprintf(fp, "\tXA:i:%d", postprocessAlgorithm)) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
	if(ColorSpace == a->space && 0 < strlen(colorError)) {
		if(0>fprintf(fp, "\tXE:Z:%s", colorError)) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}

	if(0>fprintf(fp, "\n")) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}
}

/* TODO */
void AlignedReadConvertPrintAlignedEntryToCIGAR(AlignedEntry *a,
		int32_t space,
		char *colorError,
		char *MD,
		int32_t *numEdits,
		FILE *fp)
{
	char *FnName="AlignedReadConvertPrintAlignedEntryToCIGAR";
	char read[SEQUENCE_LENGTH]="\0";
	char reference[SEQUENCE_LENGTH]="\0";
	int32_t i, MDi, MDNumMatches=0;
	int32_t prevType=0;
	int32_t numPrevType=0;
	int32_t curType=0;
	int32_t startDel, endDel, startIns, endIns, prevDel, prevIns;

	(*numEdits) = 0;

	if(0>fprintf(fp, "\t")) {
		PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
	}

	if(REVERSE == a->strand) {
		GetReverseComplimentAnyCase(a->read, read, a->length);
		GetReverseComplimentAnyCase(a->reference, reference, a->length);
		if(ColorSpace == space) {
			ReverseRead(a->colorError, colorError, a->length);
		}
	}
	else {
		strcpy(read, a->read);
		strcpy(reference, a->reference);
		if(ColorSpace == space) {
			strcpy(colorError, a->colorError);
		}
	}

	/* Move all insertions and deletions to the 5' end - cool*/
	i=0;
	prevDel = prevIns = 0;
	startDel = endDel = startIns = endIns = -1;
	while(i<a->length) {
		assert(0 == prevIns || 0 == prevDel);

		if(GAP == read[i]) {
			if(0 == prevDel) {
				startDel = i;
			}
			prevDel = 1;
			endDel = i;
			prevIns = 0;
			startIns = -1;
			endIns = -1;
			i++;
		}
		else if(GAP == reference[i]) {
			if(0 == prevIns) {
				startIns = i;
			}
			prevIns = 1;
			endIns = i;
			prevDel = 0;
			startDel = -1;
			endDel = -1;
			i++;
		}
		else {
			if(1 == prevDel) {
				assert(0 < startDel);
				assert(startDel <= endDel);
				startDel--;
				while(0 <= startDel && /* Bases remaining to examine */
						read[startDel] != GAP && /* Hit another deletion */
						reference[startDel] != GAP && /* Hit an insertion */
						reference[startDel] == reference[endDel]) { /* src ref base matches dest ref base */
					assert(GAP != reference[startDel]);
					assert(GAP != reference[endDel]);
					assert(GAP != read[startDel]);
					assert(GAP == read[endDel]);
					read[endDel] = read[startDel];
					read[startDel] = GAP;
					if(ColorSpace == space) {
						assert(GAP == colorError[endDel]); /* No color errors in the deletion */
						colorError[endDel] = colorError[startDel];
						colorError[startDel] = GAP;
					}
					startDel--;
					endDel--;
				}
				endDel++; /* We decremented when we exited the loop */
				i = endDel;
				assert(GAP != read[i]);
				assert(GAP != reference[i]);
			}
			else if(1 == prevIns) {
				assert(startIns <= endIns);
				startIns--;
				while(0 <= startIns && /* Bases remaining to examine */
						read[startIns] != GAP && /* Hit another deletion */
						reference[startIns] != GAP && /* Hit an insertion */
						read[startIns] == read[endIns]) { /* src read base matches dest read base */
					assert(GAP != read[startIns]);
					assert(GAP != read[endIns]);
					assert(GAP != reference[startIns]);
					assert(GAP == reference[endIns]);
					reference[endIns] = reference[startIns];
					reference[startIns] = GAP;
					if(ColorSpace == space) {
						assert(GAP == colorError[endIns]);
						colorError[endIns] = colorError[startIns];
						colorError[startIns] = GAP;
					}
					startIns--;
					endIns--;
				}
				endIns++; /* We decremented when we exited the loop */
				i = endIns;
				assert(GAP != read[i]);
				assert(GAP != reference[i]);
			}
			else {
				i++;
			}
			prevDel = 0;
			prevIns = 0;
			startDel = -1;
			endDel = -1;
			startIns = -1;
			endIns = -1;
		}
	}

	/* Convert to cigar format */
	prevType = 0; /* 0 - MM, 1 - I, 2 - D */
	numPrevType = 0;
	MDi = MDNumMatches = 0;
	MD[MDi]='\0';
	for(i=0;i<a->length;i++) {
		assert(0 == prevIns || 0 == prevDel);

		if(GAP == read[i]) { // Deletion
			if(0 < MDNumMatches && sprintf(MD, "%s%d", MD, MDNumMatches) < 0) {
				PrintError(FnName, "MD", "Could not create string", Exit, OutOfRange);
			}
			MDi+=(int)(1+log10(0.1+MDNumMatches));
			MDNumMatches = 0;
			curType = 2;
			(*numEdits)++;
			// ignore MD
		}
		else if(GAP == reference[i]) { // Insertion
			if(0 < MDNumMatches && sprintf(MD, "%s%d", MD, MDNumMatches) < 0) {
				PrintError(FnName, "MD", "Could not create string", Exit, OutOfRange);
			}
			MDi+=(int)(1+log10(0.1+MDNumMatches));
			MDNumMatches = 0;
			curType = 1;
			(*numEdits)++;
			if(1 != prevType) {
				MD[MDi]='^';
				MDi++;
			}
			MD[MDi]=ToUpper(read[i]);
			MDi++;
		}
		else { // Match/Mismatch
			if(ToUpper(read[i]) != ToUpper(reference[i])) {
				(*numEdits)++;
				if(0 < MDNumMatches && sprintf(MD, "%s%d", MD, MDNumMatches) < 0) {
					PrintError(FnName, "MD", "Could not create string", Exit, OutOfRange);
				}
				MDi+=(int)(1+log10(0.1+MDNumMatches));
				MDNumMatches = 0;
				MD[MDi]=ToUpper(reference[i]);
				MDi++;
			}
			else {
				MDNumMatches++;
			}
			curType = 0;
		}

		if(curType == prevType) {
			numPrevType++;
		}
		else {
			if(0 < numPrevType) {
				assert(0 <= curType && curType <= 2);
				if(0>fprintf(fp, "%d%c",
							numPrevType,
							"MID"[prevType])) {
					PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
				}
			}
			prevType = curType;
			numPrevType = 1;
		}
	}
	if(0 < numPrevType) {
		assert(0 <= prevType && prevType <= 2);
		if(0>fprintf(fp, "%d%c",
					numPrevType,
					"MID"[prevType])) {
			PrintError(FnName, NULL, "Could not write to file", Exit, WriteFileError);
		}
	}
	if(0 < MDNumMatches && sprintf(MD, "%s%d", MD, MDNumMatches) < 0) {
		PrintError(FnName, "MD", "Could not create string", Exit, OutOfRange);
	}
	MDi+=(int)(1+log10(0.1+MDNumMatches));
	MDNumMatches=0;
	MD[MDi]='\0';
}
