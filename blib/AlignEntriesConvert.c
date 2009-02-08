#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "config.h"
#include "AlignEntries.h"
#include "AlignEntry.h"
#include "BLib.h"
#include "BLibDefinitions.h"
#include "BError.h"
#include "AlignEntriesConvert.h"

/* TODO */
void AlignEntriesConvertPrintHeader(FILE *fp,
		int outputFormat) 
{
	char *FnName = "AlignEntriesConvertPrintHeader";
	switch(outputFormat) {
		case BAF:
			/* Do nothing */
			break;
		case MAF:
			if(0>fprintf(fp, "##maf version=%s scoring=%s\n",
						PACKAGE_VERSION,
						PROGRAM_NAME)) {
				PrintError(FnName,
						"header",
						"Could not write to file",
						Exit,
						WriteFileError);
			}
			break;
		case GFF:
			if(0>fprintf(fp, "##gff version=%s scoring=%s\n",
						PACKAGE_VERSION,
						PROGRAM_NAME)) {
				PrintError(FnName,
						"header",
						"Could not write to file",
						Exit,
						WriteFileError);
			}
			break;
		default:
			PrintError(FnName,
					"outputFormat",
					"Could not understand outputFormat",
					Exit,
					OutOfRange);
			break;
	}
}

/* TODO */
void AlignEntriesConvertPrintOutputFormat(AlignEntries *a, 
		RGBinary *rg,
		FILE *fp,
		int outputFormat,
		int binaryInput)
{
	char *FnName = "AlignEntriesConvertPrintOutputFormat";
	switch(outputFormat) {
		case BAF:
			AlignEntriesPrint(a, fp, binaryInput);
			break;
		case MAF:
			AlignEntriesConvertPrintMAF(a, rg, fp);
			break;
		case GFF:
			AlignEntriesConvertPrintGFF(a, fp);
			break;
		default:
			PrintError(FnName,
					"outputFormat",
					"Could not understand outputFormat",
					Exit,
					OutOfRange);
			break;
	}
}

/* TODO */
void AlignEntriesConvertPrintMAF(AlignEntries *a,
		RGBinary *rg,
		FILE *fp)
{
	int i;

	/* Get Data */
	if(0==a->pairedEnd) {
		for(i=0;i<a->numEntriesOne;i++) {
			AlignEntriesConvertPrintAlignEntryToMAF(&a->entriesOne[i], rg, a->readName, a->pairedEnd, a->space, 1, fp); 
		}
	}
	else {
		for(i=0;i<a->numEntriesOne;i++) {
			AlignEntriesConvertPrintAlignEntryToMAF(&a->entriesOne[i], rg, a->readName, a->pairedEnd, a->space, 1, fp); 
		}
		for(i=0;i<a->numEntriesTwo;i++) {
			AlignEntriesConvertPrintAlignEntryToMAF(&a->entriesTwo[i], rg, a->readName, a->pairedEnd, a->space, 2, fp); 
		}
	}

}

/* TODO */
void AlignEntriesConvertPrintAlignEntryToMAF(AlignEntry *a,
		RGBinary *rg,
		char *readName,
		int pairedEnd,
		int space,
		int readNum,
		FILE *fp)
{
	char *FnName="AlignEntriesConvertPrintAlignEntryToMAF";
	int i;
	int originalReferenceLength=0;
	int originalReadLength=0; 

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
		if(0>fprintf(fp, "a score=%lf paired-end=%d read=%d color-errors=%s contig-index=%d\n",
					a->score,
					pairedEnd,
					readNum,
					a->colorError,
					a->contig)) {
			PrintError(FnName,
					NULL,
					"Could not write to file",
					Exit,
					WriteFileError);
		}
	}
	else {
		assert(space == NTSpace);
		if(0>fprintf(fp, "a score=%lf paired-end=%d read=%d contig-index=%d\n",
					a->score,
					pairedEnd,
					readNum,
					a->contig)) {
			PrintError(FnName,
					NULL,
					"Could not write to file",
					Exit,
					WriteFileError);
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
		PrintError(FnName,
				NULL,
				"Could not write to file",
				Exit,
				WriteFileError);
	}
	/* Print the read */
	if(0>fprintf(fp, "s %s %u %d %c %d %s\n\n", /* Include a blank line */
				readName,
				0,
				originalReadLength, /* We align the full read */
				a->strand,
				originalReadLength, /* original read length */
				a->read)) {
		PrintError(FnName,
				NULL,
				"Could not write to file",
				Exit,
				WriteFileError);
	}
}

/* TODO */
void AlignEntriesConvertPrintGFF(AlignEntries *a,
		FILE *fp)
{
	int i;

	/* Get Data */
	if(0==a->pairedEnd) {
		for(i=0;i<a->numEntriesOne;i++) {
			AlignEntriesConvertPrintAlignEntryToGFF(&a->entriesOne[i], a->readName, a->pairedEnd, a->space, 1, fp); 
		}
	}
	else {
		for(i=0;i<a->numEntriesOne;i++) {
			AlignEntriesConvertPrintAlignEntryToGFF(&a->entriesOne[i], a->readName, a->pairedEnd, a->space, 1, fp); 
		}
		for(i=0;i<a->numEntriesTwo;i++) {
			AlignEntriesConvertPrintAlignEntryToGFF(&a->entriesTwo[i], a->readName, a->pairedEnd, a->space, 2, fp); 
		}
	}

}

/* TODO */
void AlignEntriesConvertPrintAlignEntryToGFF(AlignEntry *a,
		char *readName,
		int pairedEnd,
		int space,
		int readNum,
		FILE *fp)
{
	char *FnName="AlignEntriesConvertPrintAlignEntryToGFF";
	int i;
	int originalReferenceLength=0;
	int originalReadLength=0; 
	int initialized=0;
	char string[SEQUENCE_LENGTH]="\0";
	char tempString[SEQUENCE_LENGTH]="\0";
	uint8_t color;
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
	if(0>fprintf(fp, "%s\t%s\t%s\t%d\t%d\t%lf\t%c\t.\t",
				readName,
				PROGRAM_NAME,
				"read",
				a->position, /* one based */
				a->position + a->referenceLength-1,
				a->score,
				a->strand)) {
		PrintError(FnName,
				NULL,
				"Could not write to file",
				Exit,
				WriteFileError);
	}

	/* Write attributes */ 
	if(NTSpace == space) {
		/* b attribute - the base-space representation of the read */
		/* r attribute - the base-space representation of the reference */
		if(0>fprintf(fp, "b=%s;r=%s",
					a->read,
					a->reference)) {
			PrintError(FnName,
					NULL,
					"Could not write to file",
					Exit,
					WriteFileError);
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
			PrintError(FnName,
					NULL,
					"Could not write to file",
					Exit,
					WriteFileError);
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
					PrintError(FnName,
							NULL,
							"Could not convert bases to color space",
							Exit,
							OutOfRange);
				}
				string[i] = ConvertIntColorToCharColor(color);
				prevBase = a->read[i];
			}
		}
		string[a->length]='\0';
		if(0>fprintf(fp, ";g=%s", string)) {
			PrintError(FnName,
					NULL,
					"Could not write to file",
					Exit,
					WriteFileError);
		}

		/* i - ignore, since we only use one reference */

		/* p - ignore, since we do not define the mappability */

		/* q - ignore, since we do not keep track of quality values */

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
					PrintError(FnName,
							"tempString",
							"Could not write to string",
							Exit,
							OutOfRange);
				}
				strcat(string, tempString);
				initialized=1;
			}
		}
		if(1==initialized) {
			if(0>fprintf(fp, ";r=%s", string)) {
				PrintError(FnName,
						NULL,
						"Could not write to file",
						Exit,
						WriteFileError);
			}
		}

		/* s - ignore, since we don't use color rules to align our data */

		/* u - ignore, since we don't define mappability */
	}
	/* Print new-line */
	if(0>fprintf(fp, "\n")) {
		PrintError(FnName,
				NULL,
				"Could not write to file",
				Exit,
				WriteFileError);
	}
}
