#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/AlignEntries.h"
#include "bevalsim.h"

#define Name "bevalsim"

int main(int argc, char *argv[]) 
{
	char baf[MAX_FILENAME_LENGTH]="\0";

	if(argc == 2) {

		/* Get cmd line options */
		strcpy(baf, argv[1]);

		/* Check cmd line options */

		/* Run program */
		Evaluate(baf);

		/* Terminate */
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully.\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast aligned file>\n");
	}
	return 0;
}

void ReadTypeInitialize(ReadType *r)
{
	r->strand=0;
	r->chr=0;
	r->pos=0;
	r->pairedEnd=0;
	r->pairedEndLength=0;
	r->readLength=0;
	r->whichReadVariants=0;
	r->startIndel=0;
	r->indelLength=0;
	r->numSNPs=0;
	r->numErrors=0;
	r->deletionLength=0;
	r->insertionLength=0;
}

int ReadTypeReadFromBAF(ReadType *r, 
		FILE *fp)
{
	char *FnName = "ReadTypeReadFromBAF";
	AlignEntries a;
	int i;
	char r1[SEQUENCE_LENGTH]="\0";
	char r2[SEQUENCE_LENGTH]="\0";

	/* Initialize */
	AlignEntriesInitialize(&a);

	/* Read in align entries */
	if(EOF==AlignEntriesRead(&a, fp, PairedEndDoesNotMatter, SpaceDoesNotMatter)) {
		return EOF;
	}
	/* There should be only one */
	if(a.numEntriesOne > 1 ||
			(a.pairedEnd == 1 && a.numEntriesOne > 1)) {
		PrintError(FnName,
				NULL,
				"There was more than one alignment for a given read",
				Exit,
				OutOfRange);
	}

	/* Convert into read type */
	int tempPairedEnd = a.pairedEnd;
	r->space = a.colorSpace;
	/* Get the rest from read name */
	if(tempPairedEnd == 0) {
		if(EOF == sscanf(a.readName, 
					">strand=%c_chr=%d_pos=%d_pe=%d_pel=%d_rl=%d_wrv=%d_si=%d_il=%d_r1=%s",
					&r->strand,
					&r->chr,
					&r->pos,
					&r->pairedEnd,
					&r->pairedEndLength,
					&r->readLength,
					&r->whichReadVariants,
					&r->startIndel,
					&r->indelLength,
					r1)) {
			PrintError(FnName,
					a.readName,
					"Could not parse read name (0)",
					Exit,
					OutOfRange);
		}
	}
	else {
		if(EOF == sscanf(a.readName, 
					">strand=%c_chr=%d_pos=%d_pe=%d_pel=%d_rl=%d_wrv=%d_si=%d_il=%d_r1=%s_r2=%s",
					&r->strand,
					&r->chr,
					&r->pos,
					&r->pairedEnd,
					&r->pairedEndLength,
					&r->readLength,
					&r->whichReadVariants,
					&r->startIndel,
					&r->indelLength,
					r1,
					r2)) {
			PrintError(FnName,
					a.readName,
					"Could not parse read name (1)",
					Exit,
					OutOfRange);
		}
	}
	/* Parse r1 and r2 */
	assert(r->pairedEnd == tempPairedEnd);
	assert(r->readLength == (int)strlen(r1));
	assert(r->pairedEnd == 0 || r->readLength == (int)strlen(r2));
	/* HERE */
	r->numSNPs = 0;
	r->numErrors = 0;
	r->deletionLength = 0;
	r->insertionLength = 0;
	for(i=0;i<r->readLength;i++) {
		switch(r1[i]) {
			case '0':
				/* Default */
				break;
			case '1':
				/* Insertion */
				r->insertionLength++;
				break;
			case '2':
				/* SNP */
				r->numSNPs++;
				break;
			case '3':
				/* Error */
				r->numErrors++;
				break;
			case '4':
				/* InsertionAndSNP */
				r->insertionLength++;
				r->numSNPs++;
				break;
			case '5':
				/* InsertionAndError */
				r->insertionLength++;
				r->numErrors++;
				break;
			case '6':
				/* SNPAndError */
				r->numSNPs++;
				r->numErrors++;
				break;
			case '7':
				/* InsertionSNPAndError */
				r->insertionLength++;
				r->numSNPs++;
				r->numErrors++;
				break;
			default:
				PrintError(FnName,
						"r1[i]",
						"Could not understand type",
						Exit,
						OutOfRange);
		}
		if(r->pairedEnd == 1) {
			switch(r2[i]) {
				case Default:
					break;
				case Insertion:
					r->insertionLength++;
					break;
				case SNP:
					r->numSNPs++;
					break;
				case Error:
					r->numErrors++;
					break;
				case InsertionAndSNP:
					r->insertionLength++;
					r->numSNPs++;
					break;
				case InsertionAndError:
					r->insertionLength++;
					r->numErrors++;
					break;
				case SNPAndError:
					r->numSNPs++;
					r->numErrors++;
					break;
				case InsertionSNPAndError:
					r->insertionLength++;
					r->numSNPs++;
					r->numErrors++;
					break;
				default:
					PrintError(FnName,
							"r2[i]",
							"Could not understand type",
							Exit,
							OutOfRange);
			}
		}
	}
	if(r->startIndel >= 0 && r->insertionLength == 0) {
		r->deletionLength = r->indelLength;
	}

	/* Delete align entries */
	AlignEntriesFree(&a);

	return 1;
}

void Evaluate(char *baf)
{
	char *FnName="Evaluate";
	FILE *fp;
	ReadType r;

	/* Open the file */
	if(!(fp=fopen(baf, "rb"))) {
		PrintError(FnName,
				baf,
				"Could not open file for reading",
				Exit,
				OpenFileError);
	}

	ReadTypeInitialize(&r);

	while(EOF != ReadTypeReadFromBAF(&r, fp)) {
		ReadTypeInitialize(&r);
	}

	/* Close the file */
	fclose(fp);
}
