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
#define COUNT_ROTATE_NUM 100000

int main(int argc, char *argv[]) 
{
	char baf[MAX_FILENAME_LENGTH]="\0";
	char outputID[MAX_FILENAME_LENGTH]="\0";

	if(argc == 3) {

		/* Get cmd line options */
		strcpy(baf, argv[1]);
		strcpy(outputID, argv[2]);

		/* Check cmd line options */

		/* Run program */
		Evaluate(baf, outputID);

		/* Terminate */
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully.\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast aligned file>\n");
		fprintf(stderr, "\t<output id>\n");
	}
	return 0;
}

void ReadTypeInitialize(ReadType *r)
{
	r->strand=0;
	r->chr=0;
	r->pos=0;
	r->space=0;
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
	r->aChr=0;
	r->aPos=0;
	r->aStrand=0;
}

void ReadTypeCopy(ReadType *dest,
		ReadType *src)
{
	/* Only copy meta data */ 
	dest->strand=src->strand;
	dest->chr=src->chr;
	dest->pos=src->pos;
	dest->space=src->space;
	dest->pairedEnd=src->pairedEnd;
	dest->pairedEndLength=src->pairedEndLength;
	dest->readLength=src->readLength;
	dest->whichReadVariants=src->whichReadVariants;
	dest->startIndel=src->startIndel;
	dest->indelLength=src->indelLength;
	dest->numSNPs=src->numSNPs;
	dest->numErrors=src->numErrors;
	dest->deletionLength=src->deletionLength;
	dest->insertionLength=src->insertionLength;
}

int ReadTypeCompare(ReadType *a,
		ReadType *b)
{
	/* Only compare meta data */ 
	/* Nice use of if, else if, and else statements */
	if(a->space != b->space) {
		return (a->space < b->space)?-1:1;
	}
	else if(a->pairedEnd != b->pairedEnd) {
		return (a->pairedEnd < b->pairedEnd)?-1:1;
	}
	else if(a->pairedEndLength != b->pairedEndLength) {
		return (a->pairedEndLength < b->pairedEndLength)?-1:1;
	}
	else if(a->readLength != b->readLength) {
		return (a->readLength < b->readLength)?-1:1;
	}
	else if(a->whichReadVariants != b->whichReadVariants) {
		return (a->whichReadVariants < b->whichReadVariants)?-1:1;
	}
	else if(a->startIndel != b->startIndel) {
		return (a->startIndel < b->startIndel)?-1:1;
	}
	else if(a->indelLength != b->indelLength) {
		return (a->indelLength < b->indelLength)?-1:1;
	}
	else if(a->numSNPs != b->numSNPs) {
		return (a->numSNPs < b->numSNPs)?-1:1;
	}
	else if(a->numErrors != b->numErrors) {
		return (a->numErrors < b->numErrors)?-1:1;
	}
	else if(a->deletionLength != b->deletionLength) {
		return (a->deletionLength < b->deletionLength)?-1:1;
	}
	else if(a->insertionLength != b->insertionLength) {
		return (a->insertionLength < b->insertionLength)?-1:1;
	}
	else {
		return 0;
	}
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

	r->aChr = a.entriesOne->chromosome;
	r->aPos = a.entriesOne->position;
	r->aStrand = a.entriesOne->strand;
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

void StatsInitialize(Stats *s, ReadType *r)
{
	s->numReads=0;
	s->numCorrectlyAligned[0]=0;
	s->numCorrectlyAligned[1]=0;
	s->numCorrectlyAligned[2]=0;
	s->numCorrectlyAligned[3]=0;
	s->numCorrectlyAligned[4]=0;
	s->space = r->space;
	s->pairedEnd = r->pairedEnd;
	s->pairedEndLength = r->pairedEndLength;
	s->readLength = r->readLength;
	s->indelLength = r->indelLength;
	s->numSNPs = r->numSNPs;
	s->numErrors = r->numErrors;
	s->deletionLength = r->deletionLength;
	s->insertionLength = r->insertionLength;
}

void StatsPrintHeader(FILE *fp)
{
	fprintf(fp, "# COL | Description\n");
	fprintf(fp, "# 0    | number of reads\n");
	fprintf(fp, "# 1    | number of correctly aligned within 0 bases\n");
	fprintf(fp, "# 2    | number of correctly aligned within 10 bases\n");
	fprintf(fp, "# 3    | number of correctly aligned within 100 bases\n");
	fprintf(fp, "# 4    | number of correctly aligned within 1000 bases\n");
	fprintf(fp, "# 5    | number of correctly aligned within 10000 bases\n");
	fprintf(fp, "# 6    | space\n");
	fprintf(fp, "# 7    | paired end\n");
	fprintf(fp, "# 8    | paired end length\n");
	fprintf(fp, "# 9    | read length\n");
	fprintf(fp, "# 10   | indel length\n");
	fprintf(fp, "# 11   | number of snps\n");
	fprintf(fp, "# 12   | number of errors\n");
	fprintf(fp, "# 13   | deletion length\n");
	fprintf(fp, "# 14   | insertion length\n");
}

void StatsPrint(Stats *s, FILE *fp)
{
	fprintf(fp, "%d %d %d %d %d %d ",
			s->numReads,
			s->numCorrectlyAligned[0],
			s->numCorrectlyAligned[1],
			s->numCorrectlyAligned[2],
			s->numCorrectlyAligned[3],
			s->numCorrectlyAligned[4]);
	fprintf(fp, "%d %d %d %d %d %d %d %d %d\n",
			s->space,
			s->pairedEnd,
			s->pairedEndLength,
			s->readLength,
			s->indelLength,
			s->numSNPs,
			s->numErrors,
			s->deletionLength,
			s->insertionLength);
}

void StatsAdd(Stats *s, ReadType *r)
{
	int diff;
	if(r->strand == r->aStrand &&
			r->chr == r->aChr) {
		diff = (r->pos > r->aPos)?(r->pos - r->aPos):(r->aPos - r->pos);

		/* Update */
		if(diff <= 10000) {
			s->numCorrectlyAligned[4]++;
			if(diff <= 1000) {
				s->numCorrectlyAligned[3]++;
				if(diff <= 100) {
					s->numCorrectlyAligned[2]++;
					if(diff <= 10) {
						s->numCorrectlyAligned[1]++;
						if(diff <= 0) {
							s->numCorrectlyAligned[0]++;
						}
					}
				}
			}
		}
	}
	s->numReads++;
}

void Evaluate(char *baf,
		char *outputID)
{
	char *FnName="Evaluate";
	FILE *fpIn;
	FILE *fpOut;
	ReadType r, prev;
	Stats s;
	int32_t count;
	char outputFileName[MAX_FILENAME_LENGTH]="\0";

	/* Open the baf file */
	if(!(fpIn=fopen(baf, "rb"))) {
		PrintError(FnName,
				baf,
				"Could not open file for reading",
				Exit,
				OpenFileError);
	}
	/* Create output file name */
	sprintf(outputFileName, "%s.evalsim.%s.txt",
			PROGRAM_NAME,
			outputID);
	/* Open the output file */
	if(!(fpOut=fopen(outputFileName, "wb"))) {
		PrintError(FnName,
				outputFileName,
				"Could not open file for writing",
				Exit,
				WriteFileError);
	}

	ReadTypeInitialize(&prev);
	ReadTypeInitialize(&r);
	StatsInitialize(&s, &r);
	StatsPrintHeader(fpOut);

	count = 0;
	fprintf(stderr, "Currently on:\n%d", 0);
	while(EOF != ReadTypeReadFromBAF(&r, fpIn)) {
		count++;
		if(count % COUNT_ROTATE_NUM == 0) {
			fprintf(stderr, "\r%d", count);
		}

		/* Process the read */
		if(ReadTypeCompare(&r, &prev)!=0) {
			/* Print statistics */
			if(count > 1) { /* Don't print on the first go round */
				StatsPrint(&s, fpOut);
			}
			/* Copy over from r to prev*/
			ReadTypeCopy(&prev, &r);
			assert(ReadTypeCompare(&r, &prev)==0);
			/* Initialize statistics */
			StatsInitialize(&s, &r);
		}
		/*Add to statistics */
		StatsAdd(&s, &r);

		/* Reinitialize */
		ReadTypeInitialize(&r);
	}
	fprintf(stderr, "\r%d\n", count);

	/* Close the files */
	fclose(fpIn);
	fclose(fpOut);
}
