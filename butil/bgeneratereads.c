#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "../blib/BError.h"

#include "bgeneratereads.h"

#define Name "bgeneratereads"
#define MAX_COUNT 100

int main(int argc, char *argv[]) 
{
	RGBinary rg;
	char rgFileName[MAX_FILENAME_LENGTH]="\0";
	int space = 0;
	int indel = 0;
	int indelLength = 0;
	int withinIndel = 0;
	int numSNPs = 0;
	int numErrors = 0;
	int readLength = 0;
	int pairedEnd = 0;
	int pairedEndLength = 0;
	int numReads = 0;

	if(argc == 12) {

		/* Get cmd line options */
		strcpy(rgFileName, argv[1]);
		space = atoi(argv[2]);
		indel = atoi(argv[3]);
		indelLength = atoi(argv[4]);
		withinIndel = atoi(argv[5]);
		numSNPs = atoi(argv[6]);
		numErrors = atoi(argv[7]);
		readLength = atoi(argv[8]);
		pairedEnd = atoi(argv[9]);
		pairedEndLength = atoi(argv[10]);
		numReads = atoi(argv[11]);

		/* Check cmd line options */
		assert(space == 0 || space == 1);
		assert(indel == 0 || indel == 1);
		assert(indelLength > 0 || indel == 0);
		assert(withinIndel == 0 || withinIndel == 1);
		assert(numSNPs >= 0);
		if(space == 0 && numErrors > 0) {
			PrintError(Name,
					"numErrors",
					"Cannot use # of errors with nt space",
					Exit,
					OutOfRange);
		}
		assert(readLength > 0);
		assert(pairedEnd == 0 || pairedEnd == 1);
		assert(pairedEndLength > 0 || pairedEnd == 0);
		assert(numReads > 0);

		/* Get reference genome */
		RGBinaryReadBinary(&rg,
				rgFileName);

		/* Delete reference genome */
		RGBinaryDelete(&rg);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast reference genome file name (must be in nt space)>\n");
		fprintf(stderr, "\t<space 0: nt space 1: color space>\n");
		fprintf(stderr, "\t<indel 0: none 1: deletion 2: insertion>\n");
		fprintf(stderr, "\t<indel length>\n");
		fprintf(stderr, "\t<SNPs and errors within indel 0: false 1: true>\n");
		fprintf(stderr, "\t<# of SNPs>\n");
		fprintf(stderr, "\t<# of errors (color space only)>\n");
		fprintf(stderr, "\t<read length>\n");
		fprintf(stderr, "\t<paired end 0: false 1: true>\n");
		fprintf(stderr, "\t<paired end length>\n");
		fprintf(stderr, "\t<number of reads>\n");
	}
	return 0;
}

/* TODO */
void GenerateReads(RGBinary *rg,
		int space,
		int indel,
		int indelLength,
		int withinIndel,
		int numSNPs,
		int numErrors,
		int readLength,
		int pairedEnd,
		int pairedEndLength,
		int numReads)
{
	char *FnName="GenerateReads";
	Read r;
	char outFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fp=NULL;
	int i;
	int64_t rgLength = 0;

	/* Seed random number */
	srand(time(NULL));

	/* Get the reference genome length */
	for(i=0;i<rg->numChrs;i++) {
		rgLength += rg->chromosomes[i].endPos - rg->chromosomes[i].startPos + 1;
	}

	/* Create output file name */
	sprintf(outFileName, "reads.%d.%d.%d.%d.%d.%d.%d.%d.%d.%d.fa",
			space,
			indel,
			indelLength,
			withinIndel,
			numSNPs,
			numErrors,
			readLength,
			pairedEnd,
			pairedEndLength,
			numReads);

	/* Open output file */
	if(!(fp=fopen(outFileName, "wb"))) {
		PrintError(FnName,
				outFileName,
				"Could not open output file for writing.\n",
				Exit,
				OpenFileError);
	}

	/* Initialize */
	r.readLength = readLength;
	r.pairedEnd = pairedEnd;
	r.pairedEndLength = pairedEndLength;

	/* Generate the reads */
	for(i=0;i<numReads;i++) {
		/* Get the read */
		GetRandomRead(rg, 
				rgLength,
				&r);
		/* Initialize read */
		ReadInitialize(&r);
	}

	/* Close output file */
	fclose(fp);
}

void GetRandomRead(RGBinary *rg,
		int64_t rgLength,
		Read *r)
{
	char *FnName="GetRandomRead";
	int posOne, posTwo;
	int readLengthOne, readLengthTwo;
	int count = 0;

	do {
		/* Avoid infinite loop */
		count++;
		if(count > MAX_COUNT) {
			PrintError(FnName,
					"count",
					"Could not get a random read",
					Exit,
					OutOfRange);
		}

		/* Initialize read */
		ReadInitialize(r);

		/* Get the random chromosome and position */
		GetRandomChrPos(rg,
				rgLength,
				&r->chr,
				&r->pos,
				&r->strand);

		/* Get the sequence for the first read */
		RGBinaryGetSequence(rg,
				r->chr,
				r->pos,
				r->strand,
				0,
				&r->readOne,
				r->readLength,
				&readLengthOne,
				&posOne);

		/* Get the sequecne for the second read */
		if(r->pairedEnd == 1) {
			RGBinaryGetSequence(rg,
					r->chr,
					r->pos + r->readLength + r->pairedEndLength,
					r->strand,
					0,
					&r->readTwo,
					r->readLength,
					&readLengthTwo,
					&posTwo);
		}

	} while(posOne != r->pos || /* First read has the same starting position */
			readLengthOne != r->readLength || /* First read has the same length */ 
			(r->pairedEnd == 1 && posTwo != r->pos + r->readLength + r->pairedEndLength) || /* Second read has the same starting position */
			(r->pairedEnd == 1 && readLengthTwo != r->readLength)); /* Second read has the same length */
}

/* Do not change read length, paired end, or paired end length */
void ReadInitialize(Read *r)
{
	free(r->readOne);
	r->readOne = NULL;
	free(r->readTwo);
	r->readTwo = NULL;
	r->chr = 0;
	r->pos = 0;
}

/* Get the random chromosome and position */
void GetRandomChrPos(RGBinary *rg,
		int64_t rgLength,
		int *chr,
		int *pos,
		char *strand)
{
	char *FnName = "GetRandomChrPos";
	int i;
	int64_t curI;
	int64_t low, mid, high;
	int value;
	int count = 0;

	low = 1;
	high = rgLength;

	/* Flip a coin for strand */
	(*strand) = ((rand()%2)==0)?FORWARD:REVERSE;

	/* Use coin flips to find position */
	mid = (low + high)/2;
	while(low <= high) {
		mid = (low + high)/2;
		value = rand() % 2;
		if(value == 0) {
			/* lower */
			high = mid;
		}
		else {
			assert(value == 1);
			/* upper */
			low = mid;
		}
		/* To avoid an infinite loop */
		count++;
		if(count > MAX_COUNT) {
			PrintError(FnName,
					"count",
					"Could not get random chromosome and position",
					Exit,
					OutOfRange);
		}
	}

	/* Identify where it occurs */
	curI=0;
	for(i=0;i<rg->numChrs;i++) {
		curI += rg->chromosomes[i].endPos - rg->chromosomes[i].startPos + 1;
		if(mid <= curI) {
			(*chr) = rg->startChr + i;
			(*pos) = rg->chromosomes[i].startPos + (curI - mid);
			return;
		}
	}
	
	PrintError(FnName,
			"mid",
			"Mid was out of range",
			Exit,
			OutOfRange);
}
