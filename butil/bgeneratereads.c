#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "SimRead.h"
#include "bgeneratereads.h"

#define READS_ROTATE_NUM 100000
#define Name "bgeneratereads"

/* Generate synthetic reads given a number of variants and errors
 * from a reference genome. */

int main(int argc, char *argv[]) 
{
	RGBinary rg;
	char rgFileName[MAX_FILENAME_LENGTH]="\0";
	int space = 0;
	int indel = 0;
	int indelLength = 0;
	int withinInsertion = 0;
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
		withinInsertion = atoi(argv[5]);
		numSNPs = atoi(argv[6]);
		numErrors = atoi(argv[7]);
		readLength = atoi(argv[8]);
		pairedEnd = atoi(argv[9]);
		pairedEndLength = atoi(argv[10]);
		numReads = atoi(argv[11]);

		/* Check cmd line options */
		assert(space == 0 || space == 1);
		assert(indel >= 0 && indel <= 2);
		assert(indelLength > 0 || indel == 0);
		assert(withinInsertion == 0 || (indel == 2 && withinInsertion == 1));
		assert(numSNPs >= 0);
		assert(readLength > 0);
		assert(readLength < SEQUENCE_LENGTH);
		assert(pairedEnd == 0 || pairedEnd == 1);
		assert(pairedEndLength > 0 || pairedEnd == 0);
		assert(numReads > 0);

		/* Should check if we have enough bases for everything */

		/* Get reference genome */
		RGBinaryReadBinary(&rg,
				rgFileName);

		/* Generate reads */
		GenerateReads(&rg,
				space,
				indel,
				indelLength,
				withinInsertion,
				numSNPs,
				numErrors,
				readLength,
				pairedEnd,
				pairedEndLength,
				numReads);

		/* Delete reference genome */
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Deleting reference genome.\n");
		RGBinaryDelete(&rg);
		fprintf(stderr, "%s", BREAK_LINE);

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully.\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast reference genome file name (must be in nt space)>\n");
		fprintf(stderr, "\t<space 0: nt space 1: color space>\n");
		fprintf(stderr, "\t<indel 0: none 1: deletion 2: insertion>\n");
		fprintf(stderr, "\t<indel length>\n");
		fprintf(stderr, "\t<include errors within insertion 0: false 1: true>\n");
		fprintf(stderr, "\t<# of SNPs>\n");
		fprintf(stderr, "\t<# of errors>\n");
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
		int withinInsertion,
		int numSNPs,
		int numErrors,
		int readLength,
		int pairedEnd,
		int pairedEndLength,
		int numReads)
{
	char *FnName="GenerateReads";
	SimRead r;
	char outFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fp=NULL;
	int i;
	int64_t rgLength = 0;

	if(NTSpace != rg->space) {
		PrintError(FnName,
				"rg->space",
				"The reference genome must be given in nucleotide space",
				Exit,
				OutOfRange);
	}

	/* Seed random number */
	srand(time(NULL));

	/* Get the reference genome length */
	for(i=0;i<rg->numContigs;i++) {
		rgLength += rg->contigs[i].sequenceLength;
	}

	/* Create output file name */
	sprintf(outFileName, "reads.%d.%d.%d.%d.%d.%d.%d.%d.%d.%d.fa",
			space,
			indel,
			indelLength,
			withinInsertion,
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

	fprintf(stderr, "%s", BREAK_LINE);
	fprintf(stderr, "Outputting to %s.\n",
			outFileName);
	fprintf(stderr, "%s", BREAK_LINE);

	/* Initialize */
	r.readLength = readLength;
	r.pairedEnd = pairedEnd;
	r.pairedEndLength = pairedEndLength;
	r.indelLength = indelLength;
	SimReadInitialize(&r);

	/* Generate the reads */
	fprintf(stderr, "%s", BREAK_LINE);
	fprintf(stderr, "Out of %d reads, currently on:\n0", numReads);
	for(i=0;i<numReads;i++) {
		if((i+1) % READS_ROTATE_NUM==0) {
			fprintf(stderr, "\r%d",
					(i+1));
		}
		/* Get the read */
		SimReadGetRandom(rg,
				rgLength,
				&r,
				space,
				indel,
				indelLength,
				withinInsertion,
				numSNPs,
				numErrors);
		/* Output */
		r.readNum = i+1;
		SimReadPrint(&r,
				fp);

		/* Initialize read */
		SimReadDelete(&r);
	}
	fprintf(stderr, "\r%d\n%s",
			numReads,
			BREAK_LINE);

	/* Close output file */
	fclose(fp);
}
