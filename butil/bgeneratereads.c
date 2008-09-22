#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "bgeneratereads.h"

#define Name "bgeneratereads"
#define MAX_COUNT 100
#define READS_ROTATE_NUM 10000

char Colors[4] = "0123";

/* Do not change read length, paired end, or paired end length */
void ReadInitialize(Read *r)
{
	int i;
	r->readOne = NULL;
	r->readTwo = NULL;
	r->chr = 0;
	r->pos = 0;
	r->strand = 0;
	r->whichReadVariants = -1;
	r->startIndel = -1;
	for(i=0;i<SEQUENCE_LENGTH-1;i++) {
		r->readOneType[i] = Default;
		r->readTwoType[i] = Default;
	}

}

/* Do not change read length, paired end, or paired end length */
void ReadDelete(Read *r) 
{
	free(r->readOne);
	free(r->readTwo);
	ReadInitialize(r);
}

void ReadPrint(Read *r, 
		FILE *fp)
{
	int i;

	/* Read name */
	fprintf(fp, ">strand=%c_chr=%d_pos=%d_pe=%d_pel=%d_rl=%d_wrv=%d_si=%d_il=%d",
			r->strand,
			r->chr,
			r->pos,
			r->pairedEnd,
			r->pairedEndLength,
			r->readLength,
			r->whichReadVariants,
			r->startIndel,
			r->indelLength);
	fprintf(fp, "_r1=");
	for(i=0;i<r->readLength;i++) {
		fprintf(fp, "%1d", r->readOneType[i]);
	}
	if(r->pairedEnd==1) {
		fprintf(fp, "_r2=");
		for(i=0;i<r->readLength;i++) {
			fprintf(fp, "%1d", r->readTwoType[i]);
		}
	}
	fprintf(fp,"\n");

	/* Read one */
	fprintf(fp, "%s\n", r->readOne);

	/* Read two */
	if(r->pairedEnd==1) {
		fprintf(fp, "%s\n", r->readTwo);
	}
}

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
	ReadInitialize(&r);

	/* Generate the reads */
	fprintf(stderr, "%s", BREAK_LINE);
	fprintf(stderr, "Out of %d reads, currently on:\n0", numReads);
	for(i=0;i<numReads;i++) {
		if((i+1) % READS_ROTATE_NUM==0) {
			fprintf(stderr, "\r%d",
					(i+1));
		}
		/* Get the read */
		GetRandomRead(rg, 
				rgLength,
				&r);

		/* Modify the read based on indel, SNPs, errors and color space */
		if(1==ModifyRead(rg,
					&r,
					space,
					indel,
					indelLength,
					withinInsertion,
					numSNPs,
					numErrors)) {
			/* Output */
			ReadPrint(&r,
					fp);
		}
		else {
			/* Do not output, try again */
			i--;
		}


		/* Initialize read */
		ReadDelete(&r);
	}
	fprintf(stderr, "\r%d\n%s",
			numReads,
			BREAK_LINE);

	/* Close output file */
	fclose(fp);
}

void GetRandomRead(RGBinary *rg,
		int64_t rgLength,
		Read *r)
{
	char *FnName="GetRandomRead";
	int count = 0;
	int i;
	int hasNs=0;
	int readOneSuccess=0;
	int readTwoSuccess=0;

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
		if(count > 1) {
			ReadDelete(r);
		}

		/* Get the random chromosome and position */
		GetRandomChrPos(rg,
				rgLength,
				&r->chr,
				&r->pos,
				&r->strand);

		if(r->pairedEnd == 1) {
			/* Get the sequence for the first read */
			readOneSuccess = RGBinaryGetSequence(rg,
					r->chr,
					r->pos + r->readLength + r->pairedEndLength,
					r->strand,
					&r->readOne,
					r->readLength);
			/* Get the sequence for the second read */
			readTwoSuccess = RGBinaryGetSequence(rg,
					r->chr,
					r->pos,
					r->strand,
					&r->readTwo,
					r->readLength);
		}
		else {
			/* Get the sequence for the first read */
			readOneSuccess = RGBinaryGetSequence(rg,
					r->chr,
					r->pos,
					r->strand,
					&r->readOne,
					r->readLength);
		}

		/* Make sure there are no Ns */
		hasNs = 0;
		if(1==readOneSuccess) {
			for(i=0;0==hasNs && i<r->readLength;i++) {
				if(RGBinaryIsBaseN(r->readOne[i]) == 1) {
					hasNs = 1;
				}
			}
		}
		if(r->pairedEnd == 1 && 1==readTwoSuccess) {
			for(i=0;0==hasNs && i<r->readLength;i++) {
				if(RGBinaryIsBaseN(r->readTwo[i]) == 1) {
					hasNs = 1;
				}
			}
		}

	} while(
			(readOneSuccess == 0) || /* Read one was successfully read */
			(r->pairedEnd == 1 && readTwoSuccess == 0) || /* Read two was successfully read */
			(hasNs == 1) /* Either read end has an "N" */
		   );
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
	while(low < high) {
		/* HERE */
		/*
		   fprintf(stderr, "low=%lld\tmid=%lld\thigh=%lld\n",
		   (long long int)low,
		   (long long int)mid,
		   (long long int)high);
		   */
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

/* TODO */
int ModifyRead(RGBinary *rg,
		Read *r,
		int space,
		int indel,
		int indelLength,
		int withinInsertion,
		int numSNPs,
		int numErrors)
{
	/* Apply them in this order:
	 * 1. insert an indel based on indel length 
	 * 2. insert a SNP based on include within insertion
	 * 3. convert to color space (if necessary)
	 * 4. insert errors (color errors if in color space, otherwise nt errors )
	 */
	char *FnName="ModifyRead";
	int tempReadLength=r->readLength;

	/* Which read should the variants be contained within */
	r->whichReadVariants= (r->pairedEnd == 0)?0:(rand()%2); 

	/* 1. Insert an indel based on the indel length */
	switch(indel) {
		case 0:
			/* Do nothing */
			break;
		case 1:
		case 2:
			if(0==InsertIndel(rg, r, indel, indelLength)) {
				/* Could not add an indel */
				return 0;
			}
			break;
		default:
			PrintError(FnName,
					"indel",
					"indel out of range",
					Exit,
					OutOfRange);
	}

	/* 2. SNPs */
	InsertMismatches(r,
			numSNPs,
			SNP,
			0);

	if(1 == space) {
		/* 3. Convert to color space if necessary */
		tempReadLength = r->readLength;
		ConvertReadToColorSpace(&r->readOne,
				&tempReadLength);
		if(1==r->pairedEnd) {
			tempReadLength = r->readLength;
			ConvertReadToColorSpace(&r->readTwo,
					&tempReadLength);
		}

		/* 4. Insert errors if necessary */
		InsertColorErrors(r,
				numErrors,
				withinInsertion);
	}
	else {
		/* 4. Insert NT errors */
		InsertMismatches(r, 
				numErrors,
				Error,
				withinInsertion);
	}

	return 1;
}

/* TODO */
int InsertIndel(RGBinary *rg,
		Read *r,
		int indel,
		int indelLength)
{
	char *FnName="InsertIndel";
	int i;
	int start; /* starting position within the read */
	int success=1;

	/* Pick a starting position within the read to insert */
	start = rand() % (r->readLength - indelLength + 1);

	if(indel == 1) {
		/* Deletion */
		/* Remove bases */
		if(r->whichReadVariants == 0) {
			/* Free read */
			free(r->readOne);
			r->readOne = NULL;
			/* Get new read */
			success = RGBinaryGetSequence(rg,
					r->chr,
					r->pos,
					r->strand,
					&r->readOne,
					r->readLength + indelLength);
			if(success == 1) {
				/* Shift over bases */
				for(i=start;i<r->readLength-indelLength;i++) {
					r->readOne[i] = r->readOne[i+indelLength];
				}
				/* Reallocate memory */
				r->readOne = realloc(r->readOne, sizeof(char)*r->readLength);
				if(NULL==r->readOne) {
					PrintError(FnName,
							"r->readOne",
							"Could not reallocate memory",
							Exit,
							ReallocMemory);
				}
				/* Adjust position if reverse strand */
				r->pos = (r->strand==REVERSE)?(r->pos + indelLength):r->pos;
			}
		}
		else {
			/* Free read */
			free(r->readTwo);
			r->readTwo = NULL;
			/* Get new read */
			success = RGBinaryGetSequence(rg,
					r->chr,
					r->pos,
					r->strand,
					&r->readTwo,
					r->readLength + indelLength);
			if(success == 1) {
				/* Shift over bases */
				for(i=start;i<r->readLength-indelLength;i++) {
					r->readTwo[i] = r->readTwo[i+indelLength];
				}
				/* Reallocate memory */
				r->readTwo = realloc(r->readTwo, sizeof(char)*r->readLength);
				if(NULL==r->readTwo) {
					PrintError(FnName,
							"r->read",
							"Could not reallocate memory",
							Exit,
							ReallocMemory);
				}
				/* Adjust position if reverse strand */
				r->pos = (r->strand==REVERSE)?(r->pos + indelLength):r->pos;
			}
		}
	}
	else if(indel == 2) {
		/* Insertion */
		if(r->whichReadVariants == 0) {
			/* shift over all above */
			for(i = r->readLength;i >= start + indelLength;i--) {
				r->readOne[i] = r->readOne[i-indelLength];
			}
			/* insert random bases */
			for(i=start;i<start+indelLength;i++) {
				r->readOne[i] = DNA[rand()%4];
				r->readOneType[i] = Insertion;
			}
		}
		else {
			/* shift over all above */
			for(i = r->readLength;i >= start + indelLength;i--) {
				r->readTwo[i] = r->readTwo[i-indelLength];
			}
			/* insert random bases */
			for(i=start;i<start+indelLength;i++) {
				r->readTwo[i] = DNA[rand()%4];
				r->readTwoType[i] = Insertion;
			}
		}
	}
	else {
		PrintError(FnName,
				"indel",
				"indel out of range",
				Exit,
				OutOfRange);
	}

	/* Update the start of the indel */
	r->startIndel = start;

	return success;
}

void InsertMismatches(Read *r,
		int numMismatches,
		int type,
		int withinInsertion)
{
	if(r->whichReadVariants == 0) {
		InsertMismatchesHelper(r->readOne,
				r->readLength,
				r->readOneType,
				numMismatches,
				type,
				withinInsertion);
	}
	else {
		InsertMismatchesHelper(r->readTwo,
				r->readLength,
				r->readTwoType,
				numMismatches,
				type,
				withinInsertion);
	}
}

void InsertMismatchesHelper(char *read,
		int readLength,
		int *readType,
		int numMismatches,
		int type,
		int withinInsertion)
{
	char *FnName = "InsertMismatches";
	int numMismatchesLeft = numMismatches;
	int index;
	char original;
	int toAdd;

	assert(type == SNP || type == Error);

	while(numMismatchesLeft > 0) {
		/* Pick a base to modify */
		index = rand()%(readLength);
		toAdd = 0;

		switch(readType[index]) {
			case Default:
				readType[index] = (type==SNP)?SNP:Error;
				toAdd = 1;
				break;
			case Insertion:
				if(withinInsertion == 1) {
					readType[index] = (type==SNP)?InsertionAndSNP:InsertionAndError;
					toAdd = 1;
				}
				break;
			case SNP:
				if(type == Error) {
					readType[index] = SNPAndError;
					toAdd = 1;
				}
				break;
			case Error:
				/* Nothing, since we assume that SNPs were applied before errors */
				assert(type != SNP);
				break;
			case InsertionAndSNP:
				if(withinInsertion == 1 && type == Error) {
					readType[index] = InsertionSNPAndError;
					toAdd = 1;
				}
				break;
			case SNPAndError:
				/* Nothing */
				break;
			case InsertionAndError:
				/* Nothing, since we assume that SNPs were applied before errors */
				assert(type != SNP);
				break;
			case InsertionSNPAndError:
				/* Nothing */
				break;
			default:
				PrintError(FnName,
						"readType[index]",
						"Could not understand type",
						Exit,
						OutOfRange);
				break;
		}
		if(1==toAdd) {
			/* Modify base to a new base */
			for(original = read[index];
					original == read[index];
					read[index] = DNA[rand()%4]) {
			}
			numMismatchesLeft--;
		}
	}
}

void InsertColorErrors(Read *r,
		int numErrors,
		int withinInsertion)
{
	/*
	   char *FnName = "InsertColorErrors";
	   */
	int numErrorsLeft = numErrors;
	int which, index;
	char original;
	int toAdd;

	while(numErrorsLeft > 0) {
		/* Pick a read */
		which = (r->pairedEnd == 0)?0:(rand()%2);
		/* Pick a color to modify */
		index = (rand()%(r->readLength) )+ 1;
		toAdd = 0;

		/* Assumes the type can only be Default, Insertion, SNP, or InsertionAndSNP */
		if(which == 0) {
			if(withinInsertion == 1 && r->readOneType[index] == Insertion) {
				r->readOneType[index] = InsertionAndError;
				toAdd = 1;
			}
			if(withinInsertion == 1 && r->readOneType[index] == InsertionAndSNP) {
				r->readOneType[index] = InsertionSNPAndError;
				toAdd = 1;
			}
			else if(r->readOneType[index] == SNP) {
				r->readOneType[index] = SNPAndError;
				toAdd = 1;
			}
			else if(r->readOneType[index] == Default) {
				r->readOneType[index] = Error;
				toAdd = 1;
			}
			if(1==toAdd) {
				/* Modify color to a new color */
				for(original = r->readTwo[index];
						original == r->readTwo[index];
						r->readTwo[index] = Colors[rand()%4]) {
				}
				numErrorsLeft--;
			}
		}
		else {
			if(withinInsertion == 1 && r->readTwoType[index] == Insertion) {
				r->readTwoType[index] = InsertionAndError;
				toAdd = 1;
			}
			if(withinInsertion == 1 && r->readTwoType[index] == InsertionAndSNP) {
				r->readTwoType[index] = InsertionSNPAndError;
				toAdd = 1;
			}
			else if(r->readTwoType[index] == SNP) {
				r->readTwoType[index] = SNPAndError;
				toAdd = 1;
			}
			else if(r->readTwoType[index] == Default) {
				r->readTwoType[index] = Error;
				toAdd = 1;
			}
			if(1==toAdd) {
				/* Modify color to a new color */
				for(original = r->readTwo[index];
						original == r->readTwo[index];
						r->readTwo[index] = Colors[rand()%4]) {
				}
				numErrorsLeft--;
			}
		}
	}
}
