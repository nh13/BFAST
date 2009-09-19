#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <zlib.h>

#include "../BLibDefinitions.h"
#include "../BError.h"
#include "../RGIndex.h"
#include "../RGBinary.h"
#include "bexonify.h"

#define Name "bexonify"

#define FILTER_ROTATE_NUM 100000

/* Modifies an index to only include locations specified
 * by the user.  The main purpose is to allow for 
 * alignment to exons, multiple subregions or the like.
 * */

int main(int argc, char *argv[]) 
{
	/* I have to admit, this is kind of a hack, we could 
	 * just modify the data structure of the index.  Oh well.
	 * */
	if(argc == 4) {

		gzFile fp;
		char rgFileName[MAX_FILENAME_LENGTH]="\0";
		char indexFileName[MAX_FILENAME_LENGTH]="\0";
		char exonsFileName[MAX_FILENAME_LENGTH]="\0";
		char outputFileName[MAX_FILENAME_LENGTH]="\0";

		RGBinary rg;
		RGIndex index;
		Exon *exons=NULL;
		int numExons = 0;
		int space = NTSpace; // TODO

		strcpy(rgFileName, argv[1]);
		strcpy(indexFileName, argv[2]);
		strcpy(exonsFileName, argv[3]);

		/* Read in the rg binary file */
		RGBinaryReadBinary(&rg, space, rgFileName);

		/* Read the index */
		RGIndexRead(&index, indexFileName);

		/* Read in the exons */
		numExons = ReadExons(exonsFileName, &exons);

		/* Filter based on the exons */
		FilterIndexBasedOnExons(&index, &exons, numExons);

		/* Free exons */
		free(exons);
		exons=NULL;
		numExons = 0;

		/* We need to update the hash */
		/* Free hash */
		free(index.starts);
		index.starts=NULL;

		/* Fix the hash by recreating it */
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Regenerating the hash.\n");
		RGIndexCreateHash(&index, &rg);

		/* Create new file name */
		sprintf(outputFileName, "%s.index.file.%s.%s",
				PROGRAM_NAME,
				Name,
				BFAST_INDEX_FILE_EXTENSION);

		/* Print the new index */ 
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Outputting to %s.\n", outputFileName);
		if(!(fp=gzopen(outputFileName, "wb"))) {
			PrintError(Name, outputFileName, "Could not open file for writing", Exit, OpenFileError);
		}
		RGIndexPrint(fp, &index);
		gzclose(fp);

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Cleaning up.\n");
		/* Delete the index */
		RGIndexDelete(&index);
		/* Delete the rg */
		RGBinaryDelete(&rg);

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully.\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast reference genome file>\n");
		fprintf(stderr, "\t<bfast index file>\n");
		fprintf(stderr, "\t<exon list file>\n");
	}

	return 0;
}

int ReadExons(char *exonsFileName,
		Exon **exons)
{
	char *FnName = "ReadExons";
	FILE *fp;
	uint32_t contig, prevContig;
	uint32_t start, end, prevStart, prevEnd;
	int numExons = 0;

	/* Open the file */
	fprintf(stderr, "Reading in exons from %s.\n",
			exonsFileName);
	if(!(fp=fopen(exonsFileName, "rb"))) {
		PrintError(FnName, exonsFileName, "Could not open file for reading", Exit, OpenFileError);
	}

	/* Read in the exons */ 
	prevContig = 0;
	prevStart = prevEnd = 0;
	while(EOF!=fscanf(fp, "%u %u %u", &contig, &start, &end)) {

		/* Check that the exons in increasing order */
		if(contig < prevContig || 
				(contig == prevContig && start < prevStart) ||
				(contig == prevContig && start == prevStart && end <= prevEnd)) {
			PrintError(FnName, NULL, "Entries must be in increasing order with the keys=(contig, start, end)", Exit, OutOfRange);
		}

		numExons++;
		(*exons) = realloc((*exons), sizeof(Exon)*numExons);
		if(NULL == (*exons)) {
			PrintError(FnName, "(*exons)", "Could not allocate memory", Exit, MallocMemory);
		}
		assert(start <= end);
		(*exons)[numExons-1].contig = contig;
		(*exons)[numExons-1].start = start;
		(*exons)[numExons-1].end = end;
	}

	/* Close the file */
	fclose(fp);

	fprintf(stderr, "Read in %d exons.\n",
			numExons);

	return numExons;
}

void FilterIndexBasedOnExons(RGIndex *index, Exon **exons, int numExons) 
{
	char *FnName = "FilterIndexBasedOnExons";
	int64_t i, j;
	int64_t low, mid, high, found;
	int64_t indexLength=index->length;

	fprintf(stderr, "Out of %lld, currently on:\n0",
			(long long int)index->length);

	/* Go through each entry in the index */
	for(i=index->length-1;i>=0;i--) {
		if( (indexLength-i)%FILTER_ROTATE_NUM == 0) {
			fprintf(stderr, "\r%lld",
					(long long int)(indexLength-i));
		}
		/* Check if it falls within range */
		/* Binary search */
		low = 0;
		high = numExons-1;
		found = 0;
		while(low <= high && found == 0) {
			mid = (low + high)/2;
			if(index->contigType == Contig_8) {
				if(index->contigs_8[i] < (*exons)[mid].contig ||
						(index->contigs_8[i] == (*exons)[mid].contig && index->positions[i] < (*exons)[mid].start)) {
					high = mid - 1;
				}
				else if(index->contigs_8[i] > (*exons)[mid].contig ||
						(index->contigs_8[i] == (*exons)[mid].contig && index->positions[i] > (*exons)[mid].end)) {
					low = mid + 1;
				}
				else {
					found = 1;
				}
			}
			else {
				if(index->contigs_32[i] < (*exons)[mid].contig ||
						(index->contigs_32[i] == (*exons)[mid].contig && index->positions[i] < (*exons)[mid].start)) {
					high = mid - 1;
				}
				else if(index->contigs_32[i] > (*exons)[mid].contig ||
						(index->contigs_32[i] == (*exons)[mid].contig && index->positions[i] > (*exons)[mid].end)) {
					low = mid + 1;
				}
				else {
					found = 1;
				}
			}
		}
		/* If not found, remove and shift */
		if(found == 0) {
			for(j=i+1;j<index->length;j++) {
				if(index->contigType == Contig_8) {
					index->contigs_8[j-1] = index->contigs_8[j];
				}
				else {
					index->contigs_32[j-1] = index->contigs_32[j];
				}
				index->positions[j-1] = index->positions[j];
			}
			/* Decrement index length, reallocate later */
			index->length--;
		}
	}
			
			fprintf(stderr, "\r%lld\n",
					(long long int)(indexLength-i));
	fprintf(stderr, "%s", BREAK_LINE);
	fprintf(stderr, "Reallocating index.\n");

	/* Reallocate */
	if(index->contigType == Contig_8) {
		index->contigs_8 = realloc(index->contigs_8, index->length*sizeof(uint8_t));
		if(NULL==index->contigs_8) {
			PrintError(FnName, "index->contigs_8", "Could not reallocate memory", Exit, ReallocMemory);
		}
	}
	else {
		index->contigs_32 = realloc(index->contigs_32, index->length*sizeof(uint32_t));
		if(NULL==index->contigs_32) {
			PrintError(FnName, "index->contigs_32", "Could not reallocate memory", Exit, ReallocMemory);
		}
	}
	index->positions = realloc(index->positions, index->length*sizeof(uint32_t));
	if(NULL==index->positions) {
		PrintError(FnName, "index->positions", "Could not reallocate memory", Exit, ReallocMemory);
	}

	/* Update index range */
	index->startContig = (*exons)[0].contig;
	index->startPos = (*exons)[0].start;
	index->endContig = (*exons)[numExons-1].contig;
	index->endPos = (*exons)[numExons-1].end;
	
	fprintf(stderr, "Reallocating index complete.\n");
	fprintf(stderr, "%s", BREAK_LINE);

}

