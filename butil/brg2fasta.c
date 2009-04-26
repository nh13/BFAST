#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/AlignedRead.h"
#include "../blib/AlignedReadConvert.h"
#include "../blib/BError.h"
#include "../blib/BLib.h"

#define Name "bafconvert"
#define BRG2FASTA_FASTA_LINE_LENGTH 60
#define BAFCONVERT_ROTATE_NUM 100000

/* Prints the reference genome in FASTA format.
 * */

int main(int argc, char *argv[])
{
	char rgFileName[MAX_FILENAME_LENGTH]="\0";
	RGBinary rg;

	if(2 == argc) {
		strcpy(rgFileName, argv[1]);
		RGBinaryReadBinary(&rg,
				rgFileName);
		/* Unpack */
		RGBinaryUnPack(&rg);

		int32_t i, j, ctr;
		for(i=0;i<rg.numContigs;i++) {
			fprintf(stdout, ">%s\n",
					rg.contigs[i].contigName);
			for(j=ctr=0;j<rg.contigs[i].sequenceLength;j++) {
				putchar(rg.contigs[i].sequence[j]);
				ctr++;
				if(BRG2FASTA_FASTA_LINE_LENGTH < ctr) {
					putchar('\n');
					ctr=0;
				}
			}
			if(0 < ctr) {
				putchar('\n');
			}
		}

		RGBinaryDelete(&rg);

		fprintf(stderr, "Terminating successfully!\n");
	}
	else {
		fprintf(stderr, "Usage:%s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast reference genome file>\n");
	}
	return 0;
}
