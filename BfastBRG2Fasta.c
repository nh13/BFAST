#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <config.h>
#include "BLibDefinitions.h"
#include "AlignedRead.h"
#include "AlignedReadConvert.h"
#include "BError.h"
#include "BLib.h"

#define Name "bfast brg2fasta"
#define BRG2FASTA_FASTA_LINE_LENGTH 60
#define BAFCONVERT_ROTATE_NUM 100000

/* Prints the reference genome in FASTA format.
 * */

int BfastBRG2Fasta(int argc, char *argv[])
{
	char rgFileName[MAX_FILENAME_LENGTH]="\0";
	RGBinary rg;

	if(2 == argc) {
		strcpy(rgFileName, argv[1]);
		RGBinaryReadBinary(&rg,
				NTSpace, // TODO
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
		fprintf(stderr, "%s %s\n", "bfast", PACKAGE_VERSION);
		fprintf(stderr, "\nUsage:%s <bfast reference genome file>\n", Name);
		fprintf(stderr, "\nsend bugs to %s\n",
				PACKAGE_BUGREPORT);
	}
	return 0;
}
