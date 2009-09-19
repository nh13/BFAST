#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>

#include "../bfast/BLibDefinitions.h"
#include "../bfast/BError.h"
#include "../bfast/BLib.h"
#include "../bfast/RGBinary.h"
#include "brepeat.h"

#define Name "brepeat"

/* Finds all contiguous repeats in the genome specified by the index that fall within the
 * specified unit length range and minimum contiguous length.
 * */

int main(int argc, char *argv[]) 
{
	FILE *fp;
	char rgFileName[MAX_FILENAME_LENGTH]="\0";

	if(argc == 6) {
		RGBinary rg;
		int minUnitLength, maxUnitLength;
		int minLength;
		char unit[MAX_UNIT_LENGTH+1]="\0";

		/* Command line arguments */
		strcpy(rgFileName, argv[1]);
		minUnitLength = atoi(argv[2]);
		maxUnitLength = atoi(argv[3]);
		minLength = atoi(argv[4]);
		assert(minUnitLength <= maxUnitLength);
		assert(maxUnitLength <= MAX_UNIT_LENGTH);

		fprintf(stderr, "minUnitLength:%d\nmaxUnitLength:%d\nminLength:%d\noutput:%s\n",
				minUnitLength,
				maxUnitLength,
				minLength,
				argv[5]);

		if(!(fp = fopen(argv[5], "w"))) {
			PrintError("Main", argv[5], "Could not open file for writing", Exit, OpenFileError);
		}

		/* Read in the rg binary file */
		RGBinaryReadBinary(&rg, rgFileName);

		fprintf(stderr, "Currently on:\n%2d %9d", -1, -1);
		/* For each contig */
		int curContig, curContigIndex;
		for(curContig=1, curContigIndex=0;
				curContig <= rg.numContigs && curContigIndex < rg.numContigs;
				curContig++, curContigIndex++) {
			/* For each starting position */
			int curPos;
			int prevBestUnitLength=-1;
			int prevBestStart=-1;
			int prevBestLength=-1;
			int prevBestEnd=-1;
			for(curPos=1;
					curPos <= rg.contigs[curContigIndex].sequenceLength;
					curPos++) {
				if(curPos%BREPEAT_ROTATE_NUM==0) {
					fprintf(stderr, "\r%2d %9d",
							curContig,
							curPos);
				}

				if(ToLower(RGBinaryGetBase(&rg, curContig, curPos)) != 'n') {

					/* For each unit length */
					int bestEnd=-1;
					char bestUnit[MAX_UNIT_LENGTH+1]="\0";
					int bestUnitLength=-1;
					int curUnitLength;
					for(curUnitLength=minUnitLength;curUnitLength<=maxUnitLength;curUnitLength++) { /* For each unit length */
						/* Check bounds */
						if(curPos + curUnitLength - 1 <= rg.contigs[curContigIndex].sequenceLength) {
							/* Get the current unit */
							int i;
							for(i=0;i<curUnitLength;i++) {
								unit[i] = ToLower(RGBinaryGetBase(&rg,
											curContig,
											curPos+i));
							}
							unit[curUnitLength]='\0';
							int end = curPos+curUnitLength-1;

							/* extend the unit */
							int cont=1;
							int curStart;
							for(curStart=curPos + curUnitLength;cont==1 && curStart <= rg.contigs[curContigIndex].sequenceLength-curUnitLength+1;curStart+=curUnitLength) {
								/* Check that the bases match cur */
								int j;
								for(j=0;j<curUnitLength;j++) {
									if(unit[j] != ToLower(RGBinaryGetBase(&rg, curContig, curStart+j))) {
										cont = 0;
									}
								}
								if(cont == 1) {
									end = curStart+curUnitLength-1;
								}
							}

							if(end-curPos+1 >= minLength && 
									(end - curPos +1 ) > curUnitLength && 
									(bestEnd <= 0 || end-curPos+1 > (bestEnd - curPos + 1))) {
								bestEnd = end;
								strcpy(bestUnit, unit);
								bestUnitLength = curUnitLength;
							}
						}

					}
					if(bestEnd > 0 && 
							bestUnitLength < (bestEnd - curPos + 1)) {
						assert(bestEnd-curPos+1 >= minLength);

						if(curPos == prevBestStart + 1 &&
								prevBestUnitLength == bestUnitLength && 
								prevBestLength == (bestEnd-curPos+1)) {
							/* Skip */
						}
						else if(bestEnd <= prevBestEnd) {
							/* Skip */
						}
						else {
							fprintf(fp, "contig%d:%d-%d\t%d\t%d\t%s\n",
									curContig,
									curPos,
									bestEnd,
									bestEnd-curPos+1,
									bestUnitLength,
									bestUnit);
							fflush(fp);
						}
						prevBestUnitLength = bestUnitLength;
						prevBestStart = curPos;
						prevBestLength = bestEnd-curPos+1;
						prevBestEnd = bestEnd;
					}
				}
			}
		}

		fprintf(stderr, "\n%s", BREAK_LINE);
		fprintf(stderr, "Cleaning up.\n");
		/* Delete the rg */
		RGBinaryDelete(&rg);
		fclose(fp);
		fprintf(stderr, "Terminating successfully.\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "%s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast reference genome file>\n");
		fprintf(stderr, "\t<minimum unit length>\n");
		fprintf(stderr, "\t<maximum unit length>\n");
		fprintf(stderr, "\t<minimum total repeat length>\n");
		fprintf(stderr, "\t<output file>\n");
	}

	return 0;
}

