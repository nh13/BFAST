#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>

#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "../blib/RGBinary.h"
#include "brepeat.h"

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
			PrintError("Main",
					argv[5],
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}

		/* Read in the rg binary file */
		RGBinaryReadBinary(&rg, rgFileName);

		fprintf(stderr, "Currently on:\n%2d %9d", -1, -1);
		/* For each chromosome */
		int curChr, curChrIndex;
		for(curChr=rg.startChr, curChrIndex=0;
				curChr < rg.endChr && curChrIndex < rg.numChrs;
				curChr++, curChrIndex++) {
			/* For each starting position */
			int curPos;
			int prevBestUnitLength=-1;
			int prevBestStart=-1;
			int prevBestLength=-1;
			int prevBestEnd=-1;
			for(curPos=rg.chromosomes[curChrIndex].startPos;
					curPos <= rg.chromosomes[curChrIndex].endPos;
					curPos++) {
				if(curPos%BREPEAT_ROTATE_NUM==0) {
					fprintf(stderr, "\r%2d %9d",
							curChr,
							curPos);
				}

				if(ToLower(RGBinaryGetBase(&rg, curChr, curPos)) != 'n') {

					/* For each unit length */
					int bestEnd=-1;
					char bestUnit[MAX_UNIT_LENGTH+1]="\0";
					int bestUnitLength=-1;
					int curUnitLength;
					for(curUnitLength=minUnitLength;curUnitLength<=maxUnitLength;curUnitLength++) { /* For each unit length */
						/* Check bounds */
						if(curPos + curUnitLength - 1 <= rg.chromosomes[curChrIndex].endPos) {
							/* Get the current unit */
							int i;
							for(i=0;i<curUnitLength;i++) {
								unit[i] = ToLower(RGBinaryGetBase(&rg,
											curChr,
											curPos+i));
							}
							unit[curUnitLength]='\0';
							int end = curPos+curUnitLength-1;

							/* extend the unit */
							int cont=1;
							int curStart;
							for(curStart=curPos + curUnitLength;cont==1 && curStart <= rg.chromosomes[curChrIndex].endPos-curUnitLength+1;curStart+=curUnitLength) {
								/* Check that the bases match cur */
								int j;
								for(j=0;j<curUnitLength;j++) {
									if(unit[j] != ToLower(RGBinaryGetBase(&rg, curChr, curStart+j))) {
										cont = 0;
									}
								}
								if(cont == 1) {
									end = curStart+curUnitLength-1;
								}
							}

							if(end-curPos+1 >= minLength && 
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
							fprintf(fp, "chr%d:%d-%d\t%d\t%d\t%s\n",
									curChr,
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

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Cleaning up.\n");
		/* Delete the rg */
		RGBinaryDelete(&rg);
		fclose(fp);
	}
	else {
		fprintf(stderr, "Incorrect number of arguments.  See source code.\n");
	}

	return 0;
}

