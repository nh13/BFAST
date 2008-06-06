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

		if(!(fp = fopen(argv[5], "w"))) {
			PrintError("Main",
					argv[5],
					"Could not open file for writing",
					Exit,
					OpenFileError);
		}

		/* Read in the rg binary file */
		RGBinaryReadBinary(&rg, rgFileName);

		/* Compute the repeats */
		int curUnitLength;
		fprintf(stderr, "Currently on:\n%3d %2d %9d", -1, -1, -1);
		for(curUnitLength=minUnitLength;curUnitLength<=maxUnitLength;curUnitLength++) { /* For each unit length */
			/* For each chromosome */
			int curChr, curChrIndex;
			for(curChr=rg.startChr, curChrIndex=0;
					curChr < rg.endChr && curChrIndex < rg.numChrs;
					curChr++, curChrIndex++) {
				/* For each starting position */
				int curPos=rg.startPos;
				int curPosRelative=rg.startPos;
				while(curPos <= rg.endPos-curUnitLength+1) {
					if(curPos > curPosRelative + BREPEAT_ROTATE_INC) {
						fprintf(stderr, "\r%3d %2d %9d",
								curUnitLength,
								curChr,
								curPos);
						curPosRelative = curPos;
					}
					/* Get the current unit */
					int i;
					for(i=0;i<curUnitLength;i++) {
						unit[i] = ToLower(RGBinaryGetBase(&rg,
									curChr,
									curPos+i));
					}
					unit[curUnitLength]='\0';
					int start = curPos;
					int end = curPos+curUnitLength-1;

					int cont=1;
					int curStart;
					for(curStart=curPos + curUnitLength;cont==1 && curStart <= rg.endPos-curUnitLength+1;curStart+=curUnitLength) {
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

					if(end-start+1 >= minLength) {
						fprintf(fp, "chr%d:%d-%d\t%d\t%s\n",
								curChr,
								start,
								end,
								end-start+1,
								unit);
					}

					/* Update curPos */
					/* Skip over what we found */
					curPos = end+1;
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

