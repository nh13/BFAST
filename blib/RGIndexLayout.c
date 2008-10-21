#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "../blib/BError.h"
#include "../blib/BLibDefinitions.h"
#include "RGIndexLayout.h"

/* TODO */
void RGIndexLayoutRead(char *layoutFileName, RGIndexLayout *layout)
{
	char *FnName="RGIndexLayoutRead";
	int i;
	char tempMask[MAX_MASK_LENGTH]="\0";
	int32_t tempHashWidth;
	FILE *fp=NULL;

	/* Open the file */
	if((fp=fopen(layoutFileName, "r"))==0) {
		PrintError(FnName,
				NULL,
				"Could not open index layout file reading",
				Exit,
				OpenFileError);
	}

	/* Initialize the layout data structure */
	layout->numIndexes=0;
	layout->hashWidths=NULL;
	layout->masks=NULL;
	layout->widths=NULL;
	layout->keysizes=NULL;

	/* Try reading in a layout */
	while(EOF != fscanf(fp, "%d %s", &tempHashWidth, tempMask)) {
		layout->numIndexes++;
		/* Reallocate memory */
		layout->hashWidths = realloc(layout->hashWidths, sizeof(int32_t)*layout->numIndexes);
		if(NULL==layout->hashWidths) {
			PrintError(FnName,
					"layout->hashWidths",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		layout->keysizes = realloc(layout->keysizes, sizeof(int32_t)*layout->numIndexes);
		if(NULL==layout->keysizes) {
			PrintError(FnName,
					"layout->keysizes",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		layout->widths = realloc(layout->widths, sizeof(int32_t)*layout->numIndexes);
		if(NULL==layout->widths) {
			PrintError(FnName,
					"layout->widths",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		layout->masks = realloc(layout->masks, sizeof(int32_t*)*layout->numIndexes);
		if(NULL==layout->masks) {
			PrintError(FnName,
					"layout->masks",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		assert(tempHashWidth > 0);
		/* Copy over hash width and width */
		layout->hashWidths[layout->numIndexes-1] = tempHashWidth;
		layout->widths[layout->numIndexes-1] = strlen(tempMask);
		/* Allocate memory for the mask */
		layout->masks[layout->numIndexes-1] = malloc(sizeof(int32_t)*layout->widths[layout->numIndexes-1]);
		if(NULL==layout->masks[layout->numIndexes-1]) {
			PrintError(FnName,
					"layout->masks[layout->numIndexes-1]",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		/* Copy over from temp mask */
		layout->keysizes[layout->numIndexes-1]=0;
		for(i=0;i<layout->widths[layout->numIndexes-1];i++) {
			switch(tempMask[i]) {
				case '0':
					layout->masks[layout->numIndexes-1][i] = 0;
					break;
				case '1':
					layout->masks[layout->numIndexes-1][i] = 1;
					layout->keysizes[layout->numIndexes-1]++;
					break;
				default:
					PrintError(FnName,
							NULL,
							"Could not understand mask",
							Exit,
							OutOfRange);
			}
		}
		/* Check mask begins and ends with a one */
		if(layout->masks[layout->numIndexes-1][0] == 0) {
			PrintError(FnName,
					NULL,
					"Layout must begin with a one",
					Exit,
					OutOfRange);
		}
		if(layout->masks[layout->numIndexes-1][layout->widths[layout->numIndexes-1]-1] == 0) {
			PrintError(FnName,
					NULL,
					"Layout must begin with a one",
					Exit,
					OutOfRange);
		}
	}

	/* Close the file */
	fclose(fp);
}

void RGIndexLayoutDelete(RGIndexLayout *layout)
{
	int i;

	free(layout->hashWidths);
	layout->hashWidths = NULL;
	for(i=0;i<layout->numIndexes;i++) {
		free(layout->masks[i]);
		layout->masks[i] = NULL;
	}
	free(layout->masks);
	layout->masks=NULL;
	free(layout->widths);
	layout->widths=NULL;
	free(layout->keysizes);
	layout->keysizes=NULL;
	/* Reinitialize the number of indexes */
	layout->numIndexes = 0;
}
