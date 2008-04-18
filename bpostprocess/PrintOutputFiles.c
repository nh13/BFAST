#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/AlignEntry.h"
#include "Definitions.h"
#include "PrintOutputFiles.h"

/* TODO */
void ConvertAndPrint(void *entries,
		int numEntries,
		int inputFormat,
		int curChr,
		char *outputID,
		char *outputDir,
		int outputFormat)
{
	char *FnName="ConvertAndPrint";
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
	FILE *fp;

	/* Create the file name */
	switch(outputFormat) {
		case WigFile:
			sprintf(outputFileName, "%sblatter.%s.%d.%s",
					outputDir,
					outputID,
					curChr,
					".wig");
			break;
		case BedFile:
			sprintf(outputFileName, "%sblatter.%s.%d.%s",
					outputDir,
					outputID,
					curChr,
					".bed");
			break;
		default:
			PrintError(FnName,
					"outputFormat",
					"Output format not yet supported",
					Exit,
					OutOfRange);
	}

	/* Open output file */
	if(0==(fp=fopen(outputFileName, "w"))) {
		PrintError(FnName,
				outputFileName,
				"Could not open output file for writing",
				Exit,
				OpenFileError);
	}

	switch(inputFormat) {
		case BAlignFile:
			switch(outputFormat) {
				case WigFile:
					ConvertAndPrintAlignEntryToWig((AlignEntry*)entries,
							numEntries,
							fp);
					break;
				case BedFile:
					/* Do not break, not supported */
				default:
					PrintError(FnName,
							"outputFormat",
							"Output format not yet supported",
							Exit,
							OutOfRange);
			}
			break;
		default:
			PrintError(FnName,
					"inputFormat",
					"Input format not yet supported",
					Exit,
					OutOfRange);
			break;
	}

	/* Close output file */
	fclose(fp);
}

/* TODO */
void ConvertAndPrintAlignEntryToWig(AlignEntry *entries,
		int numEntries,
		FILE *fp)
{
}
