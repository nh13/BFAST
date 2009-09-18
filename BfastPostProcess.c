#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <config.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#include "BError.h"
#include "BLibDefinitions.h"
#include "BLib.h"
#include "InputOutputToFiles.h"
#include "RunAligner.h"
#include "BfastPostProcess.h"

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC, OPTIONAL_GROUP_NAME}.
   */
enum { 
	DescInputFilesTitle, DescRGFileName, DescInputFileName, 
	DescAlgoTitle, DescAlgorithm, DescQueueLength, 
	DescOutputTitle, DescOutputID, DescOutputDir, DescOutputFormat, DescTiming,
	DescMiscTitle, DescParameters, DescHelp
};

static struct argp_option options[] = {
	{0, 0, 0, 0, "=========== Input Files =============================================================", 1},
	{"brgFileName", 'r', "brgFileName", 0, "Specifies the file name of the reference genome file (not required for BAF output)", 1},
	{"alignFileName", 'i', "alignFileName", 0, "Specifies the input file from the balign program", 1},
	{0, 0, 0, 0, "=========== Algorithm Options =======================================================", 2},
	{"algorithm", 'a', "algorithm", 0, "Specifies the algorithm to choose the alignment for each end of the read after filtering:"
		"\n\t\t\t0: Specifies no filtering will occur"
			"\n\t\t\t1: Specifies that all alignments that pass the filters will be outputted"
			"\n\t\t\t2: Specifies to only consider reads that have been aligned uniquely"
			"\n\t\t\t3: Specifies to choose uniquely the alignment with the best score"
			"\n\t\t\t4: Specifies to choose all alignments with the best score",
		2},
	{"queueLength", 'Q', "queueLength", 0, "Specifies the number of reads to cache", 2},
	{0, 0, 0, 0, "=========== Output Options ==========================================================", 6},
	{"outputID", 'o', "outputID", 0, "Specifies the ID tag to identify the output files", 6},
	{"outputDir", 'd', "outputDir", 0, "Specifies the output directory for the output files", 6},
	{"outputFormat", 'O', "outputFormat", 0, "Specifies the output format 0: BAF 1: MAF 2: GFF 3: SAM", 6},
	{"timing", 't', 0, OPTION_NO_USAGE, "Specifies to output timing information", 6},
	{0, 0, 0, 0, "=========== Miscellaneous Options ===================================================", 7},
	{"Parameters", 'p', 0, OPTION_NO_USAGE, "Print program parameters", 7},
	{"Help", 'h', 0, OPTION_NO_USAGE, "Display usage summary", 7},
	{0, 0, 0, 0, 0, 0}
};

static char OptionString[]=
"a:d:i:o:r:O:hpt";

	int
BfastPostProcess(int argc, char **argv)
{
	struct arguments arguments;
	time_t startTime = time(NULL);
	time_t endTime;
	RGBinary rg;

	if(argc>1) {
		/* Set argument defaults. (overriden if user specifies them)  */ 
		BfastPostProcessAssignDefaultValues(&arguments);

		/* Parse command line args */
		if(BfastPostProcessGetOptParse(argc, argv, OptionString, &arguments)==0)
		{
			switch(arguments.programMode) {
				case ExecuteGetOptHelp:
					BfastPostProcessGetOptHelp();
					break;
				case ExecutePrintProgramParameters:
					BfastPostProcessPrintProgramParameters(stderr, &arguments);
					break;
				case ExecuteProgram:
					if(BfastPostProcessValidateInputs(&arguments)) {
						fprintf(stderr, "Input arguments look good!\n");
						fprintf(stderr, BREAK_LINE);
					}
					else {
						PrintError("PrintError",
								NULL,
								"validating command-line inputs",
								Exit,
								InputArguments);
					}
					BfastPostProcessPrintProgramParameters(stderr, &arguments);
					/* Execute program */
					if(BAF != arguments.outputFormat) {
						/* Read binary */
						RGBinaryReadBinary(&rg,
								arguments.brgFileName);
					}
					ReadInputFilterAndOutput(&rg,
							arguments.alignFileName,
							arguments.algorithm,
							arguments.queueLength,
							arguments.outputID,
							arguments.outputDir,
							arguments.outputFormat);
					if(BAF != arguments.outputFormat) {
						/* Free rg binary */
						RGBinaryDelete(&rg);
					}

					if(arguments.timing == 1) {
						/* Get the time information */
						endTime = time(NULL);
						int seconds = endTime - startTime;
						int hours = seconds/3600;
						seconds -= hours*3600;
						int minutes = seconds/60;
						seconds -= minutes*60;
						fprintf(stderr, "Total time elapsed: %d hours, %d minutes and %d seconds.\n",
								hours,
								minutes,
								seconds
							   );
					}
					fprintf(stderr, "Terminating successfully!\n");
					fprintf(stderr, "%s", BREAK_LINE);
					break;
				default:
					PrintError("PrintError",
							"programMode",
							"Could not determine program mode",
							Exit,
							OutOfRange);
			}

		}
		else {
			PrintError("PrintError",
					NULL,
					"Could not parse command line argumnets",
					Exit,
					InputArguments);
		}
		/* Free program parameters */
		BfastPostProcessFreeProgramParameters(&arguments);
	}
	else {
		BfastPostProcessGetOptHelp();
	}
	return 0;
}

/* TODO */
int BfastPostProcessValidateInputs(struct arguments *args) {

	char *FnName="BfastPostProcessValidateInputs";

	fprintf(stderr, BREAK_LINE);
	fprintf(stderr, "Checking input parameters supplied by the user ...\n");

	if(args->brgFileName!=0) {
		fprintf(stderr, "Validating brgFileName %s. \n",
				args->brgFileName);
		if(ValidateFileName(args->brgFileName)==0)
			PrintError(FnName, "brgFileName", "Command line argument", Exit, IllegalFileName);
	}

	if(args->alignFileName!=0) {
		fprintf(stderr, "Validating alignFileName %s. \n",
				args->alignFileName);
		if(ValidateFileName(args->alignFileName)==0)
			PrintError(FnName, "alignFileName", "Command line argument", Exit, IllegalFileName);
	}

	if(args->algorithm < MIN_FILTER || 
			args->algorithm > MAX_FILTER) {
		PrintError(FnName, "algorithm", "Command line argument", Exit, OutOfRange);
	}

	if(args->queueLength<=0) {
		PrintError(FnName, "queueLength", "Command line argument", Exit, OutOfRange);
	}

	if(args->outputID!=0) {
		fprintf(stderr, "Validating outputID %s. \n",
				args->outputID);
		if(ValidateFileName(args->outputID)==0)
			PrintError(FnName, "outputID", "Command line argument", Exit, IllegalFileName);
	}

	if(args->outputDir!=0) {
		fprintf(stderr, "Validating outputDir path %s. \n",
				args->outputDir);
		if(ValidateFileName(args->outputDir)==0)
			PrintError(FnName, "outputDir", "Command line argument", Exit, IllegalFileName);
	}

	if(!(args->outputFormat == BAF ||
				args->outputFormat == MAF ||
				args->outputFormat == GFF ||
				args->outputFormat == SAM)) {
		PrintError(FnName, "outputFormat", "Command line argument", Exit, OutOfRange);
	}

	assert(args->timing == 0 || args->timing == 1);

	return 1;
}

/* TODO */
	void 
BfastPostProcessAssignDefaultValues(struct arguments *args)
{
	/* Assign default values */

	args->programMode = ExecuteProgram;

	args->brgFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->brgFileName!=0);
	strcpy(args->brgFileName, DEFAULT_FILENAME);

	args->alignFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->alignFileName!=0);
	strcpy(args->alignFileName, DEFAULT_FILENAME);

	args->algorithm=0;
	args->queueLength=DEFAULT_QUEUE_LENGTH;

	args->outputID = 
		(char*)malloc(sizeof(DEFAULT_OUTPUT_ID));
	assert(args->outputID!=0);
	strcpy(args->outputID, DEFAULT_OUTPUT_ID);

	args->outputDir =
		(char*)malloc(sizeof(DEFAULT_OUTPUT_DIR));
	assert(args->outputDir!=0);
	strcpy(args->outputDir, DEFAULT_OUTPUT_DIR);

	args->outputFormat=BAF;

	args->timing = 0;

	return;
}

/* TODO */
	void 
BfastPostProcessPrintProgramParameters(FILE* fp, struct arguments *args)
{
	char programmode[3][64] = {"ExecuteGetOptHelp", "ExecuteProgram", "ExecutePrintProgramParameters"};
	char algorithm[5][64] = {"No Filtering", "Filtering Only", "Unique", "Best Score", "Best Score All"};
	char outputType[8][32] = {"BRG", "BIF", "BMF", "BAF", "MAF", "GFF", "SAM", "LastFileType"};
	fprintf(fp, BREAK_LINE);
	fprintf(fp, "Printing Program Parameters:\n");
	fprintf(fp, "programMode:\t\t%d\t[%s]\n", args->programMode, programmode[args->programMode]);
	fprintf(fp, "brgFileName:\t\t%s\n", args->brgFileName);
	fprintf(fp, "alignFileName:\t\t%s\n", args->alignFileName);
	fprintf(fp, "algorithm:\t\t%d\t[%s]\n", args->algorithm, algorithm[args->algorithm]);
	fprintf(fp, "queueLength:\t\t%d\n", args->queueLength);
	fprintf(fp, "outputID:\t\t%s\n", args->outputID);
	fprintf(fp, "outputDir:\t\t%s\n", args->outputDir);
	fprintf(fp, "outputFormat:\t\t%s\n", outputType[args->outputFormat]);
	fprintf(fp, "timing:\t\t\t%d\n", args->timing);
	fprintf(fp, BREAK_LINE);
	return;
}

/* TODO */
void BfastPostProcessFreeProgramParameters(struct arguments *args)
{
	free(args->brgFileName);
	args->brgFileName=NULL;
	free(args->alignFileName);
	args->alignFileName=NULL;
	free(args->outputID);
	args->outputID=NULL;
	free(args->outputDir);
	args->outputDir=NULL;
}

/* TODO */
void
BfastPostProcessGetOptHelp() {

	struct argp_option *a=options;
	fprintf(stderr, "%s %s\n", "bfast ", PACKAGE_VERSION);
	fprintf(stderr, "\nUsage: bfast postprocess [options]\n");
	while((*a).group>0) {
		switch((*a).key) {
			case 0:
				fprintf(stderr, "\n%s\n", (*a).doc); break;
			default:
				if((*a).arg != 0) {
					fprintf(stderr, "-%c\t%12s\t%s\n", (*a).key, (*a).arg, (*a).doc); 
				}
				else {
					fprintf(stderr, "-%c\t%12s\t%s\n", (*a).key, "", (*a).doc); 
				}
				break;
		}
		a++;
	}
	fprintf(stderr, "\nsend bugs to %s\n", 
			PACKAGE_BUGREPORT);
	return;
}

/* TODO */
	int
BfastPostProcessGetOptParse(int argc, char** argv, char OptionString[], struct arguments* arguments) 
{
	char key;
	int OptErr=0;
	while((OptErr==0) && ((key = getopt (argc, argv, OptionString)) != -1)) {
		/*
		   fprintf(stderr, "Key is %c and OptErr = %d\n", key, OptErr);
		   */
		switch (key) {
			case 'a':
				arguments->algorithm = atoi(optarg);break;
			case 'd':
				StringCopyAndReallocate(&arguments->outputDir, optarg);
				break;
			case 'h':
				arguments->programMode=ExecuteGetOptHelp;break;
			case 'i':
				StringCopyAndReallocate(&arguments->alignFileName, optarg);
				break;
			case 'o':
				StringCopyAndReallocate(&arguments->outputID, optarg);
				break;
			case 'p':
				arguments->programMode=ExecutePrintProgramParameters;break;
			case 'r':
				StringCopyAndReallocate(&arguments->brgFileName, optarg);break;
			case 't':
				arguments->timing = 1;break;
			case 'O':
				switch(atoi(optarg)) {
					case 0:
						arguments->outputFormat = BAF;
						break;
					case 1:
						arguments->outputFormat = MAF;
						break;
					case 2:
						arguments->outputFormat = GFF;
						break;
					case 3:
						arguments->outputFormat = SAM;
						break;
					default:
						arguments->outputFormat = -1;
						/* Deal with this when we validate the input parameters */
						break;
				}
				break;
			case 'Q':
				arguments->queueLength=atoi(optarg);break;
			default:
				OptErr=1;
		} /* while */
	} /* switch */
	return OptErr;
}
