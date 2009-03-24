#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>

/* 
   Based on GNU/Linux and glibc or not, argp.h may or may not be available.
   If it is not, fall back to getopt. Also see #ifdefs in the ParseInput.h file. 
   */
#ifdef HAVE_ARGP_H
#include <argp.h>
#define OPTARG arg
#elif defined HAVE_UNISTD_H
#include <unistd.h>
#define OPTARG optarg
#else
#include "GetOpt.h"
#define OPTARG optarg
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h> /* For u_int etc. */
#endif

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h> /* For Mac OS X with resource.h */
#endif

#ifdef HAVE_RESOURCE_H
#include <sys/resource.h>
#endif

#include <time.h>
#include "../blib/BError.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BLib.h"
#include "Definitions.h"
#include "InputOutputToFiles.h"
#include "ParseInput.h"
const char *argp_program_version =
"bpostprocess "
PACKAGE_VERSION
"\n"
"Copyright 2008-2009";

const char *argp_program_bug_address =
PACKAGE_BUGREPORT;

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC, OPTIONAL_GROUP_NAME}.
   */
enum { 
	DescInputFilesTitle, DescRGFileName, DescInputFileName, 
	DescAlgoTitle, DescAlgorithm, 
	DescGenFiltTitle, DescStartContig, DescStartPos, DescEndContig, DescEndPos, DescMinScore, escMaxMismatches, DescMaxColorErrors, 
	DescPairedEndTitle, DescMinDistancePaired, DescMaxDistancePaired, DescContigAbPaired, DescInversionsPaired, DescUnpaired,
	DescOutputTitle, DescOutputID, DescOutputDir, DescOutputFormat, DescTiming,
	DescMiscTitle, DescParameters, DescHelp
};

/* 
   The prototype for argp_option comes fron argp.h. If argp.h
   absent, then ParseInput.h declares it 
   */
static struct argp_option options[] = {
	{0, 0, 0, 0, "=========== Input Files =============================================================", 1},
	{"rgFileName", 'r', "rgFileName", 0, "Specifies the file name of the reference genome file (not required for BAF output)", 1},
	{"alignFileName", 'i', "alignFileName", 0, "Specifies the input file from the balign program", 1},
	/*
	   {"binaryInput", 'b', 0, OPTION_NO_USAGE, "Specifies that the input files will be in binary format", 1},
	   */
	{0, 0, 0, 0, "=========== Algorithm Options =======================================================", 2},
	{"algorithm", 'a', "algorithm", 0, "Specifies the algorithm to choose the alignment for each end of the read after filtering:"
		"\n\t\t\t0: Specifies no filtering will occur"
		"\n\t\t\t1: Specifies that all alignments that pass the filters will be outputted"
		"\n\t\t\t2: Specifies to only consider reads that have been aligned uniquely"
		"\n\t\t\t3: Specifies to choose the alignment with the best score",
		2},
	{0, 0, 0, 0, "=========== General Filter Options ==================================================", 3},
	{"startContig", 's', "startContig", 0, "Specifies the start contig for filtering", 3},
	{"startPos", 'S', "startPos", 0, "Specifies the end position for filtering", 3},
	{"endContig", 'e', "endContig", 0, "Specifies the end contig for filtering", 3},
	{"endPos", 'E', "endPos", 0, "Specifies the end postion for filtering", 3},
	{"minScore", 'm', "minScore", 0, "Specifies the minimum score to consider for a given end of a read", 3},
	{"maxMismatches", 'j', "maxMismatches", 0, "Specifies the maximum number of mismatches to consider for a given end of a read", 3},
	{"maxColorErrors", 'k', "maxColorErrors", 0, "Specifies the maximum number of color errors to consider for a given end of a read", 3},
	{0, 0, 0, 0, "=========== Paired End Filter Options ===============================================", 4},
	{"minDistancePaired", 'X', "minDistancePaired", 0, "Specifies the minimum allowable distance between the paired ends for filtering", 4},
	{"maxDistancePaired", 'Y', "maxDistancePaired", 0, "Specifies the maximum allowable distance between the paired ends for filtering", 4},
	{"contigAbPaired", 'C', 0, OPTION_NO_USAGE, "Specifies to output separately those paired reads that do not fall within the\n"
		"\t\t\tspecified distance but are on the same strand (paired end only)", 2},
	{"inversionsPaired", 'I', 0, OPTION_NO_USAGE, "Specifies to output separately those paired reads that do not fall within the\n"
		"\t\t\tspecified distance but are on the opposite strands (paired end only)", 2},
	{"unpaired", 'U', 0, OPTION_NO_USAGE, "Specifies to output separately those paired reads for which one end could not be unambiguously chosen", 2},
	{0, 0, 0, 0, "=========== Output Options ==========================================================", 5},
	{"outputID", 'o', "outputID", 0, "Specifies the ID tag to identify the output files", 5},
	{"outputDir", 'd', "outputDir", 0, "Specifies the output directory for the output files", 5},
	{"outputFormat", 'O', "outputFormat", 0, "Specifies the output format 0: BAF 1: MAF 2: GFF 3: SAM", 5},
	{"timing", 't', 0, OPTION_NO_USAGE, "Specifies to output timing information", 5},
	{0, 0, 0, 0, "=========== Miscellaneous Options ===================================================", 6},
	{"Parameters", 'p', 0, OPTION_NO_USAGE, "Print program parameters", 6},
	{"Help", 'h', 0, OPTION_NO_USAGE, "Display usage summary", 6},
	{0, 0, 0, 0, 0, 0}
};

#ifdef HAVE_ARGP_H
/*
   ARGS_DOC. Field 3 in ARGP.
   A description of the non-option command-line arguments that we accept.
   Not complete yet. So empty string
   */
static char args_doc[] = "";
/*
   DOC.  Field 4 in ARGP.  Program documentation.
   */
static char doc[] = "";
/*
   The ARGP structure itself.
   */
static struct argp argp = {options, parse_opt, args_doc, doc};
#else
/* argp.h support not available! Fall back to getopt */
static char OptionString[]=
"a:d:e:i:j:k:m:o:r:s:z:E:I:O:P:S:T:X:Y:hptCIU";
#endif

enum {ExecuteGetOptHelp, ExecuteProgram, ExecutePrintProgramParameters};

/*
   The main function. All command-line options parsed using argp_parse
   or getopt whichever available
   */
	int
main (int argc, char **argv)
{
	struct arguments arguments;
	time_t startTime = time(NULL);
	time_t endTime;
	RGBinary rg;

	if(argc>1) {
		/* Set argument defaults. (overriden if user specifies them)  */ 
		AssignDefaultValues(&arguments);

		/* Parse command line args */
#ifdef HAVE_ARGP_H
		if(argp_parse(&argp, argc, argv, 0, 0, &arguments)==0)
#else
			if(getopt_parse(argc, argv, OptionString, &arguments)==0)
#endif 
			{
				switch(arguments.programMode) {
					case ExecuteGetOptHelp:
						GetOptHelp();
						break;
					case ExecutePrintProgramParameters:
						PrintProgramParameters(stderr, &arguments);
						break;
					case ExecuteProgram:
						if(ValidateInputs(&arguments)) {
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
						PrintProgramParameters(stderr, &arguments);
						/* Execute program */
						if(BAF != arguments.outputFormat) {
							/* Read binary */
							RGBinaryReadBinary(&rg,
									arguments.rgFileName);
						}
						ReadInputFilterAndOutput(&rg,
								arguments.alignFileName,
								arguments.binaryInput,
								arguments.startContig,
								arguments.startPos,
								arguments.endContig,
								arguments.endPos,
								arguments.algorithm,
								arguments.minScore,
								arguments.maxMismatches,
								arguments.maxColorErrors,
								arguments.minDistancePaired,
								arguments.maxDistancePaired,
								arguments.useDistancePaired,
								arguments.contigAbPaired,
								arguments.inversionsPaired,
								arguments.unpaired,
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
		FreeProgramParameters(&arguments);
	}
	else {
		GetOptHelp();
#ifdef HAVE_ARGP_H
		/* fprintf(stderr, "Type \"%s --help\" to see usage\n", argv[0]); */
#else
		/*     fprintf(stderr, "Type \"%s -h\" to see usage\n", argv[0]); */
#endif
	}
	return 0;
}

/* TODO */
int ValidateInputs(struct arguments *args) {

	char *FnName="ValidateInputs";

	fprintf(stderr, BREAK_LINE);
	fprintf(stderr, "Checking input parameters supplied by the user ...\n");

	if(args->rgFileName!=0) {
		fprintf(stderr, "Validating rgFileName %s. \n",
				args->rgFileName);
		if(ValidateFileName(args->rgFileName)==0)
			PrintError(FnName, "rgFileName", "Command line argument", Exit, IllegalFileName);
	}

	if(args->alignFileName!=0) {
		fprintf(stderr, "Validating alignFileName %s. \n",
				args->alignFileName);
		if(ValidateFileName(args->alignFileName)==0)
			PrintError(FnName, "alignFileName", "Command line argument", Exit, IllegalFileName);
	}

	assert(args->binaryInput == TextInput || args->binaryInput == BinaryInput);

	if(args->startContig < 0) {
		PrintError(FnName, "startContig", "Command line argument", Exit, OutOfRange);
	}

	if(args->startPos < 0) {
		PrintError(FnName, "startPos", "Command line argument", Exit, OutOfRange);
	}

	if(args->endContig < 0) {
		PrintError(FnName, "endContig", "Command line argument", Exit, OutOfRange);
	}

	if(args->endPos < 0) {
		PrintError(FnName, "endPos", "Command line argument", Exit, OutOfRange);
	}

	if(args->algorithm < MIN_FILTER || 
			args->algorithm > MAX_FILTER) {
		PrintError(FnName, "algorithm", "Command line argument", Exit, OutOfRange);
	}
	
	if(args->maxMismatches < 0) {
		PrintError(FnName, "maxMismatches < 0", "Command line argument", Exit, OutOfRange);
	}

	if(args->maxColorErrors < 0) {
		PrintError(FnName, "maxColorErrors < 0", "Command line argument", Exit, OutOfRange);
	}

	/* Check that the min distance is less than or equal to the max distance */
	if(args->minDistancePaired > args->maxDistancePaired) {
		PrintError(FnName, "minDistancePaired > maxDistancePaired", "Command line argument", Exit, OutOfRange);
	}

	/* This must hold internally */
	assert(0 == args->useDistancePaired || 1 == args->useDistancePaired);
	assert(args->contigAbPaired == 0 || args->contigAbPaired == 1);
	assert(args->inversionsPaired == 0 || args->inversionsPaired == 1);
	assert(args->unpaired == 0 || args->unpaired == 1);

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
AssignDefaultValues(struct arguments *args)
{
	/* Assign default values */

	args->programMode = ExecuteProgram;

	args->rgFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->rgFileName!=0);
	strcpy(args->rgFileName, DEFAULT_FILENAME);

	args->alignFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->alignFileName!=0);
	strcpy(args->alignFileName, DEFAULT_FILENAME);

	args->binaryInput = BALIGN_DEFAULT_OUTPUT;

	args->startContig=0;
	args->startPos=0;
	args->endContig=0;
	args->endPos=0;

	args->algorithm=0;
	args->minScore=INT_MIN;
	args->maxMismatches=INT_MAX;
	args->maxColorErrors=INT_MAX;

	args->minDistancePaired=INT_MIN;
	args->maxDistancePaired=INT_MAX;
	args->useDistancePaired=0;
	args->contigAbPaired=0;
	args->inversionsPaired=0;
	args->unpaired=0;

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
PrintProgramParameters(FILE* fp, struct arguments *args)
{
	char programmode[3][64] = {"ExecuteGetOptHelp", "ExecuteProgram", "ExecutePrintProgramParameters"};
	char algorithm[4][64] = {"No Filtering", "Filtering Only", "Unique", "Best Score"};
	fprintf(fp, BREAK_LINE);
	fprintf(fp, "Printing Program Parameters:\n");
	fprintf(fp, "programMode:\t\t%d\t[%s]\n", args->programMode, programmode[args->programMode]);
	fprintf(fp, "rgFileName:\t\t%s\n", args->rgFileName);
	fprintf(fp, "alignFileName:\t\t%s\n", args->alignFileName);
	/*
	   fprintf(fp, "binaryInput:\t\t%d\n", args->binaryInput);
	   */
	fprintf(fp, "algorithm:\t\t%d\t[%s]\n", args->algorithm, algorithm[args->algorithm]);
	fprintf(fp, "startContig:\t\t%d\n", args->startContig);
	fprintf(fp, "startPos:\t\t%d\n", args->startPos);
	fprintf(fp, "endContig:\t\t%d\n", args->endContig);
	fprintf(fp, "endPos:\t\t\t%d\n", args->endPos);
	fprintf(fp, "minScore:\t\t%d\n", args->minScore);
	fprintf(fp, "maxMismatches:\t\t%d\n", args->maxMismatches);
	fprintf(fp, "maxColorErrors:\t\t%d\n", args->maxColorErrors);
	fprintf(fp, "minDistancePaired:\t%d\n", args->minDistancePaired);
	fprintf(fp, "maxDistancePaired:\t%d\n", args->maxDistancePaired);
	fprintf(fp, "contigAbPaired:\t\t%d\n", args->contigAbPaired);
	fprintf(fp, "inversionsPaired:\t%d\n", args->inversionsPaired);
	fprintf(fp, "unpaired:\t\t%d\n", args->unpaired);
	fprintf(fp, "outputID:\t\t%s\n", args->outputID);
	fprintf(fp, "outputDir:\t\t%s\n", args->outputDir);
	fprintf(fp, "outputFormat:\t\t%d\n", args->outputFormat);
	fprintf(fp, "timing:\t\t\t%d\n", args->timing);
	fprintf(fp, BREAK_LINE);
	return;
}

/* TODO */
void FreeProgramParameters(struct arguments *args)
{
	free(args->rgFileName);
	args->rgFileName=NULL;
	free(args->alignFileName);
	args->alignFileName=NULL;
	free(args->outputID);
	args->outputID=NULL;
	free(args->outputDir);
	args->outputDir=NULL;
}

/* TODO */
void
GetOptHelp() {

	struct argp_option *a=options;
	fprintf(stderr, "%s\n", argp_program_version);
	fprintf(stderr, "\nUsage: bpostprocess [options]\n");
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
	fprintf(stderr, "\nsend bugs to %s\n", argp_program_bug_address);
	return;
}

/* TODO */
#ifdef HAVE_ARGP_H
	static error_t
parse_opt (int key, char *arg, struct argp_state *state) 
{
	struct arguments *arguments = state->input;
#else
	int
		getopt_parse(int argc, char** argv, char OptionString[], struct arguments* arguments) 
		{
			char key;
			int OptErr=0;
			while((OptErr==0) && ((key = getopt (argc, argv, OptionString)) != -1)) {
				/*
				   fprintf(stderr, "Key is %c and OptErr = %d\n", key, OptErr);
				   */
#endif
				switch (key) {
					case 'a':
						arguments->algorithm = atoi(OPTARG);break;
						/*
						   case 'b':
						   arguments->binaryInput = 1;break;
						   */
					case 'd':
						StringCopyAndReallocate(&arguments->outputDir, OPTARG);
						break;
					case 'e':
						arguments->endContig=atoi(OPTARG);break;
					case 'h':
						arguments->programMode=ExecuteGetOptHelp;break;
					case 'i':
						StringCopyAndReallocate(&arguments->alignFileName, OPTARG);
						break;
					case 'j':
						arguments->maxMismatches=atoi(OPTARG);break;
					case 'k':
						arguments->maxColorErrors=atoi(OPTARG);break;
					case 'm':
						arguments->minScore = atoi(OPTARG);break;
					case 'o':
						StringCopyAndReallocate(&arguments->outputID, OPTARG);
						break;
					case 'p':
						arguments->programMode=ExecutePrintProgramParameters;break;
					case 'r':
						StringCopyAndReallocate(&arguments->rgFileName, OPTARG);break;
					case 's':
						arguments->startContig=atoi(OPTARG);break;
					case 't':
						arguments->timing = 1;break;
					case 'C':
						arguments->contigAbPaired = 1;break;
					case 'E':
						arguments->endPos=atoi(OPTARG);break;
					case 'I':
						arguments->inversionsPaired = 1;break;
					case 'O':
						switch(atoi(OPTARG)) {
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
					case 'S':
						arguments->startPos=atoi(OPTARG);break;
					case 'U':
						arguments->unpaired = 1;break;
					case 'X':
						arguments->useDistancePaired = 1;
						arguments->minDistancePaired = atoi(OPTARG);break;
					case 'Y':
						arguments->useDistancePaired = 1;
						arguments->maxDistancePaired = atoi(OPTARG);break;
					default:
#ifdef HAVE_ARGP_H
						return ARGP_ERR_UNKNOWN;
				} /* switch */
				return 0;
#else
				OptErr=1;
			} /* while */
		} /* switch */
	return OptErr;
#endif
}
