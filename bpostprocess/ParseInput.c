#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
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

#include "../blib/BError.h"
#include "../blib/BLibDefinitions.h"
#include "Definitions.h"
#include "ParseInput.h"

const char *argp_program_version =
"bpostprocess version 0.1.1\n"
"Copyright 2008.";

const char *argp_program_bug_address =
"Nils Homer <nhomer@cs.ucla.edu>";

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC, OPTIONAL_GROUP_NAME}.
   */
enum { 
	DescInputFilesTitle, DescInputFileName, DescBinaryInput,
	DescAlgoTitle, DescAlgorithmReads, DescMinScoreReads, DescStartChr, DescStartPos, DescEndChr, DescEndPos, 
	DescPairedEndTitle, DescPairedEnd, DescAlgorithmPaired, DescMinScoreReadsPaired, DescMinDistancePaired, DescMaxDistancePaired, DescMeanDistancePaired, DescChrAbPaired, DescInversionsPaired,
	DescOutputTitle, DescOutputID, DescOutputDir, DescTmpDir, DescOutputFormat, DescTiming,
	DescMiscTitle, DescParameters, DescHelp
};

/* 
   The prototype for argp_option comes fron argp.h. If argp.h
   absent, then ParseInput.h declares it 
   */
static struct argp_option options[] = {
	{0, 0, 0, 0, "=========== Input Files =============================================================", 1},
	{"inputFileName", 'i', "inputFileName", 0, "Specifies the input file", 1},
	{"binaryInput", 'b', 0, OPTION_NO_USAGE, "Specifies that the input files will be in binary format", 1},
	{0, 0, 0, 0, "=========== Algorithm Options =======================================================", 2},
	{"algorithmReads", 'a', "algorithmReads", 0, "Algorithm to determine alignments:"
		"\n\t\t\t\t0: Specifies to only consider reads that have been aligned uniquely"
			"\n\t\t\t\t1: Specifies to choose the alignment with the best score", 2},
	{"minScoreReads", 'm', "minScoreReads", 0, "Specifies the minimum score to consider", 2},
	{"startChr", 's', "startChr", 0, "Specifies the start chromosome", 2},
	{"startPos", 'S', "startPos", 0, "Specifies the end position", 2},
	{"endChr", 'e', "endChr", 0, "Specifies the end chromosome", 2},
	{"endPos", 'E', "endPos", 0, "Specifies the end postion", 2},
	{0, 0, 0, 0, "=========== Paired End Options ======================================================", 3},
	{"pairedEnd", '2', 0, OPTION_NO_USAGE, "Specifies that paired end data is to be expected", 3},
	{"algorithmReadsPaired", 'A', "algorithmReadsPaired", 0, "Algorithm to determine paired alignments:"
		"\n\t\t\t\t\t0: Specifies to only consider paired reads that have both been aligned uniquely"
			"\n\t\t\t\t\t1: Specifies for both ends of the pair to choose the alignment with the best score"
			"\n\t\t\t\t\t2: Specifies to prefer pairs of reads that are closest to the mean distance and are the unique pair with that idstance  (must specify -Z)", 3},
			"\n\t\t\t\t\t2: Specifies to prefer pairs of reads that are closest to the mean distance and have the best score with that distance (must specify -Z)", 3},
	{"minScoreReadsPaired", 'M', "minScoreReadsPaired", 0, "Specifies the minimum score to consider for the combination of the two paired reads", 3},
	{"minDistancePaired", 'X', "minDistancePaired", 0, "Specifies the minimum allowable distance between the paired ends", 3},
	{"maxDistancePaired", 'Y', "maxDistancePaired", 0, "Specifies the maximum allowable distance between the paired ends", 3},
	{"meanDistancePaired", 'Z', "meanDistancePaired", 0, "Specifies the mean distance between the paired ends (for use with -A 2)", 3},
	{"chrAbPaired", 'C', 0, OPTION_NO_USAGE, "Specifies to output separately those paired reads that do not fall within the specified distance but are on the same strand", 3},
	{"inversionsPaired", 'I', 0, OPTION_NO_USAGE, "Specifies to output separately those paired reads that do not fall within the specified distance but are on opposite strands", 3},
	{0, 0, 0, 0, "=========== Output Options ==========================================================", 4},
	{"outputID", 'o', "outputID", 0, "Specifies the ID tag to identify the output files", 4},
	{"outputDir", 'd', "outputDir", 0, "Specifies the output directory for the output files", 4},
	{"tmpDir", 'T', "tmpDir", 0, "Specifies the directory in which to store temporary files", 4},
	{"outputFormat", 'O', "outputFormat", 0, "Specifies the output format 1:wig 2:bed 3:bed & wig ", 4},
	{"timing", 't', 0, OPTION_NO_USAGE, "Specifies to output timing information", 4},
	{0, 0, 0, 0, "=========== Miscellaneous Options ===================================================", 5},
	{"Parameters", 'p', 0, OPTION_NO_USAGE, "Print program parameters", 5},
	{"Help", 'h', 0, OPTION_NO_USAGE, "Display usage summary", 5},
	{0, 0, 0, 0, 0, 0}
};
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

#ifdef HAVE_ARGP_H
/*
   The ARGP structure itself.
   */
static struct argp argp = {options, parse_opt, args_doc, doc};
#else
/* argp.h support not available! Fall back to getopt */
static char OptionString[]=
"a:d:e:i:m:o:s:A:E:I:M:O:S:T:X:Y:Z:2bhptCI";
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

	if(args->inputFileName!=0) {
		fprintf(stderr, "Validating inputFileName %s. \n",
				args->inputFileName);
		if(ValidateFileName(args->inputFileName)==0)
			PrintError(FnName, "inputFileName", "Command line argument", Exit, IllegalFileName);
	}

	/* binary input not currently supported */
	if(args->binaryInput == 1) {
		PrintError("ValidateInputs",
				"binaryInput",
				"Binary input not supported",
				Exit,
				InputArguments);
	}

	if(args->algorithmReads < 0 || 
			args->algorithmReads > 1) {
		PrintError(FnName, "algorithmReads", "Command line argument", Exit, OutOfRange);
	}

	/* This should hold internally */
	assert(args->pairedEnd == 0 || args->pairedEnd == 1);

	if(args->startChr <= 0) {
		PrintError(FnName, "startChr", "Command line argument", Exit, OutOfRange);
	}

	if(args->startPos <= 0) {
		PrintError(FnName, "startPos", "Command line argument", Exit, OutOfRange);
	}

	if(args->endChr <= 0) {
		PrintError(FnName, "endChr", "Command line argument", Exit, OutOfRange);
	}

	if(args->endPos <= 0) {
		PrintError(FnName, "endPos", "Command line argument", Exit, OutOfRange);
	}

	/* This should internally hold */
	if(args->algorithmReadsPaired < 0 || 
			args->algorithmReadsPaired > 2) {
		PrintError(FnName, "algorithmReads", "Command line argument", Exit, OutOfRange);
	}

	/* Check that the min distance is less than or equal to the max distance */
	if(args->minDistancePaired > args->maxDistancePaired) {
		PrintError(FnName, "minDistancePaired > maxDistancePaired", "Command line argument", Exit, OutOfRange);
	}

	/* This must hold internally */
	assert(args->chrAbPaired == 0 || args->chrAbPaired == 1);
	assert(args->inversionsPaired == 0 || args->inversionsPaired == 1);

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

	if(args->tmpDir!=0) {
		fprintf(stderr, "Validating tmpDir path %s. \n",
				args->tmpDir);
		if(ValidateFileName(args->tmpDir)==0)
			PrintError(FnName, "tmpDir", "Command line argument", Exit, IllegalFileName);
	}

	if(args->outputFormat < 0) {
		PrintError(FnName, "outputFormat", "Command line argument", Exit, OutOfRange);
	}

	assert(args->timing == 0 || args->timing == 1);

	return 1;
}

/* TODO */
	int 
ValidateFileName(char *Name) 
{
	/* 
	   Checking that strings are good: FileName = [a-zA-Z_0-9][a-zA-Z0-9-.]+
	   FileName can start with only [a-zA-Z_0-9]
	   */

	char *ptr=Name;
	int counter=0;
	/*   fprintf(stderr, "Validating FileName %s with length %d\n", ptr, strlen(Name));  */

	assert(ptr!=0);

	while(*ptr) {
		if((isalnum(*ptr) || (*ptr=='_') || (*ptr=='+') || 
					((*ptr=='.') /* && (counter>0)*/) || /* FileNames can't start  with . or - */
					((*ptr=='/')) || /* Make sure that we can navigate through folders */
					((*ptr=='-') && (counter>0)))) {
			ptr++;
			counter++;
		}
		else return 0;
	}
	return 1;
}

/* TODO */
	void 
AssignDefaultValues(struct arguments *args)
{
	/* Assign default values */

	args->programMode = ExecuteProgram;

	args->inputFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->inputFileName!=0);
	strcpy(args->inputFileName, DEFAULT_FILENAME);

	args->binaryInput = 0;

	args->algorithmReads=0;
	args->minScoreReads=INT_MIN;

	args->startChr=0;
	args->startPos=0;
	args->endChr=0;
	args->endPos=0;

	args->pairedEnd=0;
	args->algorithmReadsPaired=0;
	args->minScoreReadsPaired=INT_MIN;
	args->minDistancePaired=INT_MIN;
	args->maxDistancePaired=INT_MAX;
	args->meanDistancePaired=0;
	args->chrAbPaired=0;
	args->inversionsPaired=0;

	args->outputID = 
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->outputID!=0);
	strcpy(args->outputID, DEFAULT_FILENAME);

	args->outputDir =
		(char*)malloc(sizeof(DEFAULT_OUTPUT_DIR));
	assert(args->outputDir!=0);
	strcpy(args->outputDir, DEFAULT_OUTPUT_DIR);

	args->tmpDir =
		(char*)malloc(sizeof(DEFAULT_OUTPUT_DIR));
	assert(args->tmpDir!=0);
	strcpy(args->tmpDir, DEFAULT_OUTPUT_DIR);

	args->outputFormat=0;

	args->timing = 0;

	return;
}

/* TODO */
	void 
PrintProgramParameters(FILE* fp, struct arguments *args)
{
	char programmode[3][64] = {"ExecuteGetOptHelp", "ExecuteProgram", "ExecutePrintProgramParameters"};
	char algorithm[2][64] = {"Unique", "Best Score"};
	fprintf(fp, BREAK_LINE);
	fprintf(fp, "Printing Program Parameters:\n");
	fprintf(fp, "programMode:\t\t%d\t[%s]\n", args->programMode, programmode[args->programMode]);
	fprintf(fp, "inputFileName:\t\t%s\n", args->inputFileName);
	fprintf(fp, "binaryInput:\t\t%d\n", args->binaryInput);
	fprintf(fp, "algorithmReads:\t\t%d\t[%s]\n", args->algorithmReads, algorithm[args->algorithmReads]);
	fprintf(fp, "minScoreReads:\t\t%d\n", args->minScoreReads);
	fprintf(fp, "startChr:\t\t%d\n", args->startChr);
	fprintf(fp, "startPos:\t\t%d\n", args->startPos);
	fprintf(fp, "endChr:\t\t\t%d\n", args->endChr);
	fprintf(fp, "endPos:\t\t\t%d\n", args->endPos);
	fprintf(fp, "pairedEnd:\t\t%d\n", args->pairedEnd);
	fprintf(fp, "algorithmReadsPaired:\t%d\t[%s]\n", args->algorithmReadsPaired, algorithm[args->algorithmReadsPaired]);
	fprintf(fp, "minScoreReadsPaired:\t%d\n", args->minScoreReadsPaired);
	fprintf(fp, "minDistancePaired:\t%d\n", args->minDistancePaired);
	fprintf(fp, "maxDistancePaired:\t%d\n", args->maxDistancePaired);
	fprintf(fp, "meanDistancePaired:\t%d\n", args->meanDistancePaired);
	fprintf(fp, "chrAbPaired:\t\t%d\n", args->chrAbPaired);
	fprintf(fp, "inversionsPaired:\t%d\n", args->inversionsPaired);
	fprintf(fp, "outputID:\t\t%s\n", args->outputID);
	fprintf(fp, "outputDir:\t%s\n", args->outputDir);
	fprintf(fp, "tmpDir:\t\t\t%s\n", args->tmpDir);
	fprintf(fp, "outputFormat:\t\t%d\n", args->outputFormat);
	fprintf(fp, "timing:\t\t\t%d\n", args->timing);
	fprintf(fp, BREAK_LINE);
	return;
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
				fprintf(stderr, "-%c\t%12s\t%s\n", (*a).key, (*a).arg, (*a).doc); break;
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
					case '2':
						arguments->pairedEnd = 1;break;
					case 'a':
						arguments->algorithmReads = atoi(OPTARG);break;
					case 'b':
						arguments->binaryInput = 1;break;
					case 'd':
						if(arguments->outputDir) free(arguments->outputDir);
						arguments->outputDir = OPTARG;
						/* set the tmp directory to the output director */
						if(strcmp(arguments->tmpDir, DEFAULT_FILENAME)==0) {
							free(arguments->tmpDir);
							arguments->tmpDir = malloc(sizeof(char)*(strlen(arguments->outputDir)+1));
							strcpy(arguments->tmpDir, arguments->outputDir);
						}
						break;
					case 'e':
						arguments->endChr=atoi(OPTARG);break;
					case 'h':
						arguments->programMode=ExecuteGetOptHelp;break;
					case 'i':
						if(arguments->inputFileName) free(arguments->inputFileName);
						arguments->inputFileName = OPTARG;break;
					case 'm':
						arguments->minScoreReads = atoi(OPTARG);break;
					case 'o':
						if(arguments->outputID) free(arguments->outputID);
						arguments->outputID = OPTARG;break;
					case 'p':
						arguments->programMode=ExecutePrintProgramParameters;break;
					case 's':
						arguments->startChr=atoi(OPTARG);break;
					case 't':
						arguments->timing = 1;break;
					case 'A':
						arguments->algorithmReadsPaired = atoi(OPTARG);break;
					case 'C':
						arguments->chrAbPaired = 1;break;
					case 'E':
						arguments->endPos=atoi(OPTARG);break;
					case 'I':
						arguments->inversionsPaired = 1;break;
					case 'M':
						arguments->minScoreReadsPaired = atoi(OPTARG);break;
					case 'O':
						arguments->outputFormat=atoi(OPTARG);break;
					case 'S':
						arguments->startPos=atoi(OPTARG);break;
					case 'T':
						if(arguments->tmpDir) free(arguments->tmpDir);
						arguments->tmpDir = OPTARG;break;
					case 'X':
						arguments->minDistancePaired = atoi(OPTARG);break;
					case 'Y':
						arguments->maxDistancePaired = atoi(OPTARG);break;
					case 'Z':
						arguments->meanDistancePaired = atoi(OPTARG);break;
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
