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
"Copyright 2007.";

const char *argp_program_bug_address =
"Nils Homer <nhomer@cs.ucla.edu>";

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC, OPTIONAL_GROUP_NAME}.
   */
enum { 
	DescInputFilesTitle, DescInputFileName, DescInputFormat,
	DescAlgoTitle, DescUniqueMatches, DescBestScore, DescMinScore, DescStartChr, DescStartPos, DescEndChr, DescEndPos,
	DescOutputTitle, DescOutputFileName, DescOutputFormat, DescTiming,
	DescMiscTitle, DescParameters, DescHelp
};

/* 
   The prototype for argp_option comes fron argp.h. If argp.h
   absent, then ParseInput.h declares it 
   */
static struct argp_option options[] = {
	{0, 0, 0, 0, "=========== Input Files =============================================================", 1},
	{"inputFileName", 'i', "inputFileName", 0, "Specifies the input file", 1},
	{"inputFormat", 'I', "inputFormat", 0, "Specifies the input format", 1},
	{0, 0, 0, 0, "=========== Algorithm Options =======================================================", 2},
	{"uniqueMatches", 'u', 0, OPTION_NO_USAGE, "Specifies to only consider unique matches", 2},
	{"bestScore", 'b', 0, OPTION_NO_USAGE, "Specifies to choose the best score from the possible matches", 2},
	{"minScore", 'm', "minScore", 0, "Specifies the minimum score to consider", 2},
	{"startChr", 's', "startChr", 0, "Specifies the start chromosome", 2},
	{"startPos", 'S', "startPos", 0, "Specifies the end position", 2},
	{"endChr", 'e', "endChr", 0, "Specifies the end chromosome", 2},
	{"endPos", 'E', "endPos", 0, "Specifies the end postion", 2},
	{0, 0, 0, 0, "=========== Output Options ==========================================================", 3},
	{"outputFileName", 'o', "outputFileName", 0, "Specifies the name to identify the output files", 3},
	{"outputFormat", 'O', "outputFormat", 0, "Specifies the output format", 3},
	{"timing", 't', 0, OPTION_NO_USAGE, "Specifies to output timing information", 3},
	{0, 0, 0, 0, "=========== Miscellaneous Options ===================================================", 4},
	{"Parameters", 'p', 0, OPTION_NO_USAGE, "Print program parameters", 4},
	{"Help", 'h', 0, OPTION_NO_USAGE, "Display usage summary", 4},
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
static char doc[] ="This program was created by Nils Homer and is not intended for distribution.";

#ifdef HAVE_ARGP_H
/*
   The ARGP structure itself.
   */
static struct argp argp = {options, parse_opt, args_doc, doc};
#else
/* argp.h support not available! Fall back to getopt */
static char OptionString[]=
"d:e:i:m:o:s:E:I:O:S:bhptu";
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

	if(args->inputFormat < 0) {
		PrintError(FnName, "inputFormat", "Command line argument", Exit, OutOfRange);
	}

	/* This should hold internally */
	assert(args->uniqueMatches == 0 || args->uniqueMatches == 1);
	assert(args->bestScore == 0 || args->bestScore == 1);

	/* Check that we do not use unique matches and best score together */
	if(args->uniqueMatches == 1 && args->bestScore == 1) {
		PrintError(FnName, "uniqueMatches and bestScore", "Command line argument: cannot use both command line arguments", Exit, OutOfRange);
	}

	/* Check that at either unique matches or best score is used */
	if(args->uniqueMatches == 0 && args->bestScore == 0) {
		PrintError(FnName, "uniqueMatches and bestScore", "Command line argument: must use one of the two command line arguments", Exit, OutOfRange);
	}

	if(args->startChr < 0) {
		PrintError(FnName, "startChr", "Command line argument", Exit, OutOfRange);
	}

	if(args->startPos < 0) {
		PrintError(FnName, "startPos", "Command line argument", Exit, OutOfRange);
	}

	if(args->endChr < 0) {
		PrintError(FnName, "endChr", "Command line argument", Exit, OutOfRange);
	}

	if(args->endPos < 0) {
		PrintError(FnName, "endPos", "Command line argument", Exit, OutOfRange);
	}

	if(args->outputFileName!=0) {
		fprintf(stderr, "Validating outputFileName %s. \n",
				args->outputFileName);
		if(ValidateFileName(args->outputFileName)==0)
			PrintError(FnName, "outputFileName", "Command line argument", Exit, IllegalFileName);
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

	args->inputFormat=0;

	args->uniqueMatches=0;
	args->bestScore=0;
	args->minScore=INT_MIN;

	args->startChr=0;
	args->startPos=0;
	args->endChr=0;
	args->endPos=0;

	args->outputFileName = 
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->outputFileName!=0);
	strcpy(args->outputFileName, DEFAULT_FILENAME);

	args->outputFormat=0;

	args->timing = 0;

	return;
}

/* TODO */
	void 
PrintProgramParameters(FILE* fp, struct arguments *args)
{
	char programmode[3][64] = {"ExecuteGetOptHelp", "ExecuteProgram", "ExecutePrintProgramParameters"};
	fprintf(fp, BREAK_LINE);
	fprintf(fp, "Printing Program Parameters:\n");
	fprintf(fp, "programMode:\t\t\t\t%d\t[%s]\n", args->programMode, programmode[args->programMode]);
	fprintf(fp, "inputFileName:\t\t\t\t%s\n", args->inputFileName);
	fprintf(fp, "inputFormat:\t\t\t\t%d\n", args->inputFormat);
	fprintf(fp, "uniqueMatches:\t\t\t\t%d\n", args->uniqueMatches);
	fprintf(fp, "bestScore:\t\t\t\t%d\n", args->bestScore);
	fprintf(fp, "minScore:\t\t\t\t%d\n", args->minScore);
	fprintf(fp, "startChr:\t\t\t\t%d\n", args->startChr);
	fprintf(fp, "startPos:\t\t\t\t%d\n", args->startPos);
	fprintf(fp, "endChr:\t\t\t\t\t%d\n", args->endChr);
	fprintf(fp, "endPos:\t\t\t\t\t%d\n", args->endPos);
	fprintf(fp, "outputFileName:\t\t\t\t%s\n", args->outputFileName);
	fprintf(fp, "outputFormat:\t\t\t\t%d\n", args->outputFormat);
	fprintf(fp, "timing:\t\t\t\t\t%d\n", args->timing);
	fprintf(fp, BREAK_LINE);
	return;
}

/* TODO */
void
GetOptHelp() {

	struct argp_option *a=options;
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
	return;
}

/* TODO */
void
PrintGetOptHelp() {

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
	fprintf(stderr, "\n%s\n", argp_program_bug_address);
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
					case 'b':
						arguments->bestScore = 1;break;
					case 'e':
						arguments->endChr=atoi(OPTARG);break;
					case 'h':
						arguments->programMode=ExecuteGetOptHelp;break;
					case 'i':
						if(arguments->inputFileName) free(arguments->inputFileName);
						arguments->inputFileName = OPTARG;break;
					case 'm':
						arguments->minScore = atoi(OPTARG);break;
					case 'o':
						if(arguments->outputFileName) free(arguments->outputFileName);
						arguments->outputFileName = OPTARG;break;
					case 'p':
						arguments->programMode=ExecutePrintProgramParameters;break;
					case 's':
						arguments->startChr=atoi(OPTARG);break;
					case 't':
						arguments->timing = 1;break;
					case 'u':
						arguments->uniqueMatches = 1;break;
					case 'E':
						arguments->endPos=atoi(OPTARG);break;
					case 'I':
						arguments->inputFormat=atoi(OPTARG);break;
					case 'O':
						arguments->outputFormat=atoi(OPTARG);break;
					case 'S':
						arguments->startPos=atoi(OPTARG);break;
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
