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

#include "Definitions.h"
#include "FindMatches.h"
#include "ParseInput.h"

const char *argp_program_version =
"bmatches version 0.1.1\n"
"Copyright 2007.";

const char *argp_program_bug_address =
"Nils Homer <nhomer@cs.ucla.edu>";

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC, OPTIONAL_GROUP_NAME}.
   */
enum { 
	DescInputFilesTitle, DescBlatterTreesFileName, DescReadsFileName, 
	DescAlgoTitle, DescStartReadNum, DescEndReadNum, DescNumMismatches, DescNumInsertions, DescNumDeletions, DescPairedEnd,
	DescOutputTitle, DescOutputID, DescOutputDir,
	DescMiscTitle, DescParameters, DescHelp
};

/* 
   The prototype for argp_option comes fron argp.h. If argp.h
   absent, then ParseInput.h declares it 
   */
static struct argp_option options[] = {
	{0, 0, 0, 0, "=========== Input Files =============================================================", 1},
	{"blatterTreesFileName", 'b', "blatterTreesFileName", 0, "Specifies the file name holding the list of btf files", 1},
	{"readsFileName", 'r', "readsFileName", 0, "Specifies the file name for the reads", 1}, 
	{0, 0, 0, 0, "=========== Algorithm Options: (Unless specified, default value = 1) ================", 2},
	{"startReadNum", 's', "startReadNum", 0, "Specifies the read to begin with (skip the first startReadNum-1 lines)", 2},
	{"endReadNum", 'e', "endReadNum", 0, "Specifies the last read to use (inclusive)", 2},
	{"numMismatches", 'm', "numMistmatches", 0, "Specifies the number of mismatches to allow when searching for candidates", 2},
	{"numInsertions", 'i', "numInsertions", 0, "Specifies the number of insertions to allow when searching for candidates", 2},
	{"numDeletions", 'a', "numDeletions", 0, "Specifies the number of deletions to allow when searching for candidates", 2},
	{"pairedEnd", '2', "pairedEnd", 0, "Specifies that paired end data is to be expected", 2},
	{0, 0, 0, 0, "=========== Output Options ==========================================================", 3},
	{"outputID", 'o', "outputID", 0, "Specifies the name to identify the output files", 3},
	{"outputDir", 'd', "outputDir", 0, "Specifies the output directory for the output files", 3},
	{0, 0, 0, 0, "=========== Miscellaneous Options ===================================================", 4},
	{"Parameters", 'p', "Parameters", OPTION_NO_USAGE, "Print program parameters", 4},
	{"Help", 'h', "Help", OPTION_NO_USAGE, "Display usage summary", 4},
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
"a:b:d:e:i:m:o:r:s:2hp";
#endif

enum {ExecuteGetOptHelp, ExecuteProgram, ExecutePrintProgramParameters};

/*
   The main function. All command-line options parsed using argp_parse
   or getopt whichever available
   */
	int
main (int argc, char **argv)
{
	char outputFileName[MAX_FILENAME_LENGTH];

	struct arguments arguments;
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
							fprintf(stderr, "**** Input arguments look good!\n");
							fprintf(stderr, BREAK_LINE);
						}
						else {
							fprintf(stderr, "PrintError validating command-line inputs. Terminating!\n");
							exit(1);
						}
						PrintProgramParameters(stderr, &arguments);
						/* Execute Program */

						/* Create output file name */
						sprintf(outputFileName, "%sblatter.matches.file.%s.%d.%d.%d.%d.%d.%d.%s",
								arguments.outputDir,
								arguments.outputID,
								arguments.startReadNum,
								arguments.endReadNum,
								arguments.numMismatches,
								arguments.numInsertions,
								arguments.numDeletions,
								arguments.pairedEnd,
								BLATTER_MATCHES_FILE_EXTENSION);

						/* Run Matches */
						RunMatches(outputFileName,
								arguments.readsFileName,
								arguments.blatterTreesFileName,
								0,
								1,
								arguments.startReadNum,
								arguments.endReadNum,
								arguments.numMismatches,
								arguments.numInsertions,
								arguments.numDeletions,
								arguments.pairedEnd);

						fprintf(stderr, "Terminating successfully!\n");
						break;
					default:
						fprintf(stderr, "PrintError determining program mode. Terminating!\n");
						exit(1);
				}
			}
			else {
				fprintf(stderr, "PrintError parsing command line arguments!\n");
				exit(1);
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

	if(args->readsFileName!=0) {
		fprintf(stderr, "Validating readsFileName %s. \n",
				args->readsFileName);
		if(ValidateFileName(args->readsFileName)==0)
			PrintError(FnName, "readsFileName", "Command line argument", Exit, IllegalFileName);
	}

	if(args->numMismatches <= 0) {
		PrintError(FnName, "numMismatches", "Command line argument", Exit, OutOfRange);
	}

	if(args->numInsertions <= 0) {
		PrintError(FnName, "numInsertions", "Command line argument", Exit, OutOfRange);
	}

	if(args->numDeletions <= 0) {
		PrintError(FnName, "numDeletions", "Command line argument", Exit, OutOfRange);
	}

	if(args->outputID!=0) {
		fprintf(stderr, "Validating outputID %s. \n",
				args->outputID);
		if(ValidateFileName(args->outputID)==0)
			PrintError(FnName, "outputID", "Command line argument", Exit, IllegalFileName);
	}

	if(args->outputDir!=0) {
		fprintf(stderr, "Validating outputDir %s. \n", 
				args->outputDir);
		if(ValidateFileName(args->outputDir)==0) 
			PrintError(FnName, "outputDir", "Command line argument", Exit, IllegalFileName);
	}

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

	args->readsFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->readsFileName!=0);
	strcpy(args->readsFileName, DEFAULT_FILENAME);

	args->startReadNum = -1;
	args->endReadNum = -1;
	args->numMismatches = 0;
	args->numInsertions = 0;
	args->numDeletions = 0;
	args->pairedEnd = 0;

	args->outputID =
		(char*)malloc(sizeof(DEFAULT_OUTPUT_ID));
	assert(args->outputID!=0);
	strcpy(args->outputID, DEFAULT_OUTPUT_ID);

	args->outputDir = 
		(char*)malloc(sizeof(DEFAULT_OUTPUT_DIR));
	assert(args->outputDir!=0);
	strcpy(args->outputDir, DEFAULT_OUTPUT_DIR);

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
	fprintf(fp, "blatterTreesFileName\t\t\t%s\n", args->blatterTreesFileName);
	fprintf(fp, "readsFileName:\t\t\t\t%s\n", args->readsFileName);
	fprintf(fp, "startReadNum:\t\t\t\t%d\n", args->startReadNum);
	fprintf(fp, "endReadNum:\t\t\t\t%d\n", args->endReadNum);
	fprintf(fp, "numMismatches:\t\t\t\t%d\n", args->numMismatches);
	fprintf(fp, "numInsertions:\t\t\t\t%d\n", args->numInsertions);
	fprintf(fp, "numDeletions:\t\t\t\t%d\n", args->numDeletions);
	fprintf(fp, "pairedEnd:\t\t\t\t%d\n", args->pairedEnd);
	fprintf(fp, "outputID:\t\t\t\t%s\n", args->outputID);
	fprintf(fp, "outputDir:\t\t\t\t%s\n", args->outputDir);
	fprintf(fp, BREAK_LINE);
	return;
}

/* TODO */
void
GetOptHelp() {

	struct argp_option *a=options;
	fprintf(stderr, "\nUsage: bmatches [options]\n");
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
	fprintf(stderr, "\nUsage: bmatches [options]\n");
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
					case '2':
						arguments->pairedEnd = 1;break;
					case 'a':
						arguments->numDeletions = atoi(OPTARG);break;
					case 'b':
						if(arguments->blatterTreesFileName) free(arguments->blatterTreesFileName);
						arguments->blatterTreesFileName = OPTARG;break;
					case 'd':
						if(arguments->outputDir) free(arguments->outputDir);
						arguments->outputDir = OPTARG;break;
					case 'e':
						arguments->endReadNum = atoi(OPTARG);break;
					case 'h':
						arguments->programMode=ExecuteGetOptHelp; break;
					case 'i':
						arguments->numInsertions=atoi(OPTARG);break;
					case 'm':
						arguments->numMismatches=atoi(OPTARG);break;
					case 'r':
						if(arguments->readsFileName) free(arguments->readsFileName);
						arguments->readsFileName = OPTARG;break;
					case 's':
						arguments->startReadNum = atoi(OPTARG);break;
					case 'o':
						if(arguments->outputID) free(arguments->outputID);
						arguments->outputID = OPTARG;break;
					case 'p':
						arguments->programMode=ExecutePrintProgramParameters;break;
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
