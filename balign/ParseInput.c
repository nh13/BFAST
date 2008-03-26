#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <limits.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
#include "ParseInput.h"

const char *argp_program_version =
"balign version 0.1.1\n"
"Copyright 2007.";

const char *argp_program_bug_address =
"Nils Homer <nhomer@cs.ucla.edu>";

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC, OPTIONAL_GROUP_NAME}.
   */
enum { 
	DescInputFilesTitle, 
	DescAlgoTitle, 
	DescOutputTitle, DescOutputFileName, 
	DescMiscTitle, DescHelp
};

/* 
   The prototype for argp_option comes fron argp.h. If argp.h
   absent, then ParseInput.h declares it 
   */
static struct argp_option options[] = {
	{0, 0, 0, 0, "=========== Input Files =============================================================", 1},
	{0, 0, 0, 0, "=========== Algorithm Options: (Unless specified, default value = 1) ================", 2},
	{0, 0, 0, 0, "=========== Output Options ==========================================================", 3},
	{"OutputFileName", 'o', "OutputFileName", 0, "Specifies the output file name", 3},
	{0, 0, 0, 0, "=========== Miscellaneous Options ===================================================", 4},
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
"o:h";
#endif

enum {ExecuteGetOptHelp, ExecuteProgram};

/*
   The main function. All command-line options parsed using argp_parse
   or getopt whichever available
   */

	int
main (int argc, char **argv)
{
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
				switch(arguments.ProgramMode) {
					case ExecuteGetOptHelp:
						PrintProgramParameters(stderr, &arguments);
						break;
					case ExecuteProgram:
						if(ValidateInputs(&arguments)) {
							fprintf(stderr, "**** Input arguments look good! *****\n");
							fprintf(stderr, BREAK_LINE);
						}
						else {
							fprintf(stderr, "PrintError validating command-line inputs. Terminating!\n");
							exit(1);
						}
						PrintProgramParameters(stderr, &arguments);
						/* Execute Program */
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

int ValidateInputs(struct arguments *CommandLineArg) {

	char *FnName="ValidateInputs";

	fprintf(stderr, BREAK_LINE);
	fprintf(stderr, "Checking input parameters supplied by the user ...\n");

	if((*CommandLineArg).OutputFileName!=0) {
		fprintf(stderr, "Validating OutputFileName path %s. \n",
				(*CommandLineArg).OutputFileName);
		if(ValidateFileName((*CommandLineArg).OutputFileName)==0)
			PrintError(FnName, "OutputFileName", "Command line argument", Exit, IllegalFileName);
	}

	return 1;
}

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

	void 
AssignDefaultValues(struct arguments *args)
{
	/* Assign default values */

	(*args).ProgramMode = ExecuteProgram;

	(*args).OutputFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert((*args).OutputFileName!=0);
	strcpy((*args).OutputFileName, DEFAULT_FILENAME);

	return;
}

	void 
PrintProgramParameters(FILE* fp, struct arguments *args)
{
	char programmode[2][64] = {"ExecuteGetOptHelp", "ExecuteProgram"};
	fprintf(fp, BREAK_LINE);
	fprintf(fp, "Printing Program Parameters:\n");
	fprintf(fp, "ProgramMode:\t\t\t\t%d\t[%s]\n", (*args).ProgramMode, programmode[(*args).ProgramMode]);
	fprintf(fp, "OutputFileName:\t\t\t\t%s\n", (*args).OutputFileName);
	fprintf(fp, BREAK_LINE);
	return;
}

void
GetOptHelp() {

	struct argp_option *a=options;
	fprintf(stderr, "\nUsage: danalyze [options]\n");
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

void
PrintGetOptHelp() {

	struct argp_option *a=options;
	fprintf(stderr, "%s\n", argp_program_version);
	fprintf(stderr, "\nUsage: danalyze [options]\n");
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
					case 'h':
						arguments->ProgramMode=ExecuteGetOptHelp; break;
					case 'o':
						if(arguments->OutputFileName) free(arguments->OutputFileName);
						arguments->OutputFileName = OPTARG;break;
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
