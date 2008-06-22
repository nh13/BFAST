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

#include <time.h>

#include "Definitions.h"
#include "../blib/BError.h"
#include "FindMatches.h"
#include "ParseInput.h"

const char *argp_program_version =
"bmatches version 0.1.1\n"
"Copyright 2008.";

const char *argp_program_bug_address =
"Nils Homer <nhomer@cs.ucla.edu>";

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC, OPTIONAL_GROUP_NAME}.
   */
enum { 
	DescInputFilesTitle, DescRGFileName, DescBfastMainIndexesFileName, DescBfastSecondaryIndexesFileName, DescReadsFileName, DescOffsetsFileName, 
	DescAlgoTitle, DescStartReadNum, DescEndReadNum, DescNumMismatches, DescNumInsertions, DescNumDeletions, DescNumGapInsertions, DescNumGapDeletions, DescPairedEnd, /*DescMaxMatches, */DescNumThreads, 
	DescOutputTitle, DescOutputID, DescOutputDir, DescTmpDir, DescTiming,
	DescMiscTitle, DescParameters, DescHelp
};

/* 
   The prototype for argp_option comes fron argp.h. If argp.h
   absent, then ParseInput.h declares it 
   */
static struct argp_option options[] = {
	{0, 0, 0, 0, "=========== Input Files =============================================================", 1},
	{"rgFileName", 'r', "rgFileName", 0, "Specifies the file name of the reference genome file", 1},
	{"bfastMainIndexesFileName", 'i', "bfastMainIndexesFileName", 0, "Specifies the file name holding the list of main bif files", 1},
	{"bfastSecondaryIndexesFileName", 'I', "bfastSecondaryIndexesFileName", 0, "Specifies the file name holding the list of bif files", 1},
	{"readsFileName", 'R', "readsFileName", 0, "Specifies the file name for the reads", 1}, 
	{"offsetsFileName", 'O', "offsetsFileName", 0, "Specifies the offsets", 1},
	/*
	   {"binaryInput", 'b', 0, OPTION_NO_USAGE, "Specifies that the bfast input files will be in binary format", 1},
	   */
	{0, 0, 0, 0, "=========== Algorithm Options: (Unless specified, default value = 0) ================", 2},
	{"startReadNum", 's', "startReadNum", 0, "Specifies the read to begin with (skip the first startReadNum-1 lines)", 2},
	{"endReadNum", 'e', "endReadNum", 0, "Specifies the last read to use (inclusive)", 2},
	{"numMismatches", 'x', "numMismatches", 0, "Specifies the number of mismatches to allow when searching for candidates", 2},
	{"numInsertions", 'y', "numInsertions", 0, "Specifies the number of insertions to allow when searching for candidates", 2},
	{"numDeletions", 'z', "numDeletions", 0, "Specifies the number of deletions to allow when searching for candidates", 2},
	{"numGapInsertions", 'Y', "numGapInsertions", 0, "Specifies the number of insertions allowed in the gap between pairs", 2},
	{"numGapDeletions", 'Z', "numGapDeletions", 0, "Specifies the number of gap deletions allowd in the gap between paris", 2},
	{"pairedEnd", '2', 0, OPTION_NO_USAGE, "Specifies that paired end data is to be expected", 2},
	/*
	   {"maxMatches", 'm', "maxMatches", 0, "Specifies the maximum number of matches to consider", 2},
	   */
	{"numThreads", 'n', "numThreads", 0, "Specifies the number of threads to use (Default 1", 2},
	{0, 0, 0, 0, "=========== Output Options ==========================================================", 3},
	{"outputID", 'o', "outputID", 0, "Specifies the name to identify the output files", 3},
	{"outputDir", 'd', "outputDir", 0, "Specifies the output directory for the output files", 3},
	{"tmpDir", 'T', "tmpDir", 0, "Specifies the directory in which to store temporary files", 3},
	/*
	   {"binaryOutput", 'B', 0, OPTION_NO_USAGE, "Specifies that the output should be in binary format", 3},
	   */
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
static char doc[] = "";

#ifdef HAVE_ARGP_H
/*
   The ARGP structure itself.
   */
static struct argp argp = {options, parse_opt, args_doc, doc};
#else
/* argp.h support not available! Fall back to getopt */
static char OptionString[]=
"d:e:i:m:n:o:r:s:x:y:z:I:O:R:T:Y:Z:2hpt";
#endif

enum {ExecuteGetOptHelp, ExecuteProgram, ExecutePrintProgramParameters};

/*
   The main function. All command-line options parsed using argp_parse
   or getopt whichever available
   */
	int
main (int argc, char **argv)
{
	char outputFileName[MAX_FILENAME_LENGTH]="\0";

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
							fprintf(stderr, "**** Input arguments look good!\n");
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
						/* Execute Program */

						/* Create output file name */
						sprintf(outputFileName, "%s%s.matches.file.%s.%d.%d.%d.%d.%d.%d.%d.%d.%s",
								arguments.outputDir,
								PROGRAM_NAME,
								arguments.outputID,
								arguments.startReadNum,
								arguments.endReadNum,
								arguments.numMismatches,
								arguments.numInsertions,
								arguments.numDeletions,
								arguments.numGapInsertions,
								arguments.numGapDeletions,
								arguments.pairedEnd,
								BFAST_MATCHES_FILE_EXTENSION);

						/* Run Matches */
						FindMatches(outputFileName,
								arguments.binaryOutput,
								arguments.rgFileName,
								arguments.bfastMainIndexesFileName,
								arguments.bfastSecondaryIndexesFileName,
								arguments.readsFileName,
								arguments.offsetsFileName,
								arguments.binaryInput,
								arguments.startReadNum,
								arguments.endReadNum,
								arguments.numMismatches,
								arguments.numInsertions,
								arguments.numDeletions,
								arguments.numGapInsertions,
								arguments.numGapDeletions,
								arguments.pairedEnd,
								arguments.maxMatches,
								arguments.numThreads,
								arguments.tmpDir,
								arguments.timing);

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

	if(args->rgFileName!=0) {
		fprintf(stderr, "Validating rgFileName %s. \n",
				args->rgFileName);
		if(ValidateFileName(args->rgFileName)==0)
			PrintError(FnName, "rgFileName", "Command line argument", Exit, IllegalFileName);
	}

	if(args->bfastMainIndexesFileName!=0) {
		fprintf(stderr, "Validating bfastMainIndexesFileName %s. \n",
				args->bfastMainIndexesFileName);
		if(ValidateFileName(args->bfastMainIndexesFileName)==0)
			PrintError(FnName, "bfastMainIndexesFileName", "Command line argument", Exit, IllegalFileName);
	}

	if(args->bfastSecondaryIndexesFileName!=0) {
		fprintf(stderr, "Validating bfastSecondaryIndexesFileName %s. \n",
				args->bfastSecondaryIndexesFileName);
		if(ValidateFileName(args->bfastSecondaryIndexesFileName)==0)
			PrintError(FnName, "bfastSecondaryIndexesFileName", "Command line argument", Exit, IllegalFileName);
	}

	if(args->readsFileName!=0) {
		fprintf(stderr, "Validating readsFileName %s. \n",
				args->readsFileName);
		if(ValidateFileName(args->readsFileName)==0)
			PrintError(FnName, "readsFileName", "Command line argument", Exit, IllegalFileName);
	}

	if(args->offsetsFileName!=0) {
		fprintf(stderr, "Validating offsetsFileName %s. \n",
				args->offsetsFileName);
		if(ValidateFileName(args->offsetsFileName)==0)
			PrintError(FnName, "offsetsFileName", "Command line argument", Exit, IllegalFileName);
	}

	assert(args->binaryOutput == 0 || args->binaryOutput== 1);

	if(args->numMismatches < 0) {
		PrintError(FnName, "numMismatches", "Command line argument", Exit, OutOfRange);
	}

	if(args->numInsertions < 0) {
		PrintError(FnName, "numInsertions", "Command line argument", Exit, OutOfRange);
	}

	if(args->numDeletions < 0) {
		PrintError(FnName, "numDeletions", "Command line argument", Exit, OutOfRange);
	}

	if(args->numGapInsertions < 0) {
		PrintError(FnName, "numGapInsertions", "Command line argument", Exit, OutOfRange);
	}

	if(args->numGapDeletions < 0) {
		PrintError(FnName, "numGapDeletions", "Command line argument", Exit, OutOfRange);
	}

	if(args->pairedEnd < 0 || args->pairedEnd > 1) {
		PrintError(FnName, "pairedEnd", "Command line argument", Exit, OutOfRange);
	}

	if(args->maxMatches <= 0) {
		PrintError(FnName, "maxMatches", "Command line argument", Exit, OutOfRange);
	}

	if(args->numThreads<=0) {
		PrintError(FnName, "numThreads", "Command line argument", Exit, OutOfRange);
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

	if(args->tmpDir!=0) {
		fprintf(stderr, "Validating tmpDir path %s. \n",
				args->tmpDir);
		if(ValidateFileName(args->tmpDir)==0)
			PrintError(FnName, "tmpDir", "Command line argument", Exit, IllegalFileName);
	}

	/* If this does not hold, we have done something wrong internally */
	assert(args->timing == 0 || args->timing == 1);
	assert(args->binaryOutput == 0 || args->binaryOutput == 1);

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

	args->rgFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->rgFileName!=0);
	strcpy(args->rgFileName, DEFAULT_FILENAME);

	args->bfastMainIndexesFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->bfastMainIndexesFileName!=0);
	strcpy(args->bfastMainIndexesFileName, DEFAULT_FILENAME);

	args->bfastSecondaryIndexesFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->bfastSecondaryIndexesFileName!=0);
	strcpy(args->bfastSecondaryIndexesFileName, DEFAULT_FILENAME);

	args->readsFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->readsFileName!=0);
	strcpy(args->readsFileName, DEFAULT_FILENAME);

	args->offsetsFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->offsetsFileName!=0);
	strcpy(args->offsetsFileName, DEFAULT_FILENAME);

	args->binaryInput = 1;

	args->startReadNum = -1;
	args->endReadNum = -1;
	args->numMismatches = 0;
	args->numInsertions = 0;
	args->numDeletions = 0;
	args->numGapInsertions = 0;
	args->numGapDeletions = 0;
	args->pairedEnd = 0;
	args->maxMatches = INT_MAX;
	args->numThreads = 1;

	args->outputID =
		(char*)malloc(sizeof(DEFAULT_OUTPUT_ID));
	assert(args->outputID!=0);
	strcpy(args->outputID, DEFAULT_OUTPUT_ID);

	args->outputDir = 
		(char*)malloc(sizeof(DEFAULT_OUTPUT_DIR));
	assert(args->outputDir!=0);
	strcpy(args->outputDir, DEFAULT_OUTPUT_DIR);

	args->tmpDir =
		(char*)malloc(sizeof(DEFAULT_OUTPUT_DIR));
	assert(args->tmpDir!=0);
	strcpy(args->tmpDir, DEFAULT_OUTPUT_DIR);

	args->binaryOutput = 1;

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
	fprintf(fp, "rgFileName:\t\t\t\t%s\n", args->rgFileName);
	fprintf(fp, "bfastMainIndexesFileName\t\t%s\n", args->bfastMainIndexesFileName);
	fprintf(fp, "bfastSecondaryIndexesFileName\t\t\%s\n", args->bfastSecondaryIndexesFileName);
	fprintf(fp, "readsFileName:\t\t\t\t%s\n", args->readsFileName);
	fprintf(fp, "offsetsFileName:\t\t\t%s\n", args->offsetsFileName);
	/*
	   fprintf(fp, "binaryInput:\t\t\t\t%d\n", args->binaryInput);
	   */
	fprintf(fp, "startReadNum:\t\t\t\t%d\n", args->startReadNum);
	fprintf(fp, "endReadNum:\t\t\t\t%d\n", args->endReadNum);
	fprintf(fp, "numMismatches:\t\t\t\t%d\n", args->numMismatches);
	fprintf(fp, "numInsertions:\t\t\t\t%d\n", args->numInsertions);
	fprintf(fp, "numDeletions:\t\t\t\t%d\n", args->numDeletions);
	fprintf(fp, "numGapInsertions:\t\t\t%d\n", args->numGapInsertions);
	fprintf(fp, "numGapDeletions:\t\t\t%d\n", args->numGapDeletions);
	fprintf(fp, "pairedEnd:\t\t\t\t%d\n", args->pairedEnd);
	/*
	   fprintf(fp, "maxMatches:\t\t\t\t%d\n", args->maxMatches);
	   */
	fprintf(fp, "numThreads:\t\t\t\t%d\n", args->numThreads);
	fprintf(fp, "outputID:\t\t\t\t%s\n", args->outputID);
	fprintf(fp, "outputDir:\t\t\t\t%s\n", args->outputDir);
	fprintf(fp, "tmpDir:\t\t\t\t\t%s\n", args->tmpDir);
	/*
	   fprintf(fp, "binaryOutput:\t\t\t\t%d\n", args->binaryOutput);
	   */
	fprintf(fp, "timing:\t\t\t\t\t%d\n", args->timing);
	fprintf(fp, BREAK_LINE);
	return;
}

/* TODO */
void
GetOptHelp() {

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
						/*
						   case 'b':
						   arguments->binaryInput = 1;break;
						   */
					case 'd':
						if(arguments->outputDir) free(arguments->outputDir);
						arguments->outputDir = OPTARG;
						/* set the tmp directory to the output director */
						if(strcmp(arguments->tmpDir, DEFAULT_FILENAME)==0) {
							free(arguments->tmpDir);
							arguments->tmpDir = OPTARG;
						} 
						break;
					case 'e':
						arguments->endReadNum = atoi(OPTARG);break;
					case 'h':
						arguments->programMode=ExecuteGetOptHelp; break;
					case 'i':
						if(arguments->bfastMainIndexesFileName) free(arguments->bfastMainIndexesFileName);
						arguments->bfastMainIndexesFileName = OPTARG;break;
					case 'n':
						arguments->numThreads=atoi(OPTARG);break;
						/*
						   case 'm':
						   arguments->maxMatches=atoi(OPTARG);break;
						   */
					case 'o':
						if(arguments->outputID) free(arguments->outputID);
						arguments->outputID = OPTARG;break;
					case 'p':
						arguments->programMode=ExecutePrintProgramParameters;break;
					case 'r':
						if(arguments->rgFileName) free(arguments->rgFileName);
						arguments->rgFileName = OPTARG;break;
					case 's':
						arguments->startReadNum = atoi(OPTARG);break;
					case 't':
						arguments->timing = 1;break;
					case 'x':
						arguments->numMismatches=atoi(OPTARG);break;
					case 'y':
						arguments->numInsertions=atoi(OPTARG);break;
					case 'z':
						arguments->numDeletions = atoi(OPTARG);break;
						/*
						   case 'B':
						   arguments->binaryOutput = 1;break;
						   */
					case 'I':
						if(arguments->bfastSecondaryIndexesFileName) free(arguments->bfastSecondaryIndexesFileName);
						arguments->bfastSecondaryIndexesFileName = OPTARG;break;
					case 'O':
						if(arguments->offsetsFileName) free(arguments->offsetsFileName);
						arguments->offsetsFileName = OPTARG;break;
					case 'R':
						if(arguments->readsFileName) free(arguments->readsFileName);
						arguments->readsFileName = OPTARG;break;
					case 'T':
						if(arguments->tmpDir == arguments->outputDir) {
						}
						if(arguments->tmpDir) { 
							free(arguments->tmpDir);
						}
						arguments->tmpDir = OPTARG;
						break;
					case 'Y':
						arguments->numGapInsertions = atoi(OPTARG);break;
					case 'Z':
						arguments->numGapDeletions = atoi(OPTARG);break;
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
