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

#include "../blib/BLibDefinitions.h"
#include "../blib/RGBinary.h"
#include "../blib/RGIndexLayout.h"
#include "../blib/RGIndexExons.h"
#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "Definitions.h"
#include "GenerateIndex.h"
#include "ParseInput.h"

const char *argp_program_version =
"bpreprocess "
PACKAGE_VERSION
"\n"
"Copyright 2008";

const char *argp_program_bug_address =
PACKAGE_BUGREPORT;

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC, OPTIONAL_GROUP_NAME}.
   */
enum { 
	DescInputFilesTitle, DescRGFileName, DescIndexLayoutFileName,  
	DescAlgoTitle, DescAlgorithm, DescSpace, DescNumThreads, 
	DescIndexSpecificTitle, DescRepeatMasker, DescStartContig, DescStartPos, DescEndContig, DescEndPos, DescExonFileName, 
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
	/*
	   {"binaryInput", 'b', 0, OPTION_NO_USAGE, "Specifies that the reference genome will be in binary format", 1},
	   */
	{0, 0, 0, 0, "=========== Algorithm Options: (Unless specified, default value = 0) ================", 2},
	{"algorithm", 'a', "algorithm", 0, "Specifies the program mode 0: create a 4-bit file from the reference contigs 1: create an index", 2},
	{"space", 'A', "space", 0, "0: NT space 1: Color space", 2},
	{0, 0, 0, 0, "=========== Index Specific Options: (Unless specified, default value = 0) ================", 3},
	{"indexLayoutFileName", 'i', "indexLayOutFileName", 0, "Specifies the file name of the file that contains tile sizes and gaps", 1},
	{"repeatMasker", 'R', 0, OPTION_NO_USAGE, "Specifies that lower case bases will be ignored (default: off).", 3},
	{"startContig", 's', "startContig", 0, "Specifies the start contig", 3},
	{"startPos", 'S', "startPos", 0, "Specifies the end position", 3},
	{"endContig", 'e', "endContig", 0, "Specifies the end contig", 3},
	{"endPos", 'E', "endPos", 0, "Specifies the end postion", 3},
	{"exonsFileName", 'x', "exonsFileName", 0, "Specifies the file name that specifies the exon-like ranges to include in the index", 3},
	{"numThreads", 'n', "numThreads", 0, "Specifies the number of threads to use (Default 1)", 2},
	{0, 0, 0, 0, "=========== Output Options ==========================================================", 4},
	{"outputID", 'o', "outputID", 0, "Specifies the name to identify the output files", 4},
	{"outputDir", 'd', "outputDir", 0, "Specifies the output directory for the output files", 4},
	{"tmpDir", 'T', "tmpDir", 0, "Specifies the directory in which to store temporary files", 4},
	/*
	   {"binaryOuput", 'B', 0, OPTION_NO_USAGE, "Specifies that we write the files as binary files", 4},
	   */	
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
"a:d:e:i:n:o:r:s:x:A:E:S:T:hptR";
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
	RGBinary rg;
	RGIndexLayout rgLayout;
	RGIndexExons exons;
	char outputFileName[MAX_FILENAME_LENGTH]="\0";
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

						switch(arguments.algorithm) {
							case 0:
								/* Read fasta files */
								RGBinaryRead(arguments.rgFileName, 
										&rg,
										arguments.space);
								sprintf(outputFileName, "%s%s.rg.file.%s.%d.%s",
										arguments.outputDir,
										PROGRAM_NAME,
										arguments.outputID,
										rg.space,
										BFAST_RG_FILE_EXTENSION);
								/* Write binary */
								assert(arguments.binaryOutput == BinaryOutput);
								RGBinaryWriteBinary(&rg,
										outputFileName);
								break;
							case 1:
								/* Read binary */
								RGBinaryReadBinary(&rg,
										arguments.rgFileName);
								/* Read in the RGIndex layout */
								RGIndexLayoutRead(arguments.indexLayoutFileName, &rgLayout);
								/* Read exons, if necessary */
								if(arguments.useExons == UseExons) {
									RGIndexExonsRead(arguments.exonsFileName,
											&exons);
								}

								/* Generate the indexes */
								GenerateIndex(&rg,
										&rgLayout,
										arguments.space,
										arguments.startContig,
										arguments.startPos,
										arguments.endContig,
										arguments.endPos,
										arguments.useExons,
										&exons,
										arguments.repeatMasker,
										arguments.numThreads,
										arguments.outputID,
										arguments.outputDir,
										arguments.tmpDir,
										arguments.binaryOutput);

								/* Free the RGIndex layout */
								RGIndexLayoutDelete(&rgLayout);
								/* Free exons, if necessary */
								if(arguments.useExons == UseExons) {
									RGIndexExonsDelete(&exons);
								}
								else {
									/* Free exons file name if we did not use it */
									free(arguments.exonsFileName);
									arguments.exonsFileName=NULL;
								}
								break;
							default:
								break;
						}

						/* Free the Reference Genome */
						RGBinaryDelete(&rg);

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

	assert(args->binaryInput == TextInput || args->binaryInput == BinaryInput);

	if(args->algorithm < ALGORITHM_MIN || args->algorithm > ALGORITHM_MAX) {
		PrintError(FnName, "algorithm", "Command line argument", Exit, OutOfRange);
	}

	if(args->space != NTSpace && args->space != ColorSpace) {
		PrintError(FnName, "space", "Command line argument", Exit, OutOfRange);
	}

	if(args->indexLayoutFileName!=0) {
		fprintf(stderr, "Validating indexLayoutFileName %s. \n",
				args->indexLayoutFileName);
		if(ValidateFileName(args->indexLayoutFileName)==0)
			PrintError(FnName, "indexLayoutFileName", "Command line argument", Exit, IllegalFileName);
	}

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

	if(args->exonsFileName!=0) {
		fprintf(stderr, "Validating exonsFileName %s. \n",
				args->exonsFileName);
		if(ValidateFileName(args->exonsFileName)==0)
			PrintError(FnName, "exonsFileName", "Command line argument", Exit, IllegalFileName);
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
	assert(args->useExons == IgnoreExons || args->useExons == UseExons);
	assert(args->timing == 0 || args->timing == 1);
	assert(args->binaryOutput == TextOutput || args->binaryOutput == BinaryOutput);
	assert(args->repeatMasker == 0 || args->repeatMasker == 1);

	/* Cross-check arguments */
	if(args->startContig > args->endContig) {
		PrintError(FnName, "startContig > endContig", "Command line argument", Exit, OutOfRange);
	}
	if(args->startContig == args->endContig && args->startPos > args->endPos) {
		PrintError(FnName, "endPos < startPos with startContig == endContig", "Command line argument", Exit, OutOfRange);
	}
	if(args->useExons == UseExons && args->startContig > 0) {
		PrintError(FnName, "Cannot use -s with -x", "Command line argument", Exit, OutOfRange);
	}
	if(args->useExons == UseExons && args->startPos > 0) {
		PrintError(FnName, "Cannot use -S with -x", "Command line argument", Exit, OutOfRange);
	}
	if(args->useExons == UseExons && args->endContig < INT_MAX) {
		PrintError(FnName, "Cannot use -e with -x", "Command line argument", Exit, OutOfRange);
	}
	if(args->useExons == UseExons && args->endPos < INT_MAX) {
		PrintError(FnName, "Cannot use -E with -x", "Command line argument", Exit, OutOfRange);
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

	args->binaryInput = BPREPROCESS_DEFAULT_OUTPUT;

	args->algorithm = 0;
	args->space = NTSpace;

	args->indexLayoutFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->indexLayoutFileName!=0);
	strcpy(args->indexLayoutFileName, DEFAULT_FILENAME);

	args->repeatMasker=0;
	args->startContig=0;
	args->startPos=0;
	args->endContig=INT_MAX;
	args->endPos=INT_MAX;

	args->exonsFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->exonsFileName!=0);
	strcpy(args->exonsFileName, DEFAULT_FILENAME);
	args->useExons=IgnoreExons;

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

	args->binaryOutput = BPREPROCESS_DEFAULT_OUTPUT;
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
	/*
	   fprintf(fp, "binaryInput:\t\t\t\t%d\n", args->binaryInput);
	   */
	fprintf(fp, "algorithm:\t\t\t\t%d\n", args->algorithm);
	fprintf(fp, "space:\t\t\t\t\t%d\n", args->space);
	fprintf(fp, "indexLayoutFileName:\t\t\t%s\n", args->indexLayoutFileName);
	fprintf(fp, "repeatMasker:\t\t\t\t%d\n", args->repeatMasker);
	fprintf(fp, "startContig:\t\t\t\t%d\n", args->startContig);
	fprintf(fp, "startPos:\t\t\t\t%d\n", args->startPos);
	fprintf(fp, "endContig:\t\t\t\t%d\n", args->endContig);
	fprintf(fp, "endPos:\t\t\t\t\t%d\n", args->endPos);
	fprintf(fp, "exonsFileName:\t\t\t\t\%s\n", args->exonsFileName);
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
void FreeProgramParameters(struct arguments *args) 
{
	free(args->rgFileName);
	args->rgFileName=NULL;
	free(args->indexLayoutFileName);
	args->indexLayoutFileName=NULL;
	free(args->exonsFileName);
	args->exonsFileName=NULL;
	free(args->outputID);
	args->outputID=NULL;
	free(args->outputDir);
	args->outputDir=NULL;
	free(args->tmpDir);
	args->tmpDir=NULL;
}

/* TODO */
void
GetOptHelp() {

	struct argp_option *a=options;
	fprintf(stderr, "%s\n", argp_program_version);
	fprintf(stderr, "\nUsage: bpreprocess [options]\n");
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
						arguments->algorithm=atoi(OPTARG);break;
						/*
						   case 'b':
						   arguments->binaryInput=1;break;
						   */
					case 'd':
						StringCopyAndReallocate(&arguments->outputDir, OPTARG);
						/* set the tmp directory to the output director */
						if(strcmp(arguments->tmpDir, DEFAULT_FILENAME)==0) {
							StringCopyAndReallocate(&arguments->tmpDir, OPTARG);
						}
						break;
					case 'e':
						arguments->endContig=atoi(OPTARG);break;
					case 'h':
						arguments->programMode=ExecuteGetOptHelp;break;
					case 'i':
						StringCopyAndReallocate(&arguments->indexLayoutFileName, OPTARG);
						break;
					case 'n':
						arguments->numThreads=atoi(OPTARG); break;
					case 'o':
						StringCopyAndReallocate(&arguments->outputID, OPTARG);
						break;
					case 'p':
						arguments->programMode=ExecutePrintProgramParameters;break;
					case 'r':
						StringCopyAndReallocate(&arguments->rgFileName, OPTARG);
						break;
					case 's':
						arguments->startContig=atoi(OPTARG);break;
					case 't':
						arguments->timing = 1;break;
						/*
						   case 'B':
						   arguments->binaryOutput=1;break;
						   */
					case 'x':
						arguments->useExons=UseExons;
						StringCopyAndReallocate(&arguments->exonsFileName, OPTARG);
						break;
					case 'A':
						arguments->space=atoi(OPTARG);break;
					case 'E':
						arguments->endPos=atoi(OPTARG);break;
					case 'R':
						arguments->repeatMasker=1;break;
					case 'S':
						arguments->startPos=atoi(OPTARG);break;
					case 'T':
						StringCopyAndReallocate(&arguments->tmpDir, OPTARG);
						break;
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
