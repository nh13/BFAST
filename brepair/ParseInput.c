#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
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

#include <time.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/RGBinary.h"
#include "../blib/BLib.h"
#include "../balign/Definitions.h"
#include "RunRepair.h"
#include "ParseInput.h"

const char *argp_program_version =
"brepair "
PACKAGE_VERSION
"\n"
"Copyright 2008-2009.";

const char *argp_program_bug_address =
PACKAGE_BUGREPORT;

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC, OPTIONAL_GROUP_NAME}.
   */
enum { 
	DescInputFilesTitle, DescRGFileName, DescUnpairedFileName, DescScoringMatrixFileName, 
	DescAlgoTitle, DescAlignmentType, DescSpace, DescNumThreads,
	DescPairedEndOptionsTitle, DescMinPairedEndDistance, DescMaxPairedEndDistance, DescMirroringType, DescStrandedness, 
	DescOutputTitle, DescOutputID, DescOutputDir, DescTmpDir, DescTiming, 
	DescMiscTitle, DescHelp
};

/* 
   The prototype for argp_option comes fron argp.h. If argp.h
   absent, then ParseInput.h declares it 
   */
static struct argp_option options[] = {
	{0, 0, 0, 0, "=========== Input Files =============================================================", 1},
	{"rgFileName", 'r', "rgFileName", 0, "Specifies the file name of the reference genome file", 1},
	{"unpairedFileName", 'u', "unpairedFileName", 0, "Specifies the bfast unpaired file", 1},
	{"scoringMatrixFileName", 'x', "scoringMatrixFileName", 0, "Specifies the file name storing the scoring matrix", 1},
	/*
	   {"binaryInput", 'b', 0, OPTION_NO_USAGE, "Specifies that the input unpairedes files will be in binary format", 1},
	   */
	{0, 0, 0, 0, "=========== Algorithm Options: (Unless specified, default value = 0) ================", 2},
	{"alignmentType", 'a', "alignmentType", 0, "0: Full alignment 1: misunpairedes only", 2},
	{"space", 'A', "space", 0, "0: NT space 1: Color space", 2},
	{"numThreads", 'n', "numThreads", 0, "Specifies the number of threads to use (Default 1)", 2},
	{0, 0, 0, 0, "=========== Paired End Options ======================================================", 3},
	{"minPairedEndDistance", 'm', "minPairedEndDistance", 0, "Specifies the minimum distance from the current to consider, where the distance is specified according to the mirroring type", 3},
	{"maxPairedEndDistance", 'M', "maxPairedEndDistance", 0, "Specifies the maximum distance from the current to consider, where the distance is specified according to the mirroring type", 3},
	 {"mirroringType", 'L', "mirroringType", 0, "1: specifies that we assume that the first end is before the second end (5'-3')"
		 "\n\t\t\t2: specifies that we assume that the second end is before the first end (5'->3)"         
			 "\n\t\t\t3: specifies that we mirror CALs in both directions", 3},
	{"strandedness", 'S', "strandedness", 0, "0: Mirror on the same strand 1: Mirror on opposite strands 2: Mirror both strands", 3},
	{0, 0, 0, 0, "=========== Output Options ==========================================================", 4},
	{"outputID", 'o', "outputID", 0, "Specifies the name to identify the output files", 4},
	{"outputDir", 'd', "outputDir", 0, "Specifies the output directory for the output files", 4},
	{"tmpDir", 'T', "tmpDir", 0, "Specifies the directory in which to store temporary files", 4},
	/*
	   {"binaryOutput", 'B', 0, OPTION_NO_USAGE, "Specifies that the output aligned file will be in binary format", 4},
	   */
	{"timing", 't', 0, OPTION_NO_USAGE, "Specifies to output timing information", 4},
	{0, 0, 0, 0, "=========== Miscellaneous Options ===================================================", 5},
	{"Parameters", 'p', 0, OPTION_NO_USAGE, "Print program parameters", 5},
	{"Help", 'h', 0, OPTION_NO_USAGE, "Display usage summary", 5},
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
"a:d:m:n:o:r:u:x:A:L:M:S:T:hpt";
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
	time_t startTotalTime = time(NULL);
	time_t endTotalTime;
	time_t startTime;
	time_t endTime;
	int totalReferenceGenomeTime = 0; /* Total time to load and delete the reference genome */
	int totalAlignTime = 0; /* Total time to align the reads */
	int totalFileHandlingTime = 0; /* Total time partitioning and merging unpairedes and alignments respectively */
	int seconds, minutes, hours;
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
							fprintf(stderr, "**** Input arguments look good! *****\n");
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
						startTime = time(NULL);
						RGBinaryReadBinary(&rg,
								arguments.rgFileName);
						/* Unpack */
						/*
						   RGBinaryUnPack(&rg);
						   */
						endTime = time(NULL);
						totalReferenceGenomeTime = endTime - startTime;
						/* Run the aligner */
						RunRepair(&rg,
								arguments.unpairedFileName,
								arguments.scoringMatrixFileName,
								arguments.alignmentType,
								arguments.space,
								arguments.binaryInput,
								arguments.numThreads,
								arguments.minPairedEndDistance,
								arguments.maxPairedEndDistance,
								arguments.mirroringType,
								arguments.strandedness,
								arguments.outputID,
								arguments.outputDir,
								arguments.tmpDir,
								arguments.binaryOutput,
								&totalAlignTime,
								&totalFileHandlingTime);
						/* Free the Reference Genome */
						RGBinaryDelete(&rg);

						if(arguments.timing == 1) {
							/* Output loading reference genome time */
							seconds = totalReferenceGenomeTime;
							hours = seconds/3600;
							seconds -= hours*3600;
							minutes = seconds/60;
							seconds -= minutes*60;
							fprintf(stderr, "Reference Genome loading time took: %d hours, %d minutes and %d seconds.\n",
									hours,
									minutes,
									seconds
								   );

							/* Output aligning time */
							seconds = totalAlignTime;
							hours = seconds/3600;
							seconds -= hours*3600;
							minutes = seconds/60;
							seconds -= minutes*60;
							fprintf(stderr, "Align time took: %d hours, %d minutes and %d seconds.\n",
									hours,
									minutes,
									seconds
								   );

							/* Output file handling time */
							seconds = totalFileHandlingTime;
							hours = seconds/3600;
							seconds -= hours*3600;
							minutes = seconds/60;
							seconds -= minutes*60;
							fprintf(stderr, "File handling time took: %d hours, %d minutes and %d seconds.\n",
									hours,
									minutes,
									seconds
								   );

							/* Output total time */
							endTotalTime = time(NULL);
							seconds = endTotalTime - startTotalTime;
							hours = seconds/3600;
							seconds -= hours*3600;
							minutes = seconds/60;
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

	if(args->unpairedFileName!=0) {
		fprintf(stderr, "Validating unpairedFileName path %s. \n",
				args->unpairedFileName);
		if(ValidateFileName(args->unpairedFileName)==0)
			PrintError(FnName, "unpairedFileName", "Command line argument", Exit, IllegalFileName);
	}

	if(args->scoringMatrixFileName!=0) {
		fprintf(stderr, "Validating scoringMatrixFileName path %s. \n",
				args->scoringMatrixFileName);
		if(ValidateFileName(args->scoringMatrixFileName)==0)
			PrintError(FnName, "scoringMatrixFileName", "Command line argument", Exit, IllegalFileName);
	}
	
	if(args->alignmentType != FullAlignment && args->alignmentType != MismatchesOnly) {
		PrintError(FnName, "alignmentType", "Command line argument", Exit, OutOfRange);
	}

	if(args->space != NTSpace && args->space != ColorSpace) {
		PrintError(FnName, "space", "Command line argument", Exit, OutOfRange);
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

	/* If this does not hold, we have done something wrong internally */
	assert(args->timing == 0 || args->timing == 1);
	assert(args->binaryInput == TextInput || args->binaryInput == BinaryInput);
	assert(args->binaryOutput == TextOutput || args->binaryOutput == BinaryOutput);
	assert(MirrorForward <= args->mirroringType && args->mirroringType <= MirrorBoth);
	assert(StrandSame <= args->strandedness && args->strandedness <= StrandBoth);
	if(args->maxPairedEndDistance < args->minPairedEndDistance) {
		PrintError(FnName, "args->maxPairedEndDistance < args->minPairedEndDistance", "Command line argument", Exit, OutOfRange);
	}

	return 1;
}

	void 
AssignDefaultValues(struct arguments *args)
{
	/* Assign default values */

	args->programMode = ExecuteProgram;

	args->rgFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->rgFileName!=0);
	strcpy(args->rgFileName, DEFAULT_FILENAME);

	args->unpairedFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->unpairedFileName!=0);
	strcpy(args->unpairedFileName, DEFAULT_FILENAME);

	args->scoringMatrixFileName =
		(char*)malloc(sizeof(DEFAULT_FILENAME));
	assert(args->scoringMatrixFileName!=0);
	strcpy(args->scoringMatrixFileName, DEFAULT_FILENAME);
	
	args->binaryInput = BMATCHES_DEFAULT_OUTPUT;
	args->binaryOutput = BALIGN_DEFAULT_OUTPUT;

	args->alignmentType = FullAlignment;
	args->space = NTSpace;
	args->numThreads = 1;
	args->minPairedEndDistance = 0;
	args->maxPairedEndDistance = 0;
	args->mirroringType = NoMirroring;
	args->strandedness = ForwardStrand;

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

	args->timing = 0;

	return;
}

	void 
PrintProgramParameters(FILE* fp, struct arguments *args)
{
	char programmode[3][64] = {"ExecuteGetOptHelp", "ExecuteProgram", "ExecutePrintProgramParameters"};
	fprintf(fp, BREAK_LINE);
	fprintf(fp, "Printing Program Parameters:\n");
	fprintf(fp, "programMode:\t\t\t\t%d\t[%s]\n", args->programMode, programmode[args->programMode]);
	fprintf(fp, "rgFileName:\t\t\t\t%s\n", args->rgFileName);
	fprintf(fp, "unpairedFileName:\t\t\t%s\n", args->unpairedFileName);
	fprintf(fp, "scoringMatrixFileName:\t\t\t%s\n", args->scoringMatrixFileName);
	/*
	   fprintf(fp, "binaryInput:\t\t\t\t%d\n", args->binaryInput);
	   */
	fprintf(fp, "alignmentType:\t\t\t\t%d\n", args->alignmentType);
	fprintf(fp, "space:\t\t\t\t\t%d\n", args->space);
	fprintf(fp, "numThreads:\t\t\t\t%d\n", args->numThreads);
	fprintf(fp, "minPairedEndDistance:\t\t\t%d\n", args->minPairedEndDistance);
	fprintf(fp, "maxPairedEndDistance:\t\t\t%d\n", args->maxPairedEndDistance);
	fprintf(fp, "mirroringType:\t\t\t\t%d\n", args->mirroringType);
	fprintf(fp, "strandedness:\t\t\t\t%d\n", args->strandedness);
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
	free(args->unpairedFileName);
	args->unpairedFileName=NULL;
	free(args->scoringMatrixFileName);
	args->scoringMatrixFileName=NULL;
	free(args->outputID);
	args->outputID=NULL;
	free(args->outputDir);
	args->outputDir=NULL;
	free(args->tmpDir);
	args->tmpDir=NULL;
}

void
GetOptHelp() {

	struct argp_option *a=options;
	fprintf(stderr, "%s\n", argp_program_version);
	fprintf(stderr, "\nUsage: brepair [options]\n");
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
	fprintf(stderr, "\n send bugs to %s\n", argp_program_bug_address);
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
					case 'a':
						arguments->alignmentType=atoi(OPTARG);break;
						/*
						   case 'b':
						   arguments->binaryInput = 1;break;
						   */
					case 'd':
						StringCopyAndReallocate(&arguments->outputDir, OPTARG);
						/* set the tmp directory to the output director */
						if(strcmp(arguments->tmpDir, DEFAULT_FILENAME)==0) {
							StringCopyAndReallocate(&arguments->tmpDir, OPTARG);
						}
						break;
					case 'h':
						arguments->programMode=ExecuteGetOptHelp; break;
					case 'm':
						arguments->minPairedEndDistance=atoi(OPTARG);break;
					case 'n':
						arguments->numThreads=atoi(OPTARG); break;
					case 'o':
						StringCopyAndReallocate(&arguments->outputID, OPTARG);
						break;
					case 'p':
						arguments->programMode=ExecutePrintProgramParameters; break;
					case 'r':
						StringCopyAndReallocate(&arguments->rgFileName, OPTARG);
						break;
					case 't':
						arguments->timing = 1;break;
					case 'u':
						StringCopyAndReallocate(&arguments->unpairedFileName, OPTARG);
						break;
					case 'x':
						StringCopyAndReallocate(&arguments->scoringMatrixFileName, OPTARG);
						break;
					case 'A':
						arguments->space=atoi(OPTARG);break;
						/*
						   case 'B':
						   arguments->binaryOutput = 1;break;
						   */
					case 'L':
						arguments->mirroringType=atoi(OPTARG);break;
					case 'M':
						arguments->maxPairedEndDistance=atoi(OPTARG);break;
					case 'S':
						arguments->strandedness=atoi(OPTARG);break;
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
