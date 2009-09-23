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
#include "BLib.h"
#include "FindMatches.h"
#include "BfastMatch.h"

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC, OPTIONAL_GROUP_NAME}.
   */
enum { 
	DescInputFilesTitle, DescFastaFileName, DescMainIndexes, DescSecondaryIndexes, DescReadsFileName, DescOffsets, 
	DescAlgoTitle, DescSpace, DescStartReadNum, DescEndReadNum, 
	DescKeySize, DescMaxKeyMatches, DescMaxTotalMatches, DescWhichStrand, DescNumThreads, DescQueueLength, 
	DescOutputTitle, DescTmpDir, DescTiming,
	DescMiscTitle, DescParameters, DescHelp
};

static struct argp_option options[] = {
	{0, 0, 0, 0, "=========== Input Files =============================================================", 1},
	{"fastaFileName", 'f', "fastaFileName", 0, "Specifies the file name of the FASTA reference genome", 1},
	{"mainIndexes", 'i', "mainIndexes", 0, "Specifies the index numbers for the main bif files (comma separated)", 1},
	{"secondaryIndexes", 'I', "secondaryIndexes", 0, "Specifies the index numbers for the secondary bif files (comma separated)", 1},
	{"readsFileName", 'r', "readsFileName", 0, "Specifies the file name for the reads", 1}, 
	{"offsets", 'o', "offsets", 0, "Specifies the offsets", 1},
	{0, 0, 0, 0, "=========== Algorithm Options: (Unless specified, default value = 0) ================", 2},
	{"space", 'A', "space", 0, "0: NT space 1: Color space", 2},
	{"startReadNum", 's', "startReadNum", 0, "Specifies the read to begin with (skip the first startReadNum-1 reads)", 2},
	{"endReadNum", 'e', "endReadNum", 0, "Specifies the last read to use (inclusive)", 2},
	{"keySize", 'k', "keySize", 0, "Specifies to truncate all indexes to have the given key size (must be greater than the hash width)", 2},
	{"maxKeyMatches", 'K', "maxKeyMatches", 0, "Specifies the maximum number of matches to allow before a key is ignored", 2},
	{"maxNumMatches", 'M', "maxNumMatches", 0, "Specifies the maximum total number of matches to consider before the read is discarded", 2},
	{"whichStrand", 'w', "whichStrand", 0, "0: consider both strands 1: forward strand only 2: reverse strand only", 2},
	{"numThreads", 'n', "numThreads", 0, "Specifies the number of threads to use (Default 1)", 2},
	{"queueLength", 'Q', "queueLength", 0, "Specifies the number of reads to cache", 2},
	{0, 0, 0, 0, "=========== Output Options ==========================================================", 3},
	{"tmpDir", 'T', "tmpDir", 0, "Specifies the directory in which to store temporary files", 3},
	{"timing", 't', 0, OPTION_NO_USAGE, "Specifies to output timing information", 3},
	{0, 0, 0, 0, "=========== Miscellaneous Options ===================================================", 4},
	{"Parameters", 'p', 0, OPTION_NO_USAGE, "Print program parameters", 4},
	{"Help", 'h', 0, OPTION_NO_USAGE, "Display usage summary", 4},
	{0, 0, 0, 0, 0, 0}
};

static char OptionString[]=
"e:f:i:k:m:n:o:r:s:w:A:I:K:M:Q:T:hpt";

	int
BfastMatch(int argc, char **argv)
{
	struct arguments arguments;
	time_t startTime = time(NULL);
	time_t endTime;

	if(argc>1) {
		/* Set argument defaults. (overriden if user specifies them)  */ 
		BfastMatchAssignDefaultValues(&arguments);

		/* Parse command line args */
		if(BfastMatchGetOptParse(argc, argv, OptionString, &arguments)==0)
		{
			switch(arguments.programMode) {
				case ExecuteGetOptHelp:
					BfastMatchGetOptHelp();
					break;
				case ExecutePrintProgramParameters:
					BfastMatchPrintProgramParameters(stderr, &arguments);
					break;
				case ExecuteProgram:
					if(BfastMatchValidateInputs(&arguments)) {
						fprintf(stderr, "**** Input arguments look good!\n");
						fprintf(stderr, BREAK_LINE);
					}
					else {
						PrintError("PrintError", NULL, "validating command-line inputs", Exit, InputArguments);

					}
					BfastMatchPrintProgramParameters(stderr, &arguments);
					/* Execute Program */

					/* Run Matches */
					FindMatches(
							arguments.fastaFileName,
							arguments.mainIndexes,
							arguments.secondaryIndexes,
							arguments.readsFileName,
							arguments.offsets,
							arguments.space,
							arguments.startReadNum,
							arguments.endReadNum,
							arguments.keySize,
							arguments.maxKeyMatches,
							arguments.maxNumMatches,
							arguments.whichStrand,
							arguments.numThreads,
							arguments.queueLength,
							arguments.tmpDir,
							arguments.timing);

					if(arguments.timing == 1) {
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
					PrintError("PrintError", "programMode", "Could not determine program mode", Exit, OutOfRange);
			}
		}
		else {
			PrintError("PrintError", NULL, "Could not parse command line arguments", Exit, InputArguments);
		}
		/* Free program parameters */
		BfastMatchFreeProgramParameters(&arguments);
	}
	else {
		BfastMatchGetOptHelp();
	}

	return 0;
}

/* TODO */
int BfastMatchValidateInputs(struct arguments *args) {

	char *FnName="BfastMatchValidateInputs";

	fprintf(stderr, BREAK_LINE);
	fprintf(stderr, "Checking input parameters supplied by the user ...\n");

	if(args->fastaFileName!=0) {
		fprintf(stderr, "Validating fastaFileName %s. \n",
				args->fastaFileName);
		if(ValidateFileName(args->fastaFileName)==0)
			PrintError(FnName, "fastaFileName", "Command line argument", Exit, IllegalFileName);	
	}	
	else {		
		PrintError(FnName, "fastaFileName", "Required command line argument", Exit, IllegalFileName);	
	}

	if(args->space != NTSpace && args->space != ColorSpace) {
		PrintError(FnName, "space", "Command line argument", Exit, OutOfRange);	
	}	
	if(args->keySize < 0) {		
		PrintError(FnName, "keySize", "Command line argument", Exit, OutOfRange);
	}

	if(args->maxKeyMatches < 0) {
		PrintError(FnName, "maxKeyMatches", "Command line argument", Exit, OutOfRange);	
	}	
	if(args->maxNumMatches < 0) {		
		PrintError(FnName, "maxNumMatches", "Command line argument", Exit, OutOfRange);
	}

	if(!(args->whichStrand == BothStrands || 
				args->whichStrand == ForwardStrand || 
				args->whichStrand == ReverseStrand)) {
		PrintError(FnName, "whichStrand", "Command line argument", Exit, OutOfRange);	
	}	
	if(args->numThreads<=0) {		
		PrintError(FnName, "numThreads", "Command line argument", Exit, OutOfRange);
	} 

	if(args->queueLength<=0) {
		PrintError(FnName, "queueLength", "Command line argument", Exit, OutOfRange);	
	} 	
	if(args->tmpDir!=0) {		
		fprintf(stderr, "Validating tmpDir path %s. \n", 
				args->tmpDir);
		if(ValidateFileName(args->tmpDir)==0)
			PrintError(FnName, "tmpDir", "Command line argument", Exit, IllegalFileName);	
	}	
	/* If this does not hold, we have done something wrong internally */	
	assert(args->timing == 0 || args->timing == 1);

	return 1;
}

/* TODO */
	void 
BfastMatchAssignDefaultValues(struct arguments *args)
{
	/* Assign default values */

	args->programMode = ExecuteProgram;

	args->fastaFileName = NULL;

	args->mainIndexes = NULL;
	args->secondaryIndexes = NULL;
	args->readsFileName = NULL;
	args->offsets = NULL;

	args->space = NTSpace;

	args->startReadNum = -1;
	args->endReadNum = -1;
	args->keySize = 0;
	args->maxKeyMatches = MAX_KEY_MATCHES;
	args->maxNumMatches = MAX_NUM_MATCHES;
	args->whichStrand = BothStrands;
	args->numThreads = 1;
	args->queueLength = DEFAULT_MATCHES_QUEUE_LENGTH;

	args->tmpDir =
		(char*)malloc(sizeof(DEFAULT_OUTPUT_DIR));
	assert(args->tmpDir!=0);
	strcpy(args->tmpDir, DEFAULT_OUTPUT_DIR);

	args->timing = 0;

	return;
}

/* TODO */
	void 
BfastMatchPrintProgramParameters(FILE* fp, struct arguments *args)
{
	char programmode[3][64] = {"ExecuteGetOptHelp", "ExecuteProgram", "ExecutePrintProgramParameters"};
	char whichStrand[3][64] = {"BothStrands", "ForwardStrand", "ReverseStrand"};
	fprintf(fp, BREAK_LINE);
	fprintf(fp, "Printing Program Parameters:\n");
	fprintf(fp, "programMode:\t\t\t\t%d\t[%s]\n", args->programMode, programmode[args->programMode]);
	fprintf(fp, "fastaFileName:\t\t\t\t%s\n", args->fastaFileName);
	fprintf(fp, "mainIndexes\t\t\t\t%s\n", args->mainIndexes);
	fprintf(fp, "secondaryIndexes\t\t\t%s\n", args->secondaryIndexes);
	fprintf(fp, "readsFileName:\t\t\t\t%s\n", args->readsFileName);
	fprintf(fp, "offsets:\t\t\t\t%s\n", args->offsets);
	fprintf(fp, "space:\t\t\t\t\t%d\n", args->space);
	fprintf(fp, "startReadNum:\t\t\t\t%d\n", args->startReadNum);
	fprintf(fp, "endReadNum:\t\t\t\t%d\n", args->endReadNum);
	fprintf(fp, "keySize:\t\t\t\t%d\n", args->keySize);
	fprintf(fp, "maxKeyMatches:\t\t\t\t%d\n", args->maxKeyMatches);
	fprintf(fp, "maxNumMatches:\t\t\t\t%d\n", args->maxNumMatches);
	fprintf(fp, "whichStrand:\t\t\t\t%d\t[%s]\n", args->whichStrand, whichStrand[args->whichStrand]);
	fprintf(fp, "numThreads:\t\t\t\t%d\n", args->numThreads);
	fprintf(fp, "queueLength:\t\t\t\t%d\n", args->queueLength);
	fprintf(fp, "tmpDir:\t\t\t\t\t%s\n", args->tmpDir);
	fprintf(fp, "timing:\t\t\t\t\t%d\n", args->timing);
	fprintf(fp, BREAK_LINE);
	return;
}

/* TODO */
void BfastMatchFreeProgramParameters(struct arguments *args)
{
	free(args->fastaFileName);
	args->fastaFileName=NULL;
	free(args->mainIndexes);
	args->mainIndexes=NULL;
	free(args->secondaryIndexes);
	args->secondaryIndexes=NULL;
	free(args->readsFileName);
	args->readsFileName=NULL;
	free(args->offsets);
	args->offsets=NULL;
	free(args->tmpDir);
	args->tmpDir=NULL;
}

/* TODO */
void
BfastMatchGetOptHelp() {

	struct argp_option *a=options;
	fprintf(stderr, "%s %s\n", "bfast ", PACKAGE_VERSION);
	fprintf(stderr, "\nUsage: bfast match [options]\n");
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
BfastMatchGetOptParse(int argc, char** argv, char OptionString[], struct arguments* arguments) 
{
	char key;
	int OptErr=0;
	while((OptErr==0) && ((key = getopt (argc, argv, OptionString)) != -1)) {
		/*
		   fprintf(stderr, "Key is %c and OptErr = %d\n", key, OptErr);
		   */
		switch (key) {
			case 'e':
				arguments->endReadNum = atoi(optarg);break;
			case 'f':
				arguments->fastaFileName = strdup(optarg);break;
			case 'h':
				arguments->programMode=ExecuteGetOptHelp; break;
			case 'i':
				arguments->mainIndexes=strdup(optarg); break;
			case 'k':
				arguments->keySize = atoi(optarg);break;
			case 'n':
				arguments->numThreads=atoi(optarg);break;
			case 'o':
				arguments->offsets=strdup(optarg); break;
			case 'p':
				arguments->programMode=ExecutePrintProgramParameters;break;
			case 'r':
				arguments->readsFileName=strdup(optarg);break;
			case 's':
				arguments->startReadNum = atoi(optarg);break;
			case 't':
				arguments->timing = 1;break;
			case 'w':
				arguments->whichStrand = atoi(optarg);break;
			case 'A':
				arguments->space=atoi(optarg);break;
			case 'I':
				arguments->secondaryIndexes=strdup(optarg); break;
			case 'K':
				arguments->maxKeyMatches=atoi(optarg);break;
			case 'M':
				arguments->maxNumMatches=atoi(optarg);break;
			case 'Q':
				arguments->queueLength=atoi(optarg);break;
			case 'T':
				StringCopyAndReallocate(&arguments->tmpDir, optarg); break;
			default:
				OptErr=1;
		} /* while */
	} /* switch */
	return OptErr;
}
