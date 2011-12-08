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
#include "RunPostProcess.h"
#include "BfastPostProcess.h"

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC, OPTIONAL_GROUP_NAME}.
   */
enum { 
	DescInputFilesTitle, DescFastaFileName, DescInputFileName, 
	DescAlgoTitle, DescAlgorithm, DescSpace, DescStrandedness, DescPositioning, DescPairing, DescAvgMismatchQuality, 
	DescScoringMatrixFileName, DescRandomBest, DescMinimumMappingQuality, DescMinimumNormalizedScore,  
	DescNumThreads, DescQueueLength, 
	DescOutputTitle, DescOutputFormat, DescOutputID, DescRGFileName, DescBaseQualityType, DescTiming,
	DescMiscTitle, DescParameters, DescHelp
};

static struct argp_option options[] = {
	{0, 0, 0, 0, "=========== Input Files =============================================================", 1},
	{"fastaFileName", 'f', "fastaFileName", 0, "Specifies the file name of the FASTA reference genome", 1},
	{"alignFileName", 'i', "alignFileName", 0, "Specifies the input file from the balign program", 1},
	{0, 0, 0, 0, "=========== Algorithm Options =======================================================", 2},
	{"algorithm", 'a', "algorithm", 0, "Specifies the algorithm to choose the alignment for each end" 
		"\n\t\t\t  of the read:"
			"\n\t\t\t  0: No filtering will occur."
			"\n\t\t\t  1: All alignments that pass the filters will be outputted"
			"\n\t\t\t  2: Only consider reads that have been aligned uniquely"
			"\n\t\t\t  3: Choose uniquely the alignment with the best score"
			"\n\t\t\t  4: Choose all alignments with the best score",
		2},
	{"space", 'A', "space", 0, "0: NT space 1: Color space", 2},
	{"strandedness", 'S', 0, OPTION_NO_USAGE, "Specifies the pairing strandedness:"
			"\n\t\t\t  0: The reads should be mapped onto the same strand"
			"\n\t\t\t  1: The reads should be mapped onto the opposite strand",
                        2},
	{"positioning", 'P', 0, OPTION_NO_USAGE, "Specifies the pairing positioning:"
			"\n\t\t\t  0: The first read should be upstream of the second read (sequencing strand)"
			"\n\t\t\t  1: The second read should be upstream of the first read (sequencing strand)"
			"\n\t\t\t  2: There is no positioning",
                        2},
	{"pairing", 'Y', 0, OPTION_NO_USAGE, "Specifies the pairing options (overrides -S and -P):"
			"\n\t\t\t  0: paired ends"
			"\n\t\t\t  1: mate pairs"
                        "\n\t\t\t  2: no pairing", 
                        2},
	{"avgMismatchQuality", 'q', "avgMismatchQuality", 0, "Specifies the average mismatch quality", 2},
	{"scoringMatrixFileName", 'x', "scoringMatrixFileName", 0, "Specifies the file name storing the scoring matrix", 1},
	{"randomBest", 'z', 0, OPTION_NO_USAGE, "Specifies to output a random best scoring alignment (with -a 3)", 2},
	{"minMappingQuality", 'm', "minMappingQuality", 0, "Specifies to remove low mapping quality alignments", 2},
	{"minNormalizedScore", 'M', "minNormalizedScore", 0, "Specifies to remove low (alignment) scoring alignments", 2},
	{"insertSizeAvg", 'v', "insertSizeAvg", 0, "Specifies the mean insert size to use when pairing", 2}, 
	{"insertSizeStdDev", 's', "insertSizeStdDev", 0, "Specifies the standard deviation of the insert size to use when pairing", 2}, 
	{"numThreads", 'n', "numThreads", 0, "Specifies the number of threads to use (Default 1)", 2},
	{"queueLength", 'Q', "queueLength", 0, "Specifies the number of reads to cache", 2},
	{0, 0, 0, 0, "=========== Output Options ==========================================================", 3},
	{"outputFormat", 'O', "outputFormat", 0, "Specifies the output format 0: BAF 1: SAM", 3},
	{"outputID", 'o', "outputID", 0, "Specifies output ID to prepend to the read name (SAM only)", 3},
	{"RGFileName", 'r', "RGFileName", 0, "Specifies to add the RG in the specified file to the SAM"
		"\n\t\t\t  header and updates the RG tag (and LB/PU tags if present) in"
			"\n\t\t\t  the reads (SAM only)", 3},
        {"baseQualityType", 'b', "baseQualityType", 0, "Specifies the base quality type for SOLiD reads:"
			"\n\t\t\t  0: MAQ-style"
			"\n\t\t\t  1: Minimum (min(color 1, color 2))"
			"\n\t\t\t  2: Maximum (max(color 1, color 2))"
			"\n\t\t\t  3: Nullify (if either are an error, base quality is zero)",
                        3},
	{"timing", 't', 0, OPTION_NO_USAGE, "Specifies to output timing information", 3},
	{0, 0, 0, 0, "=========== Miscellaneous Options ===================================================", 4},
	{"Parameters", 'p', 0, OPTION_NO_USAGE, "Print program parameters", 4},
	{"Help", 'h', 0, OPTION_NO_USAGE, "Display usage summary", 4},
	{0, 0, 0, 0, 0, 0}
};

static char OptionString[]=
"a:b:i:f:m:n:o:q:r:s:v:x:A:M:O:P:S:Y:Q:hptzRU";

	int
BfastPostProcess(int argc, char **argv)
{
	struct arguments arguments;
	time_t startTime = time(NULL);
	time_t endTime;
	RGBinary rg;
	char *readGroup=NULL;

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
						if(0 <= VERBOSE) {
							fprintf(stderr, "Input arguments look good!\n");
							fprintf(stderr, BREAK_LINE);
						}
					}
					else {
						PrintError("PrintError", NULL, "validating command-line inputs", Exit, InputArguments);
					}
					BfastPostProcessPrintProgramParameters(stderr, &arguments);
					/* Execute program */
					if(BAF != arguments.outputFormat) {
						/* Read binary */
						RGBinaryReadBinary(&rg,
								NTSpace,
								arguments.fastaFileName);
						if(SAM == arguments.outputFormat && NULL != arguments.RGFileName) {
							readGroup = ReadInReadGroup(arguments.RGFileName);
						}
					}
					ReadInputFilterAndOutput(&rg,
							arguments.alignFileName,
							arguments.algorithm,
							arguments.space,
							arguments.strandedness,
							arguments.positioning,
							(3 == arguments.pairing) ? 1 : 0,
							arguments.avgMismatchQuality,
							arguments.scoringMatrixFileName,
							arguments.randomBest,
							arguments.minMappingQuality,
							arguments.minNormalizedScore,
							arguments.insertSizeSpecified,
							arguments.insertSizeAvg,
							arguments.insertSizeStdDev,
							arguments.numThreads,
							arguments.queueLength,
							arguments.outputFormat,
							arguments.outputID,
							readGroup,
                                                        arguments.baseQualityType,
							stdout);
					if(BAF != arguments.outputFormat) {
						/* Free rg binary */
						RGBinaryDelete(&rg);
						if(SAM == arguments.outputFormat && NULL != arguments.RGFileName) {
							free(readGroup);
						}
					}
					if(arguments.timing == 1) {
						/* Get the time information */
						endTime = time(NULL);
						int seconds = endTime - startTime;
						int hours = seconds/3600;
						seconds -= hours*3600;
						int minutes = seconds/60;
						seconds -= minutes*60;
						if(0 <= VERBOSE) {
							fprintf(stderr, "Total time elapsed: %d hours, %d minutes and %d seconds.\n",
									hours,
									minutes,
									seconds
								   );
						}
					}
					if(0 <= VERBOSE) {
						fprintf(stderr, "Terminating successfully!\n");
						fprintf(stderr, "%s", BREAK_LINE);
					}
					break;
				default:
					PrintError("PrintError", "programMode", "Could not determine program mode", Exit, OutOfRange);
			}

		}
		else {
			PrintError("PrintError", NULL, "Could not parse command line arguments", Exit, InputArguments);
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

	/* Check if we are piping */
	if(NULL == args->alignFileName) {
		VERBOSE = -1;
	}

	if(0 <= VERBOSE) {
		fprintf(stderr, BREAK_LINE);
		fprintf(stderr, "Checking input parameters supplied by the user ...\n");
	}

	if(args->fastaFileName!=0) {
		if(0 <= VERBOSE) {
			fprintf(stderr, "Validating fastaFileName %s. \n",
					args->fastaFileName);
		}
		if(ValidateFileName(args->fastaFileName)==0)
			PrintError(FnName, "fastaFileName", "Command line argument", Exit, IllegalFileName);	
	}	
	else if(args->outputFormat != BAF) {
		PrintError(FnName, "fastaFileName", "Required command line argument", Exit, IllegalFileName);
	}

	if(args->alignFileName!=0) {		
		if(0 <= VERBOSE) {
			fprintf(stderr, "Validating alignFileName %s. \n", 
					args->alignFileName);
		}
		if(ValidateFileName(args->alignFileName)==0)
			PrintError(FnName, "alignFileName", "Command line argument", Exit, IllegalFileName);	
	}	

	if(args->algorithm < MIN_FILTER || 			
			args->algorithm > MAX_FILTER) {
		PrintError(FnName, "algorithm", "Command line argument", Exit, OutOfRange);	
	}	

	if(args->space != NTSpace && args->space != ColorSpace) {
		PrintError(FnName, "space", "Command line argument", Exit, OutOfRange);
	}

	if(args->scoringMatrixFileName!=0) {
		if(0 <= VERBOSE) {
			fprintf(stderr, "Validating scoringMatrixFileName path %s. \n",
					args->scoringMatrixFileName);
		}
		if(ValidateFileName(args->scoringMatrixFileName)==0)
			PrintError(FnName, "scoringMatrixFileName", "Command line argument", Exit, IllegalFileName);
	}
	if(args->avgMismatchQuality <= 0) {
		PrintError(FnName, "avgMismatchQuality", "Command line argument", Exit, OutOfRange);
	}

	if(args->numThreads <= 0) {
		PrintError(FnName, "numThreads", "Command line argument", Exit, OutOfRange);
	}

	if(args->queueLength<=0) {		
		PrintError(FnName, "queueLength", "Command line argument", Exit, OutOfRange);
	}

	if(!(args->outputFormat == BAF ||				
				args->outputFormat == SAM)) {
		PrintError(FnName, "outputFormat", "Command line argument", Exit, OutOfRange);	
	}	
	assert(args->timing == 0 || args->timing == 1);
	if(args->pairing < 0 && 3 < args->pairing) {
		PrintError(FnName, "pairing", "Command line argument", Exit, OutOfRange);	
	}
	if(0 < args->pairing) {
		if(1 == args->pairing || 2 ==  args->pairing) {
			args->strandedness = (1 == args->pairing) ? 1 : 0; 
			args->positioning = (1 == args->pairing) ? 0 : 1; 
		}
		else {
			args->strandedness = -1;
			args->positioning = -1;
		}
	}
	else {
		if(args->strandedness != 0 && args->strandedness != 1) {
			PrintError(FnName, "strandedness", "Command line argument", Exit, OutOfRange);	
		}
		if(args->positioning < 0 || 2 < args->positioning) {
			PrintError(FnName, "positioning", "Command line argument", Exit, OutOfRange);	
		}
	}
	assert(args->randomBest == 0 || args->randomBest == 1);

        if(args->baseQualityType < 0 || 3 < args->baseQualityType) {
		PrintError(FnName, "baseQualityType", "Command line argument", Exit, OutOfRange);	
        }

	if(SAM != args->outputFormat && NULL != args->RGFileName) {
		PrintError(FnName, "RGFileName", "Command line argument can only be used when outputting to SAM format", Exit, OutOfRange);
	}

	if (1 == args->insertSizeSpecified) {
		if (args->insertSizeStdDev <= 0.0) {
			PrintError(FnName, "insertSizeStdDev", "When specifying insertSizeAvg, you must also specify an insertSizeStdDev > 0.", Exit, OutOfRange);
		}
		if (args->insertSizeAvg <= 0.0) {
			PrintError(FnName, "insertSizeAvg", "insertSizeAvg <= 0.", Warn, OutOfRange);
		}
	}

	return 1;
}

/* TODO */
	void 
BfastPostProcessAssignDefaultValues(struct arguments *args)
{
	/* Assign default values */

	args->programMode = ExecuteProgram;

	args->fastaFileName = NULL;
	args->alignFileName = NULL;

	args->algorithm=BestScore;
	args->space = NTSpace;
	args->strandedness=-1;
	args->positioning=-1;
        args->pairing=3;
	args->scoringMatrixFileName=NULL;
	args->randomBest=0;
	args->minMappingQuality=INT_MIN;
	args->minNormalizedScore=INT_MIN;
	args->insertSizeSpecified=0;
	args->insertSizeAvg=0.0;
	args->insertSizeStdDev=0.0;
	args->avgMismatchQuality=AVG_MISMATCH_QUALITY;
	args->numThreads=1;
	args->queueLength=DEFAULT_POSTPROCESS_QUEUE_LENGTH;

	args->outputFormat=SAM;
	args->outputID=NULL;
	args->RGFileName=NULL;
        args->baseQualityType=0;

	args->timing = 0;

	return;
}

/* TODO */
	void 
BfastPostProcessPrintProgramParameters(FILE* fp, struct arguments *args)
{
	char algorithm[5][64] = {"[No Filtering]", "[Filtering Only]", "[Unique]", "[Best Score]", "[Best Score All]"};
	char outputType[8][32] = {"[BRG]", "[BIF]", "[BMF]", "[BAF]", "[SAM]", "[LastFileType]"};
	char baseQualityType[4][32] = {"[MAQ-style]", "[Min]", "[Max]", "[Nullify]"};
        char strandedness[2][32] = {"[Same strand]", "[Opposite strand]"};
        char positioning[3][32] = {"[Read one first]", "[Read two first]", "[No Positioning]"};
        char pairing[4][32] = {"[Not Using]", "[Paired End]", "[Mate Pair]", "[No Pairing]"};
	if(0 <= VERBOSE) {
		fprintf(fp, BREAK_LINE);
		fprintf(fp, "Printing Program Parameters:\n");
		fprintf(fp, "programMode:\t\t\t%s\n", PROGRAMMODE(args->programMode));
		fprintf(fp, "fastaFileName:\t\t\t%s\n", FILEREQUIRED(args->fastaFileName));
		fprintf(fp, "alignFileName:\t\t\t%s\n", FILESTDIN(args->alignFileName));
		fprintf(fp, "algorithm:\t\t\t%s\n", algorithm[args->algorithm]);
		fprintf(fp, "space:\t\t\t\t%s\n", SPACE(args->space));
		fprintf(fp, "strandedness:\t\t\t%s\n", BOOLREQUIRED(args->strandedness, strandedness));
		fprintf(fp, "positioning:\t\t\t%s\n", BOOLREQUIRED(args->positioning, positioning));
		fprintf(fp, "pairing:\t\t\t%s\n", pairing[args->pairing]);
		fprintf(fp, "avgMismatchQuality:\t\t%d\n", args->avgMismatchQuality);
		fprintf(fp, "scoringMatrixFileName:\t\t%s\n", FILEUSING(args->scoringMatrixFileName));
		fprintf(fp, "randomBest:\t\t\t%s\n", INTUSING(args->randomBest));
		fprintf(fp, "minMappingQuality:\t\t%d\n", args->minMappingQuality);
		fprintf(fp, "minNormalizedScore:\t\t%d\n", args->minNormalizedScore);
		if (0 == args->insertSizeSpecified) {
			fprintf(fp, "insertSizeAvg:\t\t\t%lf\n", args->insertSizeAvg);
			fprintf(fp, "insertSizeStdDev:\t\t%lf\n", args->insertSizeStdDev);
		}
                else { 
			fprintf(fp, "insertSizeAvg:\t\t\t%s\n", INTUSING(0));
			fprintf(fp, "insertSizeStdDev:\t\t%s\n", INTUSING(0));
                }
		fprintf(fp, "numThreads:\t\t\t%d\n", args->numThreads);
		fprintf(fp, "queueLength:\t\t\t%d\n", args->queueLength);
		fprintf(fp, "outputFormat:\t\t\t%s\n", outputType[args->outputFormat]);
		fprintf(fp, "outputID:\t\t\t%s\n", FILEUSING(args->outputID));
		fprintf(fp, "RGFileName:\t\t\t%s\n", FILEUSING(args->RGFileName));
		fprintf(fp, "baseQualityType:\t\t\t%s\n", baseQualityType[args->baseQualityType]);
		fprintf(fp, "timing:\t\t\t\t%s\n", INTUSING(args->timing));
		fprintf(fp, BREAK_LINE);
	}
	return;
}

/* TODO */
void BfastPostProcessFreeProgramParameters(struct arguments *args)
{
	free(args->fastaFileName);
	args->fastaFileName=NULL;
	free(args->alignFileName);
	args->alignFileName=NULL;
	free(args->outputID);
	args->outputID=NULL;
	free(args->RGFileName);
	args->RGFileName=NULL;
	free(args->scoringMatrixFileName);
	args->scoringMatrixFileName=NULL;
}

/* TODO */
void
BfastPostProcessGetOptHelp() {

	struct argp_option *a=options;
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
	int key;
	int OptErr=0;
	while((OptErr==0) && ((key = getopt (argc, argv, OptionString)) != -1)) {
		/*
		   fprintf(stderr, "Key is %c and OptErr = %d\n", key, OptErr);
		   */
		switch (key) {
			case 'a':
				arguments->algorithm = atoi(optarg);break;
                        case 'b':
                                arguments->baseQualityType = atoi(optarg);break;
			case 'f':
				arguments->fastaFileName=strdup(optarg);break;
			case 'h':
				arguments->programMode=ExecuteGetOptHelp;break;
			case 'i':
				arguments->alignFileName=strdup(optarg);break;
				break;
			case 'm':
				arguments->minMappingQuality=atoi(optarg);break;
			case 'n':
				arguments->numThreads=atoi(optarg); break;
			case 'o':
				arguments->outputID=strdup(optarg);break;
			case 'p':
				arguments->programMode=ExecutePrintProgramParameters;break;
			case 'q':
				arguments->avgMismatchQuality = atoi(optarg); break;
			case 'r':
				arguments->RGFileName=strdup(optarg);break;
			case 's':
				arguments->insertSizeStdDev=atof(optarg);break;
			case 'v':
				arguments->insertSizeSpecified=1;
				arguments->insertSizeAvg=atof(optarg);break;
			case 't':
				arguments->timing = 1; break;
			case 'x':
				StringCopyAndReallocate(&arguments->scoringMatrixFileName, optarg);
				break;
			case 'z':
				arguments->randomBest = 1; break;
			case 'A':
				arguments->space=atoi(optarg);break;
			case 'M':
				arguments->minNormalizedScore=atoi(optarg);break;
			case 'O':
				switch(atoi(optarg)) {
					case 0:
						arguments->outputFormat = BAF;
						break;
					case 1:
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
			case 'S':
				arguments->pairing = 0; 
				arguments->strandedness = atoi(optarg); break;
			case 'P':
				arguments->pairing = 0; 
				arguments->positioning = atoi(optarg); break;
			case 'Y':
				arguments->pairing = atoi(optarg)+1; break;
			default:
				OptErr=1;
		} /* while */
	} /* switch */
	return OptErr;
}
