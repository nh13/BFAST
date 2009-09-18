#ifndef BFASTLOCALALIGN_H_
#define BFASTLOCALALIGN_H_

/* This structure is used by main to communicate with parse_opt. */
struct arguments
{
	char *args[1];							/* No arguments to this function */
	char *brgFileName;                   	/* -r */
	char *matchFileName;					/* -m */
	char *scoringMatrixFileName;			/* -x */
	int binaryInput;                        /* -b - not used */
	int alignmentType;						/* -a */
	int bestOnly;							/* -b */
	int space;								/* -A */
	int startReadNum;                       /* -s */
	int endReadNum;                         /* -e */
	int offsetLength;						/* -O */
	int maxNumMatches;						/* -M */
	int avgMismatchQuality;					/* -q */
	int numThreads;                         /* -n */
	int queueLength;                        /* -Q */
	int usePairedEndLength;					/* -l - companion to pairedEndLength */
	int mirroringType;						/* -L */
	int forceMirroring;						/* -f */
	int pairedEndLength;					/* -l */
	char *outputID;							/* -o */
	char *outputDir;						/* -d */
	char *tmpDir;							/* -T */
	int binaryOutput;                       /* -B - not used */
	int timing;                             /* -t */
	int programMode;						/* -h */ 
};

/* Local functions */
int BfastLocalAlignValidateInputs(struct arguments*);
void BfastLocalAlignAssignDefaultValues(struct arguments*);
void BfastLocalAlignPrintProgramParameters(FILE*, struct arguments*);
void BfastLocalAlignFreeProgramParameters(struct arguments *args);
void BfastLocalAlignPrintGetOptHelp();
void BfastLocalAlignGetOptHelp();
void BfastLocalAlignPrintGetOptHelp();
struct argp_option {
	char *name; /* Arg name */
	int key;
	char *arg; /* arg symbol */
	int flags;
	char *doc; /* short info about the arg */
	int group;
};
int BfastLocalAlignGetOptParse(int, char**, char*, struct arguments*);
#endif
