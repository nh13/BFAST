#ifndef BFASTMATCH_H_
#define BFASTMATCH_H_

#include <stdio.h>

/* This structure is used by main to communicate with parse_opt. */
struct arguments
{
	char *args[1];							/* No arguments to this function */
	char *brgFileName;                   	/* -r */
	char *bfastMainIndexesFileName;			/* -i */
	char *bfastSecondaryIndexesFileName;	/* -I */
	char *readsFileName;					/* -R */
	char *offsetsFileName;					/* -O */
	int space;								/* -A */
	int binaryInput;						/* -b - not used */
	int startReadNum;						/* -s */
	int endReadNum;							/* -e */
	int numMismatches;						/* -x - not used */
	int numInsertions;						/* -y - not used */
	int numDeletions;						/* -z - not used */
	int numGapInsertions;					/* -Y - not used */
	int numGapDeletions;					/* -Z - not used */
	int keySize;							/* -k */
	int maxKeyMatches;						/* -K */
	int maxNumMatches;						/* -M */
	int whichStrand;						/* -w */
	int numThreads;							/* -n */
	int queueLength;						/* -Q */
	char *outputID;							/* -o */
	char *outputDir;						/* -d */
	char *tmpDir;							/* -T */
	int binaryOutput;						/* -B - not used */
	int timing;								/* -t */
	int programMode;						/* -h */ 
};

/* Local functions */
int BfastMatchValidateInputs(struct arguments*);
void BfastMatchAssignDefaultValues(struct arguments*);
void BfastMatchPrintProgramParameters(FILE*, struct arguments*);
void BfastMatchFreeProgramParameters(struct arguments *args);
void BfastMatchPrintGetOptHelp();
void BfastMatchGetOptHelp();
void BfastMatchPrintGetOptHelp();
struct argp_option {
	char *name; /* Arg name */
	int key;
	char *arg; /* arg symbol */
	int flags; 
	char *doc; /* short info about the arg */
	int group;
};
int BfastMatchGetOptParse(int, char**, char*, struct arguments*); 
#endif
