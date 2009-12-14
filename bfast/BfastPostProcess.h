#ifndef BFASTPOSTPROCESS_H_
#define BFASTPOSTPROCESS_H_

#define MIN_FILTER 0
#define MAX_FILTER 4
#define DEFAULT_QUEUE_LENGTH 10000

/* This structure is used by main to communicate with parse_opt. */
struct arguments
{
	char *args[1];							/* No arguments to this function */
	char *fastaFileName;					/* -f */
	char *alignFileName;					/* -i */
	int algorithm;							/* -a */
	int pairedEndInfer;						/* -P */
	int queueLength;						/* -Q */
	char *unmappedFileName;					/* -i */
	int outputFormat;						/* -O */
	char *outputID;							/* -o */
	char *readGroup;						/* -r */
	int timing;                             /* -t */
	int programMode;						/* -h */ 
};

/* Local functions */
int BfastPostProcessValidateInputs(struct arguments*);
void BfastPostProcessAssignDefaultValues(struct arguments*);
void BfastPostProcessPrintProgramParameters(FILE*, struct arguments*);
void BfastPostProcessFreeProgramParameters(struct arguments *args);
void BfastPostProcessPrintGetOptHelp();
void BfastPostProcessGetOptHelp();
void BfastPostProcessPrintGetOptHelp();
struct argp_option {
	char *name; /* Arg name */
	int key;
	char *arg; /* arg symbol */
	int flags;
	char *doc; /* short info about the arg */
	int group;
};
int BfastPostProcessGetOptParse(int, char**, char*, struct arguments*); 
#endif
