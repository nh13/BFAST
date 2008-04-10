#include <stdio.h>

/* This structure is used by main to communicate with parse_opt. */
struct arguments
{
	char *args[1];							/* No arguments to this function */
	char *rgListFileName;                   /* -r */
	char *matchesFileName;					/* -m */
	char *scoringMatrixFileName;			/* -x */
	int algorithm;							/* -a */
	int startChr;                           /* -s */
	int startPos;                           /* -S */
	int endChr;                             /* -e */
	int endPos;                             /* -E */
	int offsetLength;						/* -O */
	int maxNumMatches;						/* -M */
	int pairedEnd;                          /* -2 */
	int numThreads;                         /* -n */
	char *outputID;							/* -o */
	char *outputDir;						/* -d */
	int timing;                             /* -t */
	int programMode;						/* -h */ 
};

/* Local functions */
static int ValidateInputs(struct arguments*);
int ValidateFileName(char*);
void AssignDefaultValues(struct arguments*);
void PrintProgramParameters(FILE*, struct arguments*);
void PrintGetOptHelp();
void GetOptHelp();
void PrintGetOptHelp();
/*
   PARSER. Field 2 in ARGP.
   Order of parameters: KEY, ARG, STATE.
   KEY  -- An integer specifying which option this is (taken
   from the KEY field ), or
   a special key specifying something else; the only
   special keys we use here are ARGP_KEY_ARG, meaning
   a non-option argument, and ARGP_KEY_END, meaning
   that all arguments have been parsed

   ARG  -- For an option KEY, the string value of its
   argument, or NULL if it has none

   STATE-- A pointer to a struct argp_state (defined in argp.h), containing
   various useful information about the parsing state; used here
   are the INPUT field, which reflects the INPUT argument to
   argp_parse, and the ARG_NUM field, which is the number of the
   current non-option argument being parsed

   It should return either 0, meaning success, ARGP_ERR_UNKNOWN, meaning the
   given KEY wasn't recognized, or an errno value indicating some other
   error.
   */
#ifdef HAVE_ARGP_H
static error_t parse_opt (int, char*, struct argp_state*);
#else
struct argp_option {
	char *name; /* Arg name */
	int key;
	char *arg; /* arg symbol */
	int flags; 
	char *doc; /* short info about the arg */
	int group;
};
#define OPTION_ARG_OPTIONAL 0
#define OPTION_NO_USAGE 0
int getopt_parse(int, char**, char*, struct arguments*); 
#endif

/* From Error handling routine */
extern void PrintError(char*, char*, char*, int, int);
