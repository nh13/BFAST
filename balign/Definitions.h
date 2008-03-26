#define DEFAULT_FILENAME "Default.txt"
#define BREAK_LINE "************************************************************\n"


/* For BError.c  */
enum {Exit, Warn, LastActionType};
enum {
	Dummy,
	OutofRange, /* e.g. command line args */
	IllegalFileName,   /*  KeepAdding */
	LastErrorType
};
