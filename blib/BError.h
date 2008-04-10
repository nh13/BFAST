#ifndef BERROR_H_
#define BERROR_H_
/* Action */
enum {Exit, Warn, LastActionType};

/* Type */
enum {
	Dummy,
	OutOfRange, /* e.g. command line args */
	InputArguments,
	IllegalFileName,   
	OpenFileError,
	EndOfFile,
	ReallocMemory,
	MallocMemory,
	ThreadError,
	LastErrorType
};

void PrintError(char*, char*, char*, int, int);

#endif
