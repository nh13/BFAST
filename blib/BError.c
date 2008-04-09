#include <stdio.h>
#include <stdlib.h>
#include "BError.h"
#include "BLibDefinitions.h"

static char ErrorString[][20]=
{ "\0", "OutOfRange", "InputArguments", "IllegalFileName", "OpenFileError", "EndOfFile", "ReallocMemory", "MallocMemory"}; 
static char ActionType[][20]={"Fatal Error", "Warning"};

	void
PrintError(char* FunctionName, char *VariableName, 
		char* Message, int Action, int type)
{

	fprintf(stderr, "\nIn function \"%s\": %s[%s]. ", 
			FunctionName, ActionType[Action], ErrorString[type]);

	/* Only print variable name if is available */
	if(VariableName) {
		fprintf(stderr, "Variable/Value: %s\n", VariableName);
	}
	/* Only print message name if is available */
	if(Message) { 
		fprintf(stderr, "Message: %s\n", Message);
	}

	switch(Action) {
		case Exit: 
			fprintf(stderr, " ***** Exiting due to errors *****\n"); 
			exit(EXIT_FAILURE); 
			break; /* Not necessary actually! */
		case Warn:
			return; break;
		default:
			fprintf(stderr, "Trouble!!!\n");
	}
}
