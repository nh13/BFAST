#include <stdio.h>
#include <stdlib.h>
#include "BError.h"
#include "BLibDefinitions.h"

static char ErrorString[][20]=
  { "\0", "OutOfRange", "IllegalFileName"}; 
static char ActionType[][20]={"Fatal Error", "Warning"};

void
PrintError(char* FunctionName, char *VariableName, 
	   char* Message, int type, int Action) 
{
  
  fprintf(stderr, "\nIn function \"%s\": %s[%s]. ", 
	  FunctionName, ActionType[Action], ErrorString[type]);

  /* 
     Based on the type of error, variable may or may not
     have a "printable" string
  */
  
  if(VariableName)
    fprintf(stderr, "VariableName: %s ", VariableName);

    fprintf(stderr, "\n");
  
  switch(Action) {
  case Exit: 
    fprintf(stderr, " ***** Exiting due to errors *****\n"); 
    exit(0); 
    break; /* Not necessary actually! */
  case Warn:
    return; break;
  default:
    fprintf(stderr, "Trouble!!!\n");
  }
}
