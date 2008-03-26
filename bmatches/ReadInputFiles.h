#include <stdio.h>
#include "../blib/SRTree.h"
#include "../blib/RGTree.h"

void ReadSequences(char*, SRNode*, int);
void ReadSequencesToTempFile(char*, FILE**, int, int, int);
int ReadNextSequence(FILE*, char**, char**, char**, int);
void ReadRGTree(char*, RGTree*);
int ReadRGTreeFileNames(char*, char***, int**, int*);
