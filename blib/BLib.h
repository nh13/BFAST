#ifndef BLIB_H_
#define BLIB_H_

char ToLower(char);
char ToUpper(char);
void GetReverseComplimentAnyCase(char*, char*, int);
char GetReverseComplimentAnyCaseBase(char);
int ValidateBasePair(char);
int IsAPowerOfTwo(unsigned int);
char TransformFromIUPAC(char);
void CheckRGIndexes(char**, int, char**, int, int, int*, int*, int*, int*);

#endif
