#ifndef SRTREE_H_
#define SRTREE_H_

#include <stdio.h>

/* TODO */
typedef struct {
	void *next[4];
} SRNode;

enum {SRT_KILOBYTES, SRT_MEGABYTES, SRT_GIGABYTES};

int SRTreeInsert(SRNode*, char*, int);
void SRTreeDelete(SRNode**, int, int);
double SRTreeGetSize(SRNode*, int, int);
double SRTreeGetSizeHelper(SRNode*, int, int, int);
int SRTreeReadFromFile(SRNode*, FILE*);
void SequenceToLower(char*, int);

#endif
