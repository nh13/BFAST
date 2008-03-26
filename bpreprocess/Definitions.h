#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

typedef struct {
	char *sequence;
	int startPos;
	int endPos;
	int chromosome;
} RGChr;

typedef struct {
	RGChr *chromosomes;
	int numChrs;
	int startChr;
	int startPos;
	int endChr;
	int endPos;
} RGList;

/* For ReadInputFile.c */
/* The amount of memory to additionally allocate when reallocating memory */ 
#define RIF_REALLOCATE_INCREMENT 65536
#define RIF_ROTATE_NUM 1000000

/* For GenerateTree.c */
#define GT_ROTATE_NUM 1000000

#endif
