#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

/* TODO */
typedef struct {
	char *tmpFileName;
	FILE *tmpFP;
	int numEntries;
} TmpFP;

enum {BAlignFile, WigFile, BedFile, BedAndWigFile, LastFileType};

#endif
