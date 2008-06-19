#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#define ROTATE_BINNING 100000

/* When splitting the input into temp files, have one temp
 * file for every 1 million bases 
 * */
#define REGION_LENGTH 1000000

/* TODO */
typedef struct {
	FILE **files;
	char **fileNames;
	int numFiles;
} ChrFiles;

enum {BAlignFile, WigFile, BedFile, LastFileType};

#endif
