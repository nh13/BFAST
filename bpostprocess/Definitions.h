#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#define ROTATE_BINNING 100000

/* When splitting the input into temp files, have one temp
 * file for every 10 million bases 
 * */
#define DEFAULT_REGION_LENGTH 10000000

/* TODO */
typedef struct {
	FILE **files;
	char **fileNames;
	int numFiles;
} ChrFiles;

/* TODO */
typedef struct {
	int startChr;
	int endChr;
	FILE **chrFiles;
	char **chrFileNames;
} RGFiles;

enum {BAlignFile, WigFile, BedFile, LastFileType};

#endif
