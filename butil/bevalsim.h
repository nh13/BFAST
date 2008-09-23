#ifndef BEVALSIM_H_
#define BEVALSIM_H_

typedef struct {
	/* Meta data */
	char strand;
	int chr;
	int pos;
	int space;
	int pairedEnd;
	int pairedEndLength;
	int readLength;
	int whichReadVariants;
	int startIndel;
	int indelLength;
	/* Inferred data */
	int numSNPs;
	int numErrors;
	int deletionLength;
	int insertionLength;
} ReadType;

void ReadTypeInitialize(ReadType*);
void ReadTypeReadFromRAF(ReadType*, FILE*);

void Evaluate(char*);

#endif
