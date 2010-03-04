#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <config.h>
#include <pthread.h>
#include <unistd.h>

#include "../bfast/BError.h"
#include "../bfast/RGMatches.h"
#include "../bfast/RGBinary.h"

#define Name "bmffix"

int bmffix(char *matchFileName, FILE *fpOut, RGBinary *rg)
{
	char *FnName="bmffix";
	RGMatches m;
	gzFile matchFP=NULL;
	gzFile outputFP=NULL;

	RGMatchesInitialize(&m);
	if((matchFP=gzopen(matchFileName, "rb"))==0) {
		PrintError(FnName, matchFileName, "Could not open file for reading", Exit, OpenFileError);
	}
	if((outputFP=gzdopen(fileno(fpOut), "wb"))==0) {
		PrintError(FnName, "stdout", "Could not open stdout file for writing", Exit, OpenFileError);
	}

	while(EOF != RGMatchesRead(matchFP, &m)) {
		RGMatchesFixConstraints(&m, rg);
		RGMatchesPrint(outputFP, &m);
		RGMatchesFree(&m);
	}

	gzclose(matchFP);
	gzclose(outputFP);

	return 0;
}

int PrintUsage()
{
	fprintf(stderr, "%s %s\n", "bfast", PACKAGE_VERSION);
	fprintf(stderr, "\nUsage:%s [options] <bmf files>\n", Name);
	fprintf(stderr, "\t-f\tFILE\tSpecifies the file name of the FASTA reference genome\n");
	    fprintf(stderr, "\t-A\tINT\t0: NT space 1: Color space\n");
	fprintf(stderr, "\t-h\t\tprints this help message\n");
	fprintf(stderr, "\nsend bugs to %s\n",
			PACKAGE_BUGREPORT);

	return 1;
}

int main(int argc, char *argv[])
{
	int c, i, space=NTSpace;
	char *fastaFileName=NULL;
	RGBinary rg;

	while((c = getopt(argc, argv, "f:hA:")) >= 0) {
		switch(c) {
			case 'f': fastaFileName=strdup(optarg); break;
			case 'A': space=atoi(optarg); break;
			case 'h': return PrintUsage();
			default: fprintf(stderr, "Unrecognized option: -%c\n", c); return 1;
		}
	}

	if(1 == argc || argc == optind) {
		return PrintUsage();
	}
	if(NULL == fastaFileName) {
		PrintError(Name, "fastaFileName", "Command line option", Exit, InputArguments);
	}

	RGBinaryReadBinary(&rg, space, fastaFileName);

	for(i=optind;i<argc;i++) {
		bmffix(argv[i], stdout, &rg); 
	}

	RGBinaryDelete(&rg);


	return 0;
}
