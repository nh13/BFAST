#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <zlib.h>

#include "../blib/BError.h"
#include "../blib/BLib.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/RGIndexAccuracy.h"
#include "../blib/RGMatches.h"
#include "../blib/AlignedRead.h"
#include "../blib/ScoringMatrix.h"
#include "../balign/RunAligner.h"
#include "SimRead.h"
#include "balignsim.h"

#define round(x)(int)(x<0?ceil((x)-0.5):floor((x)+0.5))
#define READS_ROTATE_NUM 10000
#define Name "balignsim"

/* Generate synthetic reads given a number of variants and errors
 * from a reference genome and tests the various local alignment 
 * algorithms. */

int main(int argc, char *argv[]) 
{
	RGBinary rg;
	char rgFileName[MAX_FILENAME_LENGTH]="\0";
	char scoringMatrixFileName[MAX_FILENAME_LENGTH]="\0";
	int alignmentType=FullAlignment;
	int space = 0;
	int indel = 0;
	int indelLength = 0;
	int withinInsertion = 0;
	int numSNPs = 0;
	int numErrors = 0;
	int readLength = 0;
	int numReads = 0;
	int numThreads = 1;
	char tmpDir[MAX_FILENAME_LENGTH]="./";

	char *inputFile=NULL;
	char c;

	while((c = getopt(argc, argv, "f:n:r:x:T:")) >= 0) {
		switch(c) {
			case 'f': inputFile = strdup(optarg); break;
			case 'r': strcpy(rgFileName, optarg); break;
			case 'x': strcpy(scoringMatrixFileName, optarg); break;
			case 'n': numThreads=atoi(optarg); break;
			case 'T': strcpy(tmpDir, optarg); break;
			default: break;
		}
	}

	if(NULL != inputFile) { // hidden !
		/* Get reference genome */
		RGBinaryReadBinary(&rg,
				rgFileName);

		FILE *fpIn=NULL;
		if(!(fpIn = fopen(inputFile, "r"))) {
			PrintError(Name,
					inputFile,
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}
		while(0 == feof(fpIn)) {
			if(fscanf(fpIn, "%d %d %d %d %d %d %d %d %d",
						&alignmentType,
						&space,
						&indel,
						&indelLength,
						&withinInsertion,
						&numSNPs,
						&numErrors,
						&readLength,
						&numReads) < 0) {
				break;
			}
			Run(&rg,
					scoringMatrixFileName,
					alignmentType,
					space,
					indel,
					indelLength,
					withinInsertion,
					numSNPs,
					numErrors,
					readLength,
					numReads,
					numThreads,
					tmpDir);
		}
		fclose(fpIn);
		free(inputFile);
	}
	else if(argc == 14) {

		/* Get cmd line options */
		strcpy(rgFileName, argv[1]);
		strcpy(scoringMatrixFileName, argv[2]);
		alignmentType = atoi(argv[3]);
		space = atoi(argv[4]);
		indel = atoi(argv[5]);
		indelLength = atoi(argv[6]);
		withinInsertion = atoi(argv[7]);
		numSNPs = atoi(argv[8]);
		numErrors = atoi(argv[9]);
		readLength = atoi(argv[10]);
		numReads = atoi(argv[11]);
		numThreads = atoi(argv[12]);
		strcpy(tmpDir, argv[13]);

		/* Check cmd line options */
		assert(FullAlignment == alignmentType || MismatchesOnly == alignmentType);
		assert(NTSpace == space || ColorSpace == space);
		assert(indel >= 0 && indel <= 2);
		assert(indelLength > 0 || indel == 0);
		assert(withinInsertion == 0 || (indel == 2 && withinInsertion == 1));
		assert(numSNPs >= 0);
		assert(readLength > 0);
		assert(readLength < SEQUENCE_LENGTH);
		assert(numReads > 0);
		assert(numThreads > 0);

		/* Should check if we have enough bases for everything */

		/* Get reference genome */
		RGBinaryReadBinary(&rg,
				rgFileName);

		/* Run Simulation */
		Run(&rg,
				scoringMatrixFileName,
				alignmentType,
				space,
				indel,
				indelLength,
				withinInsertion,
				numSNPs,
				numErrors,
				readLength,
				numReads,
				numThreads,
				tmpDir);

		/* Delete reference genome */
		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Deleting reference genome.\n");
		RGBinaryDelete(&rg);
		fprintf(stderr, "%s", BREAK_LINE);

		fprintf(stderr, "%s", BREAK_LINE);
		fprintf(stderr, "Terminating successfully.\n");
		fprintf(stderr, "%s", BREAK_LINE);
	}
	else {
		fprintf(stderr, "Usage: %s [OPTIONS]\n", Name);
		fprintf(stderr, "\t<bfast reference genome file name (must be in nt space)>\n");
		fprintf(stderr, "\t<scoring matrix file name>\n");
		fprintf(stderr, "\t<alignmentType 0: Full alignment 1: mismatches only>\n");
		fprintf(stderr, "\t<space 0: nt space 1: color space>\n");
		fprintf(stderr, "\t<indel 0: none 1: deletion 2: insertion>\n");
		fprintf(stderr, "\t<indel length>\n");
		fprintf(stderr, "\t<include errors within insertion 0: false 1: true>\n");
		fprintf(stderr, "\t<# of SNPs>\n");
		fprintf(stderr, "\t<# of errors>\n");
		fprintf(stderr, "\t<read length>\n");
		fprintf(stderr, "\t<number of reads>\n");
		fprintf(stderr, "\t<number of threads\n");
		fprintf(stderr, "\t<tmp file directory>\n");
	}
	return 0;
}

/* TODO */
void Run(RGBinary *rg,
		char *scoringMatrixFileName,
		int alignmentType,
		int space,
		int indel,
		int indelLength,
		int withinInsertion,
		int numSNPs,
		int numErrors,
		int readLength,
		int numReads,
		int numThreads,
		char *tmpDir)
{
	char *FnName="Run";
	int i, j;
	int64_t rgLength=0;
	gzFile matchesFP=NULL;
	char *matchesFileName=NULL;
	gzFile alignFP=NULL;
	char *alignFileName=NULL;
	gzFile notAlignedFP=NULL;
	char *notAlignedFileName=NULL;
	int32_t totalAlignTime = 0;
	int32_t totalFileHandlingTime = 0;
	ScoringMatrix sm;
	double mismatchScore;
	SimRead r;
	RGMatches m;
	AlignedRead a;
	int32_t score, prev, score_m, score_mm, score_cm, score_ce, wasInsertion;
	int32_t numScoreLessThan, numScoreEqual, numScoreGreaterThan;
	int insertionLength = (2==indel)?indelLength:0;
	int deletionLength = (1==indel)?indelLength:0;
	char string[4096]="\0";
	int ret=0;
	char *s=NULL;

	if(ColorSpace == space &&
			1 == withinInsertion && 
			0 < indelLength && 
			2 == indel) {
		PrintError(Name,
				"withinInsertion",
				"Incosistent results will occurs.  Try not using withinInsertion == 1.",
				Warn,
				OutOfRange);
	}


	score = prev = score_m = score_mm = score_cm = score_ce = 0;

	/* Check rg to make sure it is in NT Space */
	if(rg->space != NTSpace) {
		PrintError(FnName,
				"rg->space",
				"The reference genome must be in NT space",
				Exit,
				OutOfRange);
	}


	/* ********************************************************
	 * 1.
	 * Generate reads and create artificial matches file.
	 ********************************************************
	 */

	fprintf(stderr, "%s", BREAK_LINE);
	fprintf(stderr, "Generating reads and creating artificial matches.\n");

	/* Get the reference genome length */
	for(i=0;i<rg->numContigs;i++) {
		rgLength += rg->contigs[i].sequenceLength;
	}

	/* Open tmp files */
	matchesFP = OpenTmpGZFile(tmpDir, &matchesFileName);
	alignFP = OpenTmpGZFile(tmpDir, &alignFileName);
	notAlignedFP = OpenTmpGZFile(tmpDir, &notAlignedFileName);

	/* Get scoring matrix */
	ScoringMatrixInitialize(&sm);
	ScoringMatrixRead(scoringMatrixFileName,
			&sm,
			space);
	if(space == NTSpace) {
		mismatchScore = sm.ntMatch - sm.ntMismatch;
	}
	else {
		mismatchScore = sm.colorMatch - sm.colorMismatch;
	}
	/* In these sims we want the scoring matrix to have certain constraints:
	 * All scores for matches must be the same and all scores for mismatches 
	 * must be the same */

	score_m = sm.ntMatch;
	score_mm = sm.ntMismatch;
	if(ColorSpace == space) {
		score_cm = sm.colorMatch;
		score_ce = sm.colorMismatch;
	}

	/* Seed random number */
	srand(time(NULL));

	/* Create RGMatches */
	RGMatchesInitialize(&m);
	SimReadInitialize(&r);
	fprintf(stderr, "Currently on:\n0");
	for(i=0;i<numReads;i++) {
		if((i+1) % READS_ROTATE_NUM==0) {
			fprintf(stderr, "\r%d",
					(i+1));
		}

		/* Get the read */
		SimReadGetRandom(rg,
				rgLength,
				&r,
				space,
				indel,
				indelLength,
				withinInsertion,
				numSNPs,
				numErrors,
				readLength,
				1,
				0);
		r.readNum = i+1;

		/* Convert into RGMatches */
		RGMatchesInitialize(&m);
		RGMatchesReallocate(&m, 1);
		/* Get score for proper alignment and store in read name */
		score = 0;
		prev = Default;
		wasInsertion=0;
		if(NTSpace == space) {
			for(j=0;j<r.readLength;j++) {
				switch(r.readOneType[j]) {
					case Insertion:
					case InsertionAndSNP:
					case InsertionAndError:
					case InsertionSNPAndError:
						if(Insertion == prev) {
							score += sm.gapExtensionPenalty;  
						}
						else {
							score += sm.gapOpenPenalty;  
						}
						wasInsertion=1;
						break;
					case SNP:
					case Error:
					case SNPAndError:
						score += score_mm;
						break;
					case Default:
						score += score_m;
						break;
					default:
						fprintf(stderr, "r.readOneType[%d]=%d\n", j, r.readOneType[j]);
						PrintError(FnName,
								"r.readOneType[j]",
								"Could not recognize type",
								Exit,
								OutOfRange);
				}
				prev = r.readOneType[j];
			}
		} 
		else {
			for(j=0;j<r.readLength;j++) {
				switch(r.readOneType[j]) {
					case Insertion:
					case InsertionAndSNP:
					case InsertionAndError:
					case InsertionSNPAndError:
						if(Insertion == prev) {
							score += sm.gapExtensionPenalty;  
						}
						else {
							score += sm.gapOpenPenalty;  
						}
						wasInsertion=1;
						break;
					case SNP:
						score += score_mm + score_cm;
						break;
					case Error:
						score += score_m + score_ce;
						break;
					case SNPAndError:
						score += score_mm + score_ce;
						break;
					case Default:
						score += score_m + score_cm;
						break;
					default:
						fprintf(stderr, "r.readOneType[%d]=%d\n", j, r.readOneType[j]);
						PrintError(FnName,
								"r.readOneType[j]",
								"Could not recognize type",
								Exit,
								OutOfRange);
				}
				prev = r.readOneType[j];
			}
		}
		if(0 < indelLength && 0 == wasInsertion) {
			/* Add in deletion */
			score += sm.gapOpenPenalty;
			score += (r.indelLength-1)*sm.gapExtensionPenalty;
		}
		m.readName = SimReadGetName(&r);
		sprintf(string, "_score=%d", score);
		strcat((char*)m.readName, string);
		/*
		   assert(NULL==m.readName);
		   m.readName = malloc(sizeof(int8_t)*(SEQUENCE_LENGTH+1));
		   if(NULL == m.readName) {
		   PrintError(FnName,
		   "m.readName",
		   "Could not allocate memory",
		   Exit,
		   MallocMemory);
		   }
		   assert(0 <= sprintf((char*)m.readName, ">%d", score));
		   */
		m.readNameLength = strlen((char*)m.readName);
		m.ends[0].numEntries = 1;
		m.ends[0].readLength = (int)strlen(r.readOne);
		assert(r.readLength > 0);
		m.ends[0].read = malloc(sizeof(int8_t)*(m.ends[0].readLength+1));
		if(NULL==m.ends[0].read) {
			PrintError(FnName,
					"m.ends[0].read",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		assert(m.ends[0].readLength > 0);
		strcpy((char*)m.ends[0].read, r.readOne); 
		if(ColorSpace == space) {
			m.ends[0].qualLength = m.ends[0].readLength-1;
		}
		else {
			m.ends[0].qualLength = m.ends[0].readLength;
		}
		m.ends[0].qual = malloc(sizeof(char)*(m.ends[0].qualLength + 1));
		if(NULL==m.ends[0].qual) {
			PrintError(FnName,
					"m.ends[0].qual",
					"Could not allocate memory",
					Exit,
					MallocMemory);
		}
		for(j=0;j<m.ends[0].qualLength;j++) {
			m.ends[0].qual[j] = 'I';
		}
		m.ends[0].maxReached = 0;
		m.ends[0].numEntries = 1;
		m.ends[0].contigs = malloc(sizeof(uint32_t));
		assert(NULL != m.ends[0].contigs);
		m.ends[0].positions = malloc(sizeof(int32_t));
		assert(NULL != m.ends[0].positions);
		m.ends[0].strands= malloc(sizeof(int8_t));
		assert(NULL != m.ends[0].strands);
		m.ends[0].contigs[0] = r.contig;
		m.ends[0].positions[0] = r.pos;
		m.ends[0].strands[0] = r.strand;

		/* Output */
		RGMatchesPrint(matchesFP,
				&m);

		/* Clean up */
		SimReadDelete(&r);
		RGMatchesFree(&m);
	}
	fprintf(stderr, "\r%d\n", numReads);
	fprintf(stderr, "%s", BREAK_LINE);

	/* Re-initialize */
	CloseTmpGZFile(&matchesFP,
			&matchesFileName,
			0);
	if(!(matchesFP=gzopen(matchesFileName, "rb"))) {
		PrintError(FnName,
				matchesFileName,
				"Could not re-open file for reading",
				Exit,
				OpenFileError);
	}

	/* Run "RunDynamicProgramming" from balign */
	fprintf(stderr, "%s", BREAK_LINE);
	fprintf(stderr, "Running local alignment.\n");
	RunDynamicProgramming(matchesFP,
			rg,
			scoringMatrixFileName,
			alignmentType,
			AllAlignments,
			space,
			0,
			INT_MAX,
			readLength, 
			INT_MAX,
			AVG_MISMATCH_QUALITY,
			numThreads,
			100000,
			0,
			0,
			0,
			0,
			tmpDir,
			alignFP,
			notAlignedFP,
			&totalAlignTime,
			&totalFileHandlingTime);
	fprintf(stderr, "%s", BREAK_LINE);

	/* Re-initialize */
	CloseTmpGZFile(&alignFP,
			&alignFileName,
			0);
	if(!(alignFP=gzopen(alignFileName, "rb"))) {
		PrintError(FnName,
				alignFileName,
				"Could not re-open file for reading",
				Exit,
				OpenFileError);
	}

	/* Read in output and sum up accuracy */
	fprintf(stderr, "%s", BREAK_LINE);
	fprintf(stderr, "Summing up totals.\n");
	AlignedReadInitialize(&a);
	numScoreLessThan = numScoreEqual = numScoreGreaterThan = 0;
	while(EOF != AlignedReadRead(&a,
				alignFP)) {
		/* Get substring */
		s = strstr(a.readName, "score=");
		if(NULL == s) {
			PrintError(FnName,
					"a.readName",
					"Could not find \"score=\"",
					Exit,
					OutOfRange);
		}
		/* Extract score */
		ret = sscanf(s, "score=%d", &score);
		if(ret != 1) {
			fprintf(stderr, "ret=%d\nscore=%d\n", ret, score);
			fprintf(stderr, "a.readName=%s\n", a.readName);
			PrintError(FnName,
					"a.readName",
					"Could not parse read name",
					Exit,
					OutOfRange);
		}

		if(round(a.ends[0].entries[0].score) < score) {
			numScoreLessThan++;
			if(FullAlignment == alignmentType) {
				fprintf(stderr, "a.readName=%s\n", a.readName);
				fprintf(stderr, "found=%d\nexpected=%d\n",
						round(a.ends[0].entries[0].score),
						score);
				AlignedReadPrintText(&a, stderr);
				PrintError(FnName,
						"numScoreLessThan",
						"The alignment score should not be less than expected",
						Exit,
						OutOfRange);
			}
		}
		else if(score < round(a.ends[0].entries[0].score)) {
			numScoreGreaterThan++;
			/*
			   PrintError(FnName,
			   "numScoreGreaterThan",
			   "The alignment score was greater than expected",
			   Exit,
			   OutOfRange);
			   */
		}
		else {
			numScoreEqual++;
		}
		/* Free */
		AlignedReadFree(&a);
	}
	fprintf(stderr, "%s", BREAK_LINE);

	/* Close matches file */
	CloseTmpGZFile(&matchesFP, &matchesFileName, 1);
	CloseTmpGZFile(&alignFP, &alignFileName, 1);
	CloseTmpGZFile(&notAlignedFP, &notAlignedFileName, 1);

	fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
			numReads,
			numScoreLessThan,
			numScoreEqual,
			numScoreGreaterThan,
			numReads - numScoreLessThan - numScoreEqual - numScoreGreaterThan,
			numSNPs,
			numErrors,
			deletionLength,
			insertionLength,
			totalAlignTime
		   );
}
