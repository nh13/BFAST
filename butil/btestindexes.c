#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#define BREAK_LINE "****************************************************************\n"
#define GENERATE_HELPER_ROTATE_NUM 1000000
#define RANDOM_ROTATE_NUM 100
#define ROTATE_NUM 10000

#define NUM_MISMATCHES_START 1
#define NUM_MISMATCHES_END 6 
#define NUM_INSERTIONS_START 1
#define NUM_INSERTIONS_END 1
#define NUM_DELETIONS_START 1
#define NUM_DELETIONS_END 1

enum {ReadFromFile, BruteForce, Sample, ProgramParameters}; 
char Algorithm[4][2048] = {"Read From File", "Brute Force", "Sample", "Print Program Parameters"}; 
enum {NO_EVENT, MISMATCH, INSERTION, DELETION};

typedef struct {
	int64_t numIndexes;
	char **indexes;
	int *indexLengths;
} Indexes;

typedef struct {
	int algorithm;
	char indexesFileName[2048];
	int readLength;
	int indexLength;
	int numIndexes;
	int numEventsToSample;
	int numIndexesToSample;
	int numMismatchesStart;
	int numMismatchesEnd;
	int insertionLengthStart;
	int insertionLengthEnd;
	int numErrors;
} arguments;

/*********************************************************************/
/* 							Function Definitions 					  */
/*********************************************************************/
int64_t NChooseR(int64_t N, int64_t R);
void PrintIndexes(Indexes *indexes, FILE *fp);
void InitializeIndexes(Indexes *indexes); 
void AllocateIndexes(Indexes *indexes,
		int numIndexes,
		int indexLength);
void DeleteIndexes(Indexes *indexes);
void GenerateRandomIndex(char *index,
		int *length,
		int readLength,
		int indexLength);
void GenerateRandomIndexes(Indexes *indexes, 
		int readLength,
		int indexLength,
		int numMainIndexes,
		int numIndexes);
void GenerateIndexesHelper(Indexes *indexes,
		FILE *fp,
		int readLength,
		int indexLength,
		char *curIndex,
		int curIndexPos,
		int numOnes);
void GenerateIndexes(Indexes *indexes,
		FILE *fp,
		int readLength,
		int indexLength);
int ReadIndexes(Indexes *indexes, FILE *fp, char *indexesFileName, int readLength);
int CheckReadAgainstIndexes(Indexes *indexes,
		int *read,
		int readLength);
void GenerateRandomInsertionRead(int *read,
		int readLength,
		int insertionLength,
		int numErrors);
void RunInsertion(Indexes *indexes,
		int *read,
		int readLength,
		int insertionLength,
		int numErrors,
		int *totalGenerated,
		int *totalMatched);
void GenerateRandomMismatchRead(int *read,
		int readLength,
		int numMismatches);
void RunMismatches(Indexes *indexes,
		int *read,
		int readLength,
		int readPos,
		int numMismatchesLeft,
		int *totalGenerated,
		int *totalMatched);
void RunSampling(Indexes *mainIndexes,
		FILE *fp,
		int readLength,
		int indexLength,
		int numIndexes,
		int numIndexesToSample,
		int numEventsToSample,
		int numMismatchesStart,
		int numMismatchesEnd,
		int insertionLengthStart,
		int insertionLengthEnd,
		int numDeletionsStart,
		int numDeletionsEnd,
		int numErrors);
void FindBestIndexes(Indexes *indexes,
		FILE *fp,
		int readLength,
		int indexLength,
		int numIndexes,
		int numMismatchesStart,
		int numMismatchesEnd,
		int insertionLengthStart,
		int insertionLengthEnd,
		int numDeletionsStart,
		int numDeletionsEnd,
		int numErrors);
void ComputeAccuracyForEachIndex(Indexes *indexes,
		FILE *fp,
		int readLength,
		int numMismatchesStart,
		int numMismatchesEnd,
		int insertionLengthStart,
		int insertionLengthEnd,
		int numDeletionsStart,
		int numDeletionsEnd,
		int numErrors);
void PrintUsage();
void PrintProgramParmeters(arguments *args);
void AssignDefaultValues(arguments *args); 
void ValidateArguments(arguments *args);
void ParseCommandLineArguments(int argc, char *argv[], arguments *args); 

/*********************************************************************/
/* 							Implementation 							 */
/*********************************************************************/
int64_t NChooseR(int64_t N, int64_t R) 
{
	int64_t i, j, sum;
	int64_t **table=NULL;

	if(N-R < R) {
		R = N - R;
	}

	/* Allocate memory for the dynamic programming table */
	table = malloc(sizeof(int64_t*)*(N+1));
	if(NULL == table) {
		fprintf(stderr, "Error.  Could not allocate memory for table.  Terminating!\n");
		exit(1);
	}
	for(i=0;i<=N;i++) {
		table[i]= malloc(sizeof(int64_t)*(R+1));
		if(NULL == table[i]) {
			fprintf(stderr, "Error.  Could not allocate memory for table[i].  Terminating!\n");
			exit(1);
		}
	}

	/* Initialize all to zero */
	for(i=0;i<=N;i++) {
		for(j=0;j<=R;j++) {
			table[i][j] = 0;
		}
	}

	/* Initialize: N choose 0 is always 1 */
	for(i=0;i<=N;i++) {
		table[i][0] = 1;
	}
	/* Initialize: i choose i is always 1 */
	for(i=0;i<=N && i<=R;i++) {
		table[i][i] = 1;
	}

	for(j=1;j<=R;j++) {
		for(i=j+1;i<=N;i++) {
			table[i][j] = table[i-1][j-1] + table[i-1][j];
		}
	}

	/*
	   fprintf(stderr, "%lld\t%lld\n",
	   (long long int)N,
	   (long long int)R);
	   for(i=0;i<=N;i++) {
	   for(j=0;j<=R;j++) {
	   fprintf(stderr, "%lld ",
	   (long long int)table[i][j]);
	   }
	   fprintf(stderr, "\n");
	   }
	   fprintf(stderr, "\n");
	   */

	/* Get the binomial coefficient */
	sum = table[N][R];

	/* Free memory */
	for(i=0;i<=N;i++) {
		free(table[i]);
		table[i] = NULL;
	}
	free(table);
	table = NULL;

	return sum;
}

void PrintIndexes(Indexes *indexes, FILE *fp)
{
	int i, j;

	fprintf(fp, "%lld\n", (long long int)indexes->numIndexes);
	for(i=0;i<indexes->numIndexes;i++) {
		fprintf(fp, "%d ", indexes->indexLengths[i]);
		for(j=0;j<indexes->indexLengths[i];j++) {
			fprintf(fp, "%d ", indexes->indexes[i][j]);
		}
		fprintf(fp, "\n");
	}
	fflush(fp);
}

void InitializeIndexes(Indexes *indexes) 
{
	indexes->numIndexes = 0;
	indexes->indexes = NULL;
	indexes->indexLengths = NULL;
}

void AllocateIndexes(Indexes *indexes,
		int numIndexes,
		int indexLength)
{
	int i;

	indexes->numIndexes = numIndexes;
	assert(indexes->numIndexes > 0);
	indexes->indexes = malloc(sizeof(char*)*indexes->numIndexes);
	if(NULL == indexes->indexes) {
		fprintf(stderr, "Error.  Could not allocate memory for indexes->indexes.  Terminating!\n");
		exit(1);
	}
	indexes->indexLengths = malloc(sizeof(int)*indexes->numIndexes);
	if(NULL == indexes->indexLengths) {
		fprintf(stderr, "Error.  Could not allocate memory for indexes->indexLengths.  Terminating!\n");
		exit(1);
	}
	for(i=0;i<indexes->numIndexes;i++) {
		indexes->indexLengths[i] = indexLength;
		assert(indexes->indexLengths[i] > 0);
		indexes->indexes[i] = malloc(sizeof(char)*indexes->indexLengths[i]);
		if(NULL == indexes->indexes) {
			fprintf(stderr, "Error.  Could not allocate memory for indexes->indexes.  Terminating!\n");
			exit(1);
		}
	}
}

void DeleteIndexes(Indexes *indexes)
{
	int i;
	for(i=0;i<indexes->numIndexes;i++) {
		free(indexes->indexes[i]);
		indexes->indexes[i] = NULL;
	}
	free(indexes->indexes);
	indexes->indexes = NULL;
	free(indexes->indexLengths);
	indexes->indexLengths = NULL;
	indexes->numIndexes = 0;
}

void GenerateRandomIndex(char *index,
		int *length,
		int readLength,
		int indexLength)
{
	int i, j, k;
	int numLeft;
	int *bins=NULL;
	int numBins = indexLength;

	/* Allocate memory for the bins */
	bins = malloc(sizeof(int)*numBins);
	if(NULL == bins) {
		fprintf(stderr, "Error.  Could not allocate memory for bins.  Terminating!\n");
		exit(1);
	}
	for(i=0;i<numBins;i++) {
		bins[i] = 0;
	}

	/* Choose a number of zeros to insert into the bins */
	numLeft = rand()%(readLength - indexLength);
	assert(numLeft >=0 && numLeft <= readLength - indexLength);

	/* Insert into bins */
	while(numLeft > 0) {
		/* choose a bin between 1 and indexLength-1 */
		i = (rand()%numBins); /* Note: this is not truly inform, but a good approximation */
		assert(i>=0 && i<=numBins-1);
		bins[i]++;
		numLeft--;
	}

	/* Generate index based off the bins */
	/* First base is always a 1 */ 
	for(i=0, j=1, index[0] = 1;i<indexLength-1;i++, j++) {
		/* Insert zero based on the bin size */
		assert(i>=0 && i<=numBins-1);
		for(k=0;k<bins[i];k++, j++) {
			index[j] = 0;
		}
		/* Add a one */
		index[j] = 1;
	}
	(*length) = j;

	/* Free memory */
	free(bins);
	bins=NULL;
}

void GenerateRandomIndexes(Indexes *indexes, 
		int readLength,
		int indexLength,
		int numMainIndexes,
		int numIndexes)
{
	int i;

	assert(numIndexes + numMainIndexes == indexes->numIndexes);
	/* For each index */
	for(i=numMainIndexes;i<indexes->numIndexes;i++) {
		GenerateRandomIndex(indexes->indexes[i],
				&indexes->indexLengths[i],
				readLength,
				indexLength);
	}
}

void GenerateIndexesHelper(Indexes *indexes,
		FILE *fp,
		int readLength,
		int indexLength,
		char *curIndex,
		int curIndexPos,
		int numOnes)
{
	int i;

	if(indexes->numIndexes%GENERATE_HELPER_ROTATE_NUM == 0) {
		fprintf(fp, "\r%lld",
				(long long int)indexes->numIndexes);
	}
	/*
	   fprintf(stderr, "readLength:%d\tindexLength:%d\tcurIndexPos:%d\tnumOnes:%d\n",
	   readLength,
	   indexLength,
	   curIndexPos,
	   numOnes);
	   */

	if(curIndexPos > readLength) {
		return;
	}
	else if(numOnes == indexLength) {
		/* Allocate memory */
		indexes->numIndexes++;
		indexes->indexes = realloc(indexes->indexes, sizeof(char*)*indexes->numIndexes);
		if(NULL == indexes->indexes) {
			fprintf(stderr, "Error.  Could not reallocate memory for indexes->indexes.  Terminating!\n");
			exit(1);
		}
		indexes->indexLengths = realloc(indexes->indexLengths, sizeof(int)*indexes->numIndexes);
		if(NULL == indexes->indexLengths) {
			fprintf(stderr, "Error.  Could not reallocate memory for indexes->indexLengths.  Terminating!\n");
			exit(1);
		}
		indexes->indexLengths[indexes->numIndexes-1] = curIndexPos;
		assert(indexes->indexLengths[indexes->numIndexes-1] > 0);
		indexes->indexes[indexes->numIndexes-1] = malloc(sizeof(char)*indexes->indexLengths[indexes->numIndexes-1]);
		/* Copy over */
		for(i=0;i<indexes->indexLengths[indexes->numIndexes-1];i++) {
			indexes->indexes[indexes->numIndexes-1][i] = curIndex[i];
		}
	}
	else {
		/* Do not add a zero as the first entry for the index */
		if(numOnes > 0 &&
				readLength - (curIndexPos + 1) >= indexLength - numOnes) {
			/* Try zero*/
			curIndex[curIndexPos] = 0; 
			GenerateIndexesHelper(indexes,
					fp,
					readLength,
					indexLength,
					curIndex,
					curIndexPos+1,
					numOnes);
		}
		/* Try one */
		curIndex[curIndexPos] = 1; 
		GenerateIndexesHelper(indexes,
				fp,
				readLength,
				indexLength,
				curIndex,
				curIndexPos+1,
				numOnes+1);
		/* Print percentage complete */
	}
}

void GenerateIndexes(Indexes *indexes,
		FILE *fp,
		int readLength,
		int indexLength)
{
	char *curIndex;

	curIndex = malloc(sizeof(char)*readLength);
	if(NULL == curIndex) {
		fprintf(stderr, "Error.  Could not allocate memory for curIndex.  Terminating!\n");
		exit(1);
	}

	fprintf(fp, "Generating %lld indexes...\n0",
			(long long int)NChooseR(readLength-1, indexLength-1));
	fflush(fp);
	GenerateIndexesHelper(indexes,
			fp,
			readLength,
			indexLength,
			curIndex,
			0,
			0);
	fprintf(fp, "\r%lld\nGenerated %lld indexes.\n",
			(long long int)indexes->numIndexes,
			(long long int)indexes->numIndexes);
	fflush(fp);

	free(curIndex);
	curIndex=NULL;
}

int ReadIndexes(Indexes *indexes, FILE *fp, char *indexesFileName, int readLength)
{
	int i, j, tempInt;

	if(0 > fscanf(fp, "%d", &tempInt)) {
		return EOF;
	}
	indexes->numIndexes = tempInt;
	assert(indexes->numIndexes > 0);
	indexes->indexes = malloc(sizeof(char*)*indexes->numIndexes);
	if(NULL == indexes->indexes) {
		fprintf(stderr, "Error.  Could not allocate memory for indexes->indexes.  Terminating!\n");
		exit(1);
	}
	indexes->indexLengths = malloc(sizeof(int)*indexes->numIndexes);
	if(NULL == indexes->indexLengths) {
		fprintf(stderr, "Error.  Could not allocate memory for indexes->indexLengths.  Terminating!\n");
		exit(1);
	}
	for(i=0;i<indexes->numIndexes;i++) {
		if(0 > fscanf(fp, "%d", &indexes->indexLengths[i])) {
			fprintf(stderr, "Error.  Could not read in the number of indexes from %s.  Terminating!\n",
					indexesFileName);
			exit(1);
		}
		assert(indexes->indexLengths[i] > 0);
		if(indexes->indexLengths[i] > readLength) {
			fprintf(stderr, "Error.  Length for index %d was too big [%d>%d].  Terminating!\n",
					i+1,
					indexes->indexLengths[i],
					readLength);
			exit(1);
		}
		indexes->indexes[i] = malloc(sizeof(char)*indexes->indexLengths[i]);
		if(NULL == indexes->indexes) {
			fprintf(stderr, "Error.  Could not allocate memory for indexes->indexes.  Terminating!\n");
			exit(1);
		}
		for(j=0;j<indexes->indexLengths[i];j++) {
			if(0 > fscanf(fp, "%d", &tempInt)) {
				fprintf(stderr, "Error.  Could not read in value for index %d and entry %d from %s.  Terminating!\n",
						i,
						j,
						indexesFileName);
				exit(1);
			}
			indexes->indexes[i][j] = tempInt;
		}
	}
	return 1;
}

int CheckReadAgainstIndexes(Indexes *indexes,
		int *read,
		int readLength)
{
	int i;
	int offset;
	int curIndexPos, curReadPos;
	int pos;
	int found = 0;

	/*
	   fprintf(stderr, "%s", BREAK_LINE);
	   fprintf(stderr, "index:\n");
	   PrintIndexes(indexes, stderr);
	   fprintf(stderr, "read:\n");
	   for(i=0;i<readLength;i++) {
	   fprintf(stderr, "%d ", read[i]);
	   }
	   fprintf(stderr, "\n");
	   */

	for(i=0;i<indexes->numIndexes;i++) {
		/* Check all offsets of the read */
		for(offset=0;offset<=readLength-indexes->indexLengths[i];offset++) {
			/* Check against index */
			curReadPos = offset;
			curIndexPos = 0;
			pos = 0;
			found = 1;

			/*
			   fprintf(stderr, "i:%d\toffset:%d\n", i, offset);
			   */

			while(found == 1 && 
					curReadPos < readLength &&
					curIndexPos < indexes->indexLengths[i]) {
				/*
				   fprintf(stderr, "1:\tcurIndexPos:%d\tpos:%d\tread[%d]:%d\n",
				   curIndexPos,
				   pos,
				   curReadPos,
				   read[curReadPos]);
				   */
				/* Update effective position */
				while(curReadPos < readLength &&
						(INSERTION == read[curReadPos] || DELETION == read[curReadPos])) {
					switch(read[curReadPos]) {
						case INSERTION:
							pos++;
							/* Do not update position */
							break;
						case DELETION:
							/* Jump ahead one */
							pos+=2;
							break;
						default:
							fprintf(stderr, "Could not understand read [%d].  Terminating!\n", read[curReadPos]);
							exit(1);
							break;
					}
					curReadPos++;
				}

				/*
				   fprintf(stderr, "2:\tcurIndexPos:%d\tpos:%d\tread[%d]:%d\n",
				   curIndexPos,
				   pos,
				   curReadPos,
				   read[curReadPos]);
				   */
				/* Check against index */
				if((1 == indexes->indexes[i][curIndexPos]) &&
						((curReadPos < readLength && read[curReadPos] == MISMATCH) ||
						 (curReadPos < readLength && read[curReadPos] == INSERTION))){
					found = 0;
				}
				else if(pos != curIndexPos) {
					found = 0;
				}
				/* Update positions */
				pos++;
				curIndexPos++;
				curReadPos++;
			}
			if(found == 1) {
				/*
				   fprintf(stderr, "HERE found\n");
				   */
				return 1;
			}
			assert(found == 0);
		}
	}

	return 0;
}
void GenerateRandomInsertionRead(int *read,
		int readLength,
		int insertionLength,
		int numErrors)
{
	int i;
	int startPos;

	/* Initialize */
	for(i=0;i<readLength;i++) {
		read[i] = NO_EVENT;
	}
	/* Insert errors anywhere (including the insertion).  Errors are modeled
	 * as mismatches */
	GenerateRandomMismatchRead(read,
			readLength,
			numErrors);
	/* Choose a random starting position in the read */
	startPos = rand()%(readLength - insertionLength + 1);
	assert(startPos >= 0 && startPos <= readLength - insertionLength);
	/* Insertion */
	for(i=startPos;i<startPos + insertionLength;i++) {
		read[i] = INSERTION;
	}
}

void RunInsertion(Indexes *indexes,
		int *read,
		int readLength,
		int insertionLength,
		int numErrors,
		int *totalGenerated,
		int *totalMatched)
{
	int i, j;

	/* Go through all possible starting positions for the insertion */ 
	for(i=0;i<=readLength - insertionLength;i++) {
		/* Initialize read */
		for(j=0;j<readLength;j++) {
			read[j] = NO_EVENT;
		}
		/* Insert at the position */
		for(j=i;j<i+insertionLength;j++) {
			read[j] = INSERTION;
		}
		/* Run errors as mismatches */
		RunMismatches(indexes,
				read,
				readLength,
				0,
				numErrors,
				totalGenerated,
				totalMatched);
	}
}

void GenerateRandomMismatchRead(int *read,
		int readLength,
		int numMismatches)
{
	int i;
	int numMismatchesLeft;
	int theDefault, theNew;

	/* Optimization: you figure it out */
	if(readLength - numMismatches < numMismatches) {
		numMismatches = readLength - numMismatches;
		theDefault = MISMATCH;
		theNew = NO_EVENT;
	}
	else {
		theDefault = NO_EVENT;
		theNew = MISMATCH;
	}

	/* Initialize read */
	for(i=0;i<readLength;i++) {
		read[i] = theDefault;
	}

	/* Pick bases */
	numMismatchesLeft=numMismatches;

	while(numMismatchesLeft > 0) {
		/* Get a number from 0 to readLength-1 */
		i = rand()%readLength; 
		/* Only modify if it was not previously modified */
		if(read[i] == theDefault) {
			read[i] = theNew;
			numMismatchesLeft--;
		}
	}
}

void RunMismatches(Indexes *indexes,
		int *read,
		int readLength,
		int readPos,
		int numMismatchesLeft,
		int *totalGenerated,
		int *totalMatched)
{
	int prevType;

	if(readPos > readLength) {
		return;
	}
	assert(readPos >= 0 && readPos <= readLength);

	if(0 < numMismatchesLeft) {
		if(readPos < readLength) {
			prevType = read[readPos];

			/* Try adding a mismatch at the current base */
			read[readPos] = MISMATCH;
			RunMismatches(indexes,
					read,
					readLength,
					readPos+1,
					numMismatchesLeft-1,
					totalGenerated,
					totalMatched);

			/* Try not adding a mismatch at the current base */
			read[readPos] = prevType;
			RunMismatches(indexes,
					read,
					readLength,
					readPos+1,
					numMismatchesLeft,
					totalGenerated,
					totalMatched);
			assert(read[readPos] == prevType);
		}
	}
	else {
		assert(0 == numMismatchesLeft);
		/* Increment the number of reads generated */
		(*totalGenerated)++;
		/* Check if the read is matched using the indexes */
		if(1==CheckReadAgainstIndexes(indexes, 
					read, 
					readLength)) {
			(*totalMatched)++;
		}
	}
}

void RunSampling(Indexes *mainIndexes,
		FILE *fp,
		int readLength,
		int indexLength,
		int numIndexes,
		int numIndexesToSample,
		int numEventsToSample,
		int numMismatchesStart,
		int numMismatchesEnd,
		int insertionLengthStart,
		int insertionLengthEnd,
		int numDeletionsStart,
		int numDeletionsEnd,
		int numErrors)
{
	int i, j, k, l, foundOptimal;
	int *curRead=NULL;
	Indexes curIndexes;
	Indexes bestIndexes;
	int curMatched;
	int bestMatched;
	int numEvents;
	int totalMatched, totalGenerated;
	int exact;
	curRead = malloc(sizeof(int)*readLength);
	if(NULL == curRead) {
		fprintf(stderr, "Error.  Could not allocate memory for curRead.  Terminating!\n");
		exit(1);
	}
	for(i=0;i<readLength;i++) {
		curRead[i] = -1;
	}

	InitializeIndexes(&curIndexes);
	AllocateIndexes(&curIndexes,
			numIndexes+mainIndexes->numIndexes,
			readLength);
	InitializeIndexes(&bestIndexes);
	AllocateIndexes(&bestIndexes,
			numIndexes+mainIndexes->numIndexes,
			readLength);

	/* Copy over main indexes */
	for(i=0;i<mainIndexes->numIndexes;i++) {
		curIndexes.indexLengths[i] = mainIndexes->indexLengths[i];
		for(j=0;j<mainIndexes->indexLengths[i];j++) {
			curIndexes.indexes[i][j] =  mainIndexes->indexes[i][j];
		}
	}

	for(i=numMismatchesStart; 0 < i &&  i<=numMismatchesEnd;i++) {
		fprintf(fp, "%s", BREAK_LINE);
		fprintf(fp, "Currently processing %d mismatches and %d errors.\nOut of %lld, currently on:\n0",
				i,
				numErrors,
				(long long int)numIndexesToSample);
		fflush(fp);

		/* Initialize */
		bestMatched = -1;
		foundOptimal = 0;

		/* If the total number of mismatches possible is less than we wish
		 * to sample then just brute-force */
		exact = (NChooseR(readLength, i+numErrors) <= numEventsToSample)?1:0;

		/* Sample indexes */
		numEvents = 0;
		for(j=0;foundOptimal != 1 && j<numIndexesToSample;j++) {
			if(j%RANDOM_ROTATE_NUM == 0) {
				fprintf(fp, "\r%10d",
						j);
				fflush(fp);
			}
			/* Generate random indexes */
			GenerateRandomIndexes(&curIndexes,
					readLength,
					indexLength,
					mainIndexes->numIndexes,
					numIndexes);
			curMatched = 0;
			numEvents = 0;
			/* Initialize read */
			for(k=0;k<readLength;k++) {
				curRead[k] = NO_EVENT;
			}
			/* Sample mismatches */
			if(0==exact) {
				for(k=0;k<numEventsToSample;k++) {
					/* Generate random read with j mismatches */ 
					GenerateRandomMismatchRead(curRead,
							readLength,
							i+numErrors);

					/* Test the read against the indexes */
					switch(CheckReadAgainstIndexes(&curIndexes, curRead, readLength)) {
						case 0:
							/* do nothing */
							break;
						case 1:
							/* increment */
							curMatched++;
							break;
						default:
							/* do nothing */
							break;
					}
					numEvents++;
				}
			}
			else {
				RunMismatches(&curIndexes,
						curRead,
						readLength,
						0,
						i+numErrors,
						&numEvents,
						&curMatched);
			}

			/* Check the results against the best so far */
			if(curMatched > bestMatched) {
				bestMatched = curMatched;
				bestIndexes.numIndexes = curIndexes.numIndexes;

				for(k=0;k<curIndexes.numIndexes;k++) {
					bestIndexes.indexLengths[k] = curIndexes.indexLengths[k];
					for(l=0;l<curIndexes.indexLengths[k];l++) {
						bestIndexes.indexes[k][l] = curIndexes.indexes[k][l];
					}
				}
				if(bestMatched == numEvents) {
					foundOptimal = 1;
				}
			}
		}
		fprintf(fp, "\r%10d\n",
				(int)numIndexesToSample);
		/* Print asymptotic */
		fprintf(fp, "Best index found for %d mismatches and %d errors (%d out of %d [%3.3lf]):\n",
				i,
				numErrors,
				bestMatched,
				numEvents,
				(bestMatched*100.0)/numEvents);
		fflush(fp);
		totalGenerated=0;
		totalMatched=0;
		/* Initialize read */
		for(j=0;j<readLength;j++) {
			curRead[j] = NO_EVENT;
		}
		RunMismatches(&bestIndexes,
				curRead,
				readLength,
				0,
				i+numErrors,
				&totalGenerated,
				&totalMatched);
		/* Print truth for the best one */
		fprintf(fp, "Mismatches[%d,%d]: %d out of %d [%3.3lf].\n",
				i,
				numErrors,
				totalMatched,
				totalGenerated,
				(totalMatched*100.0)/totalGenerated);
		fflush(fp);
		PrintIndexes(&bestIndexes, fp);
	}

	for(i=insertionLengthStart; 0 < insertionLengthStart && i<=insertionLengthEnd;i++) {
		fprintf(fp, "%s", BREAK_LINE);
		fprintf(fp, "Currently processing an insertion of length %d and %d errors.\nOut of %lld, currently on:\n0",
				i,
				numErrors,
				(long long int)numIndexesToSample);
		fflush(fp);

		/* Initialize */
		bestMatched = -1;
		foundOptimal = 0;

		/* If the total number of insertion locations possible is less than we wish
		 * to sample then just brute-force.  Remember to factor in errors. */
		if(numErrors > 0) {
			exact = (NChooseR(readLength, numErrors)*(readLength - i + 1) <= numEventsToSample)?1:0;
		}
		else {
			exact = ((readLength - i + 1) <= numEventsToSample)?1:0;
		}

		/* Sample indexes */
		numEvents = 0;
		for(j=0;foundOptimal != 1 && j<numIndexesToSample;j++) {
			if(j%RANDOM_ROTATE_NUM == 0) {
				fprintf(fp, "\r%10d",
						j);
				fflush(fp);
			}
			/* Generate random indexes */
			GenerateRandomIndexes(&curIndexes,
					readLength,
					indexLength,
					mainIndexes->numIndexes,
					numIndexes);
			curMatched = 0;
			numEvents = 0;
			/* Initialize read */
			for(k=0;k<readLength;k++) {
				curRead[k] = NO_EVENT;
			}
			/* Sample insertions */
			if(0==exact) {
				for(k=0;k<numEventsToSample;k++) {
					/* Generate random read with j mismatches */ 
					GenerateRandomInsertionRead(curRead,
							readLength,
							i,
							numErrors);

					/* Test the read against the indexes */
					switch(CheckReadAgainstIndexes(&curIndexes, curRead, readLength)) {
						case 0:
							/* do nothing */
							break;
						case 1:
							/* increment */
							curMatched++;
							break;
						default:
							/* do nothing */
							break;
					}
					numEvents++;
				}
			}
			else {
				RunInsertion(&curIndexes,
						curRead,
						readLength,
						i,
						numErrors,
						&numEvents,
						&curMatched);
			}

			/* Check the results against the best so far */
			if(curMatched > bestMatched) {
				bestMatched = curMatched;
				bestIndexes.numIndexes = curIndexes.numIndexes;

				for(k=0;k<curIndexes.numIndexes;k++) {
					bestIndexes.indexLengths[k] = curIndexes.indexLengths[k];
					for(l=0;l<curIndexes.indexLengths[k];l++) {
						bestIndexes.indexes[k][l] = curIndexes.indexes[k][l];
					}
				}
				if(bestMatched == numEvents) {
					foundOptimal = 1;
				}
			}
		}
		fprintf(fp, "\r%10d\n",
				(int)numIndexesToSample);
		/* Print asymptotic */
		fprintf(fp, "Best index found insertion of length %d and %d errors (%d out of %d [%3.3lf]):\n",
				i,
				numErrors,
				bestMatched,
				numEvents,
				(bestMatched*100.0)/numEvents);
		fflush(fp);
		totalGenerated=0;
		totalMatched=0;
		RunInsertion(&bestIndexes,
				curRead,
				readLength,
				i,
				numErrors,
				&totalGenerated,
				&totalMatched);
		/* Print truth for the best one */
		fprintf(fp, "Insertion[%d,%d]: %d out of %d [%3.3lf].\n",
				i,
				numErrors,
				totalMatched,
				totalGenerated,
				(totalMatched*100.0)/totalGenerated);
		fflush(fp);
		PrintIndexes(&bestIndexes, fp);
	}

	/* Free memory */
	free(curRead);
	curRead=NULL;
	DeleteIndexes(&curIndexes);
	DeleteIndexes(&bestIndexes);

}

void FindBestIndexes(Indexes *indexes,
		FILE *fp,
		int readLength,
		int indexLength,
		int numIndexes,
		int numMismatchesStart,
		int numMismatchesEnd,
		int insertionLengthStart,
		int insertionLengthEnd,
		int numDeletionsStart,
		int numDeletionsEnd,
		int numErrors)
{
	int i, j, k, foundOptimal;
	int64_t counter;
	int *read=NULL;
	int *curRead=NULL;
	int totalGenerated;
	Indexes curIndexes;
	Indexes bestIndexes;
	int curMatched;
	int bestMatched;
	int *whichIndexes=NULL;

	read = malloc(sizeof(int)*readLength);
	if(NULL == read) {
		fprintf(stderr, "Error.  Could not allocate memory for read.  Terminating!\n");
		exit(1);
	}
	curRead = malloc(sizeof(int)*readLength);
	if(NULL == curRead) {
		fprintf(stderr, "Error.  Could not allocate memory for curRead.  Terminating!\n");
		exit(1);
	}
	for(i=0;i<readLength;i++) {
		read[i] = NO_EVENT;
		curRead[i] = -1;
	}
	whichIndexes = malloc(sizeof(int)*numIndexes);
	if(NULL == whichIndexes) {
		fprintf(stderr, "Error.  Could not allocate memory for whichIndexes.  Terminating!\n");
		exit(1);
	}

	InitializeIndexes(&curIndexes);
	AllocateIndexes(&curIndexes,
			numIndexes,
			readLength);
	InitializeIndexes(&bestIndexes);
	AllocateIndexes(&bestIndexes,
			numIndexes,
			readLength);

	/* For each mismatch */
	for(i=numMismatchesStart;0 < numMismatchesStart && i<=numMismatchesEnd;i++) {
		fprintf(fp, "%s", BREAK_LINE);
		fprintf(fp, "Currently processing %d mismatches and %d errors.\nOut of %lld distinct combinations of %d indexes, currently on:\n0",
				i,
				numErrors,
				(long long int)NChooseR(indexes->numIndexes, numIndexes),
				numIndexes);
		fflush(fp);
		/*
		 *
		 fprintf(fp, "Currently processing %d mismatches.\nOut of %lld distinct combinations of %d indexes, currently on:\n0",
		 i,
		 (long long int)NChooseR(indexes->numIndexes, numIndexes),
		 numIndexes);
		 */
		bestMatched = 0;
		foundOptimal = 0;
		counter = 0;
		/* Initialize start indexes */
		for(j=0;j<numIndexes;j++) {
			whichIndexes[j] = j;
		}
		while(whichIndexes[0] <= indexes->numIndexes - numIndexes &&
				foundOptimal == 0) {
			counter++;
			if(counter%ROTATE_NUM == 0) {
				fprintf(fp, "\r%lld",
						(long long int)counter);
				fflush(fp);
			}
			/*
			   fprintf(stderr, "\r");
			   for(j=0;j<numIndexes;j++) {
			   fprintf(stderr, "%lld ",
			   (long long int)whichIndexes[j]);
			   assert(whichIndexes[j] < indexes->numIndexes);
			   }
			   fprintf(stderr, "\n");
			   */

			/* Copy over index */
			curIndexes.numIndexes = numIndexes;
			for(j=0;j<numIndexes;j++) {
				curIndexes.indexLengths[j] = indexes->indexLengths[whichIndexes[j]];
				for(k=0;k<indexes->indexLengths[whichIndexes[j]];k++) {
					curIndexes.indexes[j][k] = indexes->indexes[whichIndexes[j]][k];
				}
			}

			/* Try all indexes */
			totalGenerated = 0;
			curMatched = 0;
		/* Initialize read */
		for(j=0;j<readLength;j++) {
			read[j] = NO_EVENT;
		}
			RunMismatches(&curIndexes,
					read,
					readLength,
					0,
					i+numErrors,
					&totalGenerated,
					&curMatched);

			/* Check the results against the best so far */
			if(curMatched > bestMatched) {
				bestMatched = curMatched;
				bestIndexes.numIndexes = numIndexes;
				for(j=0;j<numIndexes;j++) {
					bestIndexes.indexLengths[j] = indexes->indexLengths[whichIndexes[j]];
					for(k=0;k<indexes->indexLengths[whichIndexes[j]];k++) {
						bestIndexes.indexes[j][k] = indexes->indexes[whichIndexes[j]][k];
					}
				}
				if(bestMatched == totalGenerated) {
					foundOptimal = 1;
				}
			}

			/* Update */
			whichIndexes[numIndexes-1]++;
			j=numIndexes-1;
			/*
			   fprintf(stderr, "j:%d\t[%d,%d,%d]\n",
			   j,
			   (whichIndexes[j] == indexes->numIndexes),
			   (whichIndexes[0] <= indexes->numIndexes - numIndexes),
			   (0 < j));
			   */
			while(whichIndexes[j] == indexes->numIndexes &&
					whichIndexes[0] <= indexes->numIndexes - numIndexes &&
					0 < j) {
				assert(whichIndexes[j] <= indexes->numIndexes);
				whichIndexes[j-1]++;
				whichIndexes[j] = whichIndexes[j-1] + 1;
				if(whichIndexes[j] < numIndexes ||
						whichIndexes[j-1] == indexes->numIndexes) {
					j--;
				}
				/*
				   fprintf(stderr, "which:\t");
				   for(k=0;k<numIndexes;k++) {
				   fprintf(stderr, "%d ", whichIndexes[k]);
				   }
				   fprintf(stderr, "\n");
				   fprintf(stderr, "j:%d\t[%d,%d,%d]\n",
				   j,
				   (whichIndexes[j] == indexes->numIndexes),
				   (whichIndexes[0] <= indexes->numIndexes - numIndexes),
				   (0 < j));
				   */
			}
			for(k=j+1;k<numIndexes;k++) {
				whichIndexes[k] = whichIndexes[k-1]+1;
			}
		}
		fprintf(fp, "\nBest index found for %d mismatches and %d errors (%d out of %d [%3.3lf]):\n",
				i,
				numErrors,
				bestMatched,
				totalGenerated,
				(bestMatched*100.0)/totalGenerated);
		fflush(fp);
		PrintIndexes(&bestIndexes, fp);
	}
	if(insertionLengthStart > 0) {
		fprintf(fp, "%s", BREAK_LINE);
		fprintf(fp, "Currently processing insertions of length %d and %d errors.\nOut of %lld distinct combinations of %d indexes, currently on:\n0",
				i,
				numErrors,
				(long long int)NChooseR(indexes->numIndexes, numIndexes),
				numIndexes);
		fflush(fp);
		bestMatched = 0;
		foundOptimal = 0;
		counter = 0;
		/* Initialize start indexes */
		for(j=0;j<numIndexes;j++) {
			whichIndexes[j] = j;
		}
		while(whichIndexes[0] <= indexes->numIndexes - numIndexes &&
				foundOptimal == 0) {
			counter++;
			if(counter%ROTATE_NUM == 0) {
				fprintf(fp, "\r%lld",
						(long long int)counter);
				fflush(fp);
			}
			/*
			   fprintf(stderr, "\r");
			   for(j=0;j<numIndexes;j++) {
			   fprintf(stderr, "%lld ",
			   (long long int)whichIndexes[j]);
			   assert(whichIndexes[j] < indexes->numIndexes);
			   }
			   fprintf(stderr, "\n");
			   */

			/* Copy over index */
			curIndexes.numIndexes = numIndexes;
			for(j=0;j<numIndexes;j++) {
				curIndexes.indexLengths[j] = indexes->indexLengths[whichIndexes[j]];
				for(k=0;k<indexes->indexLengths[whichIndexes[j]];k++) {
					curIndexes.indexes[j][k] = indexes->indexes[whichIndexes[j]][k];
				}
			}

			/* Try all indexes */
			totalGenerated = 0;
			curMatched = 0;
			RunInsertion(&curIndexes,
					read,
					readLength,
					i,
					numErrors,
					&totalGenerated,
					&curMatched);

			/* Check the results against the best so far */
			if(curMatched > bestMatched) {
				bestMatched = curMatched;
				bestIndexes.numIndexes = numIndexes;
				for(j=0;j<numIndexes;j++) {
					bestIndexes.indexLengths[j] = indexes->indexLengths[whichIndexes[j]];
					for(k=0;k<indexes->indexLengths[whichIndexes[j]];k++) {
						bestIndexes.indexes[j][k] = indexes->indexes[whichIndexes[j]][k];
					}
				}
				if(bestMatched == totalGenerated) {
					foundOptimal = 1;
				}
			}

			/* Update */
			whichIndexes[numIndexes-1]++;
			j=numIndexes-1;
			/*
			   fprintf(stderr, "j:%d\t[%d,%d,%d]\n",
			   j,
			   (whichIndexes[j] == indexes->numIndexes),
			   (whichIndexes[0] <= indexes->numIndexes - numIndexes),
			   (0 < j));
			   */
			while(whichIndexes[j] == indexes->numIndexes &&
					whichIndexes[0] <= indexes->numIndexes - numIndexes &&
					0 < j) {
				assert(whichIndexes[j] <= indexes->numIndexes);
				whichIndexes[j-1]++;
				whichIndexes[j] = whichIndexes[j-1] + 1;
				if(whichIndexes[j] < numIndexes ||
						whichIndexes[j-1] == indexes->numIndexes) {
					j--;
				}
				/*
				   fprintf(stderr, "which:\t");
				   for(k=0;k<numIndexes;k++) {
				   fprintf(stderr, "%d ", whichIndexes[k]);
				   }
				   fprintf(stderr, "\n");
				   fprintf(stderr, "j:%d\t[%d,%d,%d]\n",
				   j,
				   (whichIndexes[j] == indexes->numIndexes),
				   (whichIndexes[0] <= indexes->numIndexes - numIndexes),
				   (0 < j));
				   */
			}
			for(k=j+1;k<numIndexes;k++) {
				whichIndexes[k] = whichIndexes[k-1]+1;
			}
		}
		fprintf(fp, "\nBest index found for insertion of length %d and %d errors (%d out of %d [%3.3lf]):\n",
				i,
				numErrors,
				bestMatched,
				totalGenerated,
				(bestMatched*100.0)/totalGenerated);
		fflush(fp);
		PrintIndexes(&bestIndexes, fp);
	}
	/*
	   for(i=insertionLengthStart;i<=insertionLengthEnd;i++) {
	   }
	   for(i=numDeletionsStart;i<=numDeletionsEnd;i++) {
	   }
	   */

	free(read);
	read = NULL;
	free(curRead);
	read = NULL;
	free(whichIndexes);
	whichIndexes=NULL;
	DeleteIndexes(&curIndexes);
	DeleteIndexes(&bestIndexes);
}

void ComputeAccuracyForEachIndex(Indexes *indexes,
		FILE *fp,
		int readLength,
		int numMismatchesStart,
		int numMismatchesEnd,
		int insertionLengthStart,
		int insertionLengthEnd,
		int numDeletionsStart,
		int numDeletionsEnd,
		int numErrors)
{
	int i, j;
	int *read=NULL;
	int totalGenerated;
	int totalMatched;

	read = malloc(sizeof(int)*readLength);
	if(NULL == read) {
		fprintf(stderr, "Error.  Could not allocate memory for read.  Terminating!\n");
		exit(1);
	}
	for(i=0;i<readLength;i++) {
		read[i] = NO_EVENT;
	}

	fprintf(fp, "%s", BREAK_LINE);
	fflush(fp);
	PrintIndexes(indexes, fp);
	for(i=numMismatchesStart;0 < numMismatchesStart && i<=numMismatchesEnd;i++) {
		totalGenerated=0;
		totalMatched=0;
		/* Initialize read */
		for(j=0;j<readLength;j++) {
			read[j] = NO_EVENT;
		}
		RunMismatches(indexes,
				read,
				readLength,
				0,
				i+numErrors,
				&totalGenerated,
				&totalMatched);
		fprintf(fp, "Mismatches[%d,%d]: %d out of %d [%3.3lf].\n",
				i,
				numErrors,
				totalMatched,
				totalGenerated,
				(totalMatched*100.0)/totalGenerated);
		fflush(fp);
	}

	for(i=insertionLengthStart;0 < insertionLengthStart && i<=insertionLengthEnd;i++) {
		totalGenerated=0;
		totalMatched=0;
		RunInsertion(indexes,
				read,
				readLength,
				i,
				numErrors,
				&totalGenerated,
				&totalMatched);
		fprintf(fp, "Insertion[%d,%d]: %d out of %d [%3.3lf].\n",
				i,
				numErrors,
				totalMatched,
				totalGenerated,
				(totalMatched*100.0)/totalGenerated);
		fflush(fp);
	}
	/*
	   for(i=numDeletionsStart;i<=numDeletionsEnd;i++) {
	   }
	   */

	free(read);
	read = NULL;
}

void PrintUsage()
{
	fprintf(stderr, "Usage: test.indexes [OPTIONS]...\n");
	fprintf(stderr, "******************************* Algorithm Options (no defaults) *******************************\n");
	fprintf(stderr, "\t-a\tINT\talgorithm\n\t\t\t\t0: read from file\n\t\t\t\t1: brute-force search all indexes\n\t\t\t\t2: sample indexes and events\n");
	fprintf(stderr, "\t-f\tSTRING\tinput file name (for -a 0)\n");
	fprintf(stderr, "\t-r\tINT\tread length (for all) \n");
	fprintf(stderr, "\t-l\tINT\tindex length (for -a 1 and -a 2)\n");
	fprintf(stderr, "\t-n\tINT\tnumber of indexes in a set of indexes (for -a 1 and -a 2)\n");
	fprintf(stderr, "\t-s\tINT\tnumber of indexes to sample (for -a 2)\n");
	fprintf(stderr, "\t-S\tINT\tnumber of events to sample (for -a 2)\n");
	fprintf(stderr, "******************************* Event Options (default =0 ) ***********************************\n");
	fprintf(stderr, "\t-m\tINT\tminimum number of mismatches\n");
	fprintf(stderr, "\t-M\tINT\tmaximum number of mismatches\n");
	fprintf(stderr, "\t-i\tINT\tminimum insertion length\n");
	fprintf(stderr, "\t-I\tINT\tmaximum insertion length\n");
	fprintf(stderr, "\t-e\tINT\tnumber of errors\n");
	fprintf(stderr, "******************************* Miscellaneous Options  ****************************************\n");
	fprintf(stderr, "\t-p\tNULL\tprints the program parameters\n");
	fprintf(stderr, "\t-h\tNULL\tprints this message\n");
}

void PrintProgramParameters(arguments *args)
{
	/* Print program parameters */
	fprintf(stdout, "%s", BREAK_LINE);
	fprintf(stdout, "Printing program parameters:\n");
	fprintf(stdout, "algorithm:\t\t\t%d\t[%s]\n", args->algorithm, Algorithm[args->algorithm]); 
	fprintf(stdout, "input file name:\t\t%s\n", args->indexesFileName);
	fprintf(stdout, "read length:\t\t\t%d\n", args->readLength);
	fprintf(stdout, "index length:\t\t\t%d\n", args->indexLength);
	fprintf(stdout, "number of indexes in a set:\t%d\n", args->numIndexes);
	fprintf(stdout, "number of indexes to sample:\t%d\n", args->numIndexesToSample);
	fprintf(stdout, "number of events to sample:\t%d\n", args->numEventsToSample);
	fprintf(stdout, "number of mismatches:\t\tfrom %d to %d\n", args->numMismatchesStart, args->numMismatchesEnd);
	fprintf(stdout, "insertion length:\t\tfrom %d to %d\n", args->insertionLengthStart, args->insertionLengthEnd);
	fprintf(stdout, "number of errors:\t\t%d\n", args->numErrors);
	fprintf(stdout, "%s", BREAK_LINE);
}

void AssignDefaultValues(arguments *args) 
{
	args->algorithm=0;
	strcpy(args->indexesFileName, "\0");
	args->readLength=0;
	args->indexLength=0;
	args->numIndexes=0;
	args->numEventsToSample=0;
	args->numIndexesToSample=0;
	args->numMismatchesStart=0;
	args->numMismatchesEnd=0;
	args->insertionLengthStart=0;
	args->insertionLengthEnd=0;
	args->numErrors=0;
}

void ValidateArguments(arguments *args)
{
	/* Check command line arguments */
	if(args->readLength < 0 ||
			args->numIndexes < 0 ||
			args->numMismatchesStart < 0 ||
			args->numMismatchesEnd < args->numMismatchesStart ||
			args->insertionLengthStart < 0 ||
			args->insertionLengthEnd < args->insertionLengthStart ||
			args->numErrors < 0) {
		fprintf(stderr, "*** Error.  Input parameters not valid (all must be positive with read length greater than index length).  Teriminating! ***\n");
		exit(1);
	}
	switch(args->algorithm) {
		case ReadFromFile:
			if(strlen(args->indexesFileName) == 0) {
				fprintf(stderr, "*** Error.  Please specify an input file name.  Terminating! ***\n");
				exit(1);
			}
			break;
		case BruteForce:
			if(args->indexLength <= 0 ||
					args->numIndexes <= 0) {
				fprintf(stderr, "*** Error.  Input parameters not valid (all must be positive with read length greater than index length).  Teriminating! ***\n");
				exit(1);
			}
			break;
		case Sample:
			if(args->indexLength <= 0 ||
					args->numIndexes <= 0 ||
					args->numIndexesToSample <= 0 ||
					args->numEventsToSample <= 0) {
				fprintf(stderr, "*** Error.  Input parameters not valid (all must be positive with read length greater than index length).  Teriminating! ***\n");
				exit(1);
			}
			break;
		default:
			fprintf(stderr, "*** Error.  Algorithm option %d is out of range.  Terminating! ***\n",
					args->algorithm);
			exit(1);
	}
}

void ParseCommandLineArguments(int argc, char *argv[], arguments *args) 
{
	int i;
	if(argc==1) {
		PrintUsage();
		exit(1);
	}
	for(i=1;i<argc;i+=2) {
		if(argv[i][0] != '-' ||
				strlen(argv[i]) != 2) {
			fprintf(stderr, "*** Error.  Could not understand command line option %s.  Terminating! ***\n",
					argv[i]);
			exit(1);
		}
		switch(argv[i][1]) {
			case 'a':
				args->algorithm = atoi(argv[i+1]);
				break;
			case 'e':
				args->numErrors = atoi(argv[i+1]);
				break;
			case 'f':
				strcpy(args->indexesFileName, argv[i+1]);
				break;
			case 'i':
				args->insertionLengthStart = atoi(argv[i+1]);
				break;
			case 'I':
				args->insertionLengthEnd = atoi(argv[i+1]);
				break;
			case 'h':
				PrintUsage();
				exit(1);
				break;
			case 'l':
				args->indexLength = atoi(argv[i+1]);
				break;
			case 'm':
				args->numMismatchesStart = atoi(argv[i+1]);
				break;
			case 'M':
				args->numMismatchesEnd = atoi(argv[i+1]);
				break;
			case 'n':
				args->numIndexes = atoi(argv[i+1]);
				break;
			case 'p':
				args->algorithm = ProgramParameters;
				break;
			case 'r':
				args->readLength = atoi(argv[i+1]);
				break;
			case 's':
				args->numIndexesToSample = atoi(argv[i+1]);
				break;
			case 'S':
				args->numEventsToSample = atoi(argv[i+1]);
				break;
			default:
				fprintf(stderr, "*** Error.  Could not understand command line option %s.  Terminating! ***\n",
						argv[i]);
				exit(1);
				break;

		}
	}
}

int main(int argc, char *argv[])
{
	/* Command line arguments */
	arguments args;
	/* Local variables */
	FILE *fp=NULL;
	Indexes indexes;

	/* Assign default values */
	AssignDefaultValues(&args);

	/* Parse command line arguments */
	ParseCommandLineArguments(argc, argv, &args);

	/* Validate command line arguments */
	ValidateArguments(&args);

	/* Print program parameters */
	PrintProgramParameters(&args);

	/* Initialize indexes */
	InitializeIndexes(&indexes);
	switch(args.algorithm) {
		case ReadFromFile:
			/* Read in indexes */
			if(!(fp = fopen(args.indexesFileName, "r"))) {
				fprintf(stderr, "Error.  Could not open %s for reading.  Terminating!\n",
						args.indexesFileName);
				exit(1);
			}
			while(EOF!=ReadIndexes(&indexes, fp, args.indexesFileName, args.readLength)) {

				/* Compute results */
				ComputeAccuracyForEachIndex(&indexes,
						stdout,
						args.readLength,
						args.numMismatchesStart,
						args.numMismatchesEnd,
						args.insertionLengthStart,
						args.insertionLengthEnd,
						NUM_DELETIONS_START,
						NUM_DELETIONS_END,
						args.numErrors);
				/* Delete indexes */
				DeleteIndexes(&indexes);
			}
			fclose(fp);
			break;
		case BruteForce:
			/* Generate Indexes */
			GenerateIndexes(&indexes,
					fp,
					args.readLength,
					args.indexLength);
			/* Find the best indexes */
			FindBestIndexes(&indexes,
					stdout,
					args.readLength,
					args.indexLength,
					args.numIndexes,
					args.numMismatchesStart,
					args.numMismatchesEnd,
					args.insertionLengthStart,
					args.insertionLengthEnd,
					NUM_DELETIONS_START,
					NUM_DELETIONS_END,
					args.numErrors);
			/* Delete indexes */
			DeleteIndexes(&indexes);
			break;
		case Sample:
			if(strlen(args.indexesFileName) > 0) {
				/* Read in indexes */
				if(!(fp = fopen(args.indexesFileName, "r"))) {
					fprintf(stderr, "Error.  Could not open %s for reading.  Terminating!\n",
							args.indexesFileName);
					exit(1);
				}
				assert(EOF!=ReadIndexes(&indexes, fp, args.indexesFileName, args.readLength));
				fclose(fp);
			}
			/* Initialize seed */
			srand(time(NULL));
			/* Sample */
			RunSampling(&indexes,
					stdout,
					args.readLength,
					args.indexLength,
					args.numIndexes,
					args.numIndexesToSample,
					args.numEventsToSample,
					args.numMismatchesStart,
					args.numMismatchesEnd,
					args.insertionLengthStart,
					args.insertionLengthEnd,
					NUM_DELETIONS_START,
					NUM_DELETIONS_END,
					args.numErrors);
			DeleteIndexes(&indexes);
			break;
		case ProgramParameters:
			/* Do nothing */
			break;
		default:
			fprintf(stderr, "Error.  Could not understand program mode [%d].  Terminating!\n", args.algorithm);
			exit(1);

	}
	return 0;
}
