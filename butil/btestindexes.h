#ifndef BTESTINDEXES_H_
#define BTESTINDEXES_H_

enum {SearchForIndexes, EvaluateIndexes, ProgramParameters};
char Algorithm[3][2048] = {"Search for indexes", "Evaluate indexes", "Print Program Parameters"};
enum {NO_EVENT, MISMATCH, INSERTION, DELETION};

typedef struct {
	int32_t length;
	int32_t *profile;
} Read;

typedef struct {
	int32_t numReads;
	int32_t *correct;
	int32_t length; /* lenght of correct */
	int32_t numSNPs;
	int32_t numColorErrors;
} AccuracyProfile;

typedef struct {
	int32_t *mask;
	int32_t keySize;
	int32_t keyWidth;
} Index;

typedef struct {
	int32_t numIndexes;
	Index *indexes;
} IndexSet;

typedef struct {
	int algorithm;
	char inputFileName[MAX_FILENAME_LENGTH];
	int readLength;
	int numEventsToSample;
	int numIndexesToSample;
	int keySize;
	int maxKeyWidth;
	int maxIndexSetSize;
	int accuracyThreshold;
	int space;
	int maxNumMismatches;
	int maxInsertionLength;
	int maxNumColorErrors;
} arguments;

enum{NoIndelType, DeletionType, InsertionType};

/* Functions */
void RunSearchForIndexes(int, int, int, int, int, int, int, int, int, int);
void RunEvaluateIndexes(char*, int, int, int, int, int, int);
void RunEvaluateIndexesNTSpace(IndexSet*, int, int, int, int);
void RunEvaluateIndexesColorSpace(IndexSet*, int, int, int, int, int);
int32_t GetNumCorrect(IndexSet*, int, int, int, int, int, int, int);
/* IndexSet functions */
int IndexSetContains(IndexSet*, Index*);
int32_t IndexSetCheckRead(IndexSet*, Read*);
void IndexSetPush(IndexSet*, Index*);
void IndexSetPop(IndexSet*);
void IndexSetSeed(IndexSet*, int);
void IndexSetInitialize(IndexSet*);
void IndexSetFree(IndexSet*);
void IndexSetPrint(IndexSet*, FILE*);
void IndexSetRead(IndexSet*, char*);
/* Index functions */
int IndexCompare(Index*, Index*);
int32_t IndexCheckRead(Index*, Read*);
void IndexCopy(Index*, Index*);
void IndexGetRandom(Index*, int, int);
void IndexAllocate(Index*, int, int);
void IndexInitialize(Index*);
void IndexFree(Index*);
void IndexPrint(Index*, FILE*);
int IndexRead(Index*, FILE*);
/* Accuracy Profile functions */
int AccuracyProfileCompare(AccuracyProfile*, AccuracyProfile*, int);
void AccuracyProfileCopy(AccuracyProfile*, AccuracyProfile*);
void AccuracyProfileUpdate(IndexSet*, AccuracyProfile*, int, int, int, int, int);
void AccuracyProfileAllocate(AccuracyProfile*, int, int);
void AccuracyProfileInitialize(AccuracyProfile*);
void AccuracyProfileFree(AccuracyProfile*);
/* Read functions */
void ReadSplit(Read*, Read*, Read*, int, int);
void ReadGetRandom(Read*, int, int, int, int);
void ReadInitialize(Read*);
void ReadAllocate(Read*, int);
void ReadFree(Read*);
void ReadPrint(Read*, FILE*);
/* Command line functions */
void PrintUsage();
void PrintProgramParmeters(arguments *args);
void AssignDefaultValues(arguments *args);
void ValidateArguments(arguments *args);
void ParseCommandLineArguments(int argc, char *argv[], arguments *args);



#endif
