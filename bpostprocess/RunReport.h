#ifndef RUN_REPORT_H_
#define RUN_REPORT_H_

void RunReport(char *inputFileName,
		int binaryInput,
		int uniqueReads,
		int bestScoreReads,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int pairedEnd,
		int uniqueReadsPaired,
		int bestScoreReadsPaired,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int useMeanDistancePaired,
		int meanDistancePaired,
		int chrAbPaired,
		int inversionsPaired,
		char *outputID,
		char *outputDir,
		char *tmpDir,
		int outputFormat,
		int timing);


#endif
