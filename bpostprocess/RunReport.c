#include "RunReport.h"

/* TODO */
void RunReport(char *inputFileName,
		int binaryInput,
		int algorithmReads,
		int minScoreReads,
		int startChr,
		int startPos,
		int endChr,
		int endPos,
		int pairedEnd,
		int algorithmReads,
		int minScoreReadsPaired,
		int minDistancePaired,
		int maxDistancePaired,
		int meanDistancePaired,
		int chrAbPaired,
		int inversionsPaired,
		char *outputID,
		char *outputDir,
		char *tmpDir,
		int outputFormat,
		int timing)
{
	/* Algorithm Overview 
	 * Step 1: Filter reads.  For paired end, filter each read in the pair.
	 * Step 2: Filter paired end reads.  Only for paired end data.
	 * Step 3: Split the data up by chromosome.
	 * Step 4: Split the data up until manageable and output.
	 *
	 * Note: steps 1-3 are performed in serial for each read.
	 * */
}
