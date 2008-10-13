#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

enum {BAF, MAF, LastFileType};
#define MIN_FILTER 0
#define MAX_FILTER_SE 3
#define MAX_FILTER_PE 5
enum {NoFiltering, 	/* 0 */
	AllNotFiltered, /* 1 */
	Unique, 		/* 2 */
	BestScore, 		/* 3 */
	MeanUnique, 	/* 4 */
	MeanBestScore 	/* 5 */
};
enum {First, Second};
enum {NoneFound, Found, ContigAb, Inversion, OutsideBounds};

#endif
