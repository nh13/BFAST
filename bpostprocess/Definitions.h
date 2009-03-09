#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#define MIN_FILTER 0
#define MAX_FILTER 3
enum {NoFiltering, 		/* 0 */
	AllNotFiltered, 	/* 1 */
	Unique, 			/* 2 */
	BestScore, 			/* 3 */
};
enum {First, Second};
enum {NoneFound, Found, ContigAb, Unpaired, Inversion, OutsideBounds};

#endif
