#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#define MIN_FILTER 0
#define MAX_FILTER 4
#define DEFAULT_QUEUE_LENGTH 10000

enum {NoFiltering, 		/* 0 */
	AllNotFiltered, 	/* 1 */
	Unique, 			/* 2 */
	BestScore, 			/* 3 */
	BestScoreAll,		/* 4 */
};
enum {First, Second};
enum {NoneFound, Found};

#endif
