#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

enum {MAF, LastFileType};
#define MIN_FILTER 0
#define MAX_FILTER_SE 3
#define MAX_FILTER_PE 5
enum {NoFiltering, AllNotFiltered, Unique, BestScore, MeanUnique, MeanBestScore};
enum {First, Second};
enum {NoneFound, Found, ChrAb, Inversion, OutsideBounds};

#endif
