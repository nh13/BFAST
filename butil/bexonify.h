#ifndef BEXONIFY_H_
#define BEXONIFY_H_

typedef struct {
	uint8_t chr;
	uint32_t start;
	uint32_t end;
} Exon;

int ReadExons(char*, Exon**);
void FilterIndexBasedOnExons(RGIndex*, Exon **, int);

#endif
