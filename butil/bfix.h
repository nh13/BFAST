#ifndef BFIX_H_
#define BFIX_H_

/* Add to definitions when the data structures change. At
 * some point, we should call old data structures unsupported. */

enum{V_0_1_13, V_Last};

void ConvertRGBinary(char*, int);
void ConvertRGIndex(char*, int);

/******************************/
/* For 0.1.13 Data Structures */
typedef struct {
	/* Storage type */
	int32_t id;
	/* RG storage */
	RGBinaryContig *contigs;
	int32_t numContigs;
	/* Metadata */
	int32_t space;
} RGBinary_0_1_13;
typedef struct {
	/* Storage type */
	int32_t id;
	/* Index storage */
	uint8_t *contigs_8;
	uint32_t *contigs_32;
	int32_t *positions;
	int64_t length;
	int32_t contigType;
	/* Index range */
	int32_t startContig;
	int32_t startPos;
	int32_t endContig;
	int32_t endPos;
	/* Index layout */
	int32_t width;
	int32_t keysize;
	int32_t *mask;
	/* Index properties */
	int32_t repeatMasker;
	int32_t space;
	/* Hash storage */
	uint32_t hashWidth; /* in bases */
	int64_t hashLength;
	uint32_t *starts;
	uint32_t *ends;
} RGIndex_0_1_13;
/* 0.1.13 Read Functions */
void ConvertRGBinaryFrom_0_1_13(char*, char*);
void ConvertRGIndexFrom_0_1_13(char*, char*);
void RGBinaryReadBinary_0_1_13(RGBinary_0_1_13*, char*);
void RGIndexRead_0_1_13(FILE*, RGIndex_0_1_13*, int32_t);
void RGIndexReadHeader_0_1_13(FILE*, RGIndex_0_1_13*, int32_t);
/* End for 0.1.13 Data Structures */
/**********************************/

#endif
