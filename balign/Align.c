#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "../blib/BLib.h"
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/AlignEntry.h"
#include "AlignNTSpace.h"
#include "AlignColorSpace.h"
#include "Align.h"

/* TODO */
/* Assumes the read and reference are both forward strand */
/* If we are in color space, the read should also be in 
 * color space */
int Align(char *read,
		int readLength,
		char *reference,
		int referenceLength,
		ScoringMatrix *sm,
		AlignEntry *a,
		char strand,
		int space,
		int scoringType,
		int type)
{
	char *FnName="Align";
	int returnValue=-1;
	char reverseRead[SEQUENCE_LENGTH]="\0";
	char tmpString[SEQUENCE_LENGTH]="\0";
	char *reverseReference=NULL;

	switch(space) {
		case NTSpace:
			/* NT Space */
			switch(strand) {
				case FORWARD:
					/* Matches the forward strand */
					/* Align */
					returnValue = AlignNTSpace(read,
							readLength,
							reference,
							referenceLength,
							sm,
							a,
							type);
					break;
				case REVERSE:
					/* Reverse the read to match the forward strand  */
					GetReverseComplimentAnyCase(read, reverseRead, readLength);
					/* Align */
					returnValue = AlignNTSpace(reverseRead,
							readLength,
							reference,
							referenceLength,
							sm,
							a,
							type);

					/* We must reverse the alignment to match the REVERSE stand */
					GetReverseComplimentAnyCase(a->read, tmpString, a->length);
					strcpy(a->read, tmpString);
					GetReverseComplimentAnyCase(a->reference, tmpString, a->length);
					strcpy(a->reference, tmpString);
					break;
				default:
					PrintError(FnName,
							NULL,
							"Could not understand strand",
							Exit,
							OutOfRange);
					break;

			}
			break;
		case ColorSpace:
			/* Color Space */
			switch(strand) {
				case FORWARD:
					/* Matches the forward strand */
					/* Align */
					returnValue = AlignColorSpace(read,
							readLength,
							reference,
							referenceLength,
							sm,
							a,
							FORWARD,
							type,
							scoringType);
					break;
				case REVERSE:
					/* Matches the reverse strand */
					/* Reverse compliment the reference */
					reverseReference = malloc(sizeof(char)*(referenceLength+1));
					if(NULL==reverseReference) {
						PrintError(FnName,
								"reverseReference",
								"Could not allocate memory",
								Exit,
								MallocMemory);
					}
					GetReverseComplimentAnyCase(reference,
							reverseReference,
							referenceLength);
					/* Align */
					returnValue = AlignColorSpace(read,
							readLength,
							reverseReference,
							referenceLength,
							sm,
							a,
							REVERSE,
							type,
							scoringType);
					/* Adjust for the reverse strand */
					returnValue = referenceLength - a->referenceLength - returnValue;

					free(reverseReference);
					reverseReference=NULL;
					break;
				default:
					PrintError(FnName,
							NULL,
							"Could not understand strand",
							Exit,
							OutOfRange);
					break;

			}
			break;
		default:
			PrintError(FnName,
					"space",
					"Could not understand space",
					Exit,
					OutOfRange);
	}

	return returnValue;
}
