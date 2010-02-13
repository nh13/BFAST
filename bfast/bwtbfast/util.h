#ifndef UTIL_H_
#define UTIL_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include "kseq.h"

// CONSTANTS
#define BWTA_REF_SUFFIX ".bwta_ref"
#define BWTA_SA_SUFFIX ".bwta_sa"
#define BWTA_BWT_B_SUFFIX ".bwta_bwt_b"
#define BWTA_BWT_B_REV_SUFFIX ".bwta_bwt_b_rev"
#define BWTA_OCC_SUFFIX ".bwta_occ"
#define BWTA_OCC_REV_SUFFIX ".bwta_occ_rev"
#define STR_LEN 1024
#define MAGIC_ID 'M'+'A'+'G'+'I'+'C'
#define ALPHABET_SIZE 4
#define READ_REF_ROTATE_NUM 1000000 
#define MAX_NUM_REFS (1 << (8*((int)sizeof(uint8_t)) - 1))
#define BWTA_SAM_VERSION "0.1.2" 
#define FORWARD '+'
#define REVERSE '-'

// MACROS
#define __get_min(_X, _Y)  ((_X) < (_Y) ? (_X) : (_Y))
#define __ref_t_update(_ref, _c, _i) do { \
	if(0 == (_i % 4)) { \
		_ref[(int)(_i>>2)] = 0; \
	} \
	_ref[(int)(_i>>2)] |= (base_to_int(_c) << (2 * (_i & 3))); \
} while(0)
#define ROUND(_x) ((int)((_x) + 0.5))
#define cmp_c_p_s(_ca, _pa, _sa, _cb, _pb, _sb) ((_ca < _cb) ? -1 : ((_ca == _cb && _pa < _pb) ? -1 : ((_ca == _cb && _pa == _pb && _sa < _sb) ? -1 : ((_ca == _cb && _pa == _pb && _sa == _sb) ? 0 : 1))))
#define __log2(_x) ((_x <= 0) ? 0 : (31 - __builtin_clz(_x)))

// FUNCTIONS
void *my_realloc(void *ptr, size_t size, char *fn_name);
void *my_malloc(size_t size, char *fn_name);
void *my_calloc(size_t num, size_t size, char *fn_name);
#endif
