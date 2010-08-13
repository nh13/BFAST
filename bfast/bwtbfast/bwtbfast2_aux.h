#ifndef BWTBFAST2_AUX_H_
#define BWTBFAST2_AUX_H_

#include "bwt.h"
#include "bwtaln.h"
#include "bwtbfast_aux.h"

// TODO: # of bytes in the mask should checked against the read length
#define BWTBFAST2_AUX_MASK_BYTES 16

typedef struct { // recursion stack
	bwtint_t k, l; // (k,l) is the SA region of [i,n-1]
	int32_t mask_i:16, next_i:16; 
	uint8_t strand;
	uint8_t n_mm;
	int32_t offset:16, bases_used:16;
	uint8_t mask[BWTBFAST2_AUX_MASK_BYTES]; // mask bytes
} bfast2_entry_t;

typedef struct {
	int n_entries, m_entries;
	bfast2_entry_t *entries;
} bfast2_stack1_t;

typedef struct {
	// one stack per mm level: 0, 1, 2, 3, ...
	int n_stacks, best, n_entries;
	bfast2_stack1_t *stacks;
} bfast2_stack_t;

#ifdef __cplusplus
extern "C" {
#endif

	bfast2_stack_t *bfast2_init_stack(int max_mm);
	void bfast2_destroy_stack(bfast2_stack_t *stack);

	void bfast2_match(bfast_rg_match_t *match, bwt_t *bwt, bntseq_t *bns, int32_t space, int32_t alg, int32_t seed_len, int32_t max_mm, int32_t max_seed_hits, int32_t max_hits, int32_t seed_ends_only, bfast2_stack_t *stack);

	void bfast2_rg_match_t_add(bfast_rg_match_t *match, bwt_t *bwt, bntseq_t *bns, int32_t space, int32_t seed_len, bfast2_entry_t *e);
	void bfast2_rg_match_t_sort(bfast_rg_match_t *match);
	void bfast2_rg_match_t_merge_duplicates(bfast_rg_match_t *match);


#ifdef __cplusplus
}
#endif

#endif
