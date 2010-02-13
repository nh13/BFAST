#ifndef BWTBFAST_H_
#define BWTBFAST_H_

#include "bwtbfast_aux.h"

// Constants
#define BFAST_SHELL_SORT_FACTOR 2.2
#define BWTBFAST_QUEUE_LENGTH 10000
#define BWTBFAST_THREAD_SLEEP 1

#ifdef HAVE_PTHREAD
typedef struct {
	int tid;
	bwt_t *bwt;
	bntseq_t *bns;
	int n_matches;
	bfast_rg_match_t *matches;
	int32_t space;
	int n_threads;

	int32_t max_key_matches;
	int32_t max_num_matches;
	bfast_masks_t *masks;
} bwtbfast_thread_t;
#endif

int bwtbfast_usage(int32_t space, int32_t max_key_matches, int32_t max_num_matches, int32_t num_threads, int32_t queue_length);
int bwtbfast(int argc, char *argv[]);
void bwtbfast_core(char *ref_fn, char *read_fn, bfast_masks_t *masks, int32_t space, int32_t max_key_matches, int32_t max_num_matches, int32_t num_threads, int32_t queue_length);
#ifdef HAVE_PTHREAD
void *bwtbfast_thread_worker(void *data);
#endif
void bwtbfast_core_worker(int tid, bwt_t *bwt, bntseq_t *bns, int n_matches, bfast_rg_match_t *matches, int32_t space, int n_threads, bfast_masks_t *masks, int32_t max_key_matches, int32_t max_num_matches);

void bwt_t_get_with_mask(const bwt_t *bwt, uint8_t *w_int, char *mask, int32_t len, occ_results_t *r);
void bwt_t_get_with_mask_helper(const bwt_t *bwt, uint8_t *w_int, char *mask, int32_t len, int32_t r_l_old, int32_t r_u_old, occ_results_t *r);

#endif
