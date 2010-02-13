#ifndef BWTBFAST2_H_
#define BWTBFAST2_H_

#include "bwtbfast_aux.h"

// Constants
#define BFAST_SHELL_SORT_FACTOR 2.2
#define BWTBFAST2_QUEUE_LENGTH 10000
#define BWTBFAST2_THREAD_SLEEP 1

#ifdef HAVE_PTHREAD
typedef struct {
	int tid;
	bwt_t *bwt;
	bntseq_t *bns;
	int n_matches;
	bfast_rg_match_t *matches;
	int32_t space;
	int n_threads;
	int32_t alg;
	int32_t seed_len;
	int32_t max_mm;
	int32_t max_hits;
} bwtbfast2_thread_t;
#endif

int bwtbfast2_usage(int32_t alg, int32_t space, int32_t seed_len, int32_t max_mm, int32_t num_threads, int32_t max_hits, int32_t queue_length);
int bwtbfast2(int argc, char *argv[]);
void bwtbfast2_core(char *ref_fn, char *read_fn, int32_t compression, int32_t alg, int32_t seed_len, int32_t max_mm, int32_t space, int32_t start_read_num, int32_t end_read_num, int32_t max_hits, int32_t num_threads, int32_t queue_length);
#ifdef HAVE_PTHREAD
void *bwtbfast2_thread_worker(void *data);
#endif
void bwtbfast2_core_worker(int tid, bwt_t *bwt, bntseq_t *bns, int n_matches, bfast_rg_match_t *matches, int32_t space, int n_threads, int32_t alg, int32_t seed_len, int32_t max_mm, int32_t max_hits);
#endif
