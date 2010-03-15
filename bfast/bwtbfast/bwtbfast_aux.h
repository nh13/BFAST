#ifndef BWTBFAST_AUX_H_
#define BWTBFAST_AUX_H_

// Constants
#define BWTBFAST_AUX_SHELL_SORT_FACTOR 2.2

// Macros
#define cmp_c_p_s(_ca, _pa, _sa, _cb, _pb, _sb) ((_ca < _cb) ? -1 : ((_ca == _cb && _pa < _pb) ? -1 : ((_ca == _cb && _pa == _pb && _sa < _sb) ? -1 : ((_ca == _cb && _pa == _pb && _sa == _sb) ? 0 : 1))))

#include "bwtaln.h"
#include "bntseq.h"
#include "occ_results.h"

typedef struct {
	size_t l, m;
	char *s;
} bfast_string_t;

typedef struct {
	bfast_string_t name, comment, seq, qual;
} bfast_seq_t;

typedef struct {
	bfast_string_t *masks;
	int32_t n;
} bfast_masks_t;

#ifdef __cplusplus
extern "C" {
#endif
	void bfast_masks_t_init(bfast_masks_t *masks);
	void bfast_masks_t_add(bfast_masks_t *masks, char *mask);
	void bfast_masks_t_destroy(bfast_masks_t *masks);
#ifdef __cplusplus
}
#endif

typedef struct {
	int32_t read_name_length;
	char *read_name;
	int32_t read_length;
	char *read;
	int32_t qual_length;
	char *qual;
	int32_t max_reached;
	int32_t num_entries;
	uint32_t *contigs;
	int32_t *positions;
	char *strands;
	char **masks;

	// for multi-threading only - do not print
	int32_t tid;

	// for matching only - do not print
	uint8_t *read_int;
	uint8_t *read_rc_int;
	int32_t read_int_length;
} bfast_rg_match_t;

typedef struct {
	int32_t num_ends;
	bfast_rg_match_t *ends;
} bfast_rg_matches_t;

#ifdef __cplusplus
extern "C" {
#endif
	void bfast_rg_match_t_copy_results(bfast_rg_match_t *match, bwt_t *bwt, bntseq_t *bns, bfast_string_t *mask, int32_t space, int32_t offset, occ_results_t *results_f, occ_results_t *results_r);
	void bfast_rg_match_t_init(bfast_rg_match_t *match, bfast_seq_t *seq);
	void bfast_rg_match_t_destroy(bfast_rg_match_t *match);
	void bfast_rg_match_t_append(bfast_rg_match_t *dest, bfast_rg_match_t *src);
	void bfast_rg_match_t_clear(bfast_rg_match_t *match);
	void bfast_rg_match_t_print_queue(bfast_rg_match_t *match_queue, int32_t len, gzFile fp);
	void bfast_rg_matches_t_print(bfast_rg_match_t *match_queue, int32_t from, int32_t to, gzFile fp);
	void bfast_rg_match_t_print(bfast_rg_match_t *match, gzFile fp);
	char *bfast_rg_match_create_mask(bfast_string_t *bfast_mask, int32_t space, int32_t offset, int32_t strand, int32_t read_length);
	char *bfast_rg_match_string_to_mask(bfast_string_t *bfast_mask, int32_t read_length);
	void bfast_rg_match_t_print_text(bfast_rg_match_t *match, FILE *fp); // DEBUGGING FUNCTION
	void bfast_rg_match_t_copy_from_bwa(bfast_rg_match_t *m, bwa_seq_t *seq, int32_t space);
	bfast_rg_match_t *bfast_rg_match_read(bwa_seqio_t *bs, int n_needed, int *n, int space, int trim_qual);
#ifdef __cplusplus
}
#endif

#endif

