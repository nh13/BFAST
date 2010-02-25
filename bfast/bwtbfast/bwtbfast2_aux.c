#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "BLibDefinitions.h"
#include "util.h"
#include "BError.h"
#include "bwtaln.h"
#include "bwtbfast2_aux.h"

#define STATE_M 0
#define STATE_I 1
#define STATE_D 2

bfast2_stack_t *bfast2_init_stack(int max_mm)
{
	int i;
	char *fn_name="bfast2_init_stack";
	bfast2_stack_t *stack=NULL;
	stack = (bfast2_stack_t*)my_calloc(1, sizeof(bfast2_stack_t), fn_name);
	stack->n_stacks = max_mm+1;
	stack->stacks = (bfast2_stack1_t*)my_calloc(stack->n_stacks, sizeof(bfast2_stack1_t), fn_name);
	for (i = 0; i != stack->n_stacks; ++i) {
		bfast2_stack1_t *p = stack->stacks + i;
		p->m_entries = 4;
		p->entries = (bfast2_entry_t*)my_calloc(p->m_entries, sizeof(bfast2_entry_t), fn_name);
	}
	return stack;
}

void bfast2_destroy_stack(bfast2_stack_t *stack)
{
	int i;
	for (i = 0; i != stack->n_stacks; ++i) {
		free(stack->stacks[i].entries);
	}
	free(stack->stacks);
	free(stack);
}

static void bfast2_reset_stack(bfast2_stack_t *stack, int max_mm)
{
	char *fn_name="bfast2_reset_stack";
	int i;
	for (i = 0; i != stack->n_stacks; ++i)
		stack->stacks[i].n_entries = 0;
	stack->best = 0;
	stack->n_entries = 0;
	if(stack->n_stacks <= max_mm) {
		i = stack->n_stacks;
		stack->n_stacks = max_mm+1;
		stack->stacks = (bfast2_stack1_t*)my_realloc(stack->stacks, sizeof(bfast2_stack1_t)*stack->n_stacks, fn_name);
		while(i != stack->n_stacks) {
			bfast2_stack1_t *p = stack->stacks + i;
			p->m_entries = 4;
			p->entries = (bfast2_entry_t*)my_calloc(p->m_entries, sizeof(bfast2_entry_t)*p->m_entries, fn_name);
			i++;
		}
	}
}

static inline void bfast2_seed(bfast2_stack_t *stack, 
		int32_t rl, // substring length, not full read length
		int32_t offset,
		uint8_t strand,
		bwtint_t k, 
		bwtint_t l
		)
{
	char *fn_name="bfast2_seed";
	int32_t i;
	bfast2_entry_t *p=NULL;
	bfast2_stack1_t *q=NULL;

	q = stack->stacks;
	if (q->n_entries == q->m_entries) {
		q->m_entries <<= 1;
		q->entries = (bfast2_entry_t*)my_realloc(q->entries, sizeof(bfast2_entry_t) * q->m_entries, fn_name);
	}
	p = q->entries + q->n_entries;
	p->n_mm = 0;
	p->k = k; p->l = l; 
	for(i=0;i<BWTBFAST2_AUX_MASK_BYTES;i++) { // reset all zeros
		p->mask[i] = 0;
	}
	p->next_i = rl-1; 
	p->strand = strand;
	p->mask_i = (FORWARD == p->strand) ? (rl-1) : (offset); 
	// set match on the last base
	p->mask[(int)(p->mask_i >> 3)] |= (0x01 << (p->mask_i & 7));
	p->next_i--;
	p->mask_i = (FORWARD == p->strand) ? (p->mask_i-1) : (p->mask_i+1);
	p->offset = offset;
	p->bases_used=1;

	++(q->n_entries);
	++(stack->n_entries);
}

static inline void bfast2_push(bfast2_stack_t *stack, 
		bwtint_t k, 
		bwtint_t l, 
		int is_mm,
		bfast2_entry_t *prev
		)
{
	char *fn_name="bfast2_push";
	int32_t i, j;
	bfast2_entry_t *p=NULL;
	bfast2_stack1_t *q=NULL;

	if(stack->n_stacks <= prev->n_mm + is_mm) {
		PrintError(fn_name, "stack->n_stacks <= prev->n_mm + is_mm", "Control reached unexpected point", Exit, OutOfRange);
	}

	q = stack->stacks + prev->n_mm + is_mm;
	if (q->n_entries == q->m_entries) {
		q->m_entries <<= 1;
		q->entries = (bfast2_entry_t*)my_realloc(q->entries, sizeof(bfast2_entry_t) * q->m_entries, fn_name);
	}

	// sort by # of bases used (low to high):
	p = q->entries;
	for(i=q->n_entries;0<i;i--) {
		p = q->entries + i - 1;
		if(p->bases_used <= prev->bases_used - is_mm + 1) {
			// shift over
			for(j=i;j<q->n_entries;j++) {
				q->entries[j+1] = q->entries[j];
			}
			p = q->entries + i;
			break;
		}
	}

	p->k = k; p->l = l; 
	p->n_mm = prev->n_mm + is_mm; 
	p->offset = prev->offset;

	// copy mask
	memcpy(p->mask, prev->mask, BWTBFAST2_AUX_MASK_BYTES);

	p->bases_used = prev->bases_used;
	if(1 != is_mm) {
		// set match on the last base
		p->mask[(int)(p->mask_i >> 3)] |= (0x01 << (p->mask_i & 7));
		p->bases_used++;
	}
	p->strand = prev->strand;
	p->mask_i = (FORWARD == p->strand) ? (prev->mask_i-1) : (prev->mask_i+1);
	p->next_i = prev->next_i-1;

	++(q->n_entries);
	++(stack->n_entries);
	if (stack->best < p->n_mm) stack->best = p->n_mm;
}

static inline void bfast2_pop(bfast2_stack_t *stack, bfast2_entry_t *e)
{
	bfast2_stack1_t *q;
	q = stack->stacks + stack->best;
	(*e) = q->entries[q->n_entries - 1];
	--(q->n_entries);
	--(stack->n_entries);
	if(q->n_entries == 0 && 0 < stack->n_entries) { // reset best
		int i;
		for (i = stack->best; 0 <= i; --i) {
			if (0 < stack->stacks[i].n_entries) {
				stack->best = i;
				break;
			}
		}
		stack->best = i;
		if(i < 0) {
			// should not reach here
			PrintError("bfast2_pop", NULL, "Control reached unexpected line", Exit, OutOfRange); 
		}
	} else if (stack->n_entries == 0) {
		stack->best = 0;
	}
}

static inline int int_log2(uint32_t v)
{
	int c = 0;
	if (v & 0xffff0000u) { v >>= 16; c |= 16; }
	if (v & 0xff00) { v >>= 8; c |= 8; }
	if (v & 0xf0) { v >>= 4; c |= 4; }
	if (v & 0xc) { v >>= 2; c |= 2; }
	if (v & 0x2) c |= 1;
	return c;
}

void bfast2_match(bfast_rg_match_t *match, bwt_t *bwt, bntseq_t *bns, int32_t space, int32_t alg, int32_t seed_len, int32_t max_mm, int32_t max_seed_hits, int32_t max_hits, bfast2_stack_t *stack) 
{
	char *fn_name="bfast2_match";
	int i;
	bwtint_t k, l;
	uint8_t cur_base;
	bwtint_t cntk[4], cntl[4];
	int32_t best_n_mm=-1;

	assert(match->read_int_length <= BWTBFAST2_AUX_MASK_BYTES*8);

	// initialize tmp match
	bfast_rg_match_t *tmp_match = my_calloc(1, sizeof(bfast_rg_match_t), fn_name);
	// these are need for bfast2_rg_match_t_add 
	tmp_match->read_length = match->read_length;
	tmp_match->read_int_length = match->read_int_length;

	// reset stack
	bfast2_reset_stack(stack, max_mm);

	// seed - make sure last base matches
	for(i=0;i < match->read_int_length - seed_len + 1;i++) {
		// forward
		cur_base = match->read_int[match->read_int_length-i-1];
		bwt_2occ(bwt, -1, bwt->seq_len, cur_base, &k, &l);
		k = bwt->L2[cur_base] + k + 1;
		l = bwt->L2[cur_base] + l;
		if(k <= l) {
			bfast2_seed(stack, match->read_int_length - i, i, FORWARD, k, l);
		}

		// reverse
		cur_base = match->read_rc_int[match->read_int_length-i-1];
		bwt_2occ(bwt, -1, bwt->seq_len, cur_base, &k, &l);
		k = bwt->L2[cur_base] + k + 1;
		l = bwt->L2[cur_base] + l;
		if(k <= l) {
			bfast2_seed(stack, match->read_int_length - i, i, REVERSE, k, l);
		}
	}

	while(0 < stack->n_entries) {
		bfast2_entry_t e;

		// get the best entry
		bfast2_pop(stack, &e); 

		/*
		   fprintf(stderr, "e.n_mm=%d e.next_i=%d e.k=%d e.l=%d e.offset=%d e.bases_used=%d\n",
		   e.n_mm, e.next_i, e.k, e.l, e.offset, e.bases_used);
		   */
		if(max_mm < e.n_mm) break; // too many mismatches

		// only report all matches with the same minimum # of mismatches
		if(1 == alg && 0 <= best_n_mm && best_n_mm < e.n_mm) break;
		// only report all matches with one more mismatch
		if(2 == alg && 0 <= best_n_mm && best_n_mm+1 < e.n_mm) break;

		// only continue if there are more hits to be had or we have 
		// reached the maximum # of hits and are beyond the best mm.
		if(0 < max_hits && max_hits < tmp_match->num_entries) {
				// clear and set that the max has been reached
				bfast_rg_match_t_clear(tmp_match);
				tmp_match->max_reached = 1;
				break;
		}

		if(seed_len == e.bases_used) {
			// ignore uninformative seeds
			if(max_seed_hits < 0 || e.l - e.k + 1 <= max_seed_hits) {

				if(0 < max_hits && max_hits < e.l - e.k + 1) { // too many will be added
					bfast_rg_match_t_clear(tmp_match);
					tmp_match->max_reached = 1;
					break;
				}

				// Found -> report
				bfast2_rg_match_t_add(tmp_match, bwt, bns, space, seed_len, &e);

				if(0 == alg) break; // only report the first match
				if(best_n_mm < 0) best_n_mm = e.n_mm; // store the best number of mismatches
			}
		}
		else if(0 <= e.next_i) { // bases left

			cur_base = (e.strand == FORWARD) ? match->read_int[e.next_i] : match->read_rc_int[e.next_i];

			if(0 == e.n_mm - max_mm) { // exact
				if(4 != cur_base) { // no missing bases
					bwt_2occ(bwt, e.k-1, e.l, cur_base, &k, &l);
					k = bwt->L2[cur_base] + k + 1;
					l = bwt->L2[cur_base] + l;
					if(k <= l) {
						bfast2_push(stack, k, l, 0, &e);
					}
				}
			}
			else { // try all mismatches
				assert(e.n_mm < max_mm); 
				bwt_2occ4(bwt, e.k-1, e.l, cntk, cntl);
				for(i = 0; i < 4; i++) { // go through each base
					k = bwt->L2[i] + cntk[i] + 1; 
					l = bwt->L2[i] + cntl[i];
					if(k <= l) {
						bfast2_push(stack, k, l, (cur_base == i) ? 0 : 1, &e);
					}
				}
			}
		}
	}

	// add tmp match
	bfast_rg_match_t_append(match, tmp_match);

	// destroy
	bfast_rg_match_t_destroy(tmp_match);
	free(tmp_match);
}

void bfast2_rg_match_t_add(bfast_rg_match_t *match, bwt_t *bwt, bntseq_t *bns, int32_t space, int32_t seed_len, bfast2_entry_t *e)
{
	char *fn_name="bfast2_rg_match_t_add";
	int32_t i, j, len;
	int32_t num_entries;
	bwtint_t pos;
	int seqid;
	int32_t n_skipped=0;

	if(e->l < e->k) return; // control should not reach here

	len = GETMASKNUMBYTESFROMLENGTH(match->read_length);

	// get the number of entries to be added
	num_entries = match->num_entries + e->l - e->k + 1;

	// reallocate memory
	match->contigs = my_realloc(match->contigs, sizeof(uint32_t)*num_entries, fn_name);
	match->positions = my_realloc(match->positions, sizeof(int32_t)*num_entries, fn_name);
	match->strands = my_realloc(match->strands, sizeof(char)*num_entries, fn_name);
	match->masks = my_realloc(match->masks, sizeof(char*)*num_entries, fn_name);

	// append
	// important: this assumes we will sort and merge entries later

	i=match->num_entries;
	for(j=e->k;j<=e->l;j++) {
		// use suffix array
		pos = bwt_sa(bwt, j);
		bns_coor_pac2real(bns, pos, 1, &seqid);

		// corner case checking
		{
			// check if hit spans two sequences
			if(bns->anns[seqid].len < pos + seed_len + e->n_mm) {
				n_skipped++;
				continue;
			}
		}

		//contig
		match->contigs[i] = seqid+1;
		// position
		match->positions[i] = pos - bns->anns[seqid].offset + 1;
		if(FORWARD == e->strand) match->positions[i] -= (1 + e->next_i);
		else {
			match->positions[i] += e->offset - match->read_int_length + seed_len + e->n_mm - space;
		}

		// strand
		match->strands[i] = e->strand;
		// mask
		match->masks[i] = my_calloc(len, sizeof(char), fn_name);
		memcpy(match->masks[i], e->mask, len);


		i++;

	}
	if(0 < n_skipped) {
		num_entries -= n_skipped;
		if(0 == num_entries) {
			// free memory
			free(match->contigs); match->contigs=NULL;
			free(match->positions); match->positions=NULL;
			free(match->strands); match->strands=NULL;
			free(match->masks); match->masks=NULL;
		}
		else {
			// reallocate memory
			match->contigs = my_realloc(match->contigs, sizeof(uint32_t)*num_entries, fn_name);
			match->positions = my_realloc(match->positions, sizeof(int32_t)*num_entries, fn_name);
			match->strands = my_realloc(match->strands, sizeof(char)*num_entries, fn_name);
			match->masks = my_realloc(match->masks, sizeof(char*)*num_entries, fn_name);
		}
	}
	match->num_entries = num_entries;
				
	if(0 < match->num_entries) {
		// merge duplicates
		bfast2_rg_match_t_merge_duplicates(match);
	}
}

void bfast2_rg_match_t_sort(bfast_rg_match_t *match) 
{
	// for shell-sort
	int32_t i, j, len;
	int32_t low, high, inc;
	uint32_t c_tmp;
	int32_t p_tmp;
	int8_t s_tmp;
	char *m_tmp=NULL;

	if(1 < match->num_entries) {
		// keep it generalized (low/high) in case we wish to modularize it
		low = 0;
		high = match->num_entries-1;
		len = GETMASKNUMBYTESFROMLENGTH(match->read_length);

		// Shell sort
		inc = ROUND((high - low + 1) >> 1);
		while(0 < inc) {
			for(i=inc + low;i<=high;i++) {
				c_tmp = match->contigs[i];
				p_tmp = match->positions[i];
				s_tmp = match->strands[i];
				m_tmp = match->masks[i];
				j = i;
				while(inc + low <= j &&
						cmp_c_p_s(c_tmp, p_tmp, s_tmp, match->contigs[j-inc], match->positions[j-inc], match->strands[j-inc]) < 0) {
					match->contigs[j] = match->contigs[j-inc];
					match->positions[j] = match->positions[j-inc];
					match->strands[j] = match->strands[j-inc];
					match->masks[j] = match->masks[j-inc];
					j -= inc;
				}
				match->contigs[j] = c_tmp;
				match->positions[j] = p_tmp;
				match->strands[j] = s_tmp;
				match->masks[j] = m_tmp;
			}
			inc = ROUND(inc / BWTBFAST_AUX_SHELL_SORT_FACTOR);
		}
	}
}

void bfast2_rg_match_t_merge_duplicates(bfast_rg_match_t *match)
{
	char *fn_name="bfast2_rg_match_t_merge_duplicates";
	int32_t i, j, prev_index, len;

	if(match->num_entries <= 1) return;

	len = GETMASKNUMBYTESFROMLENGTH(match->read_length);

	// sort
	bfast2_rg_match_t_sort(match);

	// remove duplicates
	prev_index=0;
	for(i=1;i<match->num_entries;i++) {
		if(0==cmp_c_p_s(match->contigs[prev_index], match->positions[prev_index], match->strands[prev_index],
					match->contigs[i], match->positions[i], match->strands[i])) {
			// union masks
			for(j=0;j<len;j++) {
				match->masks[prev_index][j] |= match->masks[i][j];
			}
			free(match->masks[i]);
			match->masks[i] = NULL;
		}
		else {
			prev_index++;
			// copy to prev_index (incremented)
			if(prev_index != i) {
				match->contigs[prev_index] = match->contigs[i];
				match->positions[prev_index] = match->positions[i];
				match->strands[prev_index] = match->strands[i];
				// copy masks too!
				free(match->masks[prev_index]);
				match->masks[prev_index] = match->masks[i];
				match->masks[i] = NULL;
			}
		}
	}
	prev_index++;
	// reallocate
	assert(prev_index <= match->num_entries);
	if(prev_index < match->num_entries) {
		match->num_entries = prev_index;
		match->contigs = my_realloc(match->contigs, sizeof(uint32_t)*match->num_entries, fn_name);
		match->positions = my_realloc(match->positions, sizeof(int32_t)*match->num_entries, fn_name);
		match->strands = my_realloc(match->strands, sizeof(int32_t)*match->num_entries, fn_name);
		match->masks = my_realloc(match->masks, sizeof(char*)*match->num_entries, fn_name);
	}
}
