#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include <limits.h>
#include <assert.h>
#include <pthread.h>
#include "BError.h"
#include "BLib.h"
#include "BLibDefinitions.h"
#include "util.h"
#include "occ_results.h"
#include "bwt.h"
#include "bntseq.h"
#include "bwtbfast_aux.h"

uint8_t bfast_nt_int_to_nt_int_comp[5] = {3, 2, 1, 0, 4};

void bfast_masks_t_init(bfast_masks_t *masks)
{
	masks->masks=NULL;
	masks->n=0;
}

void bfast_masks_t_add(bfast_masks_t *masks, char *mask)
{
	char *fn_name="bfast_masks_t_add";
	int32_t i;

	masks->n++;
	masks->masks = my_realloc(masks->masks, sizeof(bfast_string_t)*masks->n, fn_name);
	masks->masks[masks->n-1].l = strlen(mask);
	masks->masks[masks->n-1].m = 1 + strlen(mask);
	masks->masks[masks->n-1].s = my_malloc(sizeof(char)*masks->masks[masks->n-1].m, fn_name);
	strcpy(masks->masks[masks->n-1].s, mask);
	// check all is ok
	for(i=0;i<masks->masks[masks->n-1].l;i++) {
		switch(masks->masks[masks->n-1].s[i]) {
			case '0':
			case '1':
				break;
			default:
				PrintError(fn_name, mask, "Incorrect format for the mask", Exit, OutOfRange);
		}
	}
}

void bfast_masks_t_destroy(bfast_masks_t *masks)
{
	int32_t i;
	for(i=0;i<masks->n;i++) {
		free(masks->masks[i].s);
	}
	free(masks->masks);
}

void bfast_rg_match_t_copy_results(bfast_rg_match_t *match, bwt_t *bwt, bntseq_t *bns, bfast_string_t *bfast_mask, int32_t space, int32_t offset, occ_results_t *results_f, occ_results_t *results_r)
{
	char *fn_name="bfast_rg_match_t_copy_results";
	int32_t n=0, i, j, ctr;
	int32_t num_entries;
	uint8_t *c=NULL, c_tmp;
	uint32_t *p=NULL, p_tmp;
	char *s=NULL, s_tmp;
	int32_t inc, low, high, prev_index;
	int32_t start_match, start_new, end_match, end_new;
	bwtint_t pos;
	int seqid;

	if(results_f->n + results_r->n <= 0) return; 

	// get the number of entries to be added
	n=0;
	for(i=0;i<results_f->n;i++) {
		n+=results_f->r_u[i]-results_f->r_l[i]+1;
	}
	for(i=0;i<results_r->n;i++) {
		n+=results_r->r_u[i]-results_r->r_l[i]+1;
	}

	// allocate
	c = my_malloc(sizeof(uint8_t)*n, fn_name);
	p = my_malloc(sizeof(uint32_t)*n, fn_name);
	s = my_malloc(sizeof(char)*n, fn_name);

	// copy over c/p/s
	int32_t n_skipped = 0;
	ctr=0;
	for(i=0;i<results_f->n;i++) {
		for(j=results_f->r_l[i];j<=results_f->r_u[i];j++) {
			pos = bwt_sa(bwt, j);
			bns_coor_pac2real(bns, pos, 1, &seqid); // third argument ?
			// corner case checking
			{
				// check if hit spans two sequences
				if(bns->anns[seqid].len < pos + bfast_mask->l) {
					n_skipped++;
					continue;
				}
			}
			c[ctr] = seqid+1; 
			p[ctr] = pos - bns->anns[seqid].offset + 1;
			p[ctr] -= offset;
			s[ctr]=FORWARD;
			ctr++;
		}
	}
	for(i=0;i<results_r->n;i++) {
		for(j=results_r->r_l[i];j<=results_r->r_u[i];j++) {
			pos = bwt_sa(bwt, j);
			bns_coor_pac2real(bns, pos, 1, &seqid); // third argument ?
			// corner case checking
			{
				// check if hit spans two sequences
				if(bns->anns[seqid].len < pos + bfast_mask->l) {
					n_skipped++;
					continue;
				}
			}
			c[ctr] = seqid+1;
			p[ctr] = pos - bns->anns[seqid].offset + 1;
			p[ctr] += bfast_mask->l + offset - match->read_length;
			s[ctr]=REVERSE;
			ctr++;
		}
	}

	if(0 < n_skipped) {
		n-=n_skipped;
		// reallocate
		c = my_realloc(c, sizeof(uint8_t)*n, fn_name);
		p = my_realloc(p, sizeof(uint32_t)*n, fn_name);
		s = my_realloc(s, sizeof(char)*n, fn_name);
	}

	if(1 < n) {
		// sort c/p/s with shell sort
		// keep it generalized (low/high) in case we wish to modularize it
		low = 0;
		high = n-1;
		inc = ROUND((high - low + 1) >> 1);
		while(0 < inc) {
			for(i=inc + low;i<=high;i++) {
				c_tmp = c[i];
				p_tmp = p[i];
				s_tmp = s[i];
				j = i;
				while(inc + low <= j &&
						cmp_c_p_s(c_tmp, p_tmp, s_tmp, c[j-inc], p[j-inc], s[j-inc]) < 0) {
					c[j] = c[j-inc];
					p[j] = p[j-inc];
					s[j] = s[j-inc];
					j -= inc;
				}
				c[j] = c_tmp;
				p[j] = p_tmp;
				s[j] = s_tmp;
			}
			inc = ROUND(inc / BWTBFAST_AUX_SHELL_SORT_FACTOR);
		}
	}

	// add more entries
	num_entries = match->num_entries + n;
	match->contigs = my_realloc(match->contigs, sizeof(uint32_t)*num_entries, fn_name);
	match->positions = my_realloc(match->positions, sizeof(int32_t)*num_entries, fn_name);
	match->strands = my_realloc(match->strands, sizeof(char)*num_entries, fn_name);
	match->masks = my_realloc(match->masks, sizeof(char*)*num_entries, fn_name);

	// shift over
	for(i=match->num_entries-1;0<=i;i--) {
		match->contigs[i+n] = match->contigs[i];
		match->positions[i+n] = match->positions[i];
		match->strands[i+n] = match->strands[i];
		match->masks[i+n] = match->masks[i];
	}
	match->num_entries = num_entries;

	// merge match c/p/s with new c/p/s
	start_match = n;
	end_match = num_entries-1;
	start_new = 0;
	end_new = n-1;
	ctr=0;
	while((start_match <= end_match) && (start_new <= end_new)) {
		assert(ctr <= start_match);
		if(cmp_c_p_s(match->contigs[start_match], match->positions[start_match], match->strands[start_match],
					c[start_new], p[start_new], s[start_new]) <= 0) {
			match->contigs[ctr] = match->contigs[start_match];
			match->positions[ctr] = match->positions[start_match];
			match->strands[ctr] = match->strands[start_match];
			match->masks[ctr] = match->masks[start_match];
			start_match++;
		}
		else {
			match->contigs[ctr] = c[start_new];
			match->positions[ctr] = p[start_new];
			match->strands[ctr] = s[start_new];
			match->masks[ctr] = bfast_rg_match_create_mask(bfast_mask, space, offset, s[start_new], match->read_length);
			start_new++;
		}
		ctr++;
	}
	while(start_match <= end_match) {
		assert(ctr <= start_match);
		match->contigs[ctr] = match->contigs[start_match];
		match->positions[ctr] = match->positions[start_match];
		match->strands[ctr] = match->strands[start_match];
		match->masks[ctr] = match->masks[start_match];
		start_match++;
		ctr++;
	}
	while(start_new <= end_new) {
		match->contigs[ctr] = c[start_new];
		match->positions[ctr] = p[start_new];
		match->strands[ctr] = s[start_new];
		match->masks[ctr] = bfast_rg_match_create_mask(bfast_mask, space, offset, s[start_new], match->read_length);
		start_new++;
		ctr++;
	}

	// free
	free(c); c=NULL;
	free(p); p=NULL;
	free(s); s=NULL;

	// remove duplicates
	prev_index=0;
	for(i=1;i<match->num_entries;i++) {
		if(0==cmp_c_p_s(match->contigs[prev_index], match->positions[prev_index], match->strands[prev_index],
					match->contigs[i], match->positions[i], match->strands[i])) {
			// union masks
			for(j=0;j<GETMASKNUMBYTESFROMLENGTH(match->read_length);j++) {
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

void bfast_rg_match_t_init(bfast_rg_match_t *match, bfast_seq_t *seq)
{
	char *fn_name="bfast_rg_match_t_init";

	// copy seq info
	if(NULL != seq) {
		match->read_name_length = seq->name.l;
		match->read_length = seq->seq.l;
		match->qual_length = seq->qual.l;
		match->read_name = my_malloc(sizeof(char)*(1+seq->name.l), fn_name);
		match->read = my_malloc(sizeof(char)*(1+seq->seq.l), fn_name);
		match->qual = my_malloc(sizeof(char)*(1+seq->qual.l), fn_name);
		strcpy(match->read_name, seq->name.s);
		strcpy(match->read, seq->seq.s);
		strcpy(match->qual, seq->qual.s);
	}
	else {
		match->read_name_length = match->read_length = match->qual_length = 0;
		match->read_name = match->read = match->qual = NULL;
	}

	// initialize
	match->max_reached = 0;
	match->num_entries = 0;
	match->contigs = NULL;
	match->positions = NULL;
	match->strands = NULL;
	match->masks = NULL;
}

void bfast_rg_match_t_destroy(bfast_rg_match_t *match)
{
	int32_t i;

	// free
	free(match->read_name);
	free(match->read);
	free(match->qual);
	free(match->contigs);
	free(match->positions);
	free(match->strands);
	for(i=0;i<match->num_entries;i++) {
		free(match->masks[i]);
	}
	free(match->masks);

	// re-initialize
	match->read_name_length = 0;
	match->read_length = 0;
	match->qual_length = 0;
	match->read_name = NULL;
	match->read = NULL;
	match->qual = NULL;
	match->max_reached = 0;
	match->num_entries = 0;
	match->contigs = NULL;
	match->positions = NULL;
	match->strands = NULL;
	match->masks = NULL;
}

void bfast_rg_match_t_append(bfast_rg_match_t *dest, bfast_rg_match_t *src)
{
	char *fn_name="bfast_rg_match_t_append";
	int32_t num_entries, i, j, len;

	if(0 == src->num_entries) return;

	num_entries = dest->num_entries + src->num_entries;
	len = GETMASKNUMBYTESFROMLENGTH(dest->read_length);
	assert(0 < len);

	// reallocate memory
	dest->contigs = my_realloc(dest->contigs, sizeof(uint32_t)*num_entries, fn_name);
	dest->positions = my_realloc(dest->positions, sizeof(int32_t)*num_entries, fn_name);
	dest->strands = my_realloc(dest->strands, sizeof(char)*num_entries, fn_name);
	dest->masks = my_realloc(dest->masks, sizeof(char*)*num_entries, fn_name);

	// copy
	for(i=dest->num_entries;i<num_entries;i++) {
		// allocate mask
		dest->masks[i] = my_calloc(len, sizeof(char), fn_name);
		// copy
		dest->contigs[i] = src->contigs[i-dest->num_entries];
		dest->positions[i] = src->positions[i-dest->num_entries];
		dest->strands[i] = src->strands[i-dest->num_entries];
		for(j=0;j<len;j++) {
			dest->masks[i][j] = src->masks[i-dest->num_entries][j];
		}
	}

	dest->num_entries = num_entries;
}

void bfast_rg_match_t_clear(bfast_rg_match_t *match)
{
	int32_t i;
	// do not delete read & qual

	// free
	free(match->contigs);
	free(match->positions);
	free(match->strands);
	for(i=0;i<match->num_entries;i++) {
		free(match->masks[i]);
	}
	free(match->masks);

	// re-initialize
	match->max_reached = 0;
	match->num_entries = 0;
	match->contigs = NULL;
	match->positions = NULL;
	match->strands = NULL;
	match->masks = NULL;
}

int32_t bfast_rg_match_t_print_queue(bfast_rg_match_t *match_queue, int32_t len, gzFile fp, int32_t flush)
{
	int32_t i, prev_i, j;

	for(i=prev_i=0;i<len;i++) {
		if(i != prev_i && 0 != strcmp(match_queue[prev_i].read_name, match_queue[i].read_name)) {
			// print prev_i to i-1
			bfast_rg_matches_t_print(match_queue, prev_i, i-1, fp);
			// free
			for(j=prev_i;j<=i-1;j++) {
				free(match_queue[j].read_int); free(match_queue[j].read_rc_int);
				match_queue[j].read_int=match_queue[j].read_rc_int=NULL;
				bfast_rg_match_t_destroy(&match_queue[j]);
			}

			// update prev
			prev_i = i;
		}
	}
	if(1 == flush) {
		// print prev_i to len-1
		bfast_rg_matches_t_print(match_queue, prev_i, len-1, fp);
		for(j=prev_i;j<=len-1;j++) {
			free(match_queue[j].read_int); free(match_queue[j].read_rc_int);
			bfast_rg_match_t_destroy(&match_queue[j]);
		}
		return 0;
	}
	else {
		// copy to the front
		for(j=prev_i;j<=len-1;j++) {
			match_queue[j-prev_i].read_int = match_queue[j].read_int;
			match_queue[j].read_int = NULL;
			match_queue[j-prev_i].read_rc_int = match_queue[j].read_rc_int;
			match_queue[j].read_rc_int = NULL;
			match_queue[j-prev_i].read_name_length = match_queue[j].read_name_length;
			match_queue[j].read_name_length=0;
			match_queue[j-prev_i].read_name = match_queue[j].read_name;
			match_queue[j].read_name=NULL;
			match_queue[j-prev_i].read_length = match_queue[j].read_length;
			match_queue[j].read_length=0;
			match_queue[j-prev_i].read = match_queue[j].read;
			match_queue[j].read=NULL;
			match_queue[j-prev_i].qual_length = match_queue[j].qual_length;
			match_queue[j].qual_length=0;
			match_queue[j-prev_i].qual = match_queue[j].qual;
			match_queue[j].qual=NULL;
			match_queue[j-prev_i].max_reached = match_queue[j].max_reached;
			match_queue[j].max_reached=0;
			match_queue[j-prev_i].num_entries = match_queue[j].num_entries;
			match_queue[j].num_entries=0;
			match_queue[j-prev_i].contigs = match_queue[j].contigs;
			match_queue[j].contigs=NULL;
			match_queue[j-prev_i].positions = match_queue[j].positions;
			match_queue[j].positions=NULL;
			match_queue[j-prev_i].strands = match_queue[j].strands;
			match_queue[j].strands=NULL;
			match_queue[j-prev_i].masks = match_queue[j].masks;
			match_queue[j].masks=NULL;
		}
		return len-prev_i;
	}
}

void bfast_rg_matches_t_print(bfast_rg_match_t *match_queue, int32_t from, int32_t to, gzFile fp)
{
	char *fn_name="bfast_rg_matches_t_print";
	int32_t i, num_ends;

	num_ends = to - from + 1;
	if(num_ends <= 0) return;

	if(gzwrite64(fp, &match_queue[from].read_name_length, sizeof(int32_t))!=sizeof(int32_t) ||
			gzwrite64(fp, match_queue[from].read_name, sizeof(char)*match_queue[from].read_name_length)!=sizeof(char)*match_queue[from].read_name_length ||
			gzwrite64(fp, &num_ends, sizeof(int32_t)) != sizeof(int32_t)) {
		PrintError(fn_name, NULL, "Could not write m->readName_length, m->readName, and m->numEnds", Exit, WriteFileError);
	}

	for(i=from;i<=to;i++) {
		bfast_rg_match_t_print(&match_queue[i], fp);
	}
}

void bfast_rg_match_t_print(bfast_rg_match_t *match, gzFile fp)
{
	char *fn_name="bfast_rg_match_t_print";
	int32_t i;

	if(gzwrite64(fp, &match->read_length, sizeof(int32_t))!=sizeof(int32_t) ||
			gzwrite64(fp, &match->qual_length, sizeof(int32_t))!=sizeof(int32_t) ||
			gzwrite64(fp, match->read, sizeof(char)*match->read_length)!=sizeof(char)*match->read_length ||
			gzwrite64(fp, match->qual, sizeof(char)*match->qual_length)!=sizeof(char)*match->qual_length ||
			gzwrite64(fp, &match->max_reached, sizeof(int32_t))!=sizeof(int32_t) ||
			gzwrite64(fp, &match->num_entries, sizeof(int32_t))!=sizeof(int32_t)) {
		PrintError(fn_name, NULL, "Could not write match->read_length, match->qual_length, match->read, match->qual, match->maxReached, and match->num_entries", Exit, WriteFileError);
	}

	/* Print the contigs, positions, and strands */
	if(gzwrite64(fp, match->contigs, sizeof(uint32_t)*match->num_entries)!=sizeof(uint32_t)*match->num_entries ||
			gzwrite64(fp, match->positions, sizeof(int32_t)*match->num_entries)!=sizeof(int32_t)*match->num_entries ||
			gzwrite64(fp, match->strands, sizeof(char)*match->num_entries)!=sizeof(char)*match->num_entries) {
		PrintError(fn_name, NULL, "Could not write contigs, positions and strands", Exit, WriteFileError);
	}
	for(i=0;i<match->num_entries;i++) {
		if(gzwrite64(fp, match->masks[i], sizeof(char)*GETMASKNUMBYTESFROMLENGTH(match->read_length))!=sizeof(char)*GETMASKNUMBYTESFROMLENGTH(match->read_length)) {
			PrintError(fn_name, NULL, "Could not write masks[i]", Exit, WriteFileError);
		}
	}
}

char *bfast_rg_match_create_mask(bfast_string_t *bfast_mask, int32_t space, int32_t offset, int32_t strand, int32_t read_length)
{
	char *fn_name="bfast_rg_match_create_mask";
	int32_t i, pos, cur_byte, cur_byte_index;

	char *mask = my_calloc((size_t)GETMASKNUMBYTESFROMLENGTH(read_length), sizeof(char), fn_name);

	for(i=0;i<bfast_mask->l;i++) {
		if('1' == bfast_mask->s[i]) {
			pos = offset + i;
			if(1 == space) pos++; // adaptor
			cur_byte = GETMASKBYTE(pos);
			cur_byte_index = pos & 7;
			mask[cur_byte] |= (0x01 << cur_byte_index);
		}
	}

	return mask;
}

void bfast_rg_match_t_print_text(bfast_rg_match_t *match, FILE *fp)
{
	char *fn_name="bfast_rg_match_t_print_text";
	char *mask=NULL;
	int32_t i, j, cur_byte, cur_byte_index;
	uint8_t byte;
	fprintf(fp, "%s\t%s\t%d\t%d",
			match->read,
			match->qual,
			match->max_reached,
			match->num_entries);

	for(i=0;i<match->num_entries;i++) {

		mask = my_malloc(sizeof(char)*(1+match->read_length), fn_name);
		for(j=0;j<match->read_length;j++) {
			cur_byte = GETMASKBYTE(j);
			cur_byte_index = j & 7;
			byte = match->masks[i][cur_byte];
			byte = byte << (8 - 1 - cur_byte_index);
			byte = byte >> 7;
			mask[j] = (0 < byte) ? '1' : '0';
		}
		mask[match->read_length] = '\0';

		fprintf(fp, "\t%u\t%d\t%c\t%s",
				match->contigs[i],
				match->positions[i],
				match->strands[i],
				mask);
		free(mask);
	}
	fputc('\n', fp);

}

void bfast_rg_match_t_copy_from_bwa(bfast_rg_match_t *m, bwa_seq_t *seq, int32_t space)
{
	char *fn_name="bfast_rg_match_t_copy_from_bwa";
	int32_t i;
	m->read_name = strdup(seq->name);
	if(NULL == m->read_name) PrintError(fn_name, "m->read_name", "Could not allocate memory", Exit, MallocMemory);
	m->read_name_length = strlen(seq->name);
	m->read_length = seq->len;
	m->read = my_malloc(sizeof(char)*(1+m->read_length), fn_name);
	m->read_int = my_malloc(sizeof(uint8_t)*m->read_length, fn_name);
	m->read_rc_int = my_malloc(sizeof(uint8_t)*m->read_length, fn_name);
	if(0 == space) {
		m->read_int_length = m->read_length;
		for(i=0;i<m->read_length;i++) {
			m->read[i] = "ACGTN"[seq->seq[i]];
			m->read_int[i] = seq->seq[i];
			m->read_rc_int[m->read_length-i-1] = bfast_nt_int_to_nt_int_comp[seq->seq[i]]; // reverse compliment
		}
		m->read[i]='\0';
	}
	else {
		m->read_int_length = m->read_length-1; // no adapter
		m->read[0] = "ACGTN"[seq->seq[0]]; // adapter
		for(i=1;i<m->read_length;i++) {
			m->read[i] = "01234"[seq->seq[i]];
			// ignore first base
			m->read_int[i-1] = seq->seq[i];
			m->read_rc_int[m->read_length-i-1] = seq->seq[i]; // reverse
		}
		m->read[i]='\0';
	}
	m->qual = strdup((char*)seq->qual);
	if(NULL == m->qual) PrintError(fn_name, "m->qual", "Could not allocate memory", Exit, MallocMemory);
	m->qual_length = strlen((char*)seq->qual);
	m->tid = seq->tid;
}

int32_t bfast_rg_match_read(bwa_seqio_t *bs, int n_needed, bfast_rg_match_t *m, int *n, int space, int trim_qual)
{
	//char *fn_name="bfast_rg_match_t";
	bwa_seq_t *seqs=NULL;
	int32_t i, prev_n;

	// Use BWA built-in
	prev_n = (*n);
	seqs = bwa_read_seq(bs, n_needed - prev_n, n, 1-space, trim_qual);

	for(i=0;i<(*n);i++) {
		bfast_rg_match_t_init(&m[i+prev_n], NULL);
		seq_reverse(seqs[i].len, seqs[i].seq, 0); // reverse back sequence
		bfast_rg_match_t_copy_from_bwa(&m[i+prev_n], &seqs[i], space);
	}

	// Free
	bwa_free_read_seq((*n), seqs);

	(*n) += prev_n;

	if((*n) == prev_n) return 1;
	else return 0;
}
