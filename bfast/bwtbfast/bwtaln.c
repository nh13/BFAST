#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <assert.h>
#include "../aflib.h"
#include "bwtaln.h"
#include "bwtgap.h"
#include "utils.h"
#include "util.h"
#include "bwtbfast_aux.h"
#include "../BLibDefinitions.h"

#ifdef HAVE_LIBPTHREAD
#define THREAD_BLOCK_SIZE 1024
#include <pthread.h>
static pthread_mutex_t g_seq_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

gap_opt_t *gap_init_opt()
{
	gap_opt_t *o;
	o = (gap_opt_t*)calloc(1, sizeof(gap_opt_t));
	/* IMPORTANT: s_mm*10 should be about the average base error
	   rate. Voilating this requirement will break pairing! */
	o->s_mm = 3; o->s_gapo = 11; o->s_gape = 4;
	o->max_diff = -1; o->max_gapo = 1; o->max_gape = 6;
	o->indel_end_skip = 5; o->max_del_occ = 10; o->max_entries = 2000000;
	o->mode = BWA_MODE_GAPE | BWA_MODE_COMPREAD;
	o->seed_len = 32; o->max_seed_diff = 2;
	o->max_num_m = MAX_NUM_MATCHES;
	o->fnr = 0.04;
	o->n_threads = 1;
	o->max_top2 = 30;
	o->trim_qual = 0;
	return o;
}

int bwa_cal_maxdiff(int l, double err, double thres)
{
	double elambda = exp(-l * err);
	double sum, y = 1.0;
	int k, x = 1;
	for (k = 1, sum = elambda; k < 1000; ++k) {
		y *= l * err;
		x *= k;
		sum += elambda * y / x;
		if (1.0 - sum < thres) return k;
	}
	return 2;
}

// width must be filled as zero
static int bwt_cal_width(const bwt_t *rbwt, int len, const ubyte_t *str, bwt_width_t *width)
{
	bwtint_t k, l, ok, ol;
	int i, bid;
	bid = 0;
	k = 0; l = rbwt->seq_len;
	for (i = 0; i < len; ++i) {
		ubyte_t c = str[i];
		if (c < 4) {
			bwt_2occ(rbwt, k - 1, l, c, &ok, &ol);
			k = rbwt->L2[c] + ok + 1;
			l = rbwt->L2[c] + ol;
		}
		if (k > l || c > 3) { // then restart
			k = 0;
			l = rbwt->seq_len;
			++bid;
		}
		width[i].w = l - k + 1;
		width[i].bid = bid;
	}
	width[len].w = 0;
	width[len].bid = ++bid;
	return bid;
}

void bwa_cal_sa_reg_gap(int tid, bwt_t *const bwt[2], int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt)
{
	int i, max_l = 0, max_len;
	gap_stack_t *stack;
	bwt_width_t *w[2], *seed_w[2];
	const ubyte_t *seq[2];
	int32_t len;
	gap_opt_t local_opt = *opt;

	// initiate priority stack
	for (i = max_len = 0; i != n_seqs; ++i)
		if (seqs[i].len > max_len) max_len = seqs[i].len;
	if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
	if (local_opt.max_diff < local_opt.max_gapo) local_opt.max_gapo = local_opt.max_diff;
	stack = gap_init_stack(local_opt.max_diff, local_opt.max_gapo, local_opt.max_gape, &local_opt);

	seed_w[0] = (bwt_width_t*)calloc(opt->seed_len+1, sizeof(bwt_width_t));
	seed_w[1] = (bwt_width_t*)calloc(opt->seed_len+1, sizeof(bwt_width_t));
	w[0] = w[1] = 0;
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p = seqs + i;
#ifdef HAVE_LIBPTHREAD
		if (opt->n_threads > 1) {
			pthread_mutex_lock(&g_seq_lock);
			if (p->tid < 0) { // unassigned
				int j;
				for (j = i; j < n_seqs && j < i + THREAD_BLOCK_SIZE; ++j)
					seqs[j].tid = tid;
			} else if (p->tid != tid) {
				pthread_mutex_unlock(&g_seq_lock);
				continue;
			}
			pthread_mutex_unlock(&g_seq_lock);
		}
#endif
		p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;
		seq[0] = p->seq; seq[1] = p->rseq;
		len = p->len;
		if(!(BWA_MODE_COMPREAD & opt->mode)) { // adjust color space
			// ignore adaptor and first color
			len -= 2; 
			//seq[0] += 2;
			// reverse adaptor and first color ignored by length decrement
		}
		if (max_l < len) { // should not go here?
			max_l = len;
			w[0] = (bwt_width_t*)calloc(max_l + 1, sizeof(bwt_width_t));
			w[1] = (bwt_width_t*)calloc(max_l + 1, sizeof(bwt_width_t));
		}
		bwt_cal_width(bwt[0], len, seq[0], w[0]);
		bwt_cal_width(bwt[1], len, seq[1], w[1]);
		if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(len, BWA_AVG_ERR, opt->fnr);
		local_opt.seed_len = opt->seed_len < len? opt->seed_len : 0x7fffffff;
		if (len > opt->seed_len) {
			bwt_cal_width(bwt[0], opt->seed_len, seq[0] + (len - opt->seed_len), seed_w[0]);
			bwt_cal_width(bwt[1], opt->seed_len, seq[1] + (len - opt->seed_len), seed_w[1]);
		}
		// core function
		p->aln = bwt_match_gap(bwt, len, seq, w, len <= opt->seed_len? 0 : seed_w, &local_opt, &p->n_aln, stack);
		// store the alignment
		// free(p->name); we need the name 
		//free(p->seq); free(p->rseq); free(p->qual);
		//p->seq = p->rseq = p->qual = 0;
	}
	free(seed_w[0]); free(seed_w[1]);
	free(w[0]); free(w[1]);
	gap_destroy_stack(stack);
}

#ifdef HAVE_LIBPTHREAD
typedef struct {
	int tid;
	bwt_t *bwt[2];
	int n_seqs;
	bwa_seq_t *seqs;
	const gap_opt_t *opt;
} thread_aux_t;

static void *worker(void *data)
{
	thread_aux_t *d = (thread_aux_t*)data;
	bwa_cal_sa_reg_gap(d->tid, d->bwt, d->n_seqs, d->seqs, d->opt);
	return 0;
}
#endif

void bwa_aln_core(const char *prefix, const char *fn_fa, const gap_opt_t *opt)
{
	char *fn_name="bwa_aln_core";
	int i, n_seqs, tot_seqs = 0;
	bwa_seq_t *seqs;
	bwa_seqio_t *ks;
	clock_t t;
	bwt_t *bwt[2];
	gzFile fp_out=NULL;
	bntseq_t *bns=NULL;
	int32_t space = (opt->mode & BWA_MODE_COMPREAD) ? NTSpace : ColorSpace;
	char *str = (char*)calloc(strlen(prefix) + 16, 1);
	bfast_rg_match_t *matches=NULL;
	int32_t n_matches = 0;

	// initialization
	ks = bwa_seq_open(fn_fa, AFILE_NO_COMPRESSION);
	//ks = bwa_seq_open(fn_fa, AFILE_GZ_COMPRESSION);

	{ // load BWT

		strcpy(str, prefix); strcat(str, "."); strcat(str, SPACENAME(space));
		strcat(str, ".bwt");  bwt[0] = bwt_restore_bwt(str);

		strcpy(str, prefix); strcat(str, "."); strcat(str, SPACENAME(space));
		strcat(str, ".rbwt"); bwt[1] = bwt_restore_bwt(str);

		strcpy(str, prefix); strcat(str, "."); strcat(str, SPACENAME(space));
		bns = bns_restore(str);

		strcpy(str, prefix); strcat(str, "."); strcat(str, SPACENAME(space));
		strcat(str, ".sa"); bwt_restore_sa(str, bwt[0]);

		strcpy(str, prefix); strcat(str, "."); strcat(str, SPACENAME(space));
		strcat(str, ".rsa"); bwt_restore_sa(str, bwt[1]);

		free(str);
	}

	fp_out = err_xzopen_core(fn_name, "-", "wb"); 

	// core loop
	while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, 1-space, opt->trim_qual)) != 0) {
		tot_seqs += n_seqs;
		t = clock();

		fprintf(stderr, "[bwa_aln_core] calculate SA coordinate... ");

#ifdef HAVE_LIBPTHREAD
		if (opt->n_threads <= 1) { // no multi-threading at all
			bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
		} else {
			pthread_t *tid;
			pthread_attr_t attr;
			thread_aux_t *data;
			int j;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
			tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
			for (j = 0; j < opt->n_threads; ++j) {
				data[j].tid = j; data[j].bwt[0] = bwt[0]; data[j].bwt[1] = bwt[1];
				data[j].n_seqs = n_seqs; data[j].seqs = seqs; data[j].opt = opt;
				pthread_create(&tid[j], &attr, worker, data + j);
			}
			for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
			free(data); free(tid);
		}
#else
		bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
#endif

		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		t = clock();
		fprintf(stderr, "[bwa_aln_core] write to the disk... ");
		matches = my_realloc(matches, sizeof(bfast_rg_match_t)*(n_matches + n_seqs), fn_name);
		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p = seqs + i;
			// BWA to BMF
			{
				int ctr=0, n_skipped=0, seqid;
				bwtint_t j, k;
				int64_t pos;

				assert(NULL != p->seq);
				seq_reverse(p->len, p->seq, 0);
				matches[i+n_matches].read_name_length = strlen(p->name);
				matches[i+n_matches].read_name = my_malloc(sizeof(char)*(1+matches[i+n_matches].read_name_length), fn_name);
				matches[i+n_matches].read_int = matches[i+n_matches].read_rc_int = NULL;
				strcpy(matches[i+n_matches].read_name, p->name);
				matches[i+n_matches].read_length = p->len;
				matches[i+n_matches].read = my_malloc(sizeof(char)*(1+p->len), fn_name);
				for(j=0;j<p->len;j++) {
					if(0 == j || NTSpace == space) {
						matches[i+n_matches].read[j] = "ACGTN"[p->seq[j]];
					}
					else {
						matches[i+n_matches].read[j] = "01234"[p->seq[j]];
					}
				}
				matches[i+n_matches].qual_length = (NTSpace == space) ? matches[i+n_matches].read_length : matches[i+n_matches].read_length-1;
				matches[i+n_matches].qual = my_malloc(sizeof(char)*(1+matches[i+n_matches].qual_length), fn_name);
				strcpy(matches[i+n_matches].qual, (char*)p->qual);
				matches[i+n_matches].max_reached = 0;

				matches[i+n_matches].num_entries = 0;
				for(j=0;j<p->n_aln;j++) {
					matches[i+n_matches].num_entries += p->aln[j].l - p->aln[j].k + 1;
				}
				if(opt->max_num_m < matches[i+n_matches].num_entries) {
					// Ignore
					matches[i+n_matches].num_entries = 0;
					matches[i+n_matches].contigs = NULL;
					matches[i+n_matches].positions = NULL;
					matches[i+n_matches].strands = NULL;
					matches[i+n_matches].masks = NULL;
				}
				else {
					matches[i+n_matches].contigs = my_malloc(sizeof(uint32_t)*matches[i+n_matches].num_entries, fn_name);
					matches[i+n_matches].positions = my_malloc(sizeof(int32_t)*matches[i+n_matches].num_entries, fn_name);
					matches[i+n_matches].strands = my_malloc(sizeof(char)*matches[i+n_matches].num_entries, fn_name);
					matches[i+n_matches].masks = my_malloc(sizeof(char*)*matches[i+n_matches].num_entries, fn_name);

					for(j=ctr=0;j<p->n_aln;j++) {
						for(k=p->aln[j].k;k<=p->aln[j].l;k++) {
							// p->aln[j].a == 1 means it is the reverse 

							// use suffix array
							if(p->aln[j].a) { // reverse
								pos = bwt_sa(bwt[0], k);
							}
							else {
								pos = (int64_t) (bwt[1]->seq_len - (bwt_sa(bwt[1], k) + p->len));
							}
							
							// weird cases when [pos (pac_coor) >= bns->l_pac]
							if(pos >= bns->l_pac) {
								pos=0;
							}

							// Get sequence id
							bns_coor_pac2real(bns, pos, 1, &seqid);
							
							pos = pos - bns->anns[seqid].offset;

							// check if hit spans two sequences
							if(bns->anns[seqid].len < pos) {
								n_skipped++;
								continue;
							}

							matches[i+n_matches].contigs[ctr] = seqid+1; //contig
							matches[i+n_matches].positions[ctr] = pos+1; // position
							matches[i+n_matches].strands[ctr] = (0 == p->aln[j].a) ? FORWARD : REVERSE; // strand
							matches[i+n_matches].masks[ctr] = my_calloc(GETMASKNUMBYTESFROMLENGTH(matches[i+n_matches].read_length), sizeof(char), fn_name); // mask

							// adjust position
							if(ColorSpace == space) {
								if(FORWARD == matches[i+n_matches].strands[ctr]) {
									matches[i+n_matches].positions[ctr]+=1;
								}
								else {
									matches[i+n_matches].positions[ctr]--;
								}
							}
							ctr++;
						}
					}
					if(0 < n_skipped) {
						matches[i+n_matches].num_entries -= n_skipped;
						if (matches[i+n_matches].num_entries != 0) {
							matches[i+n_matches].contigs = my_realloc(matches[i+n_matches].contigs, sizeof(uint32_t)*matches[i+n_matches].num_entries, fn_name);
							matches[i+n_matches].positions = my_realloc(matches[i+n_matches].positions, sizeof(int32_t)*matches[i+n_matches].num_entries, fn_name);
							matches[i+n_matches].strands = my_realloc(matches[i+n_matches].strands, sizeof(char)*matches[i+n_matches].num_entries, fn_name);
							matches[i+n_matches].masks = my_realloc(matches[i+n_matches].masks, sizeof(char*)*matches[i+n_matches].num_entries, fn_name);
						}
					}
				}
			}
		} // Finished converting to BF the BWA alignments 
		n_matches += n_seqs;
		// print
		n_matches = bfast_rg_match_t_print_queue(matches, n_matches, fp_out, 0);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		bwa_free_read_seq(n_seqs, seqs);
		fprintf(stderr, "[bwa_aln_core] %d sequences have been processed.\n", tot_seqs);
	} // There are not more reads to process

	// flush
	n_matches = bfast_rg_match_t_print_queue(matches, n_matches, fp_out, 1);
	assert(0 == n_matches);

	gzclose(fp_out);

	// destroy
	bwt_destroy(bwt[0]); bwt_destroy(bwt[1]);
	bns_destroy(bns);
	bwa_seq_close(ks);
}

int bwa_aln(int argc, char *argv[])
{
	int c, opte = -1;
	gap_opt_t *opt;

	opt = gap_init_opt();
	while ((c = getopt(argc, argv, "n:o:e:i:d:l:k:cLR:m:t:NM:O:E:q:f:")) >= 0) {
		switch (c) {
			case 'n':
				if (strstr(optarg, ".")) opt->fnr = atof(optarg), opt->max_diff = -1;
				else opt->max_diff = atoi(optarg), opt->fnr = -1.0;
				break;
			case 'o': opt->max_gapo = atoi(optarg); break;
			case 'e': opte = atoi(optarg); break;
			case 'M': opt->s_mm = atoi(optarg); break;
			case 'O': opt->s_gapo = atoi(optarg); break;
			case 'E': opt->s_gape = atoi(optarg); break;
			case 'd': opt->max_del_occ = atoi(optarg); break;
			case 'i': opt->indel_end_skip = atoi(optarg); break;
			case 'l': opt->seed_len = atoi(optarg); break;
			case 'k': opt->max_seed_diff = atoi(optarg); break;
			case 'm': opt->max_entries = atoi(optarg); break;
			case 't': opt->n_threads = atoi(optarg); break;
			case 'L': opt->mode |= BWA_MODE_LOGGAP; break;
			case 'R': opt->max_top2 = atoi(optarg); break;
			case 'q': opt->trim_qual = atoi(optarg); break;
			case 'c': opt->mode &= ~BWA_MODE_COMPREAD; break;
			case 'N': opt->mode |= BWA_MODE_NONSTOP; opt->max_top2 = 0x7fffffff; break;
			case 'f': freopen(optarg, "wb", stdout); break;
			default: return 1;
		}
	}
	if (opte > 0) {
		opt->max_gape = opte;
		opt->mode &= ~BWA_MODE_GAPE;
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa aln [options] <ref.fa> <in.fq>\n\n");
		fprintf(stderr, "Options: -n NUM    max #diff (int) or missing prob under %.2f err rate (float) [%.2f]\n",
				BWA_AVG_ERR, opt->fnr);
		fprintf(stderr, "         -o INT    maximum number or fraction of gap opens [%d]\n", opt->max_gapo);
		fprintf(stderr, "         -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]\n");
		fprintf(stderr, "         -i INT    do not put an indel within INT bp towards the ends [%d]\n", opt->indel_end_skip);
		fprintf(stderr, "         -d INT    maximum occurrences for extending a long deletion [%d]\n", opt->max_del_occ);
		fprintf(stderr, "         -l INT    seed length [%d]\n", opt->seed_len);
		fprintf(stderr, "         -k INT    maximum differences in the seed [%d]\n", opt->max_seed_diff);
		fprintf(stderr, "         -m INT    maximum entries in the queue [%d]\n", opt->max_entries);
		fprintf(stderr, "         -t INT    number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "         -M INT    mismatch penalty [%d]\n", opt->s_mm);
		fprintf(stderr, "         -O INT    gap open penalty [%d]\n", opt->s_gapo);
		fprintf(stderr, "         -E INT    gap extension penalty [%d]\n", opt->s_gape);
		fprintf(stderr, "         -R INT    stop searching when there are >INT equally best hits [%d]\n", opt->max_top2);
		fprintf(stderr, "         -q INT    quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, opt->trim_qual);
		fprintf(stderr, "         -c        input sequences are in the color space\n");
		fprintf(stderr, "         -L        log-scaled gap penalty for long deletions\n");
		fprintf(stderr, "         -N        non-iterative mode: search for all n-difference hits (slooow)\n");
		fprintf(stderr, "         -f FILE   file to write output to instead of stdout\n");
		fprintf(stderr, "         -X INT    maximum number of matches to output [%d]\n", opt->max_num_m);
		fprintf(stderr, "\n");
		return 1;
	}
	if (opt->fnr > 0.0) {
		int i, k;
		for (i = 17, k = 0; i <= 250; ++i) {
			int l = bwa_cal_maxdiff(i, BWA_AVG_ERR, opt->fnr);
			if (l != k) fprintf(stderr, "[bwa_aln] %dbp reads: max_diff = %d\n", i, l);
			k = l;
		}
	}
	bwa_aln_core(argv[optind], argv[optind+1], opt);
	free(opt);
	return 0;
}

/* rgoya: Temporary clone of aln_path2cigar to accomodate for bwa_cigar_t,
 * __cigar_op and __cigar_len while keeping stdaln stand alone */
bwa_cigar_t *bwa_aln_path2cigar(const path_t *path, int path_len, int *n_cigar)
{
	uint32_t *cigar32;
	bwa_cigar_t *cigar;
	int i;
	cigar32 = aln_path2cigar32((path_t*) path, path_len, n_cigar);
	cigar = (bwa_cigar_t*)cigar32;
	for (i = 0; i < *n_cigar; ++i)
		cigar[i] = __cigar_create( (cigar32[i]&0xf), (cigar32[i]>>4) );
	return cigar;
}
