#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include <limits.h>
#include <assert.h>
#include <pthread.h>
#include "config.h"
#include "BError.h"
#include "BLibDefinitions.h"
#include "bwt.h"
#include "util.h"
#include "aflib.h"
#include "occ_results.h"
#include "bwtbfast_aux.h"
#include "bwtbfast.h"

#define Name "bwtbfast"
#define QUIET 0

#ifndef THREAD_BLOCK_SIZE
#define THREAD_BLOCK_SIZE 102
#endif

#ifdef HAVE_LIBPTHREAD
static pthread_mutex_t bfast_g_seq_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

int bwtbfast_usage(int32_t space, int32_t max_key_matches, int32_t max_num_matches, int32_t num_threads, int32_t queue_length)
{
	fprintf(stderr, "\nUsage:%s %s [options]\n", PACKAGE_NAME, Name);
	fprintf(stderr, "\t-f\t\tSpecifies the file name of the FASTA reference genome\n");
	fprintf(stderr, "\t-r\t\tread file name\n");
	fprintf(stderr, "\t-m\t\tbfast mask\n");
	fprintf(stderr, "\t-A\t\tspace (0: NT space 1: Color space) [%d]\n", space);
	fprintf(stderr, "\t-K\t\tmaximum number of matches per key before a key is ignored [%d]\n", max_key_matches);
	fprintf(stderr, "\t-M\t\tmaximum number of matches per read before a read is ignored [%d]\n", max_num_matches);
	fprintf(stderr, "\t-n\t\tnumber of threads [%d]\n", num_threads);
	fprintf(stderr, "\t-Q\t\tnumber of reads to process at one time [%d]\n", queue_length);
	fprintf(stderr, "\t-h\t\tprints this help message\n");
	fprintf(stderr, "\nsend bugs to %s\n",
			PACKAGE_BUGREPORT);

	return 1;
}

int bwtbfast(int argc, char *argv[])
{
	int c;
	char *fn_name="bwtbfast";
	char *ref_fn=NULL;
	char *read_fn=NULL;
	int32_t space = 0;
	int32_t max_key_matches = 8;
	int32_t max_num_matches = 384;
	int32_t num_threads = 0;
	int32_t queue_length = 0x40000;
	bfast_masks_t masks;

	bfast_masks_t_init(&masks);

	// Get parameters
	while((c = getopt(argc, argv, "f:m:n:r:A:K:M:Q:h")) >= 0) {
		switch(c) {
			case 'f': ref_fn=strdup(optarg); break;
			case 'r': read_fn=strdup(optarg); break;
			case 'm': bfast_masks_t_add(&masks, optarg); break;
			case 'n': num_threads=atoi(optarg); break;
			case 'A': space=atoi(optarg); break;
			case 'K': max_key_matches=atoi(optarg); break;
			case 'M': max_num_matches=atoi(optarg); break;
			case 'Q': queue_length=atoi(optarg); break;
			case 'h': 
					  return bwtbfast_usage(space,
							  max_key_matches,
							  max_num_matches,
							  num_threads,
							  queue_length);
			default: fprintf(stderr, "Unrecognized option: -%c\n", c); return 1;
		}
	}

	if(1 == argc || optind != argc) {
		return bwtbfast_usage(space,
				max_key_matches,
				max_num_matches,
				num_threads,
				queue_length);
	}

	// check parameters
	if(NULL == ref_fn) {
		PrintError(fn_name, "-f", "Option is required", Exit, OutOfRange);
	}
	if(NULL == read_fn) {
		PrintError(fn_name, "-r", "Option is required", Exit, OutOfRange);
	}
	if(masks.n <= 0) {
		PrintError(fn_name, "-m", "Option is required", Exit, OutOfRange);
	}
	if(space < 0 || 1 < space) {
		PrintError(fn_name, "-A", "Option is out of range", Exit, OutOfRange);
	}
	if(num_threads < 0) {
		PrintError(fn_name, "-n", "Option is out of range", Exit, OutOfRange);
	}
	if(queue_length <= 0) {
		PrintError(fn_name, "-Q", "Option is out of range", Exit, OutOfRange);
	}

	// bfast
	bwtbfast_core(ref_fn, read_fn, &masks, space, max_key_matches, max_num_matches, num_threads, queue_length);

	// free file names
	free(ref_fn);
	free(read_fn);
	bfast_masks_t_destroy(&masks);

	return 0;
}

void bwtbfast_core(char *ref_fn, char *read_fn, bfast_masks_t *masks, int32_t space, int32_t max_key_matches, int32_t max_num_matches, int32_t n_threads, int32_t queue_length)
{
	char *fn_name="bwtbfast_core";
	int i;
	int n_matches, tot_matches= 0;
	bfast_rg_match_t *matches=NULL;
	bwa_seqio_t *ks;
	clock_t t;
	bwt_t *bwt=NULL;
	bntseq_t *bns=NULL;
	gzFile fp_out=NULL;

	// initialization
	ks = bwa_seq_open(read_fn, AFILE_NO_COMPRESSION);
	fp_out = gzdopen(fileno(stdout), "wb");
	if(NULL == fp_out) PrintError(fn_name, "stdout", "Coult not open for writing", Exit, OpenFileError);

	{ 
		// load BWT
		char str[1024]="\0";
		char prefix[1024]="\0";
		strcpy(prefix, ref_fn); 
		strcat(prefix, "."); strcat(prefix, SPACENAME(space));
		strcpy(str, prefix); strcat(str, ".bwt");
		bwt = bwt_restore_bwt(str);
		// load BNS
		bns = bns_restore(prefix);
		// load SA
		strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
	}

	while ((matches = bfast_rg_match_read(ks, queue_length, &n_matches, space, 0)) != 0) {
		tot_matches += n_matches;
		t = clock();

		fprintf(stderr, "[bwtbfast_core] matching... ");
#ifdef HAVE_LIBPTHREAD
		if(n_threads <= 1) {
			// no threads
			bwtbfast_core_worker(0, bwt, bns, n_matches, matches, space, 1, masks, max_key_matches, max_num_matches);
		}
		else {
			// threads
			pthread_t *tid;
			pthread_attr_t attr;
			bwtbfast_thread_t *data;
			int j;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			data = (bwtbfast_thread_t*)calloc(n_threads, sizeof(bwtbfast_thread_t));
			tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
			for (j = 0; j < n_threads; ++j) {
				data[j].tid = j; data[j].bwt = bwt; data[j].bns = bns;
				data[j].n_matches = n_matches; data[j].matches = matches; 
				data[j].space = space;
				data[j].n_threads = n_threads;
				data[j].masks = masks; 
				data[j].max_key_matches = max_key_matches; 
				data[j].max_num_matches = max_num_matches;
				pthread_create(&tid[j], &attr, bwtbfast_thread_worker, data + j);
			}
			for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
			free(data); free(tid);
		}
#else
		// no threads
		bwtbfast_core_worker(0, bwt, bns, n_matches, matches, space, 1, masks, max_key_matches, max_num_matches);
#endif

		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		// Print
		t = clock();
		fprintf(stderr, "[bwtbfast_core] write to the disk... ");
		bfast_rg_match_t_print_queue(matches, n_matches, fp_out);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		// Free
		for(i=0;i<n_matches;i++) {
			bfast_rg_match_t_destroy(&matches[i]);
		}
		free(matches); matches = NULL;

		fprintf(stderr, "[bwtbfast_core] %d sequences have been processed.\n", tot_matches);

	}

	// destroy
	bwt_destroy(bwt);
	bns_destroy(bns);
	bwa_seq_close(ks);
	gzclose(fp_out);
}

#ifdef HAVE_LIBPTHREAD
void *bwtbfast_thread_worker(void *data)
{
	bwtbfast_thread_t *d = (bwtbfast_thread_t*)data;
	bwtbfast_core_worker(d->tid, d->bwt, d->bns, d->n_matches, d->matches, d->space, d->n_threads, d->masks, d->max_key_matches, d->max_num_matches);
	return 0;
}
#endif

void bwtbfast_core_worker(int tid, bwt_t *bwt, bntseq_t *bns, int n_matches, bfast_rg_match_t *matches, int32_t space, int n_threads, bfast_masks_t *masks, int32_t max_key_matches, int32_t max_num_matches) 
{
	//char *fn_name="bwtbfast_core_worker";
	int32_t i, j, m;
	occ_results_t results_f, results_r;

	for (i = 0; i != n_matches; ++i) {
#ifdef HAVE_LIBPTHREAD
		if(1 < n_threads) {
			pthread_mutex_lock(&bfast_g_seq_lock);
			if (matches[i].tid < 0) { // unassigned
				int j;
				for (j = i; j < n_matches && j < i + THREAD_BLOCK_SIZE; ++j)
					matches[j].tid = tid;
			} else if (matches[i].tid != tid) {
				pthread_mutex_unlock(&bfast_g_seq_lock);
				continue;
			}
			pthread_mutex_unlock(&bfast_g_seq_lock);
		}
#endif
		// bfast
		for(m=0;0 == matches[i].max_reached && m<masks->n;m++) {
			for(j=0;j<matches[i].read_int_length-masks->masks[m].l+1;j++) {
				occ_results_t_init(&results_f);
				occ_results_t_init(&results_r);

				// forward
				bwt_t_get_with_mask(bwt, 
						(matches[i].read_int + j), 
						masks->masks[m].s, 
						masks->masks[m].l, 
						&results_f);
				// reverse 
				bwt_t_get_with_mask(bwt, 
						(matches[i].read_rc_int + matches[i].read_int_length - masks->masks[m].l - j), 
						masks->masks[m].s, 
						masks->masks[m].l, 
						&results_r);

				// copy over
				if(0 < results_f.n + results_r.n && results_f.n + results_r.n <= max_key_matches) {
					// Add to BFAST rgmatches
					bfast_rg_match_t_copy_results(&matches[i], bwt, bns, &masks->masks[m], space, j, &results_f, &results_r);
				}

				// destroy
				occ_results_t_destroy(&results_r);
				occ_results_t_destroy(&results_f);

				// check if too many matches.  If so, mark and exit
				if(max_num_matches < matches[i].num_entries) {
					// mark and free
					bfast_rg_match_t_clear(&matches[i]);
					matches[i].max_reached = 1;
					break;
				}

			}
		}
		// free -> we don't need these anymore
		free(matches[i].read_int); matches[i].read_int = NULL;
		free(matches[i].read_rc_int); matches[i].read_rc_int = NULL;
	}
}

void bwt_t_get_with_mask_helper(const bwt_t *bwt, uint8_t *w_int, char *mask, int32_t len, int32_t r_l_old, int32_t r_u_old, occ_results_t *r)
{
    int32_t base, r_l_new, r_u_new;
    uint32_t occ_r_l[4], occ_r_u[4];
    if(len <= 0) {
        //fprintf(stdout, "Adding %d %d\n", r_l_old, r_u_old); // HERE
        occ_results_t_add(r, r_l_old, r_u_old);
        return;
    }
    bwt_2occ4(bwt, r_l_old-1, r_u_old, occ_r_l, occ_r_u);
    for(base=0;base<4;base++) {
        r_l_new = bwt->L2[base] + occ_r_l[base] + 1;
        r_u_new = bwt->L2[base] + occ_r_u[base];
        if(r_l_new <= r_u_new) {
            //fprintf(stdout, "base=%c %d\n", "ACGT"[base], len); // HERE
            if(base == w_int[len-1]) {
                bwt_t_get_with_mask_helper(bwt, w_int, mask, len-1, r_l_new, r_u_new, r);
            }
            else if('0' == mask[len-1]) {
                bwt_t_get_with_mask_helper(bwt, w_int, mask, len-1, r_l_new, r_u_new, r);
            }
        }
    }
}

void bwt_t_get_with_mask(const bwt_t *bwt, uint8_t *w_int, char *mask, int32_t len, occ_results_t *r)
{
    bwt_t_get_with_mask_helper(bwt, w_int, mask, len, 0, bwt->seq_len, r);
}
