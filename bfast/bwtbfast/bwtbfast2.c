#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include <limits.h>
#include <assert.h>
#include <pthread.h>
#include <config.h>
#include "BError.h"
#include "BLibDefinitions.h"
#include "util.h"
#include "occ_results.h"
#include "kseq.h"
#include "aflib.h"
#include "bwt.h"
#include "bwtbfast_aux.h"
#include "bwtbfast2_aux.h"
#include "bwtbfast2.h"

#define Name "bwtbfast2"
#define QUIET 0

#ifndef THREAD_BLOCK_SIZE
#define THREAD_BLOCK_SIZE 102
#endif

#ifdef HAVE_LIBPTHREAD
static pthread_mutex_t bfast_g_seq_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

int bwtbfast2_usage(int32_t alg, int32_t space, int32_t seed_len, int32_t max_mm, int32_t num_threads, int32_t max_hits, int32_t queue_length)
{
	fprintf(stderr, "\nUsage:%s %s [options]\n", PACKAGE_NAME, Name);
	fprintf(stderr, "\n=========== Input Files =============================================================\n");
	fprintf(stderr, "\t-f\tFILE\tSpecifies the file name of the FASTA reference genome (or prefix)\n");
	fprintf(stderr, "\t-r\tFILE\tread file name\n");
	fprintf(stderr, "\t-j\t\tSpecifies that the input reads are bz2 compressed (bzip2)\n");
	fprintf(stderr, "\t-z\t\tSpecifies that the input reads are gz compressed (gzip)\n");
	fprintf(stderr, "\n=========== Algorithm Options =======================================================\n");
	fprintf(stderr, "\t-A\tINT\tspace (0: NT space 1: Color space) [%d]\n", space);
	fprintf(stderr, "\t-s\tINT\tSpecifies the read to begin with (skip the first startReadNum-1 reads)\n");
	fprintf(stderr, "\t-e\tINT\tSpecifies the last read to use (inclusive)\n");
	fprintf(stderr, "\t-a\tINT\talgorithm mode 0: first 1: all min mm 2: all min mm+1 3: all [%d]\n", alg);
	fprintf(stderr, "\t-l\tINT\tseed length [%d]\n", seed_len);
	fprintf(stderr, "\t-m\tINT\tmaximum number of mismatches [%d]\n", max_mm);
	fprintf(stderr, "\t-M\tINT\tstop searching when there are >INT equally best hits [%d]\n", max_hits);
	fprintf(stderr, "\t-n\tINT\tnumber of threads [%d]\n", num_threads);
	fprintf(stderr, "\t-Q\tINT\tnumber of reads to process at one time [%d]\n", queue_length);
	// TODO
	//fprintf(stderr, "\n=========== Output Options ==========================================================\n");
	//fprintf(stderr, "\t-t\t\tSpecifies to output timing information\n");
	fprintf(stderr, "\n=========== Miscellaneous Options ===================================================\n");
	// TODO
	//fprintf(stderr, "\t-p\t\tPrint program parameters\n");
	fprintf(stderr, "\t-h\t\tDisplay usage summary\n");
	fprintf(stderr, "\nsend bugs to %s\n",
			PACKAGE_BUGREPORT);

	return 1;
}

int bwtbfast2(int argc, char *argv[])
{
	int c;
	char *fn_name="bwtbfast2";

	// TODO: these could all be put in a struct
	char *ref_fn=NULL;
	char *read_fn=NULL;
	int32_t compression=AFILE_NO_COMPRESSION;
	int32_t alg=1;
	int32_t space = 0;
	int32_t start_read_num = 1;
	int32_t end_read_num = INT_MAX;
	int32_t seed_len = 22;
	int32_t max_mm = 0;
	int32_t max_hits = 30;
	int32_t num_threads = 0;
	int32_t queue_length = 0x40000;

	// Get parameters
	while((c = getopt(argc, argv, "a:f:l:m:n:r:A:Q:M:hjz")) >= 0) {
		switch(c) {
			case 'a': alg=atoi(optarg); break;
			case 'f': ref_fn=strdup(optarg); break;
			case 'r': read_fn=strdup(optarg); break;
			case 'l': seed_len=atoi(optarg); break;
			case 'm': max_mm=atoi(optarg); break;
			case 'n': num_threads=atoi(optarg); break;
			case 'j': compression=AFILE_BZ2_COMPRESSION; break;
			case 'z': compression=AFILE_GZ_COMPRESSION; break;
			case 'A': space=atoi(optarg); break;
			case 's': start_read_num=atoi(optarg); break;
			case 'e': end_read_num=atoi(optarg); break;
			case 'M': max_hits=atoi(optarg); break;
			case 'Q': queue_length=atoi(optarg); break;
			case 'h': 
					  return bwtbfast2_usage(alg,
							  space,
							  seed_len,
							  max_mm,
							  num_threads,
							  max_hits,
							  queue_length);
			default: fprintf(stderr, "Unrecognized option: -%c\n", c); return 1;
		}
	}

	if(1 == argc || optind != argc) {
		return bwtbfast2_usage(alg,
				space,
				seed_len,
				max_mm,
				num_threads,
				max_hits,
				queue_length);
	}

	// check parameters
	if(NULL == ref_fn) {
		PrintError(fn_name, "-f", "Option is required", Exit, OutOfRange);
	}
	if(NULL == read_fn) {
		PrintError(fn_name, "-r", "Option is required", Exit, OutOfRange);
	}
	if(alg < 0 || 3 < alg) {
		PrintError(fn_name, "-a", "Option is out of range", Exit, OutOfRange);
	}
	if(seed_len <= 0) {
		PrintError(fn_name, "-l", "Option is out of range", Exit, OutOfRange);
	}
	if(max_mm < 0) {
		PrintError(fn_name, "-m", "Option is out of range", Exit, OutOfRange);
	}
	if(space < 0 || 1 < space) {
		PrintError(fn_name, "-A", "Option is out of range", Exit, OutOfRange);
	}
	if(start_read_num < 1 || end_read_num < start_read_num) {
		PrintError(fn_name, "-s & -e", "Options are out of range", Exit, OutOfRange);
	}
	if(num_threads < 0) {
		PrintError(fn_name, "-n", "Option is out of range", Exit, OutOfRange);
	}
	if(max_hits <= 0) {
		PrintError(fn_name, "-M", "Option is out of range", Exit, OutOfRange);
	}
	if(queue_length <= 0) {
		PrintError(fn_name, "-Q", "Option is out of range", Exit, OutOfRange);
	}

	// bfast
	bwtbfast2_core(ref_fn, read_fn, compression, alg, seed_len, max_mm, space, start_read_num, end_read_num, max_hits, num_threads, queue_length);

	// free file names
	free(ref_fn);
	free(read_fn);

	return 0;
}

void bwtbfast2_core(char *ref_fn, char *read_fn, int32_t compression, int32_t alg, int32_t seed_len, int32_t max_mm, int32_t space, int32_t start_read_num, int32_t end_read_num, int32_t max_hits, int32_t n_threads, int32_t queue_length)
{
	char *fn_name="bwtbfast2_core";
	int i;
	int n_matches, tot_matches= 0;
	int32_t n_matches_left=0;
	bfast_rg_match_t *matches=NULL;
	bwa_seqio_t *ks;
	clock_t t;
	bwt_t *bwt=NULL;
	bntseq_t *bns=NULL;
	gzFile fp_out=NULL;

	// initialization
	ks = bwa_seq_open(read_fn, compression);
	fp_out = gzdopen(fileno(stdout), "wb");
	if(NULL == fp_out) PrintError(fn_name, "stdout", "Coult not open for writing", Exit, OpenFileError);
	{ 
		// load BWT
		char str[1024];
		strcpy(str, ref_fn); strcat(str, ".bwt"); bwt = bwt_restore_bwt(str);
		// load BNS
		bns = bns_restore(ref_fn);
		// load SA
		strcpy(str, ref_fn); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
	}

	// Skip over start_num_reads
	if(1 < start_read_num) {
		fprintf(stderr, "[bwtbfast2_core] skipping over %d reads... ", start_read_num);
		t = clock();
		tot_matches = start_read_num;
		while (0 < tot_matches &&
				(matches = bfast_rg_match_read(ks, GETMIN(queue_length, tot_matches), &n_matches, space, 0)) != 0) {
			// Free
			for(i=0;i<n_matches;i++) {
				bfast_rg_match_t_destroy(&matches[i]);
			}
			free(matches); matches = NULL;

			tot_matches -= n_matches;
			n_matches_left -= n_matches;
		}
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();
	}
	n_matches_left = end_read_num - start_read_num + 1;

	tot_matches = 0;
	while (0 < n_matches_left &&
			(matches = bfast_rg_match_read(ks, GETMIN(n_matches_left, queue_length), &n_matches, space, 0)) != 0) {
		tot_matches += n_matches;
		t = clock();

		fprintf(stderr, "[bwtbfast2_core] matching... ");
			
#ifdef HAVE_LIBPTHREAD
		if(n_threads <= 1) {
			// no threads
			bwtbfast2_core_worker(0, bwt, bns, n_matches, matches, space, 1, alg, seed_len, max_mm, max_hits);
		}
		else {
			// threads
			pthread_t *tid;
			pthread_attr_t attr;
			bwtbfast2_thread_t *data;
			int j;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			data = (bwtbfast2_thread_t*)calloc(n_threads, sizeof(bwtbfast2_thread_t));
			tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
			for (j = 0; j < n_threads; ++j) {
				data[j].tid = j; data[j].bwt = bwt; data[j].bns = bns;
				data[j].n_matches = n_matches; data[j].matches = matches; 
				data[j].n_threads = n_threads;
				data[j].alg = alg; data[j].seed_len = seed_len; data[j].max_mm = max_mm;
				data[j].space = space;
				data[j].max_hits = max_hits;
				pthread_create(&tid[j], &attr, bwtbfast2_thread_worker, data + j);
			}
			for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
			free(data); free(tid);
		}
#else
		// no threads
		bwtbfast2_core_worker(0, bwt, bns, n_matches, matches, space, 1, alg, seed_len, max_mm, max_hits);
#endif

		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		// Print
		t = clock();
		fprintf(stderr, "[bwtbfast2_core] write to the disk... ");
		bfast_rg_match_t_print_queue(matches, n_matches, fp_out);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		// Free
		for(i=0;i<n_matches;i++) {
			bfast_rg_match_t_destroy(&matches[i]);
		}
		free(matches); matches = NULL;
			
		n_matches_left -= n_matches;

		fprintf(stderr, "[bwtbfast2_core] %d sequences have been processed.\n", tot_matches);

	}

	// destroy
	bwt_destroy(bwt);
	bns_destroy(bns);
	bwa_seq_close(ks);
	gzclose(fp_out);
}

#ifdef HAVE_LIBPTHREAD
void *bwtbfast2_thread_worker(void *data)
{
	bwtbfast2_thread_t *d = (bwtbfast2_thread_t*)data;
	bwtbfast2_core_worker(d->tid, d->bwt, d->bns, d->n_matches, d->matches, d->space, d->n_threads, d->alg, d->seed_len, d->max_mm, d->max_hits);
	return 0;
}
#endif

void bwtbfast2_core_worker(int tid, bwt_t *bwt, bntseq_t *bns, int n_matches, bfast_rg_match_t *matches, int32_t space, int n_threads, int32_t alg, int32_t seed_len, int32_t max_mm, int32_t max_hits) 
{
	//char *fn_name="bwtbfast2_core_worker";
	int32_t i;
	bfast2_stack_t *stack=NULL;

	stack = bfast2_init_stack(max_mm);

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
		bfast2_match(&matches[i],
				bwt,
				bns,
				space,
				alg,
				seed_len,
				max_mm,
				max_hits,
				stack);

		// remove duplicates
		bfast2_rg_match_t_merge_duplicates(&matches[i]);

		// free -> we don't need these anymore
		free(matches[i].read_int); matches[i].read_int = NULL;
		free(matches[i].read_rc_int); matches[i].read_rc_int = NULL;
	}
	bfast2_destroy_stack(stack);
}
