#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "util.h"
#include "BError.h"
#include "occ_results.h"

void occ_results_t_init(occ_results_t *r)
{
	r->r_l = r->r_u = NULL;
	r->n = 0;
}

void occ_results_t_add(occ_results_t *r, int32_t r_l, int32_t r_u)
{
	char *fn_name="occ_results_t_add";
	int32_t index_r_l=0, index_r_u=r->n-1;
	int32_t i, j, low, high;

	/*
	fprintf(stderr, "HERE PRINTING\n");
	for(i=0;i<r->n;i++) {
		fprintf(stderr, "%d %d\n", r->r_l[i], r->r_u[i]);
	}
	fprintf(stderr, "HERE PRINTING DONE\n");
	*/

	if(r_u < r_l) {
		// do nothing
	}
	else if(0 == r->n) {
		r->n++;
		r->r_l = my_realloc(r->r_l, sizeof(int32_t)*r->n, fn_name);
		r->r_u = my_realloc(r->r_u, sizeof(int32_t)*r->n, fn_name);

		r->r_l[r->n-1] = r_l;
		r->r_u[r->n-1] = r_u;
	}
	else if(r_u < r->r_l[0]) {
		if(r_u + 1 == r->r_l[0]) {
			r->r_l[0] = r_l;
		}
		else {
			// append to front
			r->n++;
			r->r_l = my_realloc(r->r_l, sizeof(int32_t)*r->n, fn_name);
			r->r_u = my_realloc(r->r_u, sizeof(int32_t)*r->n, fn_name);
			// shift over
			for(i=r->n-1;0<i;i--) {
				r->r_l[i] = r->r_l[i-1];
				r->r_u[i] = r->r_u[i-1];
			}
			// append to the front
			r->r_l[0] = r_l;
			r->r_u[0] = r_u;
		}
	}
	else if(r->r_u[r->n-1] < r_l) {
		if(r->r_u[r->n-1] + 1 == r_l) {
			r->r_u[r->n-1] = r_u;
		}
		else {
			// append to end
			r->n++;
			r->r_l = my_realloc(r->r_l, sizeof(int32_t)*r->n, fn_name);
			r->r_u = my_realloc(r->r_u, sizeof(int32_t)*r->n, fn_name);
			r->r_l[r->n-1] = r_l;
			r->r_u[r->n-1] = r_u;
		}
	}
	else {

		// get lower index
		low = 0;
		high = r->n-1;
		//fprintf(stderr, "starting low = %d high = %d\n", low, high);
		while(low <= high) {
			index_r_l = (low + high) >> 1;
			if(r_l-1 < r->r_l[index_r_l]) high = index_r_l-1;
			else if(r->r_l[index_r_l] < r_l-1) low = index_r_l+1;
			else break;
			//fprintf(stderr, "next index_r_l = %d low = %d high = %d\n", index_r_l, low, high);
		}
		while(0 <= index_r_l && r_l < r->r_l[index_r_l]) {
			index_r_l--;
		}
		//fprintf(stderr, "low = %d high = %d index_r_l = %d\n", low, high, index_r_l);
		// get upper index
		low = (0 <= index_r_l) ? ((index_r_l < r->n) ? index_r_l : r->n-1) : 0;
		high = r->n-1;
		while(low <= high) {
			index_r_u = (low + high) >> 1;
			if(r->r_u[index_r_u] < r_u + 1) low = index_r_u+1;
			else if(r_u + 1 < r->r_u[index_r_u]) high = index_r_u-1;
			else break;
		}
		while(index_r_u < r->n && r->r_u[index_r_u] < r_u) {
			index_r_u++;
		}
		//fprintf(stderr, "low = %d high = %d index_r_u = %d\n", low, high, index_r_u);

		assert(index_r_l < r->n);

		// check bounds
		if(index_r_l < 0) {
			assert(0 <= index_r_u && index_r_u < r->n);
			index_r_l = 0;
			r->r_l[0] = r_l;
		}

		if(r->r_u[index_r_l] < r_l - 1 &&
				index_r_u < r->n &&
				r_u + 1 < r->r_l[index_r_u]) {
			if(index_r_u == 1 + index_r_l) {
				//fprintf(stderr, "NEW RANGE\n");
				// New range
				r->n++;
				r->r_l = my_realloc(r->r_l, sizeof(int32_t)*r->n, fn_name);
				r->r_u = my_realloc(r->r_u, sizeof(int32_t)*r->n, fn_name);
				// shift up
				for(i=r->n-1;index_r_u<i;i--) {
					r->r_l[i] = r->r_l[i-1];
					r->r_u[i] = r->r_u[i-1];
				}
				r->r_l[index_r_u] = r_l;
				r->r_u[index_r_u] = r_u;
			}
			else {
				assert(index_r_l + 2 == index_r_u);
				assert(r_l == r->r_l[index_r_l+1]);
				assert(r_u == r->r_u[index_r_u-1]);
			}
		}
		else if(r_l - 1 <= r->r_u[index_r_l]  &&
				(r->n <= index_r_u || r_u + 1 < r->r_l[index_r_u])) {
			// r_l fits r_u does not
			//fprintf(stderr, "R_L FITS\n");
			r->r_u[index_r_l] = r_u;
		}
		else if(r->r_u[index_r_l] < r_l - 1 &&
				index_r_u < r->n &&
				r->r_l[index_r_u] <= r_u + 1) {
			// r_u fits r_l does not
			//fprintf(stderr, "R_U FITS\n");
			r->r_l[index_r_u] = r_l;
		}
		else {
			//fprintf(stderr, "BOTH FIT\n");
			// both fall within ranges
			r->r_u[index_r_l] = r->r_u[index_r_u]; 
			// merge together ranges
			for(i=index_r_l+1,j=index_r_u+1;j<r->n;j++) {
				r->r_l[i] = r->r_l[j];
				r->r_u[i] = r->r_u[j];
			}
			r->n -= (j-i);
			r->r_l = my_realloc(r->r_l, sizeof(int32_t)*r->n, fn_name);
			r->r_u = my_realloc(r->r_u, sizeof(int32_t)*r->n, fn_name);
		}
	}

	// HERE: CHECK THAT THE RANGES ARE OK
	/*
	int32_t fail=0;
	for(i=0;i<r->n;i++) {
		if(r->r_l[i] > r->r_u[i]) {
			fail=1; break;
		}
		if(i<r->n-1) {
			if(r->r_u[i] >= r->r_l[i+1]) {
				fail=1; break;
			}
		}
	}
	if(1 == fail) {
		fprintf(stderr, "HERE PRINTING\n");
		for(i=0;i<r->n;i++) {
			fprintf(stderr, "%d %d\n", r->r_l[i], r->r_u[i]);
		}
		fprintf(stderr, "HERE PRINTING DONE\n");
	}
	fprintf(stderr, "\n");
		*/
	for(i=0;i<r->n;i++) {
		assert(r->r_l[i] <= r->r_u[i]);
		if(i<r->n-1) assert(r->r_u[i] + 1 < r->r_l[i+1]);
	}
}

void occ_results_t_destroy(occ_results_t *r)
{
	free(r->r_l);
	free(r->r_u);
	occ_results_t_init(r);
}

void occ_results_t_print(occ_results_t *r, FILE *fp)
{
	int32_t i;
	fprintf(fp, "n: %d\n", r->n);
	for(i=0;i<r->n;i++) {
		fprintf(fp, "%d %d\n", r->r_l[i], r->r_u[i]);
	}
}
