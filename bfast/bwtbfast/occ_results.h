#ifndef OCC_RESULTS_H_
#define OCC_RESULTS_H_

#include <stdio.h>

typedef struct {
	int32_t *r_l;
	int32_t *r_u;
	int32_t n;
} occ_results_t;

void occ_results_t_init(occ_results_t *r);
void occ_results_t_add(occ_results_t *r, int32_t r_l, int32_t r_u);
void occ_results_t_destroy(occ_results_t *r);
void occ_results_t_print(occ_results_t *r, FILE *fp);

#endif
