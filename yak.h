#ifndef YAK_H
#define YAK_H

#define YAK_VERSION "r30"

#include <stdint.h>

#define YAK_MAX_KMER     63
#define YAK_COUNTER_BITS 8
#define YAK_N_COUNTS     (1<<YAK_COUNTER_BITS)

typedef struct {
	int32_t bf_shift, bf_n_hashes;
	int32_t k;
	int32_t b_pre;
	int32_t n_threads;
	int64_t chunk_size;
} yak_copt_t;

typedef struct {
	int32_t print_each;
	int32_t min_len;
	int32_t n_threads;
	double min_frac;
	double fpr;
	int64_t chunk_size;
} yak_qopt_t;

typedef struct {
	int64_t tot;
	double qv_raw, qv, cov, err;
	double fpr_lower, fpr_upper;
	double adj_cnt[1<<YAK_COUNTER_BITS];
} yak_qstat_t;

extern int yak_verbose;

struct bfc_ch_s;
typedef struct bfc_ch_s bfc_ch_t;

void yak_copt_init(yak_copt_t *opt);
bfc_ch_t *bfc_count(const char *fn, const yak_copt_t *opt, bfc_ch_t *ch0);

void yak_qopt_init(yak_qopt_t *opt);
void yak_qv(const yak_qopt_t *opt, const char *fn, const bfc_ch_t *ch, int64_t *cnt);
int yak_qv_solve(const int64_t *hist, const int64_t *cnt, int kmer, double fpr, yak_qstat_t *qs);

bfc_ch_t *bfc_ch_init(int k, int l_pre);
void bfc_ch_destroy(bfc_ch_t *ch);
int bfc_ch_insert(bfc_ch_t *ch, const uint64_t x[2], int forced, int old_only);
int bfc_ch_get(const bfc_ch_t *ch, const uint64_t x[2]);
int bfc_ch_get_direct(const bfc_ch_t *ch, int pre, uint64_t key);
void bfc_ch_reset(bfc_ch_t *ch);
int64_t bfc_ch_del2(bfc_ch_t *ch);
uint64_t bfc_ch_count(const bfc_ch_t *ch);
int bfc_ch_hist(const bfc_ch_t *ch, int64_t cnt[1<<YAK_COUNTER_BITS]);
int bfc_ch_dump(const bfc_ch_t *ch, const char *fn);
bfc_ch_t *bfc_ch_restore(const char *fn);
int bfc_ch_get_k(const bfc_ch_t *ch);

int bfc_ch_kmer_occ(const bfc_ch_t *ch, const uint64_t z[2]);

#endif
