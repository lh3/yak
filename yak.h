#ifndef YAK_H
#define YAK_H

#define YAKS_VERSION "0.1-r89-dirty"

#include <stdint.h>

#define YAK_MAX_KMER     31
#define YAK_COUNTER_BITS 10
#define YAK_N_COUNTS     (1<<YAK_COUNTER_BITS)
#define YAK_MAX_COUNT    ((1<<YAK_COUNTER_BITS)-1)

#define YAK_BLK_SHIFT  9 // 64 bytes, the size of a cache line
#define YAK_BLK_MASK   ((1<<(YAK_BLK_SHIFT)) - 1)

#define YAK_LOAD_ALL       1
#define YAK_LOAD_TRIOBIN1  2
#define YAK_LOAD_TRIOBIN2  3
#define YAK_LOAD_SEXCHR1   4
#define YAK_LOAD_SEXCHR2   5
#define YAK_LOAD_SEXCHR3   6

#define YAK_MAGIC "YAK\2"

typedef struct {
	int32_t bf_shift, bf_n_hash;
	int32_t k;
	int32_t pre;
	int32_t n_thread;
	int64_t chunk_size;
} yak_copt_t;

typedef struct {
	int32_t print_each, print_err_kmer;
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

typedef struct {
	int n_shift, n_hashes;
	uint8_t *b;
} yak_bf_t;

struct yak_ht_t;

typedef struct {
	struct yak_ht_t *h;
	yak_bf_t *b;
} yak_ch1_t;

typedef struct {
	int k, pre, n_hash, n_shift;
	uint64_t tot;
	yak_ch1_t *h;
} yak_ch_t;

typedef struct {
	uint64_t x;
	int c;
} yak_knt_t;

extern int yak_verbose;
extern unsigned char seq_nt4_table[256];

void yak_copt_init(yak_copt_t *opt);

yak_bf_t *yak_bf_init(int n_shift, int n_hashes);
void yak_bf_destroy(yak_bf_t *b);
int yak_bf_insert(yak_bf_t *b, uint64_t hash);

yak_ch_t *yak_ch_init(int k, int pre, int n_hash, int n_shift);
void yak_ch_destroy(yak_ch_t *h);
void yak_ch_destroy_bf(yak_ch_t *h);
int yak_ch_insert_list(yak_ch_t *h, int create_new, int n, const uint64_t *a);
int yak_ch_get(const yak_ch_t *h, uint64_t x);
int yak_ch_inc(yak_ch_t *h, uint64_t x);
yak_knt_t *yak_ch_getseq(const yak_ch_t *h, int w, uint32_t *n);

void yak_ch_tighten(yak_ch_t *h);
void yak_ch_clear(yak_ch_t *h, int n_thread);
void yak_ch_hist(const yak_ch_t *h, int64_t cnt[YAK_N_COUNTS], int n_thread);
void yak_ch_shrink(yak_ch_t *h, int min, int max, int n_thread);
void yak_ch_merge(yak_ch_t *h0, yak_ch_t *h1, int min, int max, int n_thread); // merge h1 into h0 and destroy h1 afterwards
void yak_ch_setcnt(yak_ch_t *h, int cnt, int n_thread);
void yak_ch_subtract(yak_ch_t *h0, const yak_ch_t *h1, int n_thread);

int yak_ch_dump(const yak_ch_t *h, const char *fn);
yak_ch_t *yak_ch_restore(const char *fn);
yak_ch_t *yak_ch_restore_core(yak_ch_t *ch0, const char *fn, int mode, ...);

yak_ch_t *yak_count(const char *fn, const yak_copt_t *opt, yak_ch_t *h0);
void yak_recount(const char *fn, yak_ch_t *h);

void yak_qopt_init(yak_qopt_t *opt);
void yak_qv(const yak_qopt_t *opt, const char *fn, const yak_ch_t *ch, int64_t *cnt);
int yak_qv_solve(const int64_t *hist, const int64_t *cnt, int kmer, double fpr, yak_qstat_t *qs);

#endif
