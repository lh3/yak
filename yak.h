#ifndef YAK_H
#define YAK_H

#define YAK_VERSION "r10"

#include <stdint.h>

#define YAK_MAX_KMER     63

typedef struct {
	int32_t bf_shift, bf_n_hashes;
	int32_t k;
	int32_t b_pre;
	int32_t n_threads;
	int64_t chunk_size;
} yak_opt_t;

#define YAK_COUNTER_BITS 8

extern int yak_verbose;

struct bfc_ch_s;
typedef struct bfc_ch_s bfc_ch_t;

bfc_ch_t *bfc_ch_init(int k, int l_pre);
void bfc_ch_destroy(bfc_ch_t *ch);
int bfc_ch_insert(bfc_ch_t *ch, const uint64_t x[2], int forced, int old_only);
int bfc_ch_get(const bfc_ch_t *ch, const uint64_t x[2]);
void bfc_ch_reset(bfc_ch_t *ch);
int64_t bfc_ch_del2(bfc_ch_t *ch);
uint64_t bfc_ch_count(const bfc_ch_t *ch);
//int bfc_ch_hist(const bfc_ch_t *ch, uint64_t cnt[256], uint64_t high[64]);
int bfc_ch_dump(const bfc_ch_t *ch, const char *fn);
bfc_ch_t *bfc_ch_restore(const char *fn);
int bfc_ch_get_k(const bfc_ch_t *ch);

int bfc_ch_kmer_occ(const bfc_ch_t *ch, const uint64_t z[2]);

bfc_ch_t *bfc_count(const char *fn, const yak_opt_t *opt, bfc_ch_t *ch0);

#endif
