#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include "yak.h"
#include "kmer.h"
#include "khash.h"

#define YAK_MAX_COUNT ((1<<YAK_COUNTER_BITS) - 1)

#define _cnt_eq(a, b) ((a)>>YAK_COUNTER_BITS == (b)>>YAK_COUNTER_BITS)
#define _cnt_hash(a) ((a)>>YAK_COUNTER_BITS)
KHASH_INIT(cnt, uint64_t, char, 0, _cnt_hash, _cnt_eq)
typedef khash_t(cnt) cnthash_t;

struct bfc_ch_s {
	cnthash_t **h;
	int k;
	// private
	int b_pre;
};

bfc_ch_t *bfc_ch_init(int k, int b_pre)
{
	bfc_ch_t *ch;
	int i;
	if (k > 63) return 0;
	if (2 * k - b_pre < 16) return 0;
	if (64 + b_pre < 2 * k + YAK_COUNTER_BITS) return 0;
	ch = calloc(1, sizeof(bfc_ch_t));
	ch->k = k, ch->b_pre = b_pre;
	ch->h = calloc(1<<ch->b_pre, sizeof(void*));
	for (i = 0; i < 1<<ch->b_pre; ++i)
		ch->h[i] = kh_init(cnt);
	return ch;
}

void bfc_ch_destroy(bfc_ch_t *ch)
{
	int i;
	if (ch == 0) return;
	for (i = 0; i < 1<<ch->b_pre; ++i)
		kh_destroy(cnt, ch->h[i]);
	free(ch->h); free(ch);
}

static inline cnthash_t *get_subhash(const bfc_ch_t *ch, const uint64_t x[2], uint64_t *key)
{
	if (ch->k <= 32) {
		int t = ch->k * 2 - ch->b_pre;
		uint64_t z = x[0] << ch->k | x[1];
		*key = (z & ((1ULL<<t) - 1)) << YAK_COUNTER_BITS;
		return ch->h[z>>t];
	} else {
		int max_key_bits = 64 - YAK_COUNTER_BITS;
		int t = ch->k - ch->b_pre;
		int shift = t + ch->k < max_key_bits? ch->k : max_key_bits - t;
		*key = ((x[0] & ((1ULL<<t) - 1)) << shift ^ x[1]) << YAK_COUNTER_BITS;
		return ch->h[x[0]>>t];
	}
}

int bfc_ch_insert(bfc_ch_t *ch, const uint64_t x[2], int forced, int old_only)
{
	int absent;
	uint64_t key;
	cnthash_t *h;
	khint_t k;
	h = get_subhash(ch, x, &key);
	if (__sync_lock_test_and_set(&h->lock, 1)) {
		if (forced) // then wait until the hash table is unlocked by the thread using it
			while (__sync_lock_test_and_set(&h->lock, 1))
				while (h->lock); // lock
		else return -1;
	}
	if (old_only) {
		k = kh_get(cnt, h, key);
		if (k != kh_end(h) && (kh_key(h, k) & YAK_MAX_COUNT) != YAK_MAX_COUNT)
			++kh_key(h, k);
	} else {
		k = kh_put(cnt, h, key, &absent);
		if (absent) {
			kh_key(h, k) |= 1;
		} else {
			if ((kh_key(h, k) & YAK_MAX_COUNT) != YAK_MAX_COUNT) ++kh_key(h, k);
		}
	}
	__sync_lock_release(&h->lock); // unlock
	return 0;
}

int bfc_ch_get(const bfc_ch_t *ch, const uint64_t x[2])
{
	uint64_t key;
	cnthash_t *h;
	khint_t itr;
	h = get_subhash(ch, x, &key);
	itr = kh_get(cnt, h, key);
	return itr == kh_end(h)? -1 : kh_key(h, itr) & YAK_MAX_COUNT;
}

int bfc_ch_get_direct(const bfc_ch_t *ch, int pre, uint64_t key)
{
	cnthash_t *h = ch->h[pre];
	khint_t k;
	k = kh_get(cnt, h, key);
	return k == kh_end(h)? -1 : kh_key(h, k) & YAK_MAX_COUNT;
}

int bfc_ch_kmer_occ(const bfc_ch_t *ch, const uint64_t z[2])
{
	uint64_t x[2];
	bfc_kmer_hash(ch->k, z, x);
	return bfc_ch_get(ch, x);
}

uint64_t bfc_ch_count(const bfc_ch_t *ch)
{
	int i;
	uint64_t cnt = 0;
	for (i = 0; i < 1<<ch->b_pre; ++i)
		cnt += kh_size(ch->h[i]);
	return cnt;
}

void bfc_ch_reset(bfc_ch_t *ch)
{
	int i;
	uint64_t mask = (~0ULL) << YAK_COUNTER_BITS;
	for (i = 0; i < 1<<ch->b_pre; ++i) {
		cnthash_t *h = ch->h[i];
		khint_t k;
		for (k = 0; k != kh_end(h); ++k)
			if (kh_exist(h, k))
				kh_key(h, k) &= mask;
	}
}

int64_t bfc_ch_del2(bfc_ch_t *ch)
{
	int i;
	int64_t del = 0, cnt = 0;
	for (i = 0; i < 1<<ch->b_pre; ++i) {
		cnthash_t *h = ch->h[i];
		khint_t k;
		cnt += kh_size(h);
		for (k = 0; k != kh_end(h); ++k)
			if (kh_exist(h, k) && (kh_key(h, k) & YAK_MAX_COUNT) < 2) {
				kh_del(cnt, h, k);
				++del;
			}
	}
	if (yak_verbose >= 3)
		fprintf(stderr, "[M] false positive rate: %ld / %ld = %.3f%%\n",
				(long)del, (long)cnt, 100.0 * del / cnt);
	return del;
}

int bfc_ch_hist(const bfc_ch_t *ch, int64_t cnt[1<<YAK_COUNTER_BITS])
{
	int i, max_i = -1;
	uint64_t max;
	memset(cnt, 0, (1<<YAK_COUNTER_BITS) * 8);
	for (i = 0; i < 1<<ch->b_pre; ++i) {
		khint_t k;
		cnthash_t *h = ch->h[i];
		for (k = 0; k != kh_end(h); ++k)
			if (kh_exist(h, k))
				++cnt[kh_key(h, k) & YAK_MAX_COUNT];
	}
	for (i = 3, max = 0; i <= YAK_MAX_COUNT; ++i)
		if (cnt[i] > max)
			max = cnt[i], max_i = i;
	return max_i;
}

int bfc_ch_dump(const bfc_ch_t *ch, const char *fn)
{
	FILE *fp;
	uint32_t t[3];
	int i;
	if ((fp = strcmp(fn, "-")? fopen(fn, "wb") : stdout) == 0) return -1;
	t[0] = ch->k, t[1] = ch->b_pre, t[2] = YAK_COUNTER_BITS;
	fwrite(t, 4, 3, fp);
	for (i = 0; i < 1<<ch->b_pre; ++i) {
		cnthash_t *h = ch->h[i];
		khint_t k;
		t[0] = kh_n_buckets(h), t[1] = kh_size(h);
		fwrite(t, 4, 2, fp);
		for (k = 0; k < kh_end(h); ++k)
			if (kh_exist(h, k))
				fwrite(&kh_key(h, k), 8, 1, fp);
	}
	fprintf(stderr, "[M::%s] dumpped the hash table to file '%s'.\n", __func__, fn);
	fclose(fp);
	return 0;
}

bfc_ch_t *bfc_ch_restore(const char *fn)
{
	FILE *fp;
	uint32_t t[3];
	int i, j, absent;
	bfc_ch_t *ch;

	if ((fp = fopen(fn, "rb")) == 0) return 0;
	fread(t, 4, 3, fp);
	if (t[2] != YAK_COUNTER_BITS) {
		fprintf(stderr, "ERROR: saved counter bits: %d; compile-time counter bits: %d\n", t[2], YAK_COUNTER_BITS);
		fclose(fp);
		return 0;
	}
	ch = bfc_ch_init(t[0], t[1]);
	assert((int)t[1] == ch->b_pre);
	for (i = 0; i < 1<<ch->b_pre; ++i) {
		cnthash_t *h = ch->h[i];
		fread(t, 4, 2, fp);
		kh_resize(cnt, h, t[0]);
		for (j = 0; j < t[1]; ++j) {
			uint64_t key;
			fread(&key, 8, 1, fp);
			kh_put(cnt, h, key, &absent);
			assert(absent);
		}
	}
	fclose(fp);
	fprintf(stderr, "[M::%s] restored the hash table from file '%s'.\n", __func__, fn);
	return ch;
}

int bfc_ch_get_k(const bfc_ch_t *ch)
{
	return ch->k;
}

/**************************
 * Specialized operations *
 **************************/

bfc_ch_t *bfc_ch_restore_core(bfc_ch_t *ch0, const char *fn, int mode, ...)
{
	va_list ap;
	FILE *fp;
	uint32_t t[3];
	int i, j, absent, min_cnt = 0, mid_cnt = 0, mode_err = 0;
	uint64_t mask = (1ULL<<YAK_COUNTER_BITS) - 1, n_ins = 0, n_new = 0;
	bfc_ch_t *ch;

	va_start(ap, mode);
	if (mode == YAK_LOAD_ALL) { // do nothing
	} else if (mode == YAK_LOAD_TRIOBIN1 || mode == YAK_LOAD_TRIOBIN2) {
		assert(YAK_COUNTER_BITS >= 4);
		min_cnt = va_arg(ap, int);
		mid_cnt = va_arg(ap, int);
		if (ch0 == 0 && mode == YAK_LOAD_TRIOBIN2)
			mode_err = 1;
	} else mode_err = 1;
	va_end(ap);
	if (mode_err) return 0;

	if ((fp = fopen(fn, "rb")) == 0) return 0;
	fread(t, 4, 3, fp);
	if (t[2] != YAK_COUNTER_BITS) {
		fprintf(stderr, "ERROR: saved counter bits: %d; compile-time counter bits: %d\n", t[2], YAK_COUNTER_BITS);
		fclose(fp);
		return 0;
	}

	ch = ch0 == 0? bfc_ch_init(t[0], t[1]) : ch0;
	assert((int)t[0] == ch->k && (int)t[1] == ch->b_pre);
	for (i = 0; i < 1<<ch->b_pre; ++i) {
		cnthash_t *h = ch->h[i];
		fread(t, 4, 2, fp);
		if (ch0 == 0) kh_resize(cnt, h, t[0]);
		for (j = 0; j < t[1]; ++j) {
			uint64_t key;
			fread(&key, 8, 1, fp);
			if (mode == YAK_LOAD_ALL) {
				++n_ins;
				kh_put(cnt, h, key, &absent);
				if (absent) ++n_new;
			} else if (mode == YAK_LOAD_TRIOBIN1 || mode == YAK_LOAD_TRIOBIN2) {
				int cnt = key & mask, x, shift = mode == YAK_LOAD_TRIOBIN1? 0 : 2;
				if (cnt >= mid_cnt) x = 2<<shift;
				else if (cnt >= min_cnt) x = 1<<shift;
				else x = -1;
				if (x >= 0) {
					khint_t k;
					key = (key & ~mask) | x;
					++n_ins;
					k = kh_put(cnt, h, key, &absent);
					if (absent) ++n_new;
					else kh_key(h, k) = kh_key(h, k) | x;
				}
			}
		}
	}
	fclose(fp);
	fprintf(stderr, "[M::%s] inserted %ld k-mers, of which %ld are new\n", __func__, (long)n_ins, (long)n_new);
	return ch;
}
