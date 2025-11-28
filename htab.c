#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include "kthread.h"
#include "yak-priv.h"

#include "khashl.h" // hash table
#define yak_ch_eq(a, b) ((a)>>YAK_COUNTER_BITS == (b)>>YAK_COUNTER_BITS) // lower 8 bits for counts; higher bits for k-mer
#define yak_ch_hash(a) ((a)>>YAK_COUNTER_BITS)
KHASHL_SET_INIT(, yak_ht_t, yak_ht, uint64_t, yak_ch_hash, yak_ch_eq)

yak_ch_t *yak_ch_init(int k, int pre, int n_hash, int n_shift)
{
	yak_ch_t *h;
	int i;
	if (pre < YAK_COUNTER_BITS) return 0;
	CALLOC(h, 1);
	h->k = k, h->pre = pre;
	CALLOC(h->h, 1<<h->pre);
	for (i = 0; i < 1<<h->pre; ++i)
		h->h[i].h = yak_ht_init();
	if (n_hash > 0 && n_shift > h->pre) {
		h->n_hash = n_hash, h->n_shift = n_shift;
		for (i = 0; i < 1<<h->pre; ++i)
			h->h[i].b = yak_bf_init(h->n_shift - h->pre, h->n_hash);
	}
	return h;
}

void yak_ch_destroy_bf(yak_ch_t *h)
{
	int i;
	for (i = 0; i < 1<<h->pre; ++i) {
		if (h->h[i].b)
			yak_bf_destroy(h->h[i].b);
		h->h[i].b = 0;
	}
}

void yak_ch_destroy(yak_ch_t *h)
{
	int i;
	if (h == 0) return;
	yak_ch_destroy_bf(h);
	for (i = 0; i < 1<<h->pre; ++i)
		yak_ht_destroy(h->h[i].h);
	free(h->h); free(h);
}

int yak_ch_insert_list(yak_ch_t *h, int create_new, int n, const uint64_t *a)
{
	int j, mask = (1<<h->pre) - 1, n_ins = 0;
	yak_ch1_t *g;
	if (n == 0) return 0;
	g = &h->h[a[0]&mask];
	for (j = 0; j < n; ++j) {
		int ins = 1, absent;
		uint64_t x = a[j] >> h->pre;
		khint_t k;
		if ((a[j]&mask) != (a[0]&mask)) continue;
		if (create_new) {
			if (g->b)
				ins = (yak_bf_insert(g->b, x) == h->n_hash);
			if (ins) {
				k = yak_ht_put(g->h, x<<YAK_COUNTER_BITS, &absent);
				if (absent) ++n_ins;
				if ((kh_key(g->h, k)&YAK_MAX_COUNT) < YAK_MAX_COUNT)
					++kh_key(g->h, k);
			}
		} else {
			k = yak_ht_get(g->h, x<<YAK_COUNTER_BITS);
			if (k != kh_end(g->h) && (kh_key(g->h, k)&YAK_MAX_COUNT) < YAK_MAX_COUNT)
				++kh_key(g->h, k);
		}
	}
	return n_ins;
}

int yak_ch_inc(yak_ch_t *h, uint64_t x)
{
	int mask = (1<<h->pre) - 1;
	yak_ht_t *g = h->h[x&mask].h;
	khint_t k;
	k = yak_ht_get(g, x >> h->pre << YAK_COUNTER_BITS);
	if (k != kh_end(g)) {
		if ((kh_key(g, k)&YAK_MAX_COUNT) < YAK_MAX_COUNT)
			++kh_key(g, k);
		return kh_key(g, k)&YAK_MAX_COUNT;
	} else return -1;
}

int yak_ch_get(const yak_ch_t *h, uint64_t x)
{
	int mask = (1<<h->pre) - 1;
	const yak_ht_t *g = h->h[x&mask].h;
	khint_t k;
	k = yak_ht_get(g, x >> h->pre << YAK_COUNTER_BITS);
	return k == kh_end(g)? -1 : kh_key(g, k)&YAK_MAX_COUNT;
}

/*************************
 * Clear all counts to 0 *
 *************************/

static void worker_clear(void *data, long i, int tid) // callback for kt_for()
{
	yak_ch_t *h = (yak_ch_t*)data;
	yak_ht_t *g = h->h[i].h;
	khint_t k;
	uint64_t mask = ~1ULL >> YAK_COUNTER_BITS << YAK_COUNTER_BITS;
	for (k = 0; k < kh_end(g); ++k)
		if (kh_exist(g, k))
			kh_key(g, k) &= mask;
}

void yak_ch_clear(yak_ch_t *h, int n_thread)
{
	kt_for(n_thread, worker_clear, h, 1<<h->pre);
}

/*************
 * Histogram *
 *************/

typedef struct {
	uint64_t c[YAK_N_COUNTS];
} buf_cnt_t;

typedef struct {
	const yak_ch_t *h;
	buf_cnt_t *cnt;
} hist_aux_t;

static void worker_hist(void *data, long i, int tid) // callback for kt_for()
{
	hist_aux_t *a = (hist_aux_t*)data;
	uint64_t *cnt = a->cnt[tid].c;
	yak_ht_t *g = a->h->h[i].h;
	khint_t k;
	for (k = 0; k < kh_end(g); ++k)
		if (kh_exist(g, k))
			++cnt[kh_key(g, k)&YAK_MAX_COUNT];
}

void yak_ch_hist(const yak_ch_t *h, int64_t cnt[YAK_N_COUNTS], int n_thread)
{
	hist_aux_t a;
	int i, j;
	a.h = h;
	memset(cnt, 0, YAK_N_COUNTS * sizeof(uint64_t));
	CALLOC(a.cnt, n_thread);
	kt_for(n_thread, worker_hist, &a, 1<<h->pre);
	for (i = 0; i < YAK_N_COUNTS; ++i) cnt[i] = 0;
	for (j = 0; j < n_thread; ++j)
		for (i = 0; i < YAK_N_COUNTS; ++i)
			cnt[i] += a.cnt[j].c[i];
	free(a.cnt);
}

/**********
 * Shrink *
 **********/

typedef struct {
	int min, max;
	yak_ch_t *h;
} shrink_aux_t;

static void worker_shrink(void *data, long i, int tid) // callback for kt_for()
{
	shrink_aux_t *a = (shrink_aux_t*)data;
	yak_ch_t *h = a->h;
	yak_ht_t *g = h->h[i].h, *f;
	khint_t k;
	f = yak_ht_init();
	yak_ht_resize(f, kh_size(g));
	for (k = 0; k < kh_end(g); ++k) {
		if (kh_exist(g, k)) {
			int absent, c = kh_key(g, k) & YAK_MAX_COUNT;
			if (c >= a->min && c <= a->max)
				yak_ht_put(f, kh_key(g, k), &absent);
		}
	}
	yak_ht_destroy(g);
	h->h[i].h = f;
}

void yak_ch_shrink(yak_ch_t *h, int min, int max, int n_thread)
{
	int i;
	shrink_aux_t a;
	a.h = h, a.min = min;
	a.max = max >= min && max <= YAK_MAX_COUNT? max : YAK_MAX_COUNT;
	kt_for(n_thread, worker_shrink, &a, 1<<h->pre);
	for (i = 0, h->tot = 0; i < 1<<h->pre; ++i)
		h->tot += kh_size(h->h[i].h);
}

/*********
 * merge *
 *********/

typedef struct {
	int min, max;
	yak_ch_mtype_t type;
	yak_ch_t *h0, *h1;
} merge_aux_t;

static void worker_merge(void *data, long i, int tid)
{
	khint_t k, l;
	merge_aux_t *a = (merge_aux_t*)data;
	yak_ch_t *h0 = a->h0, *h1 = a->h1;
	yak_ht_t *g0 = h0->h[i].h, *g1 = h1->h[i].h, *f;
	f = yak_ht_init();
	yak_ht_resize(f, kh_size(g0));
	for (k = 0; k < kh_end(g0); ++k) {
		if (kh_exist(g0, k)) {
			int absent, c;
			l = yak_ht_get(g1, kh_key(g0, k));
			c = l != kh_end(g0)? (kh_key(g1, l) & YAK_MAX_COUNT) : 0;
			if (c > a->max) continue;
			if (a->type == YAK_MT_ISEC_MAX || c >= a->min)
				yak_ht_put(f, kh_key(g0, k), &absent);
		}
	}
	yak_ht_destroy(g0);
	yak_ht_destroy(g1);
	if (h1->h[i].b) yak_bf_destroy(h1->h[i].b);
	h0->h[i].h = f;
}

static void worker_merge_add(void *data, long i, int tid)
{
	khint_t k;
	merge_aux_t *a = (merge_aux_t*)data;
	yak_ch_t *h0 = a->h0, *h1 = a->h1;
	yak_ht_t *g0 = h0->h[i].h, *g1 = h1->h[i].h;
	for (k = 0; k < kh_end(g1); ++k) {
		if (kh_exist(g1, k)) {
			int absent, c;
			c = (kh_key(g1, k) & YAK_MAX_COUNT);
			if (c >= a->min && c <= a->max)
				yak_ht_put(g0, kh_key(g1, k), &absent);
		}
	}
	yak_ht_destroy(g1);
	if (h1->h[i].b) yak_bf_destroy(h1->h[i].b);
}

void yak_ch_merge(yak_ch_t *h0, yak_ch_t *h1, int min, int max, yak_ch_mtype_t type, int n_thread) // h1 merged into h0; h1 is destroyed afterwards
{
	merge_aux_t a;
	int i;
	a.h0 = h0, a.h1 = h1, a.type = type, a.min = min;
	a.max = type == YAK_MT_ISEC_MAX || (max >= min && max <= YAK_MAX_COUNT)? max : YAK_MAX_COUNT;
	if (type == YAK_MT_ADD)
		kt_for(n_thread, worker_merge_add, &a, 1<<h0->pre);
	else if (type == YAK_MT_ISEC_RANGE || type == YAK_MT_ISEC_MAX)
		kt_for(n_thread, worker_merge, &a, 1<<h0->pre);
	free(h1->h); free(h1);
	for (i = 0, h0->tot = 0; i < 1<<h0->pre; ++i)
		h0->tot += kh_size(h0->h[i].h);
}

/**********
 * getseq *
 **********/

yak_knt_t *yak_ch_getseq(const yak_ch_t *h, int w, uint32_t *n)
{
	yak_knt_t *a;
	uint64_t i, mask = (1ULL<<h->k*2) - 1;
	khint_t k;
	yak_ht_t *g;
	assert(h->k < 32 && w < 1<<h->pre);
	g = h->h[w].h;
	*n = kh_size(g);
	CALLOC(a, *n);
	for (i = 0, k = 0; k < kh_end(g); ++k)
		if (kh_exist(g, k))
			a[i].x = yak_hash64_inv(kh_key(g, k) >> YAK_COUNTER_BITS << h->pre | w, mask), a[i++].c = kh_key(g, k)&YAK_MAX_COUNT;
	return a;
}

/*******
 * I/O *
 *******/

int yak_ch_dump(const yak_ch_t *ch, const char *fn)
{
	FILE *fp;
	uint32_t t[3];
	int i;
	if ((fp = strcmp(fn, "-")? fopen(fn, "wb") : stdout) == 0) return -1;
	fwrite(YAK_MAGIC, 1, 4, fp);
	t[0] = ch->k, t[1] = ch->pre, t[2] = YAK_COUNTER_BITS;
	fwrite(t, 4, 3, fp);
	for (i = 0; i < 1<<ch->pre; ++i) {
		yak_ht_t *h = ch->h[i].h;
		khint_t k;
		t[0] = kh_capacity(h), t[1] = kh_size(h);
		fwrite(t, 4, 2, fp);
		for (k = 0; k < kh_end(h); ++k)
			if (kh_exist(h, k))
				fwrite(&kh_key(h, k), 8, 1, fp);
	}
	fprintf(stderr, "[M::%s] dumpped the hash table to file '%s'.\n", __func__, fn);
	fclose(fp);
	return 0;
}

yak_ch_t *yak_ch_restore_core(yak_ch_t *ch0, const char *fn, int mode, ...)
{
	va_list ap;
	FILE *fp;
	uint32_t t[3];
	char magic[4];
	int i, j, absent, min_cnt = 0, mid_cnt = 0, mode_err = 0;
	uint64_t mask = (1ULL<<YAK_COUNTER_BITS) - 1, n_ins = 0, n_new = 0;
	yak_ch_t *ch;

	va_start(ap, mode);
	if (mode == YAK_LOAD_ALL) { // do nothing
	} else if (mode == YAK_LOAD_TRIOBIN1 || mode == YAK_LOAD_TRIOBIN2) {
		assert(YAK_COUNTER_BITS >= 4);
		min_cnt = va_arg(ap, int);
		mid_cnt = va_arg(ap, int);
		if (ch0 == 0 && mode == YAK_LOAD_TRIOBIN2)
			mode_err = 1;
	} else if (mode == YAK_LOAD_SEXCHR1 || mode == YAK_LOAD_SEXCHR2 || mode == YAK_LOAD_SEXCHR3) {
		assert(YAK_COUNTER_BITS >= 3);
		if (ch0 == 0 && (mode == YAK_LOAD_SEXCHR2 || mode == YAK_LOAD_SEXCHR3))
			mode_err = 1;
	} else mode_err = 1;
	va_end(ap);
	if (mode_err) return 0;

	if ((fp = fopen(fn, "rb")) == 0) return 0;
	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, YAK_MAGIC, 4) != 0) {
		fprintf(stderr, "ERROR: wrong file magic.\n");
		fclose(fp);
		return 0;
	}
	fread(t, 4, 3, fp);
	if (t[2] != YAK_COUNTER_BITS) {
		fprintf(stderr, "ERROR: saved counter bits: %d; compile-time counter bits: %d\n", t[2], YAK_COUNTER_BITS);
		fclose(fp);
		return 0;
	}

	ch = ch0 == 0? yak_ch_init(t[0], t[1], 0, 0) : ch0;
	assert((int)t[0] == ch->k && (int)t[1] == ch->pre);
	for (i = 0; i < 1<<ch->pre; ++i) {
		yak_ht_t *h = ch->h[i].h;
		fread(t, 4, 2, fp);
		yak_ht_resize(h, t[0]);
		for (j = 0; j < t[1]; ++j) {
			uint64_t key;
			fread(&key, 8, 1, fp);
			if (mode == YAK_LOAD_ALL) {
				++n_ins;
				yak_ht_put(h, key, &absent);
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
					k = yak_ht_put(h, key, &absent);
					if (absent) ++n_new;
					else kh_key(h, k) = kh_key(h, k) | x;
				}
			} else if (mode == YAK_LOAD_SEXCHR1 || mode == YAK_LOAD_SEXCHR2 || mode == YAK_LOAD_SEXCHR3) {
				int shift = mode == YAK_LOAD_SEXCHR1? 0 : mode == YAK_LOAD_SEXCHR2? 1 : 2, x = 1<<shift;
				khint_t k;
				key = (key & ~mask) | x;
				++n_ins;
				k = yak_ht_put(h, key, &absent);
				if (absent) ++n_new;
				else kh_key(h, k) = kh_key(h, k) | x;
			}
		}
	}
	fclose(fp);
	fprintf(stderr, "[M::%s] inserted %ld k-mers, of which %ld are new\n", __func__, (long)n_ins, (long)n_new);
	return ch;
}

yak_ch_t *yak_ch_restore(const char *fn)
{
	return yak_ch_restore_core(0, fn, YAK_LOAD_ALL);
}
