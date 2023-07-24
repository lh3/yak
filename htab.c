#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <math.h>
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

int yak_ch_get(const yak_ch_t *h, uint64_t x)
{
	int mask = (1<<h->pre) - 1;
	yak_ht_t *g = h->h[x&mask].h;
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
	a.h = h, a.min = min, a.max = max;
	kt_for(n_thread, worker_shrink, &a, 1<<h->pre);
	for (i = 0, h->tot = 0; i < 1<<h->pre; ++i)
		h->tot += kh_size(h->h[i].h);
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


//SET OPERATIONS

struct setops_dat {
	char **filenames;
	yak_ch_t **hts;
	volatile int *lhts;
	int n_files;
	int n_tables;
	uint64_t *ca, *cb, *cab, *aob;
	double *mashd;
};

void tri2rowcol(long i, size_t *row, size_t *col)
{
	*row = (size_t) floor(sqrt(2.0 * i + 0.25) - 0.5);  // Triangle index to 0-based row
	*col = i-((*row+1)*((*row))/2);
	*row +=1;
}

static void worker_setops(void *data, long i, int tid) // callback for kt_for()
{
	struct setops_dat *so = data;
	khint_t k;

	size_t row, col;
	tri2rowcol(i, &row, &col);

	yak_ch_t *A = yak_ch_restore(so->filenames[row]);
	yak_ch_t *B = yak_ch_restore(so->filenames[col]);
	//if (so->lhts[row] == 0) {
	//	so->lhts[row] = -1;
	//	printf("%i load %i\n", i, row);
	//	so->hts[row] = yak_ch_restore(so->filenames[row]);
	//	so->lhts[row] = 1;
	//}
	//while (so->lhts[row] == -1) {
	//	printf("%i wait row %i %i\n", i, row, so->lhts[row] );
	//}
	//yak_ch_t *A = so->hts[row];

	//if (so->lhts[col] == 0) {
	//	so->lhts[col] = -1;
	//	printf("%i load %i\n", i, col);
	//	so->hts[col] = yak_ch_restore(so->filenames[col]);
	//	so->lhts[col] = 1;
	//}
	//while (so->lhts[col] == -1) {
	//	printf("%i wait %i %i\n", i, col, so->lhts[col] );
	//}
	//yak_ch_t *B = so->hts[col];

	uint64_t ca=0, cb=0, cab=0, aob=0;
	for (int t = 0; t<so->n_tables; t++) {
		yak_ht_t *ah = A->h[t].h;
		yak_ht_t *bh = B->h[t].h;
		for (k = 0; k < kh_end(ah); ++k) {
			if (!kh_exist(ah, k)) continue;
			uint64_t x = kh_key(ah, k);
			khint_t bk = yak_ht_get(bh, x >> A->pre << YAK_COUNTER_BITS);
			int a = kh_exist(ah, k);
			int b = kh_exist(bh, bk);
			if (a && b) {
				cab++;
				ca++;
				cb++;
				aob++;
			} else if (a) {
				aob++;
				ca++;
			}
		}
		for (k = 0; k < kh_end(bh); ++k) {
			if (!kh_exist(bh, k)) continue;
			uint64_t x = kh_key(bh, k);
			khint_t ak = yak_ht_get(ah, x >> A->pre << YAK_COUNTER_BITS);
			int a = kh_exist(ah, ak);
			int b = kh_exist(bh, k);
			if (b && !a) {
				cb++;
				aob++;
			}
		}
	}
	so->ca[i] = ca;
	so->cb[i] = cb;
	so->cab[i] = cab;
	so->aob[i] = aob;
	double j = (double)cab/(double)aob;
	so->mashd[i] = cab == aob ? 0 : -(1.0/A->k)*log(2*j/(1+j));
	yak_ch_destroy(A);
	yak_ch_destroy(B);
	fprintf(stderr, "[M::%s::%.3f*%.2f] processed %s vs %s -> %lu\n", __func__,
			yak_realtime(), yak_cputime() / yak_realtime(), 
			so->filenames[row], so->filenames[col], aob);
}

#include "ketopt.h"
int main_setops(int argc, char *argv[])
{
	struct setops_dat so;
	ketopt_t o = KETOPT_INIT;
	int n_thread = 1;
	int sz = 1<<10;

	char c;
	while ((c = ketopt(&o, argc, argv, 1, "n:t:", 0)) >= 0) {
		if (c == 't') n_thread = atoi(o.arg);
		if (c == 'n') sz = atoi(o.arg);
	}

	int ni = argc - o.ind;
	if (ni < 2 || sz > 1024 || sz < 0) {
		fprintf(stderr, "USAGE: yak setops [options] <a.yak> <b.yak>...\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "OPTIONS:\n");
		fprintf(stderr, "\t-t INT  Use INT threads\n");
		fprintf(stderr, "\t-n INT  Use only first INT/1024 tables\n");
		return 1;
	}
	
	so.filenames = argv + o.ind;
	so.n_tables = sz;

	size_t n = (ni * (ni-1)) /2;
	so.ca = calloc(n, 8);
	so.cb = calloc(n, 8);
	so.cab = calloc(n, 8);
	so.aob = calloc(n, 8);
	//so.hts = calloc(ni, 8);
	//so.lhts = calloc(ni, 8);
	so.mashd = calloc(n, sizeof(double));

	kt_for(n_thread, worker_setops, &so, n);

	fprintf(stdout, "A\tB\tn_tbls\tkmers_A\tkmers_B\tkmers_AandB\tkmers_AorB\tmashd\n");
	for (size_t i = 0; i < n; i++) {
		size_t row, col;
		tri2rowcol(i, &row, &col);
		fprintf(stdout, "%s\t%s\t%d\t%lu\t%lu\t%lu\t%lu\t%lf\n",
				so.filenames[row], so.filenames[col], sz, so.ca[i], so.cb[i], so.cab[i], so.aob[i], so.mashd[i]);
	}
	free(so.ca);
	free(so.cb);
	free(so.cab);
	free(so.aob);
	free(so.mashd);
	return EXIT_SUCCESS;
}

