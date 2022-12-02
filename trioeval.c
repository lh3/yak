#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kthread.h"
#include "ketopt.h"
#include "bseq.h"
#include "yak-priv.h"

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define CHUNK_SIZE 1000000000

typedef struct {
	int nk, c[4], d[2];
} te_cnt_t;

typedef struct {
	int max;
	uint32_t *s;
} te_buf_t;

typedef struct {
	int k, n_threads, min_n, print_err, print_frag;
	bseq_file_t *fp;
	const yak_ch_t *ch;
	te_buf_t *buf;
	int64_t n_pair, n_site, n_switch, n_err;
	int64_t n_par[2];
} te_shared_t;

typedef struct {
	int n_seq;
	te_shared_t *aux;
	bseq1_t *seq;
	te_cnt_t *cnt;
} te_step_t;

static void te_worker(void *_data, long k, int tid)
{
	te_step_t *t = (te_step_t*)_data;
	te_shared_t *aux = t->aux;
	bseq1_t *s = &t->seq[k];
	te_buf_t *b = &aux->buf[tid];
	uint64_t x[4], mask;
	int i, l, shift, last;
	int f_st, f_en, f_type, f_cnt;
	if (aux->ch->k < 32) {
		mask = (1ULL<<2*aux->ch->k) - 1;
		shift = 2 * (aux->ch->k - 1);
	} else {
		mask = (1ULL<<aux->ch->k) - 1;
		shift = aux->ch->k - 1;
	}
	if (s->l_seq > b->max) {
		b->max = s->l_seq;
		kroundup32(b->max);
		b->s = (uint32_t*)realloc(b->s, b->max * sizeof(uint32_t));
	}
	memset(b->s, 0, s->l_seq * sizeof(uint32_t));
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < s->l_seq; ++i) {
		int flag, c = seq_nt4_table[(uint8_t)s->seq[i]];
		if (c < 4) {
			if (aux->ch->k < 32) {
				x[0] = (x[0] << 2 | c) & mask;
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
			} else {
				x[0] = (x[0] << 1 | (c&1))  & mask;
				x[1] = (x[1] << 1 | (c>>1)) & mask;
				x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
				x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
			}
			if (++l >= aux->k) {
				int type = 0, c1, c2;
				uint64_t y;
				++t->cnt[k].nk;
				if (aux->ch->k < 32)
					y = yak_hash64(x[0] < x[1]? x[0] : x[1], mask);
				else
					y = yak_hash_long(x);
				flag = yak_ch_get(aux->ch, y);
				if (flag < 0) flag = 0;
				c1 = flag&3, c2 = flag>>2&3;
				if (c1 == 2 && c2 == 0) type = 1;
				else if (c2 == 2 && c1 == 0) type = 2;
				b->s[i] = type;
			}
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0;
	}
	f_type = f_st = f_en = f_cnt = 0;
	for (l = 0, i = 1, last = 0; i <= s->l_seq; ++i) {
		if (i == s->l_seq || b->s[i] != b->s[l]) { // found a streak
			if (b->s[l] > 0 && i - l >= aux->min_n) { // skip singletons
				int n = (i - l + aux->k - 1) / aux->k;
				int c = b->s[l] - 1;
				t->cnt[k].c[c << 1 | c] += n - 1;
				t->cnt[k].d[c] += n;
				if (last > 0) {
					++t->cnt[k].c[(last - 1) << 1 | c];
					if (aux->print_err && last - 1 != c)
						printf("E\t%s\t%d\t%d\t%d\n", s->name, i, last, c+1);
				}
				if (f_type != b->s[l]) {
					if (f_type > 0 && aux->print_frag)
						printf("F\t%s\t%d\t%d\t%d\t%d\n", s->name, f_type, f_st, f_en, f_cnt);
					f_type = b->s[l], f_st = l + 1 - aux->ch->k, f_cnt = 0;
				}
				++f_cnt, f_en = i + 1;
				last = b->s[l];
			}
			l = i;
		}
	}
	if (f_type > 0 && aux->print_frag)
		printf("F\t%s\t%d\t%d\t%d\t%d\n", s->name, f_type, f_st, f_en, f_cnt);
}

static void *te_pipeline(void *shared, int step, void *_data)
{
	te_shared_t *aux = (te_shared_t*)shared;
	if (step == 0) {
		te_step_t *s;
		s = (te_step_t*)calloc(1, sizeof(te_step_t));
		s->seq = bseq_read(aux->fp, CHUNK_SIZE, 0, &s->n_seq);
		s->aux = aux;
		if (s->n_seq) {
			s->cnt = (te_cnt_t*)calloc(s->n_seq, sizeof(te_cnt_t));
			fprintf(stderr, "[M::%s] read %d sequences\n", __func__, s->n_seq);
			return s;
		} else free(s);
	} else if (step == 1) {
		int i;
		te_step_t *s = (te_step_t*)_data;
		kt_for(aux->n_threads, te_worker, s, s->n_seq);
		for (i = 0; i < s->n_seq; ++i) {
			int *c = s->cnt[i].c, *d = s->cnt[i].d;
			aux->n_par[0] += d[0];
			aux->n_par[1] += d[1];
			if (d[0] + d[1] >= 2) {
				aux->n_pair += c[0] + c[1] + c[2] + c[3];
				aux->n_switch += c[1] + c[2];
				aux->n_site += d[0] + d[1];
				aux->n_err += d[0] < d[1]? d[0] : d[1];
			}
			printf("S\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", s->seq[i].name, d[0], d[1], c[0], c[1], c[2], c[3], s->seq[i].l_seq);
			free(s->seq[i].name); free(s->seq[i].seq); free(s->seq[i].qual); free(s->seq[i].comment);
		}
		free(s->seq); free(s->cnt); free(s);
	}
	return 0;
}

int main_trioeval(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int i, c, min_cnt = 2, mid_cnt = 5;
	int64_t cnt[YAK_N_COUNTS];
	yak_ch_t *ch;
	te_shared_t aux;

	memset(&aux, 0, sizeof(te_shared_t));
	aux.n_threads = 8, aux.min_n = 2, aux.print_frag = 1;
	while ((c = ketopt(&o, argc, argv, 1, "c:d:t:n:eF", 0)) >= 0) {
		if (c == 'c') min_cnt = atoi(o.arg);
		else if (c == 'd') mid_cnt = atoi(o.arg);
		else if (c == 't') aux.n_threads = atoi(o.arg);
		else if (c == 'n') aux.min_n = atoi(o.arg);
		else if (c == 'e') aux.print_err = 1;
		else if (c == 'F') aux.print_frag = 0;
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: yak trioeval [options] <pat.yak> <mat.yak> <seq.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -c INT     min occurrence [%d]\n", min_cnt);
		fprintf(stderr, "  -d INT     mid occurrence [%d]\n", mid_cnt);
		fprintf(stderr, "  -n INT     min streak [%d]\n", aux.min_n);
		fprintf(stderr, "  -t INT     number of threads [%d]\n", aux.n_threads);
		fprintf(stderr, "  -e         print error positions (out of order)\n");
		return 1;
	}

	ch = yak_ch_restore_core(0,  argv[o.ind],     YAK_LOAD_TRIOBIN1, min_cnt, mid_cnt);
	ch = yak_ch_restore_core(ch, argv[o.ind + 1], YAK_LOAD_TRIOBIN2, min_cnt, mid_cnt);
	yak_ch_hist(ch, cnt, aux.n_threads);
	fprintf(stderr, "[M::%s] %ld file1-specific k-mers and %ld file2-specific k-mers\n", __func__,
			(long)cnt[0<<2|2], (long)cnt[2<<2|0]);

	aux.k = ch->k;
	aux.fp = bseq_open(argv[o.ind+2]);
	if (aux.fp == 0) {
		fprintf(stderr, "ERROR: fail to open file '%s'\n", argv[o.ind+2]);
		exit(1);
	}
	printf("C\tS  seqName     #patKmer  #matKmer  #pat-pat  #pat-mat  #mat-pat  #mat-mat  seqLen\n");
	printf("C\tF  seqName     type      startPos  endPos    count\n");
	printf("C\tW  #switchErr  denominator  switchErrRate\n");
	printf("C\tH  #hammingErr denominator  hammingErrRate\n");
	printf("C\tN  #totPatKmer #totMatKmer  errRate\n");
	printf("C\n");
	aux.ch = ch;
	aux.buf = (te_buf_t*)calloc(aux.n_threads, sizeof(te_buf_t));
	kt_pipeline(2, te_pipeline, &aux, 2);
	bseq_close(aux.fp);
	yak_ch_destroy(ch);
	for (i = 0; i < aux.n_threads; ++i) free(aux.buf[i].s);
	printf("W\t%ld\t%ld\t%.6f\n", (long)aux.n_switch, (long)aux.n_pair, (double)aux.n_switch/aux.n_pair);
	printf("H\t%ld\t%ld\t%.6f\n", (long)aux.n_err, (long)aux.n_site, (double)aux.n_err/aux.n_site);
	printf("N\t%ld\t%ld\t%.6f\n", (long)aux.n_par[0], (long)aux.n_par[1], (double)(aux.n_par[0] < aux.n_par[1]? aux.n_par[0] : aux.n_par[1]) / (aux.n_par[0] + aux.n_par[1]));
	free(aux.buf);
	return 0;
}
