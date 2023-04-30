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

static long sc_chunk_size = 1000000000L;

typedef struct {
	int k, n_threads, hap;
	bseq_file_t *fp;
	const yak_ch_t *ch;
	uint32_t m_info, n_info;
} sc_shared_t;

typedef struct {
	int n_seq;
	sc_shared_t *aux;
	bseq1_t *seq;
} sc_step_t;

static void sc_worker(void *_data, long k, int tid)
{
	sc_step_t *t = (sc_step_t*)_data;
	sc_shared_t *aux = t->aux;
	bseq1_t *s = &t->seq[k];
	uint64_t x[4], mask;
	int i, l, shift, cnt[4];
	if (aux->ch->k < 32) {
		mask = (1ULL<<2*aux->ch->k) - 1;
		shift = 2 * (aux->ch->k - 1);
	} else {
		mask = (1ULL<<aux->ch->k) - 1;
		shift = aux->ch->k - 1;
	}
	cnt[0] = cnt[1] = cnt[2] = cnt[3] = 0;
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
				uint64_t y;
				if (aux->ch->k < 32)
					y = yak_hash64(x[0] < x[1]? x[0] : x[1], mask);
				else
					y = yak_hash_long(x);
				flag = yak_ch_get(aux->ch, y);
				if (flag < 0) flag = 0;
				++cnt[flag&3];
			}
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0;
	}
	printf("S\t%s\t%d\t%d\t0\t%d\t%d\t%d\t%d\n", s->name, s->l_seq, aux->hap, cnt[0], cnt[1], cnt[2], cnt[3]);
}

static void *sc_pipeline(void *shared, int step, void *_data)
{
	sc_shared_t *aux = (sc_shared_t*)shared;
	if (step == 0) {
		sc_step_t *s;
		s = (sc_step_t*)calloc(1, sizeof(sc_step_t));
		s->seq = bseq_read(aux->fp, sc_chunk_size, 0, &s->n_seq);
		s->aux = aux;
		if (s->n_seq) {
			fprintf(stderr, "[M::%s] read %d sequences\n", __func__, s->n_seq);
			return s;
		} else free(s);
	} else if (step == 1) {
		int i;
		sc_step_t *s = (sc_step_t*)_data;
		kt_for(aux->n_threads, sc_worker, s, s->n_seq);
		for (i = 0; i < s->n_seq; ++i) {
			free(s->seq[i].name); free(s->seq[i].seq); free(s->seq[i].qual); free(s->seq[i].comment);
		}
		free(s->seq); free(s);
	}
	return 0;
}

int main_sexchr(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int i, c;
	int64_t cnt[YAK_N_COUNTS];
	yak_ch_t *ch;
	sc_shared_t aux;

	memset(&aux, 0, sizeof(sc_shared_t));
	aux.n_threads = 8;
	while ((c = ketopt(&o, argc, argv, 1, "t:K:", 0)) >= 0) {
		if (c == 't') aux.n_threads = atoi(o.arg);
		else if (c == 'K') sc_chunk_size = mm_parse_num(o.arg);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: yak sexchr [options] <chrY.yak> <chrX.yak> <hap1.fa> <hap2.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -t INT     number of threads [%d]\n", aux.n_threads);
		fprintf(stderr, "  -K NUM     chunk size [1g]\n");
		return 1;
	}

	ch = yak_ch_restore_core(0,  argv[o.ind],     YAK_LOAD_SEXCHR1);
	ch = yak_ch_restore_core(ch, argv[o.ind + 1], YAK_LOAD_SEXCHR2);
	yak_ch_hist(ch, cnt, aux.n_threads);
	fprintf(stderr, "[M::%s] %ld file1-specific k-mers and %ld file2-specific k-mers\n", __func__,
			(long)cnt[0<<1|1], (long)cnt[1<<1|0]);

	printf("C\tS  seqName  seqLen  originalHap  0  #notsexchr  #sexchr1  #sexchr2  #sexchrboth\n");
	printf("C\n");

	aux.k = ch->k;
	aux.ch = ch;
	for (i = 1; i <= 2; ++i) {
		aux.hap = i;
		aux.fp = bseq_open(argv[o.ind+i+1]);
		if (aux.fp == 0) {
			fprintf(stderr, "ERROR: fail to open file '%s'\n", argv[o.ind+i+1]);
			exit(1);
		}
		kt_pipeline(2, sc_pipeline, &aux, 2);
		bseq_close(aux.fp);
	}
	yak_ch_destroy(ch);
	return 0;
}
