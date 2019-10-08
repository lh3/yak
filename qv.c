#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kthread.h"
#include "yak.h"
#include "kmer.h"
#include "bseq.h"
#include "sys.h"

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

typedef struct {
	int64_t c[1<<YAK_COUNTER_BITS];
	int32_t max;
	uint64_t *s;
} qv_cntbuf_t;

typedef struct {
	int k;
	const yak_qopt_t *opt;
	bseq_file_t *ks;
	const bfc_ch_t *ch;
	qv_cntbuf_t *buf;
} qv_shared_t;

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	qv_shared_t *qs;
} qv_step_t;

static int qv_count(qv_shared_t *qs, const bfc_kmer_t *x)
{
	int c;
	uint64_t y[2];
	bfc_kmer_hash(qs->k, x->x, y);
	c = bfc_ch_get(qs->ch, y);
	return c > 0? c : 0;
}

static void worker_qv(void *_data, long k, int tid)
{
	extern bfc_kmer_t bfc_kmer_null;
	qv_step_t *data = (qv_step_t*)_data;
	qv_shared_t *qs = data->qs;
	bseq1_t *s = &data->seqs[k];
	qv_cntbuf_t *b = &qs->buf[tid];
	int i, l, tot, non0;
	bfc_kmer_t x = bfc_kmer_null;

	if (s->l_seq < qs->opt->min_len) return;
	if (b->max < s->l_seq) {
		b->max = s->l_seq;
		kroundup32(b->max);
		b->s = (uint64_t*)realloc(b->s, b->max * sizeof(uint64_t));
	}
	for (i = l = 0, tot = non0 = 0; i < s->l_seq; ++i) {
		int c = seq_nt6_table[(uint8_t)s->seq[i]] - 1;
		if (c < 4) {
			bfc_kmer_append(qs->k, x.x, c);
			if (++l >= qs->k) {
				int t;
				t = qv_count(qs, &x);
				if (t > 0) ++non0;
				b->s[tot++] = (uint64_t)i<<32 | t;
			}
		} else l = 0, x = bfc_kmer_null;
	}

	if (qs->opt->print_each) {
		double qv = -1.0;
		if (tot > 0) {
			if (non0 > 0) {
				if (tot > non0) {
					qv = log((double)tot / non0) / qs->k;
					qv = 4.3429448190325175 * log(qv);
				} else qv = 99.0;
			} else qv = 0.0;
		}
		printf("S\t%s\t%d\t%d\t%d\t%.2f\n", s->name, s->l_seq, tot, non0, qv);
	}

	if (non0 < tot * qs->opt->min_frac) return;
	for (i = 0; i < tot; ++i)
		++b->c[(int32_t)b->s[i]];
}

static void *yak_qv_cb(void *shared, int step, void *_data)
{
	qv_shared_t *qs = (qv_shared_t*)shared;
	if (step == 0) {
		qv_step_t *ret;
		ret = calloc(1, sizeof(qv_step_t));
		ret->seqs = bseq_read(qs->ks, qs->opt->chunk_size, 0, &ret->n_seqs);
		ret->qs = qs;
		fprintf(stderr, "[M::%s] read %d sequences\n", __func__, ret->n_seqs);
		if (ret->seqs) return ret;
		else free(ret);
	} else if (step == 1) {
		int i;
		double rt, eff;
		qv_step_t *data = (qv_step_t*)_data;
		kt_for(qs->opt->n_threads, worker_qv, data, data->n_seqs);
		rt = yak_realtime();
		eff = yak_cputime() / (rt + 1e-6);
		fprintf(stderr, "[M::%s@%.2f*%.2f] processed %d sequences\n", __func__, rt, eff, data->n_seqs);
		for (i = 0; i < data->n_seqs; ++i) {
			bseq1_t *s = &data->seqs[i];
			free(s->seq); free(s->qual); free(s->comment); free(s->name);
		}
		free(data->seqs); free(data);
	}
	return 0;
}

void yak_qv(const yak_qopt_t *opt, const char *fn, const bfc_ch_t *ch, int64_t *cnt)
{
	qv_shared_t qs;
	int i, j, n_cnt = 1<<YAK_COUNTER_BITS;
	memset(&qs, 0, sizeof(qv_shared_t));
	qs.k = bfc_ch_get_k(ch);
	qs.opt = opt;
	qs.ch = ch;
	qs.buf = calloc(opt->n_threads, sizeof(qv_cntbuf_t));
	qs.ks = bseq_open(fn);
	kt_pipeline(2, yak_qv_cb, &qs, 2);
	bseq_close(qs.ks);
	memset(cnt, 0, n_cnt * sizeof(int64_t));
	for (i = 0; i < opt->n_threads; ++i) {
		for (j = 0; j < n_cnt; ++j)
			cnt[j] += qs.buf[i].c[j];
		free(qs.buf[i].s);
	}
	free(qs.buf);
}

void yak_qopt_init(yak_qopt_t *opt)
{
	memset(opt, 0, sizeof(yak_qopt_t));
	opt->chunk_size = 1000000000;
	opt->n_threads = 4;
	opt->min_frac = 0.5;
}
