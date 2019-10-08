#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kthread.h"
#include "yak.h"
#include "bbf.h"
#include "kmer.h"
#include "bseq.h"
#include "sys.h"

typedef struct {
	int64_t c[1<<YAK_COUNTER_BITS];
} qv_cntbuf_t;

typedef struct {
	int k, n_threads;
	int64_t chunk_size;
	bseq_file_t *ks;
	const bfc_ch_t *ch;
	qv_cntbuf_t *buf;
} qv_shared_t;

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	qv_shared_t *qs;
} qv_step_t;

static void qv_count(qv_shared_t *qs, const bfc_kmer_t *x, int tid)
{
	int c;
	uint64_t y[2], hash;
	hash = bfc_kmer_hash(qs->k, x->x, y);
	c = bfc_ch_get(qs->ch, y);
	if (c < 0) c = 0;
	++qs->buf[tid].c[c];
}

static void worker_qv(void *_data, long k, int tid)
{
	extern bfc_kmer_t bfc_kmer_null;
	qv_step_t *data = (qv_step_t*)_data;
	qv_shared_t *qs = data->qs;
	bseq1_t *s = &data->seqs[k];
	int i, l;
	bfc_kmer_t x = bfc_kmer_null;
	for (i = l = 0; i < s->l_seq; ++i) {
		int c = seq_nt6_table[(uint8_t)s->seq[i]] - 1;
		if (c < 4) {
			bfc_kmer_append(qs->k, x.x, c);
			if (++l >= qs->k) qv_count(qs, &x, tid);
		} else l = 0, x = bfc_kmer_null;
	}
}

static void *yak_qv_cb(void *shared, int step, void *_data)
{
	qv_shared_t *qs = (qv_shared_t*)shared;
	if (step == 0) {
		qv_step_t *ret;
		ret = calloc(1, sizeof(qv_step_t));
		ret->seqs = bseq_read(qs->ks, qs->chunk_size, 0, &ret->n_seqs);
		ret->qs = qs;
		fprintf(stderr, "[M::%s] read %d sequences\n", __func__, ret->n_seqs);
		if (ret->seqs) return ret;
		else free(ret);
	} else if (step == 1) {
		int i;
		double rt, eff;
		qv_step_t *data = (qv_step_t*)_data;
		kt_for(qs->n_threads, worker_qv, data, data->n_seqs);
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

void yak_qv(const char *fn, const bfc_ch_t *ch, int64_t chunk_size, int n_threads, int64_t *cnt)
{
	qv_shared_t qs;
	int i, j, n_cnt = 1<<YAK_COUNTER_BITS;
	memset(&qs, 0, sizeof(qv_shared_t));
	qs.k = bfc_ch_get_k(ch);
	qs.chunk_size = chunk_size;
	qs.n_threads = n_threads;
	qs.ch = ch;
	qs.buf = calloc(n_threads, sizeof(qv_cntbuf_t));
	qs.ks = bseq_open(fn);
	kt_pipeline(2, yak_qv_cb, &qs, 2);
	bseq_close(qs.ks);
	memset(cnt, 0, n_cnt * sizeof(int64_t));
	for (i = 0; i < n_threads; ++i)
		for (j = 0; j < n_cnt; ++j)
			cnt[j] += qs.buf[i].c[j];
	free(qs.buf);
}
