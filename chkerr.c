#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kthread.h"
#include "ketopt.h"
#include "bseq.h"
#include "yak-priv.h"

typedef struct {
	int k, n_threads, min_cnt;
	int64_t chunk_size;
	bseq_file_t *fp;
	const yak_ch_t *ch;
} ce_shared_t;

typedef struct {
	int n_seq;
	ce_shared_t *aux;
	bseq1_t *seq;
} ce_step_t;

static void te_worker(void *_data, long k, int tid)
{
	ce_step_t *t = (ce_step_t*)_data;
	ce_shared_t *aux = t->aux;
	bseq1_t *s = &t->seq[k];
	uint64_t x[4], mask;
	int i, l, last, streak, shift;
	if (aux->ch->k < 32) {
		mask = (1ULL<<2*aux->ch->k) - 1;
		shift = 2 * (aux->ch->k - 1);
	} else {
		mask = (1ULL<<aux->ch->k) - 1;
		shift = aux->ch->k - 1;
	}
	last = -1, streak = 0;
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < s->l_seq; ++i) {
		int cnt, c = seq_nt4_table[(uint8_t)s->seq[i]];
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
				cnt = yak_ch_get(aux->ch, y);
				if (cnt < aux->min_cnt) {
					if (i != last + 1) {
						if (streak > 0)
							printf("%s\t%d\t%d\t%d\n", s->name, i + 1 - aux->k - streak, i, streak);
						streak = 1;
					} else ++streak;
					last = i;
				}
			}
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0;
	}
}

static void *ce_pipeline(void *shared, int step, void *_data)
{
	ce_shared_t *aux = (ce_shared_t*)shared;
	if (step == 0) {
		ce_step_t *s;
		s = (ce_step_t*)calloc(1, sizeof(ce_step_t));
		s->seq = bseq_read(aux->fp, aux->chunk_size, 0, &s->n_seq);
		s->aux = aux;
		if (s->n_seq) {
			fprintf(stderr, "[M::%s] read %d sequences\n", __func__, s->n_seq);
			return s;
		} else free(s);
	} else if (step == 1) {
		int i;
		ce_step_t *s = (ce_step_t*)_data;
		kt_for(aux->n_threads, te_worker, s, s->n_seq);
		for (i = 0; i < s->n_seq; ++i) {
			free(s->seq[i].name); free(s->seq[i].seq); free(s->seq[i].qual);
			free(s->seq[i].comment);
		}
		free(s->seq); free(s);
	}
	return 0;
}

int main_chkerr(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c;
	yak_ch_t *ch;
	ce_shared_t aux;

	memset(&aux, 0, sizeof(ce_shared_t));
	aux.chunk_size = 1000000000;
	aux.n_threads = 8, aux.min_cnt = 3;
	while ((c = ketopt(&o, argc, argv, 1, "t:c:", 0)) >= 0) {
		if (c == 't') aux.n_threads = atoi(o.arg);
		else if (c == 'c') aux.min_cnt = atoi(o.arg);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: yak chkerr [options] <count.yak> <seq.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -t INT    number of threads [%d]\n", aux.n_threads);
		fprintf(stderr, "  -c INT    min k-mer count [%d]\n", aux.min_cnt);
		return 1;
	}

	ch = yak_ch_restore(argv[o.ind]);

	aux.k = ch->k;
	aux.fp = bseq_open(argv[o.ind+1]);
	if (aux.fp == 0) {
		fprintf(stderr, "ERROR: fail to open file '%s'\n", argv[o.ind+1]);
		exit(1);
	}
	printf("C\n");
	aux.ch = ch;
	kt_pipeline(2, ce_pipeline, &aux, 2);
	bseq_close(aux.fp);
	yak_ch_destroy(ch);
	return 0;
}
