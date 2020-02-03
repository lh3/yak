#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kthread.h"
#include "ketopt.h"
#include "kmer.h"
#include "bseq.h"
#include "yak.h"

#define CHUNK_SIZE 200000000

typedef struct {
	int c[4];
} tb_cnt_t;

typedef struct {
	int k, n_threads, print_diff;
	bseq_file_t *fp;
	const bfc_ch_t *ch;
} tb_shared_t;

typedef struct {
	int n_seq;
	tb_shared_t *aux;
	bseq1_t *seq;
	tb_cnt_t *cnt;
} tb_step_t;

static inline int tb_count(tb_shared_t *aux, const bfc_kmer_t *x)
{
	int c;
	uint64_t y[2];
	bfc_kmer_hash(aux->k, x->x, y);
	c = bfc_ch_get(aux->ch, y);
	return c < 0? 0 : c;
}

static void tb_worker(void *_data, long k, int tid)
{
	extern bfc_kmer_t bfc_kmer_null;
	tb_step_t *t = (tb_step_t*)_data;
	tb_shared_t *aux = t->aux;
	bseq1_t *s = &t->seq[k];
	bfc_kmer_t x = bfc_kmer_null;
	int i, l;
	for (i = l = 0; i < s->l_seq; ++i) {
		int flag, c = seq_nt6_table[(uint8_t)s->seq[i]] - 1;
		if (c < 4) {
			bfc_kmer_append(aux->k, x.x, c);
			if (++l >= aux->k) {
				flag = tb_count(aux, &x);
				++t->cnt[k].c[flag];
				if (aux->print_diff && (flag == 1 || flag == 2))
					printf("D\t%s\t%d\t%d\n", s->name, i, flag);
			}
		} else l = 0, x = bfc_kmer_null;
	}
}

static void *tb_pipeline(void *shared, int step, void *_data)
{
	tb_shared_t *aux = (tb_shared_t*)shared;
	if (step == 0) {
		tb_step_t *s;
		s = (tb_step_t*)calloc(1, sizeof(tb_step_t));
		s->seq = bseq_read(aux->fp, CHUNK_SIZE, 0, &s->n_seq);
		s->aux = aux;
		if (s->n_seq) {
			s->cnt = (tb_cnt_t*)calloc(s->n_seq, sizeof(tb_cnt_t));
			fprintf(stderr, "[M::%s] read %d sequences\n", __func__, s->n_seq);
			return s;
		} else free(s);
	} else if (step == 1) {
		int i;
		tb_step_t *s = (tb_step_t*)_data;
		kt_for(aux->n_threads, tb_worker, s, s->n_seq);
		for (i = 0; i < s->n_seq; ++i) {
			printf("%s\t%d\t%d\t%d\t%d\n", s->seq[i].name, s->cnt[i].c[0], s->cnt[i].c[1], s->cnt[i].c[2], s->cnt[i].c[3]);
			free(s->seq[i].name); free(s->seq[i].seq); free(s->seq[i].qual); free(s->seq[i].comment);
		}
		free(s->seq); free(s->cnt); free(s);
	}
	return 0;
}

int main_triobin(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, min_cnt = 5;
	bfc_ch_t *ch;
	tb_shared_t aux;

	memset(&aux, 0, sizeof(tb_shared_t));
	aux.n_threads = 8, aux.print_diff = 0;
	while ((c = ketopt(&o, argc, argv, 1, "c:t:p", 0)) >= 0) {
		if (c == 'c') min_cnt = atoi(o.arg);
		else if (c == 't') aux.n_threads = atoi(o.arg);
		else if (c == 'p') aux.print_diff = 1;
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: yak triobin [options] <pat.yak> <mat.yak> <seq.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -c INT     min occurrence [%d]\n", min_cnt);
		fprintf(stderr, "  -t INT     number of threads [%d]\n", aux.n_threads);
		return 1;
	}

	ch = bfc_ch_flag_restore(0, argv[o.ind], min_cnt, 1);
	ch = bfc_ch_flag_restore(ch, argv[o.ind + 1], min_cnt, 2);

	aux.k = bfc_ch_get_k(ch);
	aux.fp = bseq_open(argv[o.ind+2]);
	aux.ch = ch;
	kt_pipeline(2, tb_pipeline, &aux, 2);
	bseq_close(aux.fp);
	bfc_ch_destroy(ch);
	return 0;
}
