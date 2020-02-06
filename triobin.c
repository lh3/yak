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
	int c[16];
} tb_cnt_t;

typedef struct {
	int k, n_threads, print_diff;
	int count_thres;
	double ratio_thres;
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
				if (aux->print_diff && (flag>>2&3) != (flag&3))
					printf("D\t%s\t%d\t%d\t%d\n", s->name, i, flag&3, flag>>2&3);
			}
		} else l = 0, x = bfc_kmer_null;
	}
}

static char tb_classify(const int *c, double ratio_thres, int count_thres)
{
	char type = '0';
	if (c[0<<2|2] == c[2<<2|0]) { // equal counts
		if (c[0<<2|2] == 0) type = '0';
		else type = 'a';
	} else {
		if (c[0<<2|2] * ratio_thres >= c[2<<2|0] && c[0<<2|2] - c[2<<2|0] >= count_thres) type = 'p';
		else if (c[2<<2|0] * ratio_thres >= c[0<<2|2] && c[2<<2|0] - c[0<<2|2] >= count_thres) type = 'm';
		else type = 'a';
	}
	return type;
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
			int *c = s->cnt[i].c;
			char type;
			type = tb_classify(c, aux->ratio_thres, aux->count_thres);
			printf("%s\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", s->seq[i].name, type, c[0], c[2<<2|2], c[0<<2|2], c[2<<2|0], c[1<<2|1], c[1<<2|2], c[2<<2|1], c[0<<2|1], c[1<<2|0]);
			free(s->seq[i].name); free(s->seq[i].seq); free(s->seq[i].qual); free(s->seq[i].comment);
		}
		free(s->seq); free(s->cnt); free(s);
	}
	return 0;
}

int main_triobin(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, min_cnt = 2, mid_cnt = 5;
	bfc_ch_t *ch;
	tb_shared_t aux;

	memset(&aux, 0, sizeof(tb_shared_t));
	aux.n_threads = 8, aux.print_diff = 0;
	aux.count_thres = 5, aux.ratio_thres = 0.2;
	while ((c = ketopt(&o, argc, argv, 1, "c:d:t:pr:n:", 0)) >= 0) {
		if (c == 'c') min_cnt = atoi(o.arg);
		else if (c == 'd') mid_cnt = atoi(o.arg);
		else if (c == 't') aux.n_threads = atoi(o.arg);
		else if (c == 'p') aux.print_diff = 1;
		else if (c == 'r') aux.ratio_thres = atof(o.arg);
		else if (c == 'n') aux.count_thres = atoi(o.arg);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: yak triobin [options] <pat.yak> <mat.yak> <seq.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -c INT     min occurrence [%d]\n", min_cnt);
		fprintf(stderr, "  -d INT     mid occurrence [%d]\n", mid_cnt);
		fprintf(stderr, "  -t INT     number of threads [%d]\n", aux.n_threads);
		fprintf(stderr, "Output: ctg err strongMixed sPat sMat weakMixed wPat1 wMat1 wPat2 wMat2\n");
		return 1;
	}

	ch = bfc_ch_restore_core(0,  argv[o.ind],     YAK_LOAD_TRIOBIN1, min_cnt, mid_cnt);
	ch = bfc_ch_restore_core(ch, argv[o.ind + 1], YAK_LOAD_TRIOBIN2, min_cnt, mid_cnt);

	aux.k = bfc_ch_get_k(ch);
	aux.fp = bseq_open(argv[o.ind+2]);
	aux.ch = ch;
	kt_pipeline(2, tb_pipeline, &aux, 2);
	bseq_close(aux.fp);
	bfc_ch_destroy(ch);
	return 0;
}
