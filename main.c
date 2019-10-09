#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include "ketopt.h"
#include "yak.h"
#include "sys.h"

static inline int64_t mm_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9;
	else if (*p == 'M' || *p == 'm') x *= 1e6;
	else if (*p == 'K' || *p == 'k') x *= 1e3;
	return (int64_t)(x + .499);
}

static void usage_count(FILE *fp, yak_copt_t *o)
{
	fprintf(fp, "Usage: yak count [options] <in.fq> [in.fq]\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -p INT       prefix length [%d]\n", o->b_pre);
	fprintf(fp, "  -k INT       k-mer length [%d]\n", o->k);
	fprintf(fp, "  -b INT       set Bloom filter size to 2**INT bits; 0 to disable [%d]\n", o->bf_shift);
	fprintf(fp, "  -H INT       use INT hash functions for Bloom filter [%d]\n", o->bf_n_hashes);
	fprintf(fp, "  -t INT       number of threads [%d]\n", o->n_threads);
	fprintf(fp, "  -o FILE      write counts to FILE []\n");
	fprintf(fp, "  -K NUM       batch size [%ld]\n", (long)o->chunk_size);
	fprintf(fp, "Note: -b37 is recommended for human reads\n");
}

int main_count(int argc, char *argv[])
{
	yak_copt_t opt;
	bfc_ch_t *ch = 0;
	ketopt_t o = KETOPT_INIT;
	char *out_fn = 0;
	int c;

	yak_copt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "k:b:t:H:K:v:o:p:", 0)) >= 0) {
		if (c == 'p') opt.b_pre = atoi(o.arg);
		else if (c == 'k') opt.k = atoi(o.arg);
		else if (c == 'b') opt.bf_shift = atoi(o.arg);
		else if (c == 'H') opt.bf_n_hashes = atoi(o.arg);
		else if (c == 't') opt.n_threads = atoi(o.arg);
		else if (c == 'v') yak_verbose = atoi(o.arg);
		else if (c == 'K') opt.chunk_size = mm_parse_num(o.arg);
		else if (c == 'o') out_fn = o.arg;
	}

	if (o.ind == argc) {
		usage_count(stderr, &opt);
		return 1;
	}
	if (opt.b_pre + 64 < 2 * opt.k + YAK_COUNTER_BITS)
		fprintf(stderr, "Warning: counting hashes instead. Please increase -p or reduce -k.\n");

	ch = bfc_count(argv[o.ind], &opt, 0);
	if (o.ind + 1 < argc) {
		bfc_ch_reset(ch);
		ch = bfc_count(argv[o.ind+1], &opt, ch);
	}
	if (out_fn) bfc_ch_dump(ch, out_fn);
	bfc_ch_destroy(ch);
	return 0;
}

int main_qv(int argc, char *argv[])
{
	yak_qopt_t opt;
	bfc_ch_t *ch = 0;
	ketopt_t o = KETOPT_INIT;
	int64_t cnt[1<<YAK_COUNTER_BITS], tot, acc;
	int c, i, kmer;

	yak_qopt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "K:t:l:f:p", 0)) >= 0) {
		if (c == 'K') opt.chunk_size = mm_parse_num(o.arg);
		else if (c == 'l') opt.min_len = mm_parse_num(o.arg);
		else if (c == 'f') opt.min_frac = atof(o.arg);
		else if (c == 't') opt.n_threads = atoi(o.arg);
		else if (c == 'p') opt.print_each = 1;
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: yak qv [options] <kmer.hash> <seq.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l NUM      min sequence length [%d]\n", opt.min_len);
		fprintf(stderr, "  -f FLOAT    min k-mer fraction [%g]\n", opt.min_frac);
		fprintf(stderr, "  -p          print QV for each sequence\n");
		fprintf(stderr, "  -t INT      number of threads [%d]\n", opt.n_threads);
		fprintf(stderr, "  -K NUM      batch size [1g]\n");
		return 1;
	}

	ch = bfc_ch_restore(argv[o.ind]);
	assert(ch);
	kmer = bfc_ch_get_k(ch);
	yak_qv(&opt, argv[o.ind+1], ch, cnt);
	for (i = 0, tot = 0; i < 1<<YAK_COUNTER_BITS; ++i)
		tot += cnt[i];
	for (i = (1<<YAK_COUNTER_BITS) - 1, acc = 0; i >= 1; --i) {
		double x;
		acc += cnt[i];
		x = log((double)tot / acc) / kmer;
		x = -10.0 * log(x) / log(10);
		printf("Q\t%d\t%ld\t%.3f\n", i, (long)cnt[i], x);
	}
	printf("Q\t0\t%ld\t0\n", (long)cnt[0]);
	bfc_ch_destroy(ch);
	return 0;
}

int main(int argc, char *argv[])
{
	extern int main_inspect(int argc, char *argv[]);
	int ret = 0, i;
	yak_reset_realtime();
	if (argc == 1) {
		fprintf(stderr, "Usage: yak <command> <argument>\n");
		fprintf(stderr, "Command:\n");
		fprintf(stderr, "  count     count k-mers\n");
		fprintf(stderr, "  inspect   inspect k-mer hash tables\n");
		fprintf(stderr, "  qv        evaluate quality values\n");
		fprintf(stderr, "  version   print version number\n");
		return 1;
	}
	if (strcmp(argv[1], "count") == 0) ret = main_count(argc-1, argv+1);
	else if (strcmp(argv[1], "qv") == 0) ret = main_qv(argc-1, argv+1);
	else if (strcmp(argv[1], "inspect") == 0) ret = main_inspect(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		puts(YAK_VERSION);
		return 0;
	} else {
		fprintf(stderr, "[E::%s] unknown command\n", __func__);
		return 1;
	}
	if (ret == 0) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, YAK_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, yak_realtime(), yak_cputime(), yak_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return ret;
}
