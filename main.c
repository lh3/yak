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

int yak_verbose = 3;

void yak_opt_init(yak_opt_t *opt)
{
	memset(opt, 0, sizeof(yak_opt_t));
	opt->chunk_size = 100000000;
	opt->k = 31;
	opt->b_pre = 20;
	opt->bf_shift = 0;
	opt->bf_n_hashes = 4;
	opt->n_threads = 4;
}

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

static void usage_count(FILE *fp, yak_opt_t *o)
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
	yak_opt_t opt;
	bfc_ch_t *ch = 0;
	ketopt_t o = KETOPT_INIT;
	char *out_fn = 0;
	int c;

	yak_opt_init(&opt);
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

	ch = (bfc_ch_t*)bfc_count(argv[o.ind], &opt);
	if (out_fn) bfc_ch_dump(ch, out_fn);
	bfc_ch_destroy(ch);
	return 0;
}

int main(int argc, char *argv[])
{
	extern int main_print(int argc, char *argv[]);
	int ret = 0, i;
	yak_reset_realtime();
	if (argc == 1) {
		fprintf(stderr, "Usage: yak <command> <argument>\n");
		fprintf(stderr, "Command:\n");
		fprintf(stderr, "  count     count k-mers\n");
		fprintf(stderr, "  print     print dumped k-mer counts or hash history\n");
		return 1;
	}
	if (strcmp(argv[1], "count") == 0) ret = main_count(argc-1, argv+1);
	else if (strcmp(argv[1], "print") == 0) ret = main_print(argc-1, argv+1);
	else {
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
