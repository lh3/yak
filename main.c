#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <unistd.h>
#include <math.h>
#include <sys/stat.h>
#include "ketopt.h"
#include "yak-priv.h"

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

int main_count(int argc, char *argv[])
{
	yak_ch_t *h;
	int c;
	char *fn_out = 0;
	yak_copt_t opt;
	ketopt_t o = KETOPT_INIT;
	yak_copt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "k:p:K:t:b:H:o:", 0)) >= 0) {
		if (c == 'k') opt.k = atoi(o.arg);
		else if (c == 'p') opt.pre = atoi(o.arg);
		else if (c == 'K') opt.chunk_size = atoi(o.arg);
		else if (c == 't') opt.n_thread = atoi(o.arg);
		else if (c == 'b') opt.bf_shift = atoi(o.arg);
		else if (c == 'H') opt.bf_n_hash = mm_parse_num(o.arg);
		else if (c == 'o') fn_out = o.arg;
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: yak count [options] <in.fa> [in.fa]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", opt.k);
		fprintf(stderr, "  -p INT     prefix length [%d]\n", opt.pre);
		fprintf(stderr, "  -b INT     set Bloom filter size to 2**INT bits; 0 to disable [%d]\n", opt.bf_shift);
		fprintf(stderr, "  -H INT     use INT hash functions for Bloom filter [%d]\n", opt.bf_n_hash);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", opt.n_thread);
		fprintf(stderr, "  -o FILE    dump the count hash table to FILE []\n");
		fprintf(stderr, "  -K INT     chunk size [100m]\n");
		fprintf(stderr, "Note: -b37 is recommended for human reads\n");
		return 1;
	}
	if (opt.pre < YAK_COUNTER_BITS) {
		fprintf(stderr, "ERROR: -p should be at least %d\n", YAK_COUNTER_BITS);
		return 1;
	}
	if (opt.k >= 64) {
		fprintf(stderr, "ERROR: -k must be smaller than 64\n");
		return 1;
	} else if (opt.k >= 32) {
		fprintf(stderr, "WARNING: counts are inexact if -k is greater than 31\n");
	}
	h = yak_count(argv[o.ind], &opt, 0);
	if (opt.bf_shift > 0) {
		yak_ch_destroy_bf(h);
		yak_ch_clear(h, opt.n_thread);
		h = yak_count(argc - o.ind >= 2? argv[o.ind+1] : argv[o.ind], &opt, h);
		yak_ch_shrink(h, 2, YAK_MAX_COUNT, opt.n_thread);
		fprintf(stderr, "[M::%s] %ld distinct k-mers after shrinking\n", __func__, (long)h->tot);
	}
	if (fn_out) yak_ch_dump(h, fn_out);
	yak_ch_destroy(h);
	return 0;
}

int main_qv(int argc, char *argv[])
{
	yak_qopt_t opt;
	yak_ch_t *ch = 0;
	ketopt_t o = KETOPT_INIT;
	int64_t cnt[YAK_N_COUNTS], hist[YAK_N_COUNTS];
	int c, i, kmer;
	yak_qstat_t qs;

	yak_qopt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "K:t:l:f:pe:", 0)) >= 0) {
		if (c == 'K') opt.chunk_size = mm_parse_num(o.arg);
		else if (c == 'l') opt.min_len = mm_parse_num(o.arg);
		else if (c == 'f') opt.min_frac = atof(o.arg);
		else if (c == 't') opt.n_threads = atoi(o.arg);
		else if (c == 'p') opt.print_each = 1;
		else if (c == 'e') opt.fpr = atof(o.arg);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: yak qv [options] <kmer.hash> <seq.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l NUM      min sequence length [%d]\n", opt.min_len);
		fprintf(stderr, "  -f FLOAT    min k-mer fraction [%g]\n", opt.min_frac);
		fprintf(stderr, "  -e FLOAT    false positive rate [%g]\n", opt.fpr);
		fprintf(stderr, "  -p          print QV for each sequence\n");
		fprintf(stderr, "  -t INT      number of threads [%d]\n", opt.n_threads);
		fprintf(stderr, "  -K NUM      batch size [1g]\n");
		return 1;
	}

	ch = yak_ch_restore(argv[o.ind]);
	assert(ch);
	kmer = ch->k;
	yak_ch_hist(ch, hist, opt.n_threads);
	printf("CC\tCT  kmer_occurrence    short_read_kmer_count  raw_input_kmer_count  adjusted_input_kmer_count\n");
	printf("CC\tFR  fpr_lower_bound    fpr_upper_bound\n");
	printf("CC\tER  total_input_kmers  adjusted_error_kmers\n");
	printf("CC\tCV  coverage\n");
	printf("CC\tQV  raw_quality_value  adjusted_quality_value\n");
	printf("CC\n");
	yak_qv(&opt, argv[o.ind+1], ch, cnt);
	yak_qv_solve(hist, cnt, kmer, opt.fpr, &qs);
	for (i = (1<<YAK_COUNTER_BITS) - 1; i >= 0; --i)
		printf("CT\t%d\t%ld\t%ld\t%.3f\n", i, (long)hist[i], (long)cnt[i], qs.adj_cnt[i]);
	printf("FR\t%.3g\t%.3g\n", qs.fpr_lower, qs.fpr_upper);
	printf("ER\t%ld\t%.3f\n", (long)qs.tot, qs.err);
	printf("CV\t%.3f\n", qs.cov);
	printf("QV\t%.3f\t%.3f\n", qs.qv_raw, qs.qv);
	yak_ch_destroy(ch);
	return 0;
}

int main(int argc, char *argv[])
{
	extern int main_triobin(int argc, char *argv[]);
	extern int main_trioeval(int argc, char *argv[]);
	extern int main_inspect(int argc, char *argv[]);
	int ret = 0, i;
	yak_reset_realtime();
	if (argc == 1) {
		fprintf(stderr, "Usage: yak <command> <argument>\n");
		fprintf(stderr, "Command:\n");
		fprintf(stderr, "  count     count k-mers\n");
		fprintf(stderr, "  qv        evaluate quality values\n");
		fprintf(stderr, "  triobin   trio binning\n");
		fprintf(stderr, "  trioeval  evaluate phasing accuracy with trio\n");
		fprintf(stderr, "  inspect   k-mer hash tables\n");
		fprintf(stderr, "  version   print version number\n");
		return 1;
	}
	if (strcmp(argv[1], "count") == 0) ret = main_count(argc-1, argv+1);
	else if (strcmp(argv[1], "qv") == 0) ret = main_qv(argc-1, argv+1);
	else if (strcmp(argv[1], "triobin") == 0) ret = main_triobin(argc-1, argv+1);
	else if (strcmp(argv[1], "trioeval") == 0) ret = main_trioeval(argc-1, argv+1);
	else if (strcmp(argv[1], "inspect") == 0) ret = main_inspect(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		puts(YAKS_VERSION);
		return 0;
	} else {
		fprintf(stderr, "[E::%s] unknown command\n", __func__);
		return 1;
	}
	if (ret == 0) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, YAKS_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, yak_realtime(), yak_cputime(), yak_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return ret;
}
