#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "ketopt.h"
#include "kmer.h"
#include "yak.h"

int main_inspect(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, i, j, max_cnt = 10, kmer, l_pre, counter_bits, print_kmer = 0;
	uint32_t t[3];
	int64_t *cnt = 0, tot[YAK_N_COUNTS], hist[YAK_N_COUNTS];
	uint64_t mask;
	FILE *fp;
	double fpr = 0.00004;
	bfc_ch_t *ch = 0;
	char *fn_cmp = 0;

	while ((c = ketopt(&o, argc, argv, 1, "m:p", 0)) >= 0) {
		if (c == 'm') max_cnt = atoi(o.arg);
		else if (c == 'p') print_kmer = 1;
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: yak inspect [options] <in1.yak> [in2.yak]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -m INT    max count (effective with in2.yak) [%d]\n", max_cnt);
		fprintf(stderr, "  -p        print k-mers\n");
		fprintf(stderr, "Notes: when in2.yak is present, inspect evaluates the k-mer QV of in1.yak and\n");
		fprintf(stderr, "  the k-mer sensitivity of in2.yak.\n");
		return 1;
	}
	if (argc - o.ind >= 2) fn_cmp = argv[o.ind + 1];

	for (i = 0; i < YAK_N_COUNTS; ++i) hist[i] = tot[i] = 0;
	if (fn_cmp) {
		ch = bfc_ch_restore(fn_cmp);
		bfc_ch_hist(ch, hist);
		assert(ch);
		// cnt[in1.yak][in2.tak]
		cnt = (int64_t*)calloc(YAK_N_COUNTS * YAK_N_COUNTS, sizeof(int64_t));
	}
	fp = fopen(argv[o.ind], "rb");
	assert(fp);
	fread(t, 4, 3, fp);
	kmer = t[0], l_pre = t[1], counter_bits = t[2];
	assert(counter_bits == YAK_COUNTER_BITS);
	mask = (1ULL<<kmer) - 1;
	if (print_kmer && l_pre + 64 < 2 * kmer + counter_bits) {
		fprintf(stderr, "Warning: the hash table keeps hashes only; no k-mers are printed\n");
		print_kmer = 0;
	}

	for (i = 0; i < 1<<l_pre; ++i) {
		fread(t, 4, 2, fp);
		for (j = 0; j < t[1]; ++j) {
			uint64_t key;
			int cnt0, to_print = print_kmer;
			fread(&key, 8, 1, fp);
			cnt0 = key & ((1<<YAK_COUNTER_BITS) - 1);
			++tot[cnt0];
			if (ch) {
				int cnt1;
				cnt1 = bfc_ch_get_direct(ch, i, key);
				if (cnt1 < 0) cnt1 = 0;
				if (cnt1 > 0) to_print = 0;
				++cnt[cnt0 * YAK_N_COUNTS + cnt1];
			}
			if (to_print) {
				uint64_t h[2], y[2];
				char buf[64];
				if (kmer <= 32) {
					uint64_t z = (uint64_t)i << (kmer*2 - l_pre) | key >> counter_bits;
					h[0] = z >> kmer;
					h[1] = z & mask;
				} else {
					h[0] = (uint64_t)i << (kmer - l_pre) | key >> (counter_bits + kmer);
					h[1] = key >> counter_bits & mask;
				}
				bfc_kmer_hash_inv(kmer, h, y);
				printf("KS\t%s\t%d\n", bfc_kmer_2str(kmer, y, buf), cnt0);
			}
		}
	}
	fclose(fp);
	bfc_ch_destroy(ch);

	if (fn_cmp) {
		int64_t acc_tot = 0, *acc, acc_cnt[YAK_N_COUNTS];
		acc = (int64_t*)calloc(YAK_N_COUNTS * YAK_N_COUNTS, sizeof(int64_t));
		memcpy(acc, cnt, YAK_N_COUNTS * YAK_N_COUNTS * sizeof(int64_t));
		for (i = 0; i < YAK_N_COUNTS; ++i) acc_cnt[i] = 0;
		for (j = YAK_N_COUNTS - 2; j >= 1; --j)
			for (i = 0; i < YAK_N_COUNTS; ++i)
				acc[i * YAK_N_COUNTS + j] += acc[i * YAK_N_COUNTS + (j+1)];
		for (i = YAK_N_COUNTS - 1; i >= 0; --i) {
			acc_tot += tot[i];
			if (acc_tot == 0) continue;
			if (tot[i] == 0) continue;
			printf("SN\t%d\t%ld\t%ld", i, (long)tot[i], (long)hist[i]);
			for (j = 1; j <= max_cnt; ++j) {
				acc_cnt[j] += acc[i * YAK_N_COUNTS + j];
				printf("\t%.4f", (double)acc_cnt[j] / acc_tot);
			}
			printf("\n");
		}
		memcpy(acc, cnt, YAK_N_COUNTS * YAK_N_COUNTS * sizeof(int64_t));
		for (i = YAK_N_COUNTS - 2; i >= 0; --i)
			for (j = 0; j < YAK_N_COUNTS; ++j)
				acc[i * YAK_N_COUNTS + j] += acc[(i+1) * YAK_N_COUNTS + j];
		for (i = 1; i < YAK_N_COUNTS; ++i) {
			yak_qstat_t qs;
			if (tot[i] == 0) continue;
			yak_qv_solve(hist, &acc[i * YAK_N_COUNTS], kmer, fpr, &qs);
			printf("QV\t%d\t%ld\t%ld\t%.3f\t%.3f\n", i, (long)qs.tot, (long)acc[i * YAK_N_COUNTS], qs.qv_raw, qs.qv);
		}
		free(acc);
	} else {
		int64_t acc_tot = 0;
		for (i = YAK_N_COUNTS - 1; i >= 0; --i) {
			acc_tot += tot[i];
			if (acc_tot == 0) continue;
			printf("HS\t%d\t%ld\t%ld\t%ld\n", i, (long)hist[i], (long)tot[i], (long)acc_tot);
		}
	}
	free(cnt);
	return 0;
}
