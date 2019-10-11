#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "ketopt.h"
#include "kmer.h"
#include "yak.h"

int main_inspect(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, i, j, max_cnt = 10, k, l_pre, counter_bits, print_kmer = 0;
	uint32_t t[3];
	int64_t *cnt, tot[YAK_N_COUNTS], hist[YAK_N_COUNTS];
	uint64_t mask;
	FILE *fp;
	bfc_ch_t *ch = 0;
	char *fn_cmp = 0;

	while ((c = ketopt(&o, argc, argv, 1, "m:c:p", 0)) >= 0) {
		if (c == 'm') max_cnt = atoi(o.arg);
		else if (c == 'c') fn_cmp = o.arg;
		else if (c == 'p') print_kmer = 1;
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: yak inspect [options] <in.yak>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -c FILE   k-mer hash table for comparison []\n");
		fprintf(stderr, "  -m INT    max count (effective with -c) [%d]\n", max_cnt);
		fprintf(stderr, "  -p        print k-mers\n");
		return 1;
	}

	for (i = 0; i < YAK_N_COUNTS; ++i) hist[i] = tot[i] = 0;
	cnt = (int64_t*)calloc(YAK_N_COUNTS * (max_cnt + 1), sizeof(int64_t));
	if (fn_cmp) {
		ch = bfc_ch_restore(fn_cmp);
		bfc_ch_hist(ch, hist);
		assert(ch);
	}
	fp = fopen(argv[o.ind], "rb");
	assert(fp);
	fread(t, 4, 3, fp);
	k = t[0], l_pre = t[1], counter_bits = t[2];
	assert(counter_bits == YAK_COUNTER_BITS);
	mask = (1ULL<<k) - 1;
	if (print_kmer && l_pre + 64 < 2 * k + counter_bits) {
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
				if (cnt1 > max_cnt) cnt1 = max_cnt;
				if (cnt1 > 0) to_print = 0;
				++cnt[cnt1 * YAK_N_COUNTS + cnt0];
			}
			if (to_print) {
				uint64_t h[2], y[2];
				char buf[64];
				if (k <= 32) {
					uint64_t z = (uint64_t)i << (k*2 - l_pre) | key >> counter_bits;
					h[0] = z >> k;
					h[1] = z & mask;
				} else {
					h[0] = (uint64_t)i << (k - l_pre) | key >> (counter_bits + k);
					h[1] = key >> counter_bits & mask;
				}
				bfc_kmer_hash_inv(k, h, y);
				printf("KS\t%s\t%d\n", bfc_kmer_2str(k, y, buf), cnt0);
			}
		}
	}
	fclose(fp);
	bfc_ch_destroy(ch);

	if (fn_cmp) {
		int64_t acc_tot = 0, acc_cnt[YAK_N_COUNTS];
		for (i = 0; i < YAK_N_COUNTS; ++i) acc_cnt[i] = 0;
		for (j = max_cnt - 1; j >= 1; --j) {
			for (i = 0; i < YAK_N_COUNTS; ++i)
				cnt[j * YAK_N_COUNTS + i] += cnt[(j+1) * YAK_N_COUNTS + i];
		}
		for (i = YAK_N_COUNTS - 1; i >= 0; --i) {
			acc_tot += tot[i];
			if (acc_tot == 0) continue;
			printf("SN\t%d\t%ld\t%ld", i, (long)hist[i], (long)tot[i]);
			for (j = 1; j <= max_cnt; ++j) {
				acc_cnt[j] += cnt[j * YAK_N_COUNTS + i];
				printf("\t%.4f", (double)acc_cnt[j] / acc_tot);
			}
			printf("\n");
		}
	} else {
		int64_t acc_tot = 0;
		for (i = YAK_N_COUNTS - 1; i >= 0; --i) {
			acc_tot += tot[i];
			if (acc_tot == 0) continue;
			printf("HS\t%d\t%ld\t%ld\t%ld\n", i, (long)hist[i], (long)tot[i], (long)acc_tot);
		}
	}
	return 0;
}
