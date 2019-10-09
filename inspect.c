#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "ketopt.h"
#include "kmer.h"
#include "yak.h"

int main_inspect(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, i, j, min_cnt = 1, k, l_pre, counter_bits, print_kmer = 0;
	uint32_t t[3];
	uint64_t tot[1<<YAK_COUNTER_BITS], hit[1<<YAK_COUNTER_BITS], mask;
	FILE *fp;
	bfc_ch_t *ch = 0;
	char *fn_cmp = 0;

	while ((c = ketopt(&o, argc, argv, 1, "m:c:p", 0)) >= 0) {
		if (c == 'm') min_cnt = atoi(o.arg);
		else if (c == 'c') fn_cmp = o.arg;
		else if (c == 'p') print_kmer = 1;
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: yak inspect [options] <in.yak>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -c FILE   k-mer hash table for comparison []\n");
		fprintf(stderr, "  -m INT    min count (effective with -c) [%d]\n", min_cnt);
		fprintf(stderr, "  -p        print k-mers\n");
		return 1;
	}

	for (i = 0; i < 1<<YAK_COUNTER_BITS; ++i) tot[i] = hit[i] = 0;
	if (fn_cmp) {
		ch = bfc_ch_restore(fn_cmp);
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
				if (cnt1 >= min_cnt)
					++hit[cnt0], to_print = 0;
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
				printf("K\t%s\t%d\n", bfc_kmer_2str(k, y, buf), cnt0);
			}
		}
	}
	fclose(fp);
	bfc_ch_destroy(ch);

	if (fn_cmp) {
		uint64_t acc_tot = 0, acc_hit = 0;
		for (i = (1<<YAK_COUNTER_BITS) - 1; i >= 0; --i) {
			acc_tot += tot[i], acc_hit += hit[i];
			if (acc_tot == 0) continue;
			printf("C\t%d\t%lld\t%lld\t%lld\t%lld\t%.6f\n", i, (long long)tot[i], (long long)hit[i], (long long)acc_tot, (long long)acc_hit, (double)acc_hit / acc_tot);
		}
	} else {
		uint64_t acc_tot = 0;
		for (i = (1<<YAK_COUNTER_BITS) - 1; i >= 0; --i) {
			acc_tot += tot[i];
			if (acc_tot == 0) continue;
			printf("H\t%d\t%lld\t%lld\n", i, (long long)tot[i], (long long)acc_tot);
		}
	}
	return 0;
}
