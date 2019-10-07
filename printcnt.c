#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include "kmer.h"
#include "ketopt.h"
#include "yak.h"

int main_print(int argc, char *argv[])
{
	int c, subcnt_only = 0, hist_only = 0, l_pre, k, i, j;
	int min_cnt = 0, counter_bits = 0;
	ketopt_t o = KETOPT_INIT;
	uint32_t t[3];
	uint64_t mask, *hist;
	char buf[64];
	FILE *fp;

	while ((c = ketopt(&o, argc, argv, 1, "shm:", 0)) >= 0) {
		if (c == 's') subcnt_only = 1;
		else if (c == 'h') hist_only = 1;
		else if (c == 'm') min_cnt = atoi(o.arg);
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: yak print [options] <dump.hash>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -s       only show # elements in each sub- hash table\n");
		fprintf(stderr, "  -h       only show k-mer histogram\n");
		fprintf(stderr, "  -m INT   occ >= INT [0]\n");
		return 1;
	}

	if ((fp = fopen(argv[o.ind], "rb")) == 0) return 1;
	fread(t, 4, 3, fp);
	k = t[0], l_pre = t[1], counter_bits = t[2];
	hist = (uint64_t*)calloc(1<<counter_bits, sizeof(uint64_t));
	if (l_pre + 64 < 2 * k + counter_bits) {
		fprintf(stderr, "Warning: the hash table keeps hashes only; no k-mer counts. Print hash histogram only\n");
		hist_only = 1;
	}
	mask = (1ULL<<k) - 1;
	for (i = 0; i < 1<<l_pre; ++i) {
		uint64_t tmp;
		fread(t, 4, 2, fp);
		if (subcnt_only) printf("%d\n", t[1]);
		for (j = 0; j < t[1]; ++j) {
			int cnt;
			fread(&tmp, 8, 1, fp);
			cnt = tmp & ((1<<counter_bits) - 1);
			++hist[cnt];
			if (!subcnt_only && !hist_only && cnt >= min_cnt) {
				uint64_t h[2], y[2];
				if (k <= 32) {
					uint64_t z = (uint64_t)i << (k*2 - l_pre) | tmp >> counter_bits;
					h[0] = z >> k;
					h[1] = z & mask;
				} else {
					h[0] = (uint64_t)i << (k - l_pre) | tmp >> (counter_bits + k);
					h[1] = tmp >> counter_bits & mask;
				}
				bfc_kmer_hash_inv(k, h, y);
				printf("%s\t%d\n", bfc_kmer_2str(k, y, buf), cnt);
			}
		}
	}
	fclose(fp);

	if (hist_only)
		for (i = 0; i < 1<<counter_bits; ++i)
			printf("%d\t%llu\n", i, (unsigned long long)hist[i]);
	free(hist);
	return 0;
}
