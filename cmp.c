#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "ketopt.h"
#include "yak.h"

int main_cmp(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, min_cnt = 1, l_pre, k, i, j, counter_bits = 0;
	uint32_t t[3];
	uint64_t tot[1<<YAK_COUNTER_BITS], hit[1<<YAK_COUNTER_BITS], acc_tot, acc_hit;
	FILE *fp;
	bfc_ch_t *ch;

	while ((c = ketopt(&o, argc, argv, 1, "m:", 0)) >= 0) {
		if (c == 'm') min_cnt = atoi(o.arg);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: yak cmp [options] <base.yak> <cmp.yak>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -m INT    min count [%d]\n", min_cnt);
		return 1;
	}

	for (i = 0; i < 1<<YAK_COUNTER_BITS; ++i) tot[i] = hit[i] = 0;
	ch = bfc_ch_restore(argv[o.ind + 1]);
	assert(ch);
	fp = fopen(argv[o.ind], "rb");
	assert(fp);
	fread(t, 4, 3, fp);
	k = t[0], l_pre = t[1], counter_bits = t[2];
	assert(counter_bits == YAK_COUNTER_BITS);

	for (i = 0; i < 1<<l_pre; ++i) {
		fread(t, 4, 2, fp);
		for (j = 0; j < t[1]; ++j) {
			uint64_t key;
			int cnt1, cnt0;
			fread(&key, 8, 1, fp);
			cnt0 = key & ((1<<YAK_COUNTER_BITS) - 1);
			cnt1 = bfc_ch_get_direct(ch, i, key);
			++tot[cnt0];
			if (cnt1 >= min_cnt) ++hit[cnt0];
		}
	}
	fclose(fp);
	bfc_ch_destroy(ch);

	for (i = (1<<YAK_COUNTER_BITS) - 1, acc_tot = acc_hit = 0; i >= 0; --i) {
		acc_tot += tot[i], acc_hit += hit[i];
		if (acc_tot == 0) continue;
		printf("%d\t%lld\t%lld\t%lld\t%.6f\n", i, (long long)tot[i], (long long)acc_tot, (long long)acc_hit, (double)acc_hit / acc_tot);
	}
	return 0;
}
