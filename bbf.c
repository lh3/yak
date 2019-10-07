#include <stdlib.h>
#include <string.h>
#include "bbf.h"

bfc_bf_t *bfc_bf_init(int n_shift, int n_hashes)
{
	bfc_bf_t *b;
	void *ptr = 0;
	if (n_shift + BFC_BLK_SHIFT > 64 || n_shift < BFC_BLK_SHIFT) return 0;
	b = calloc(1, sizeof(bfc_bf_t));
	b->n_shift = n_shift;
	b->n_hashes = n_hashes;
	posix_memalign(&ptr, 1<<(BFC_BLK_SHIFT-3), 1ULL<<(n_shift-3));
	b->b = ptr;
	bzero(b->b, 1ULL<<(n_shift-3));
	return b;
}

void bfc_bf_destroy(bfc_bf_t *b)
{
	if (b == 0) return;
	free(b->b); free(b);
}

int bfc_bf_insert(bfc_bf_t *b, uint64_t hash)
{
	int x = b->n_shift - BFC_BLK_SHIFT;
	uint64_t y = hash & ((1ULL<<x) - 1);
	int h1 = hash >> x & BFC_BLK_MASK;
	int h2 = hash >> b->n_shift & BFC_BLK_MASK;
	uint8_t *p = &b->b[y<<(BFC_BLK_SHIFT-3)];
	int i, z = h1, cnt = 0;
	if ((h2&31) == 0) h2 = (h2 + 1) & BFC_BLK_MASK; // otherwise we may repeatedly use a few bits
	while (__sync_lock_test_and_set(p, 1)); // lock
	for (i = 0; i < b->n_hashes; z = (z + h2) & BFC_BLK_MASK) {
		uint8_t *q = &p[z>>3], u;
		if (p == q) continue; // don't use the first byte. It is a spin lock.
		u = 1<<(z&7);
		cnt += !!(*q & u);
		*q |= u;
		++i;
	}
	__sync_lock_release(p); // unlock
	return cnt;
}

int bfc_bf_get(const bfc_bf_t *b, uint64_t hash)
{
	int x = b->n_shift - BFC_BLK_SHIFT;
	uint64_t y = hash & ((1ULL<<x) - 1);
	int h1 = hash >> x & BFC_BLK_MASK;
	int h2 = hash >> b->n_shift & BFC_BLK_MASK;
	uint8_t *p = &b->b[y<<(BFC_BLK_SHIFT-3)];
	int i, z = h1, cnt = 0;
	if ((h2&31) == 0) h2 = (h2 + 1) & BFC_BLK_MASK; // otherwise we may repeatedly use a few bits
	for (i = 0; i < b->n_hashes; z = (z + h2) & BFC_BLK_MASK) {
		uint8_t *q = &p[z>>3];
		if (p == q) continue; // don't use the first byte. It is a spin lock.
		cnt += !!(*q & 1<<(z&7));
		++i;
	}
	return cnt;
}

static inline uint64_t popcount64(uint64_t y) // standard popcount; from wikipedia
{
	y -= ((y >> 1) & 0x5555555555555555ull);
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

uint64_t bfc_bf_load(const bfc_bf_t *b)
{
	uint64_t n = 0, i, size = b->n_shift-3-3;
	uint64_t *p = (uint64_t*)b->b;
	for (i = 0; i < size; ++i)
		n += popcount64(p[i]);
	return n;
}
