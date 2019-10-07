#ifndef BFC_BBF_H
#define BFC_BBF_H

#include <stdint.h>

#define BFC_BLK_SHIFT  9 // 64 bytes, the size of a cache line
#define BFC_BLK_MASK   ((1<<(BFC_BLK_SHIFT)) - 1)

typedef struct {
	int n_shift, n_hashes;
	uint8_t *b;
} bfc_bf_t;

bfc_bf_t *bfc_bf_init(int n_shift, int n_hashes);
void bfc_bf_destroy(bfc_bf_t *b);
int bfc_bf_insert(bfc_bf_t *b, uint64_t hash);
int bfc_bf_get(const bfc_bf_t *b, uint64_t hash);

#endif
