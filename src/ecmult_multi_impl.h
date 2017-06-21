/**********************************************************************
 * Copyright (c) 2017 Andrew Poelstra                                 *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#include "ecmult_multi.h"

/* Heap operations: parent(i) = i/2; left_child(i) = 2*i; right_child(i) = 2*i + 1 */
#define SWAP(i, j) (idx[i] ^= idx[j], idx[j] ^= idx[i], idx[i] ^= idx[j])
static void secp256k1_heap_siftdown(const secp256k1_scalar *sc, unsigned char *idx, size_t n, size_t root_idx) {
    while (2 * root_idx <= n) {
        size_t swap_idx = root_idx;
        /* If parent < lchild, swap with lchild */
        if (secp256k1_scalar_cmp_var(&sc[idx[root_idx - 1]], &sc[idx[2*root_idx - 1]]) < 0) {
            swap_idx = 2*root_idx;
        }
        /* If parent < rchild, and lchild < rchild, swap with rchild */
        if (2 * root_idx + 1 <= n && secp256k1_scalar_cmp_var(&sc[idx[swap_idx - 1]], &sc[idx[2*root_idx]]) < 0) {
            swap_idx = 2*root_idx + 1;
        }
        /* If we're swapping, do it */
        if (root_idx != swap_idx) {
            SWAP(root_idx - 1, swap_idx - 1);
            root_idx = swap_idx;
        } else {
            break;
        }
    }
}
#undef SWAP

static void secp256k1_heapify(const secp256k1_scalar *sc, unsigned char *idx, size_t n) {
    size_t i;
    for (i = n / 2; i > 0; i--) {
        secp256k1_heap_siftdown(sc, idx, n, i);
    }
}

static void secp256k1_heap_remove(const secp256k1_scalar *sc, unsigned char *idx, size_t *n) {
    VERIFY_CHECK(*n > 0);
    /* overwrite the root */
    idx[0] = idx[*n - 1];
    *n -= 1;
    /* sift the new root into place */
    secp256k1_heap_siftdown(sc, idx, *n, 1);
}

/** Multi-multiply: R = sum_i ni * Ai */
static void secp256k1_ecmult_multi(secp256k1_gej *r, secp256k1_scalar *sc, secp256k1_gej *pt, size_t n) {
    unsigned char heap_idx[SECP256K1_ECMULT_MULTI_MAX_N];
    size_t heap_n = 0;
    size_t i = 0;

    VERIFY_CHECK(n <= SECP256K1_ECMULT_MULTI_MAX_N);

    for (i = 0; i < n; i++) {
        if (!secp256k1_scalar_is_zero(&sc[i])) {
            heap_idx[heap_n++] = i;
        }
    }
    secp256k1_heapify(sc, heap_idx, heap_n);

    if (heap_n == 0) {
        secp256k1_gej_set_infinity(r);
        return;
    }

    while (heap_n > 1) {
        /* Don't remove the first-largest, we will update the second-largest first (not
         * moving it, its scalar doesn't change) and then update this one (moving it,
         * which will break our indices but that's OK as we'll do it last). */
        int max1i = heap_idx[0];
        /* Don't remove the second-largest, we won't update its key so we can work on
         * it in-place. */
        int max2i;
        if (heap_n == 2 || secp256k1_scalar_cmp_var(&sc[heap_idx[1]], &sc[heap_idx[2]]) > 0) {
            max2i = heap_idx[1];
        } else {
            max2i = heap_idx[2];
        }
        /* Observe that nX + mY = (n-m)X + m(X + Y), and if n > m this transformation
         * reduces the magnitude of the larger scalar, on average by half. So by
         * repeating this we will quickly zero out all but one exponent, which will
         * be small. */
        secp256k1_gej_add_var(&pt[max2i], &pt[max1i], &pt[max2i], NULL);  /* Y -> X + Y */
        if (!secp256k1_scalar_eq(&sc[max1i], &sc[max2i])) {
            secp256k1_scalar_numsub(&sc[max1i], &sc[max1i], &sc[max2i]);  /* n -> n - m */
            /* sift the updated root elem into place */
            secp256k1_heap_siftdown(sc, heap_idx, heap_n, 1);
        } else {
            secp256k1_heap_remove(sc, heap_idx, &heap_n);
        }
    }
    VERIFY_CHECK(heap_n == 1);
    VERIFY_CHECK(!secp256k1_scalar_is_zero(&sc[heap_idx[0]]));

    /* Now the desired result is heap_sc[0] * heap_pt[0], and for random scalars it is
     * very likely that heap_sc[0] = 1, and extremely likely heap_sc[0] < 5. (After
     * about 100k trials I saw around 200 2's and one 3.) So use a binary ladder rather
     * than any heavy machinery to finish it off. */
    secp256k1_gej_set_infinity(r);
    if (!secp256k1_gej_is_infinity(&pt[heap_idx[0]])) {
        while (!secp256k1_scalar_is_zero(&sc[heap_idx[0]])) {
            if (secp256k1_scalar_shr_int(&sc[heap_idx[0]], 1) == 1) {
                secp256k1_gej_add_var(r, r, &pt[heap_idx[0]], NULL);
            }
            secp256k1_gej_double_nonzero(&pt[heap_idx[0]], &pt[heap_idx[0]], NULL);
        }
    }
}

