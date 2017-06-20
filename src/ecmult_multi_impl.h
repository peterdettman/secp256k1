/**********************************************************************
 * Copyright (c) 2017 Andrew Poelstra                                 *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#include "ecmult_multi.h"

/* Heap operations: parent(i) = i/2; left_child(i) = 2*i; right_child(i) = 2*i + 1 */
static void secp256k1_heap_insert(const secp256k1_scalar *sc, unsigned char *idx, size_t *n, unsigned char ins_idx) {
    size_t ins = *n + 1;
    while (ins > 1 && secp256k1_scalar_cmp_var(&sc[idx[ins / 2 - 1]], &sc[ins_idx]) < 0) {
        idx[ins - 1] = idx[ins / 2 - 1];
        ins /= 2;
    }
    idx[ins - 1] = ins_idx;
    *n += 1;
}

#define SWAP(i, j) (idx[i] ^= idx[j], idx[j] ^= idx[i], idx[i] ^= idx[j])
static unsigned char secp256k1_heap_remove(const secp256k1_scalar *sc, unsigned char *idx, size_t *n) {
    unsigned char ret;
    size_t rem = 1;
    VERIFY_CHECK(*n > 0);
    /* swap-delete the root */
    ret = idx[0];
    idx[0] = idx[*n - 1];
    /* sift the new root into the correct place */
    for (;;) {
        /* if parent < lchild... */
        if (2 * rem - 1 < *n && secp256k1_scalar_cmp_var(&sc[idx[rem - 1]], &sc[idx[2 * rem - 1]]) < 0) {
            /* ...and if lchild < rchild, then parent swap with rchild */
            if (2 * rem < *n && secp256k1_scalar_cmp_var(&sc[idx[2 * rem - 1]], &sc[idx[2 * rem]]) < 0) {
                SWAP(rem - 1, 2 * rem);
                rem = 2 * rem + 1;
            /* ...and if lchild >= rchild, then parent swap with lchild */
            } else {
                SWAP(rem - 1, 2 * rem - 1);
                rem = 2 * rem;
            }
        /* if parent >= lchild... */
        } else {
            /* ...and if parent < rchild, then parent swap with rchild */
            if (2 * rem < *n && secp256k1_scalar_cmp_var(&sc[idx[rem - 1]], &sc[idx[2 * rem]]) < 0) {
                SWAP(rem - 1, 2 * rem);
                rem = 2 * rem + 1;
            /* ...and if parent >= rchild, then we're done */
            } else {
                break;
            }
        }
    }
    *n -= 1;
    return ret;
}
#undef SWAP

/** Multi-multiply: R = sum_i ni * Ai */
static void secp256k1_ecmult_multi(secp256k1_gej *r, secp256k1_scalar *sc, secp256k1_gej *pt, size_t n) {
    unsigned char heap_idx[SECP256K1_ECMULT_MULTI_MAX_N];
    size_t heap_n = 0;
    size_t i = 0;

    VERIFY_CHECK(n <= SECP256K1_ECMULT_MULTI_MAX_N);

    for (i = 0; i < n; i++) {
        if (!secp256k1_scalar_is_zero(&sc[i])) {
            secp256k1_heap_insert(sc, heap_idx, &heap_n, i);
        }
    }

    if (heap_n == 0) {
        secp256k1_gej_set_infinity(r);
        return;
    }

    while (heap_n > 1) {
        int max1i = secp256k1_heap_remove(sc, heap_idx, &heap_n);
        int max2i = heap_idx[0];  /* Don't even remove, we're just going to put it back in */
        /* Observe that nX + mY = (n-m)X + m(X + Y), and if n > m this transformation
         * reduces the magnitude of the larger scalar, on average by half. So by
         * repeating this we will quickly zero out all but one exponent, which will
         * be small. */
        secp256k1_gej_add_var(&pt[max2i], &pt[max1i], &pt[max2i], NULL);  /* Y -> X + Y */
        if (!secp256k1_scalar_eq(&sc[max1i], &sc[max2i])) {
            secp256k1_scalar_numsub(&sc[max1i], &sc[max1i], &sc[max2i]);  /* n -> n - m */
            secp256k1_heap_insert(sc, heap_idx, &heap_n, max1i);
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

