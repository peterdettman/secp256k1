/**********************************************************************
 * Copyright (c) 2017 Andrew Poelstra                                 *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#include "ecmult_multi.h"

/* Heap operations: parent(i) = i/2; left_child(i) = 2*i; right_child(i) = 2*i + 1 */
static void secp256k1_heap_insert(secp256k1_scalar *sc, secp256k1_gej *pt, size_t *n, const secp256k1_scalar *ins_sc, const secp256k1_gej *ins_pt) {
    size_t ins = *n + 1;
    while (ins > 1 && secp256k1_scalar_cmp_var(&sc[ins / 2 - 1], ins_sc) < 0) {
        sc[ins - 1] = sc[ins / 2 - 1];
        pt[ins - 1] = pt[ins / 2 - 1];
        ins /= 2;
    }
    sc[ins - 1] = *ins_sc;
    pt[ins - 1] = *ins_pt;
    *n += 1;
}

#define SWAP(i, j) do{\
    secp256k1_scalar tmps;\
    secp256k1_gej tmpj;\
    tmps = sc[i];\
    tmpj = pt[i];\
    sc[i] = sc[j];\
    pt[i] = pt[j];\
    sc[j] = tmps;\
    pt[j] = tmpj;\
} while(0)
static void secp256k1_heap_remove(secp256k1_scalar *sc, secp256k1_gej *pt, size_t *n, secp256k1_scalar *out_sc, secp256k1_gej *out_pt) {
    size_t rem = 1;
    VERIFY_CHECK(*n > 0);
    /* swap-delete the root */
    *out_sc = sc[0];
    *out_pt = pt[0];
    sc[0] = sc[*n - 1];
    pt[0] = pt[*n - 1];
    /* sift the new root into the correct place */
    for (;;) {
        /* if parent < lchild... */
        if (2 * rem - 1 < *n && secp256k1_scalar_cmp_var(&sc[rem - 1], &sc[2 * rem - 1]) < 0) {
            /* ...and if lchild < rchild, then parent swap with rchild */
            if (2 * rem < *n && secp256k1_scalar_cmp_var(&sc[2 * rem - 1], &sc[2 * rem]) < 0) {
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
            if (2 * rem < *n && secp256k1_scalar_cmp_var(&sc[rem - 1], &sc[2 * rem]) < 0) {
                SWAP(rem - 1, 2 * rem);
                rem = 2 * rem + 1;
            /* ...and if parent >= rchild, then we're done */
            } else {
                break;
            }
        }
    }
    *n -= 1;
}
#undef SWAP

/** Multi-multiply: R = sum_i ni * Ai */
static void secp256k1_ecmult_multi(secp256k1_gej *r, const secp256k1_scalar *sc, const secp256k1_gej *pt, size_t n) {
    secp256k1_scalar heap_sc[SECP256K1_ECMULT_MULTI_MAX_N];
    secp256k1_gej heap_pt[SECP256K1_ECMULT_MULTI_MAX_N];
    size_t heap_n = 0;
    size_t i = 0;

    VERIFY_CHECK(n <= SECP256K1_ECMULT_MULTI_MAX_N);

    for (i = 0; i < n; i++) {
        if (!secp256k1_scalar_is_zero(&sc[i])) {
            secp256k1_heap_insert(heap_sc, heap_pt, &heap_n, &sc[i], &pt[i]);
        }
    }

    if (heap_n == 0) {
        secp256k1_gej_set_infinity(r);
        return;
    }

    while (heap_n > 1) {
        secp256k1_scalar max_s;
        secp256k1_scalar max_s1;
        secp256k1_gej max_p;
        secp256k1_gej max_p1;

        secp256k1_heap_remove(heap_sc, heap_pt, &heap_n, &max_s, &max_p);
        secp256k1_heap_remove(heap_sc, heap_pt, &heap_n, &max_s1, &max_p1);
        /* Observe that nX + mY = (n-m)X + m(X + Y), and if n > m this transformation
         * reduces the magnitude of the larger scalar, on average by half. So by
         * repeating this we will quickly zero out all but one exponent, which will
         * be small. */
        secp256k1_gej_add_var(&max_p1, &max_p, &max_p1, NULL);  /* Y -> X + Y */
        secp256k1_heap_insert(heap_sc, heap_pt, &heap_n, &max_s1, &max_p1);
        if (!secp256k1_scalar_eq(&max_s, &max_s1)) {
            secp256k1_scalar_negate(&max_s1, &max_s1);
            secp256k1_scalar_add(&max_s, &max_s, &max_s1);    /* n -> n - m */
            secp256k1_heap_insert(heap_sc, heap_pt, &heap_n, &max_s, &max_p);
        }
    }
    VERIFY_CHECK(heap_n == 1);
    VERIFY_CHECK(!secp256k1_scalar_is_zero(&heap_sc[0]));

    /* Now the desired result is heap_sc[0] * heap_pt[0], and for random scalars it is
     * very likely that heap_sc[0] = 1, and extremely likely heap_sc[0] < 5. (After
     * about 100k trials I saw around 200 2's and one 3.) So use a binary ladder rather
     * than any heavy machinery to finish it off. */
    secp256k1_gej_set_infinity(r);
    if (!secp256k1_gej_is_infinity(&heap_pt[0])) {
        while (!secp256k1_scalar_is_zero(&heap_sc[0])) {
            if (secp256k1_scalar_shr_int(&heap_sc[0], 1) == 1) {
                secp256k1_gej_add_var(r, r, &heap_pt[0], NULL);
            }
            secp256k1_gej_double_nonzero(&heap_pt[0], &heap_pt[0], NULL);
        }
    }
}

