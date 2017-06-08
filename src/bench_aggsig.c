/**********************************************************************
 * Copyright (c) 2014-2015 Pieter Wuille                              *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#include "include/secp256k1.h"
#include "include/secp256k1_aggsig.h"
#include "util.h"
#include "bench.h"

typedef struct {
    secp256k1_context *ctx;
    unsigned char msg[32];
    unsigned char sig[64];
} bench_aggsig_t;

void bench_aggsig(void* arg) {
    /* TODO */
}

void bench_aggsig_setup(void* arg) {
    /* TODO */
}

int main(void) {
    bench_aggsig_t data;

    data.ctx = secp256k1_context_create(SECP256K1_CONTEXT_VERIFY);

    run_benchmark("aggsig_10", bench_aggsig, bench_aggsig_setup, NULL, &data, 10, 20000);

    secp256k1_context_destroy(data.ctx);
    return 0;
}
