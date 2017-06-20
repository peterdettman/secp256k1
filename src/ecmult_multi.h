/**********************************************************************
 * Copyright (c) 2017 Andrew Poelstra                                 *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _SECP256K1_ECMULT_MULTI_
#define _SECP256K1_ECMULT_MULTI_

#include "group.h"
#include "scalar.h"

#define SECP256K1_ECMULT_MULTI_MAX_N	32

/** Multi-multiply: R = sum_i ni * Ai. Will trash the sc and pt arrays. */
static void secp256k1_ecmult_multi(secp256k1_gej *r, secp256k1_scalar *sc, secp256k1_gej *pt, size_t n);

#endif
