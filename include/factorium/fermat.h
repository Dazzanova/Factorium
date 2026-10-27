/*
 * factorium/fermat.h — Fermat's difference-of-squares factorization.
 *
 * Identity: n = a^2 - b^2 = (a - b)(a + b) with a = ceil(sqrt(n)),
 * b^2 = a^2 - n. Each step increments a and updates b^2 += 2a + 1
 * (no per-step division or modular arithmetic).
 *
 * Why close factors are favorable: for n = p*q (p <= q, odd), the first
 * solution is at a = (p+q)/2, i.e. after (q - p)/2 + 1 iterations. The
 * closer p and q are, the fewer iterations; equally the farther apart the
 * factors (e.g. a small prime times a large one), the slower Fermat is.
 *
 * A budget-exhausted result (FCT_NO_FACTOR) says only "no close factor
 * pair was found within the budget" — it is NOT a primality proof.
 *
 * The inner loop may use unsigned __int128 (documented compiler extension,
 * GCC/Clang) when a^2 would not fit in uint64_t; on platforms without it
 * the search stops with FCT_WORK_LIMIT once a reaches 2^32.
 */
#ifndef FACTORIUM_FERMAT_H
#define FACTORIUM_FERMAT_H

#include "factorium/common.h"

#define FCT_FERMAT_DEFAULT_MAX_ITERATIONS (1ULL << 20)

typedef struct {
    uint64_t max_iterations; /* 0 means FCT_FERMAT_DEFAULT_MAX_ITERATIONS */
} fct_fermat_config;

fct_fermat_config fct_fermat_config_default(void);

/* Finds a nontrivial factor of n >= 4.
 *  - even n returns factor 2 immediately;
 *  - n < 4 is FCT_INVALID_INPUT;
 *  - if a solution exists within the budget, FCT_OK with
 *    factor = a - b, cofactor = a + b (factor <= sqrt(n));
 *  - otherwise FCT_NO_FACTOR / FCT_TERM_BUDGET_EXHAUSTED.
 * gcd_calls and modmul_ops are always 0 for this algorithm. */
fct_factor_result fct_fermat(uint64_t n, const fct_fermat_config *cfg);

#endif /* FACTORIUM_FERMAT_H */
