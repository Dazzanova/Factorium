/*
 * factorium/factorization.h — shared complete-factorization layer.
 *
 * Splits a composite recursively by trying, in order, on each unresolved
 * node:
 *   1. trial division by the primes <= 97 (handles even and tiny factors),
 *   2. Fermat (cheap when the two sides are close),
 *   3. Pollard Rho (general case),
 *   4. Pollard p-1 (when p-1 / q-1 are smooth).
 * Miller-Rabin decides whether a node is prime and may be emitted.
 *
 * Guarantees on FCT_OK (complete):
 *   - factors are prime (deterministic Miller-Rabin),
 *   - sorted ascending, multiplicities preserved,
 *   - product_verified means the product of factors was reconstructed with
 *     checked multiplication and equals n exactly.
 *
 * On FCT_INCOMPLETE: at least one composite cofactor could not be split
 * within the configured limits; it is listed in `unresolved_cofactors`,
 * it is NOT included in `factors`, and product_verified is false. An
 * incomplete result is never a claim that the cofactors are prime.
 */
#ifndef FACTORIUM_FACTORIZATION_H
#define FACTORIUM_FACTORIZATION_H

#include "factorium/common.h"

#define FCT_FACTORIZATION_DEFAULT_MAX_ITERATIONS (1ULL << 20)
#define FCT_FACTORIZATION_DEFAULT_MAX_RETRIES 8
#define FCT_FACTORIZATION_DEFAULT_SMOOTHNESS_BOUND 100000ULL
#define FCT_FACTORIZATION_DEFAULT_MAX_BASES 5
#define FCT_FACTORIZATION_DEFAULT_SEED 1
/* A uint64_t has at most 64 prime factors (2^64), and the recursion is at
 * most one level deeper than that. */
#define FCT_FACTORIZATION_DEFAULT_MAX_FACTORS 64
#define FCT_FACTORIZATION_DEFAULT_MAX_DEPTH 64

typedef struct {
    uint64_t max_iterations;      /* per-attempt budget for Fermat and Rho;
                                   * 0 => default 2^20                  */
    uint64_t max_retries;         /* Rho attempts per node; 0 => default 8 */
    uint64_t smoothness_bound;    /* p-1 bound; 0 => default 100000     */
    uint64_t max_bases;           /* p-1 bases; 0 => default 5          */
    uint64_t seed;                /* Rho seed; 0 => default 1           */
    uint64_t max_depth;           /* recursion cap; 0 => default 64     */
    uint64_t max_factors;         /* factor-array cap; 0 => default 64  */
} fct_factorization_config;

fct_factorization_config fct_factorization_config_default(void);

typedef struct {
    fct_status status;            /* FCT_OK, FCT_INCOMPLETE, ...        */
    uint64_t n;
    uint64_t *factors;            /* heap; sorted ascending with
                                   * multiplicity; may be NULL if none  */
    size_t factor_count;
    uint64_t *unresolved_cofactors; /* heap; composites left unsplit     */
    size_t unresolved_count;
    uint64_t product;             /* product of factors                 */
    bool product_verified;        /* true iff product == n exactly      */
    uint64_t iterations;          /* summed over all algorithm calls    */
    uint64_t gcd_calls;
    uint64_t modmul_ops;
    fct_termination reason;       /* the limiting factor, if incomplete */
} fct_factorization;

/* Factors n (n >= 2). Never fails silently: every outcome is explicit.
 * The caller must release the result with fct_factorization_free(). */
fct_factorization fct_factorize(uint64_t n, const fct_factorization_config *cfg);

void fct_factorization_free(fct_factorization *f);

#endif /* FACTORIUM_FACTORIZATION_H */
