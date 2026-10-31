/*
 * factorium/factorization.h — complete factorization of any uint64_t.
 *
 * The number is split recursively. At each composite node we try, in order:
 *   1. trial division by the primes up to 97,
 *   2. Fermat (cheap when the factors are close),
 *   3. Pollard Rho (general case),
 *   4. Pollard p-1 (when p - 1 is smooth).
 * Miller-Rabin decides when a node is prime.
 *
 * When the result is RESULT_COMPLETE the factors are prime, sorted and
 * their product was checked to equal n exactly. An incomplete result lists
 * the cofactors that could not be split; those cofactors are not claimed
 * to be prime.
 */
#ifndef FACTORIUM_FACTORIZATION_H
#define FACTORIUM_FACTORIZATION_H

#include "factorium/common.h"

#define FACTORIZE_DEFAULT_ITERATIONS (1ULL << 20)
#define FACTORIZE_DEFAULT_RETRIES 8
#define FACTORIZE_DEFAULT_BOUND 100000ULL
#define FACTORIZE_DEFAULT_BASES 5
#define FACTORIZE_DEFAULT_SEED 1
#define FACTORIZE_DEFAULT_DEPTH 64
#define FACTORIZE_DEFAULT_MAX_FACTORS 64

/* A uint64_t has at most 64 prime factors (2^64 = 2 * 2 * ...). */
#define MAX_FACTORS 64

typedef struct {
    uint64_t max_iterations;   /* Fermat and Rho budget per node    */
    uint64_t max_retries;      /* Rho attempts per node             */
    uint64_t seed;             /* Rho seed                          */
    uint64_t smoothness_bound; /* p-1 bound B                       */
    uint64_t max_bases;        /* p-1 bases per node                */
    uint64_t max_depth;        /* recursion depth cap               */
    uint64_t max_factors;      /* factor array cap (<= MAX_FACTORS) */
} FactorizeOptions;

FactorizeOptions default_options(void);

typedef struct {
    ResultCode code;
    uint64_t n;
    uint64_t factors[MAX_FACTORS];
    int factor_count;
    uint64_t unresolved[MAX_FACTORS];
    int unresolved_count;
    uint64_t product;
    bool product_verified;
    uint64_t iterations;
    uint64_t gcd_calls;
    uint64_t modmul_ops;
} Factorization;

/* Factors n >= 2. Pass NULL for the default options. */
Factorization factorize(uint64_t n, const FactorizeOptions *opt);

#endif /* FACTORIUM_FACTORIZATION_H */
