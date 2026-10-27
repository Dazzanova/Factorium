/*
 * factorium/src/pollard_p_minus_one.c — Pollard's p-1, stage 1.
 *
 * M is applied prime power by prime power: x = a^M mod n is built up as
 * x <- x^(p^k) for every prime p <= B with p^k <= B, which costs
 * O(log p^k) modular multiplications per prime instead of O(log M) once.
 */
#include "factorium/pollard_p_minus_one.h"
#include "factorium/modular_arithmetic.h"

#include <stdlib.h>
#include <string.h>

fct_pminus1_config fct_pminus1_config_default(void)
{
    fct_pminus1_config c;
    c.smoothness_bound = FCT_PM1_DEFAULT_SMOOTHNESS_BOUND;
    c.start_base = FCT_PM1_DEFAULT_START_BASE;
    c.max_bases = FCT_PM1_DEFAULT_MAX_BASES;
    c.max_modmul_ops = FCT_PM1_DEFAULT_MAX_MODMUL_OPS;
    return c;
}

/* Classic byte sieve of Eratosthenes over [0, B]; 1 marks a prime.
 * Returns NULL if allocation fails (caller only frees a non-NULL result). */
static unsigned char *sieve_primes(uint64_t B)
{
    if (B >= SIZE_MAX)
        return NULL;
    size_t len = (size_t)B + 1;
    unsigned char *s = (unsigned char *)malloc(len);
    if (s == NULL)
        return NULL;
    memset(s, 1, len);
    s[0] = 0;
    if (len > 1)
        s[1] = 0;
    for (uint64_t p = 2; p * p <= B; p++) {
        if (!s[p])
            continue;
        for (uint64_t m = p * p; m <= B; m += p)
            s[m] = 0;
    }
    return s;
}

fct_factor_result fct_pollard_p_minus_one(uint64_t n,
                                          const fct_pminus1_config *cfg)
{
    fct_factor_result res;
    res.status = FCT_NO_FACTOR;
    res.factor = 0;
    res.cofactor = 0;
    res.iterations = 0;
    res.gcd_calls = 0;
    res.modmul_ops = 0;
    res.seed = 0;
    res.c = 0;
    res.reason = FCT_TERM_BASES_EXHAUSTED;

    uint64_t B = (cfg && cfg->smoothness_bound != 0)
                     ? cfg->smoothness_bound
                     : FCT_PM1_DEFAULT_SMOOTHNESS_BOUND;

    if (n < 4 || B > FCT_PM1_MAX_SMOOTHNESS_BOUND) {
        res.status = FCT_INVALID_INPUT;
        res.reason = FCT_TERM_INVALID_INPUT;
        return res;
    }
    if ((n & 1) == 0) {
        res.status = FCT_OK;
        res.reason = FCT_TERM_FOUND;
        res.factor = 2;
        res.cofactor = n / 2;
        return res;
    }

    uint64_t start_base = (cfg && cfg->start_base != 0)
                              ? cfg->start_base
                              : FCT_PM1_DEFAULT_START_BASE;
    uint64_t max_bases = (cfg && cfg->max_bases != 0)
                             ? cfg->max_bases
                             : FCT_PM1_DEFAULT_MAX_BASES;
    uint64_t cap = (cfg && cfg->max_modmul_ops != 0)
                       ? cfg->max_modmul_ops
                       : FCT_PM1_DEFAULT_MAX_MODMUL_OPS;

    unsigned char *sieve = sieve_primes(B);
    if (sieve == NULL) {
        res.status = FCT_MEMORY;
        res.reason = FCT_TERM_NONE;
        return res;
    }

    for (uint64_t bi = 0; bi < max_bases; bi++) {
        res.iterations++;
        uint64_t a = (start_base + bi) % n;
        if (a == 0)
            continue; /* n divides the base: no information */

        /* Cheap pre-check: the base may already share a factor with n. */
        uint64_t g = fct_gcd(a, n);
        res.gcd_calls++;
        if (g > 1) { /* a < n implies g < n */
            res.status = FCT_OK;
            res.reason = FCT_TERM_FOUND;
            res.factor = g;
            res.cofactor = n / g;
            free(sieve);
            return res;
        }

        uint64_t x = a;
        for (uint64_t p = 2; p <= B; p++) {
            if (!sieve[p])
                continue;
            uint64_t e = p; /* largest prime power p^e <= B */
            while (e <= B / p)
                e *= p;
            uint64_t ops = 0;
            x = fct_mod_pow_counted(x, e, n, &ops);
            res.modmul_ops += ops;
            if (res.modmul_ops >= cap) {
                res.status = FCT_WORK_LIMIT;
                res.reason = FCT_TERM_WORK_LIMIT;
                free(sieve);
                return res;
            }
        }

        /* d = gcd(a^M - 1, n), computed without forming a^M - 1. */
        uint64_t xm1 = (x == 0) ? (n - 1) : (x - 1);
        uint64_t d = fct_gcd(xm1, n);
        res.gcd_calls++;
        if (d > 1 && d < n) {
            res.status = FCT_OK;
            res.reason = FCT_TERM_FOUND;
            res.factor = d;
            res.cofactor = n / d;
            free(sieve);
            return res;
        }
        /* d == 1: this base missed both factors.
         * d == n: a^M = 1 mod n (both factors smooth-minus-one, or n
         * prime); abandon this base — it is not a split. */
    }

    free(sieve);
    return res;
}
