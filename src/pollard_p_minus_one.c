/*
 * factorium/src/pollard_p_minus_one.c — Pollard's p-1, stage 1.
 *
 * Instead of building the huge number M, we raise x to p^k for every prime
 * power p^k <= B. The result is the same a^M mod n, with fewer
 * multiplications.
 */
#include "factorium/pollard_p_minus_one.h"
#include "factorium/modular_arithmetic.h"

#include <stdlib.h>
#include <string.h>

/* Sieve of Eratosthenes: s[i] = 1 means i is prime. NULL on failure. */
static unsigned char *make_sieve(uint64_t B)
{
    size_t len;
    unsigned char *s;
    uint64_t p;

    if (B >= (uint64_t)SIZE_MAX)
        return NULL;
    len = (size_t)B + 1;
    s = (unsigned char *)malloc(len);
    if (s == NULL)
        return NULL;
    memset(s, 1, len);
    s[0] = 0;
    if (len > 1)
        s[1] = 0;
    for (p = 2; p * p <= B; p++) {
        if (!s[p])
            continue;
        for (uint64_t m = p * p; m <= B; m += p)
            s[m] = 0;
    }
    return s;
}

FactorResult pminus1(uint64_t n, uint64_t B, uint64_t start_base,
                     uint64_t max_bases, uint64_t max_modmul_ops)
{
    FactorResult res = {0};
    unsigned char *sieve;
    uint64_t i;

    res.code = RESULT_NO_BASE;

    if (n < 4 || B > PM1_MAX_BOUND) {
        res.code = RESULT_INVALID;
        return res;
    }
    if ((n & 1) == 0) {
        res.code = RESULT_FOUND;
        res.factor = 2;
        res.cofactor = n / 2;
        return res;
    }

    sieve = make_sieve(B);
    if (sieve == NULL) {
        res.code = RESULT_MEMORY;
        return res;
    }

    for (i = 0; i < max_bases; i++) {
        uint64_t a = (start_base + i) % n;
        uint64_t x;
        uint64_t p;

        res.iterations++;

        if (a == 0)
            continue; /* the base is a multiple of n: unusable */

        /* The base may already share a factor with n. */
        x = gcd(a, n);
        res.gcd_calls++;
        if (x > 1) {
            res.code = RESULT_FOUND;
            res.factor = x;
            res.cofactor = n / x;
            free(sieve);
            return res;
        }

        /* x = a^M mod n, one prime power at a time. */
        x = a;
        for (p = 2; p <= B; p++) {
            uint64_t power;
            uint64_t ops = 0;

            if (!sieve[p])
                continue;

            power = p; /* largest p^k <= B */
            while (power <= B / p)
                power *= p;

            x = pow_mod(x, power, n, &ops);
            res.modmul_ops += ops;
            if (res.modmul_ops >= max_modmul_ops) {
                res.code = RESULT_WORK_LIMIT;
                free(sieve);
                return res;
            }
        }

        /* d = gcd(a^M - 1, n), built without forming a^M - 1. */
        {
            uint64_t xm1 = (x == 0) ? n - 1 : x - 1;
            uint64_t d = gcd(xm1, n);
            res.gcd_calls++;
            if (d > 1 && d < n) {
                res.code = RESULT_FOUND;
                res.factor = d;
                res.cofactor = n / d;
                free(sieve);
                return res;
            }
            /* d == 1: this base missed. d == n: both factors are smooth
             * minus one (a degeneracy), so the base is no use either. */
        }
    }

    free(sieve);
    return res;
}
