/*
 * factorium/src/fermat.c — Fermat's difference-of-squares factorization.
 *
 * On __int128-capable compilers (documented assumption) the running
 * b^2 is kept in unsigned __int128 so the search is correct for the full
 * 64-bit domain (already at a = 2^32 the value a^2 = 2^64 overflows
 * uint64_t). Without the extension the loop is confined to a < 2^32 and
 * stops with FCT_WORK_LIMIT.
 */
#include "factorium/fermat.h"
#include "factorium/modular_arithmetic.h"

fct_fermat_config fct_fermat_config_default(void)
{
    fct_fermat_config c;
    c.max_iterations = FCT_FERMAT_DEFAULT_MAX_ITERATIONS;
    return c;
}

#if defined(__SIZEOF_INT128__)
/* The u128 machinery is a documented compiler extension; keep -Wpedantic
 * clean around everything that mentions it. */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"

typedef unsigned __int128 fct_u128;

/* floor(sqrt(n)) for u128: same integer Newton as fct_isqrt; every
 * formed sum x + n/x is < 2^128 because x <= n/2 + 1 and, while
 * x > sqrt(n), n/x < x. */
static fct_u128 isqrt128(fct_u128 n)
{
    if (n < 2)
        return n;
    fct_u128 x = n / 2 + 1;
    for (;;) {
        fct_u128 y = (x + n / x) / 2;
        if (y >= x)
            break;
        x = y;
    }
    while (x > n / x)
        x--;
    return x;
}

/* a0 must already be ceil(sqrt n); budget >= 1. */
static fct_factor_result fermat_native(uint64_t n, uint64_t budget,
                                       uint64_t a0)
{
    fct_factor_result res;
    res.status = FCT_NO_FACTOR;
    res.reason = FCT_TERM_BUDGET_EXHAUSTED;
    res.iterations = 0;
    res.gcd_calls = 0;
    res.modmul_ops = 0;
    res.seed = 0;
    res.c = 0;

    uint64_t a = a0;
    fct_u128 b2 = (fct_u128)a * a - n;
    for (uint64_t tried = 0; tried < budget; tried++) {
        fct_u128 b = isqrt128(b2);
        if (b * b == b2) {
            res.iterations = tried + 1;
            if (b < a) { /* b < a always holds for n >= 1; defensive */
                uint64_t p = a - (uint64_t)b;
                if (p > 1) {
                    res.status = FCT_OK;
                    res.reason = FCT_TERM_FOUND;
                    res.factor = p;
                    res.cofactor = n / p;
                    return res;
                }
            }
            /* p == 1 is the trivial representation a^2 - b^2 = 1*n,
             * reached only for prime n; by contract Fermat reports no
             * nontrivial factor rather than claiming primality. */
            break;
        }
        b2 += (fct_u128)2 * a + 1; /* (a+1)^2 - n */
        a++;
    }
    res.iterations = budget;
    return res;
}

#pragma GCC diagnostic pop
#endif

fct_factor_result fct_fermat(uint64_t n, const fct_fermat_config *cfg)
{
    fct_factor_result res;
    res.status = FCT_NO_FACTOR;
    res.reason = FCT_TERM_BUDGET_EXHAUSTED;
    res.iterations = 0;
    res.gcd_calls = 0;
    res.modmul_ops = 0;
    res.seed = 0;
    res.c = 0;

    if (n < 4) {
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

    uint64_t budget = cfg ? cfg->max_iterations
                          : FCT_FERMAT_DEFAULT_MAX_ITERATIONS;
    if (budget == 0)
        return res;

    uint64_t a = fct_isqrt(n);
    if (a * a < n) /* safe: isqrt yields a <= 2^32-1 */
        a++;       /* a = ceil(sqrt n); hence a^2 >= n */

#if defined(__SIZEOF_INT128__)
    return fermat_native(n, budget, a);
#else /* no __int128: confined to a < 2^32 so a^2 fits in uint64_t */
    for (uint64_t tried = 0; tried < budget; tried++) {
        res.iterations = tried + 1;
        if (a >= (1ULL << 32)) {
            res.status = FCT_WORK_LIMIT;
            res.reason = FCT_TERM_WORK_LIMIT;
            return res;
        }
        uint64_t b2 = a * a - n; /* a < 2^32 => a^2 < 2^64, and a^2 >= n */
        uint64_t b;
        if (fct_is_square(b2, &b) && b < a) {
            uint64_t p = a - b;
            if (p > 1) {
                res.status = FCT_OK;
                res.reason = FCT_TERM_FOUND;
                res.factor = p;
                res.cofactor = n / p;
                return res;
            }
            break;
        }
        a++;
    }
    return res;
#endif
}
