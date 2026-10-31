/*
 * factorium/src/pollard_rho.c — Pollard's Rho, Floyd's version.
 */
#include "factorium/pollard_rho.h"
#include "factorium/modular_arithmetic.h"

FactorResult rho(uint64_t n, uint64_t max_iterations,
                 uint64_t max_retries, uint64_t seed)
{
    FactorResult res = {0};
    uint64_t attempt;

    res.code = RESULT_NO_FACTOR;

    if (n < 4) {
        res.code = RESULT_INVALID;
        return res;
    }
    if ((n & 1) == 0) {
        res.code = RESULT_FOUND;
        res.factor = 2;
        res.cofactor = n / 2;
        return res;
    }

    for (attempt = 0; attempt < max_retries; attempt++) {
        uint64_t x0 = (seed + attempt) % (n - 1) + 1; /* 1 .. n-1 */
        uint64_t c = attempt % n;
        uint64_t x = x0;
        uint64_t y = x0;
        uint64_t step;

        for (step = 0; step < max_iterations; step++) {
            uint64_t diff;
            uint64_t g;

            /* Tortoise takes one step, hare takes two. */
            x = add_mod(mul_mod(x, x, n), c, n);
            y = add_mod(mul_mod(y, y, n), c, n);
            y = add_mod(mul_mod(y, y, n), c, n);

            diff = (x > y) ? x - y : y - x;
            g = gcd(diff, n);
            res.iterations++;
            res.gcd_calls++;
            res.modmul_ops += 3;

            if (g == n) {
                /* The walk closed on itself: this attempt is useless. */
                break;
            }
            if (g > 1) {
                res.code = RESULT_FOUND;
                res.factor = g;
                res.cofactor = n / g;
                res.seed = x0;
                res.c = c;
                return res;
            }
        }
    }
    return res;
}
