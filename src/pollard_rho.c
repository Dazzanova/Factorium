/*
 * factorium/src/pollard_rho.c — Pollard's Rho with Floyd cycle detection.
 *
 * Per attempt: x0 = (seed + attempt) % (n - 1) + 1, c = attempt % n.
 * Both make the run bit-reproducible for a given seed and config.
 */
#include "factorium/pollard_rho.h"
#include "factorium/modular_arithmetic.h"

fct_rho_config fct_rho_config_default(void)
{
    fct_rho_config c;
    c.max_iterations_per_attempt = FCT_RHO_DEFAULT_MAX_ITERATIONS;
    c.max_retries = FCT_RHO_DEFAULT_MAX_RETRIES;
    c.seed = FCT_RHO_DEFAULT_SEED;
    return c;
}

fct_factor_result fct_pollard_rho(uint64_t n, const fct_rho_config *cfg)
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
    res.reason = FCT_TERM_BUDGET_EXHAUSTED;

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

    uint64_t iters = (cfg && cfg->max_iterations_per_attempt != 0)
                         ? cfg->max_iterations_per_attempt
                         : FCT_RHO_DEFAULT_MAX_ITERATIONS;
    uint64_t retries = (cfg && cfg->max_retries != 0)
                           ? cfg->max_retries
                           : FCT_RHO_DEFAULT_MAX_RETRIES;
    uint64_t seed = cfg ? cfg->seed : FCT_RHO_DEFAULT_SEED;

    for (uint64_t attempt = 0; attempt < retries; attempt++) {
        uint64_t x0 = (seed + attempt) % (n - 1) + 1; /* in [1, n-1] */
        uint64_t c = attempt % n;                     /* in [0, n-1] */
        uint64_t x = x0;
        uint64_t y = x0;

        for (uint64_t step = 0; step < iters; step++) {
            /* Floyd: one step for x, two for y (f(t) = t^2 + c mod n). */
            x = fct_mod_add(fct_mod_mul(x, x, n), c, n);
            y = fct_mod_add(fct_mod_mul(y, y, n), c, n);
            y = fct_mod_add(fct_mod_mul(y, y, n), c, n);

            uint64_t diff = (x > y) ? (x - y) : (y - x);
            uint64_t g = fct_gcd(diff, n);
            res.iterations++;
            res.gcd_calls++;
            res.modmul_ops += 3;

            if (g == n) {
                /* Cycle collapsed modulo both factors (includes diff == 0):
                 * this attempt cannot produce a nontrivial split. */
                break;
            }
            if (g > 1) {
                res.status = FCT_OK;
                res.reason = FCT_TERM_FOUND;
                res.factor = g;
                res.cofactor = n / g;
                res.seed = x0; /* effective initial iterate of this attempt */
                res.c = c;
                return res;
            }
        }
    }
    return res;
}
