/*
 * factorium/src/factorization.c — recursive complete factorization.
 */
#include "factorium/factorization.h"
#include "factorium/fermat.h"
#include "factorium/modular_arithmetic.h"
#include "factorium/pollard_p_minus_one.h"
#include "factorium/pollard_rho.h"
#include "factorium/primality.h"

#include <stdlib.h>
#include <string.h>

static const uint64_t SMALL_PRIMES[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
    41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97
};

typedef struct {
    FactorizeOptions opt;
    Factorization f;
    ResultCode stop;      /* why an incomplete result stopped     */
    bool factor_overflow; /* the factor array ran out of room    */
    bool internal_error;  /* an algorithm returned something odd */
} Work;

FactorizeOptions default_options(void)
{
    FactorizeOptions o;
    o.max_iterations = FACTORIZE_DEFAULT_ITERATIONS;
    o.max_retries = FACTORIZE_DEFAULT_RETRIES;
    o.seed = FACTORIZE_DEFAULT_SEED;
    o.smoothness_bound = FACTORIZE_DEFAULT_BOUND;
    o.max_bases = FACTORIZE_DEFAULT_BASES;
    o.max_depth = FACTORIZE_DEFAULT_DEPTH;
    o.max_factors = FACTORIZE_DEFAULT_MAX_FACTORS;
    return o;
}

static bool add_factor(Work *w, uint64_t p)
{
    if (w->f.factor_count >= (int)w->opt.max_factors ||
        w->f.factor_count >= MAX_FACTORS) {
        w->factor_overflow = true;
        w->stop = RESULT_FACTOR_LIMIT;
        return false;
    }
    w->f.factors[w->f.factor_count++] = p;
    return true;
}

static bool add_unresolved(Work *w, uint64_t value)
{
    if (w->f.unresolved_count >= MAX_FACTORS) {
        w->factor_overflow = true;
        w->stop = RESULT_FACTOR_LIMIT;
        return false;
    }
    w->f.unresolved[w->f.unresolved_count++] = value;
    return true;
}

static uint64_t trial_division(uint64_t n)
{
    size_t i;
    for (i = 0; i < sizeof(SMALL_PRIMES) / sizeof(SMALL_PRIMES[0]); i++) {
        if (n % SMALL_PRIMES[i] == 0)
            return SMALL_PRIMES[i];
    }
    return 0;
}

static void split(Work *w, uint64_t n, int depth)
{
    uint64_t d;
    FactorResult r;

    if (w->factor_overflow || w->internal_error)
        return;

    if (is_prime(n)) {
        add_factor(w, n);
        return;
    }

    if (depth > (int)w->opt.max_depth) {
        add_unresolved(w, n);
        w->stop = RESULT_DEPTH_LIMIT;
        return;
    }

    d = trial_division(n);

    if (d == 0) {
        r = fermat(n, w->opt.max_iterations);
        w->f.iterations += r.iterations;
        w->f.gcd_calls += r.gcd_calls;
        w->f.modmul_ops += r.modmul_ops;
        if (r.code == RESULT_FOUND)
            d = r.factor;
    }

    if (d == 0) {
        r = rho(n, w->opt.max_iterations, w->opt.max_retries, w->opt.seed);
        w->f.iterations += r.iterations;
        w->f.gcd_calls += r.gcd_calls;
        w->f.modmul_ops += r.modmul_ops;
        if (r.code == RESULT_FOUND)
            d = r.factor;
    }

    if (d == 0) {
        r = pminus1(n, w->opt.smoothness_bound, PM1_DEFAULT_START_BASE,
                    w->opt.max_bases, PM1_DEFAULT_MAX_OPS);
        w->f.iterations += r.iterations;
        w->f.gcd_calls += r.gcd_calls;
        w->f.modmul_ops += r.modmul_ops;
        if (r.code == RESULT_FOUND)
            d = r.factor;
    }

    if (d == 0) {
        add_unresolved(w, n);
        w->stop = RESULT_INCOMPLETE;
        return;
    }

    /* An algorithm must never hand back a non-split. */
    if (d <= 1 || d >= n || n % d != 0) {
        w->internal_error = true;
        return;
    }

    split(w, d, depth + 1);
    split(w, n / d, depth + 1);
}

static int compare_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a;
    uint64_t y = *(const uint64_t *)b;
    if (x < y)
        return -1;
    if (x > y)
        return 1;
    return 0;
}

Factorization factorize(uint64_t n, const FactorizeOptions *opt_in)
{
    Work w;
    uint64_t product = 1;
    bool product_ok = true;
    int i;

    memset(&w, 0, sizeof(w));
    w.f.n = n;
    w.stop = RESULT_INCOMPLETE;

    if (n < 2) {
        w.f.code = RESULT_INVALID;
        return w.f;
    }

    w.opt = default_options();
    if (opt_in != NULL) {
        if (opt_in->max_iterations != 0) w.opt.max_iterations = opt_in->max_iterations;
        if (opt_in->max_retries != 0) w.opt.max_retries = opt_in->max_retries;
        if (opt_in->seed != 0) w.opt.seed = opt_in->seed;
        if (opt_in->smoothness_bound != 0) w.opt.smoothness_bound = opt_in->smoothness_bound;
        if (opt_in->max_bases != 0) w.opt.max_bases = opt_in->max_bases;
        if (opt_in->max_depth != 0) w.opt.max_depth = opt_in->max_depth;
        if (opt_in->max_factors != 0) w.opt.max_factors = opt_in->max_factors;
    }
    if (w.opt.max_factors > MAX_FACTORS)
        w.opt.max_factors = MAX_FACTORS;

    split(&w, n, 0);

    if (w.f.factor_count > 1)
        qsort(w.f.factors, (size_t)w.f.factor_count, sizeof(uint64_t), compare_u64);

    for (i = 0; i < w.f.factor_count; i++) {
        if (!mul_checked(product, w.f.factors[i], &product)) {
            product_ok = false;
            break;
        }
    }

    w.f.product = product;
    w.f.product_verified = product_ok && product == n;

    if (w.internal_error) {
        w.f.code = RESULT_INTERNAL;
    } else if (w.factor_overflow) {
        w.f.code = RESULT_FACTOR_LIMIT;
    } else if (w.f.unresolved_count > 0) {
        w.f.code = w.stop;
    } else if (w.f.product_verified) {
        w.f.code = RESULT_COMPLETE;
    } else {
        w.f.code = RESULT_BAD_PRODUCT;
    }
    return w.f;
}
