/*
 * factorium/src/factorization.c — recursive splitting of composites.
 *
 * Order per node: Miller-Rabin (prime => emit), depth cap, trial division
 * by the primes <= 97, Fermat, Pollard Rho, Pollard p-1. The first
 * nontrivial factor splits the node; the two children recurse. Every
 * outcome — complete, incomplete, limit — is recorded explicitly.
 */
#include "factorium/factorization.h"
#include "factorium/modular_arithmetic.h"
#include "factorium/primality.h"
#include "factorium/fermat.h"
#include "factorium/pollard_rho.h"
#include "factorium/pollard_p_minus_one.h"

#include <stdlib.h>
#include <string.h>

static const uint64_t SMALL_PRIMES[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
    41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97
};

typedef struct {
    fct_factorization_config cfg;
    uint64_t *factors;
    size_t factor_count;
    size_t factor_cap;
    uint64_t *unresolved;
    size_t unresolved_count;
    size_t unresolved_cap;
    uint64_t iterations;
    uint64_t gcd_calls;
    uint64_t modmul_ops;
    fct_termination reason;
    bool factor_overflow;
    bool mem_fail;
    bool internal_error;
} split_ctx;

fct_factorization_config fct_factorization_config_default(void)
{
    fct_factorization_config c;
    c.max_iterations = FCT_FACTORIZATION_DEFAULT_MAX_ITERATIONS;
    c.max_retries = FCT_FACTORIZATION_DEFAULT_MAX_RETRIES;
    c.smoothness_bound = FCT_FACTORIZATION_DEFAULT_SMOOTHNESS_BOUND;
    c.max_bases = FCT_FACTORIZATION_DEFAULT_MAX_BASES;
    c.seed = FCT_FACTORIZATION_DEFAULT_SEED;
    c.max_depth = FCT_FACTORIZATION_DEFAULT_MAX_DEPTH;
    c.max_factors = FCT_FACTORIZATION_DEFAULT_MAX_FACTORS;
    return c;
}

static bool append_factor(split_ctx *c, uint64_t p)
{
    if (c->factor_count >= c->factor_cap) {
        c->factor_overflow = true;
        c->reason = FCT_TERM_FACTOR_LIMIT;
        return false;
    }
    c->factors[c->factor_count++] = p;
    return true;
}

static bool append_unresolved(split_ctx *c, uint64_t v)
{
    if (c->unresolved_count == c->unresolved_cap) {
        size_t ncap = c->unresolved_cap ? c->unresolved_cap * 2 : 8;
        uint64_t *grown = (uint64_t *)realloc(c->unresolved,
                                              ncap * sizeof(uint64_t));
        if (grown == NULL) {
            c->mem_fail = true;
            return false;
        }
        c->unresolved = grown;
        c->unresolved_cap = ncap;
    }
    c->unresolved[c->unresolved_count++] = v;
    return true;
}

static void accumulate(split_ctx *c, const fct_factor_result *r)
{
    c->iterations += r->iterations;
    c->gcd_calls += r->gcd_calls;
    c->modmul_ops += r->modmul_ops;
}

static void split(split_ctx *c, uint64_t v, uint64_t depth)
{
    if (c->factor_overflow || c->mem_fail || c->internal_error)
        return;

    if (fct_is_prime(v)) {
        append_factor(c, v);
        return;
    }

    if (depth > c->cfg.max_depth) {
        append_unresolved(c, v);
        c->reason = FCT_TERM_RECURSION_LIMIT;
        return;
    }

    /* 1. Trial division by the small primes. */
    uint64_t d = 0;
    for (size_t i = 0; i < sizeof(SMALL_PRIMES) / sizeof(SMALL_PRIMES[0]); i++) {
        if (v % SMALL_PRIMES[i] == 0) {
            d = SMALL_PRIMES[i];
            break;
        }
    }

    /* 2. Fermat — favorable when the two sides are close. */
    if (d == 0) {
        fct_fermat_config fc;
        fc.max_iterations = c->cfg.max_iterations;
        fct_factor_result r = fct_fermat(v, &fc);
        accumulate(c, &r);
        if (r.status == FCT_OK)
            d = r.factor;
    }

    /* 3. Pollard Rho — the general workhorse. */
    if (d == 0) {
        fct_rho_config rc;
        rc.max_iterations_per_attempt = c->cfg.max_iterations;
        rc.max_retries = c->cfg.max_retries;
        rc.seed = c->cfg.seed;
        fct_factor_result r = fct_pollard_rho(v, &rc);
        accumulate(c, &r);
        if (r.status == FCT_OK)
            d = r.factor;
    }

    /* 4. Pollard p-1 — favorable when p-1 / q-1 are smooth. */
    if (d == 0) {
        fct_pminus1_config pc;
        pc.smoothness_bound = c->cfg.smoothness_bound;
        pc.start_base = 0;   /* 0 => default 2 */
        pc.max_bases = c->cfg.max_bases;
        pc.max_modmul_ops = 0; /* 0 => default cap */
        fct_factor_result r = fct_pollard_p_minus_one(v, &pc);
        accumulate(c, &r);
        if (r.status == FCT_OK)
            d = r.factor;
    }

    if (d == 0) {
        append_unresolved(c, v);
        c->reason = FCT_TERM_BUDGET_EXHAUSTED;
        return;
    }

    /* Defensive: an algorithm must never hand back a non-split. */
    if (d <= 1 || d >= v || v % d != 0) {
        c->internal_error = true;
        c->reason = FCT_TERM_NONE;
        return;
    }

    split(c, d, depth + 1);
    split(c, v / d, depth + 1);
}

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a;
    uint64_t y = *(const uint64_t *)b;
    if (x < y)
        return -1;
    if (x > y)
        return 1;
    return 0;
}

fct_factorization fct_factorize(uint64_t n, const fct_factorization_config *cfg_in)
{
    fct_factorization res;
    memset(&res, 0, sizeof(res));
    res.status = FCT_INTERNAL;
    res.n = n;
    res.reason = FCT_TERM_NONE;

    if (n < 2) {
        res.status = FCT_INVALID_INPUT;
        res.reason = FCT_TERM_INVALID_INPUT;
        return res;
    }

    split_ctx c;
    memset(&c, 0, sizeof(c));
    c.cfg = fct_factorization_config_default();
    if (cfg_in) {
        if (cfg_in->max_iterations != 0) c.cfg.max_iterations = cfg_in->max_iterations;
        if (cfg_in->max_retries != 0) c.cfg.max_retries = cfg_in->max_retries;
        if (cfg_in->smoothness_bound != 0) c.cfg.smoothness_bound = cfg_in->smoothness_bound;
        if (cfg_in->max_bases != 0) c.cfg.max_bases = cfg_in->max_bases;
        if (cfg_in->seed != 0) c.cfg.seed = cfg_in->seed;
        if (cfg_in->max_depth != 0) c.cfg.max_depth = cfg_in->max_depth;
        if (cfg_in->max_factors != 0) c.cfg.max_factors = cfg_in->max_factors;
    }
    c.factor_cap = (size_t)c.cfg.max_factors;
    c.reason = FCT_TERM_NONE;
    c.factors = (uint64_t *)malloc(c.factor_cap * sizeof(uint64_t));
    if (c.factors == NULL) {
        res.status = FCT_MEMORY;
        return res;
    }

    split(&c, n, 0);

    if (c.factor_count > 1)
        qsort(c.factors, c.factor_count, sizeof(uint64_t), cmp_u64);

    uint64_t prod = 1;
    bool prod_ok = true;
    for (size_t i = 0; i < c.factor_count; i++) {
        if (!fct_mul_checked(prod, c.factors[i], &prod)) {
            prod_ok = false;
            break;
        }
    }

    res.factors = c.factors;
    res.factor_count = c.factor_count;
    res.unresolved_cofactors = c.unresolved;
    res.unresolved_count = c.unresolved_count;
    res.product = prod;
    res.product_verified = prod_ok && (prod == n);
    res.iterations = c.iterations;
    res.gcd_calls = c.gcd_calls;
    res.modmul_ops = c.modmul_ops;
    res.reason = c.reason;

    if (c.mem_fail) {
        res.status = FCT_MEMORY;
    } else if (c.internal_error) {
        res.status = FCT_INTERNAL;
    } else if (c.factor_overflow) {
        res.status = FCT_WORK_LIMIT;
    } else if (res.unresolved_count > 0) {
        res.status = FCT_INCOMPLETE;
    } else if (res.product_verified) {
        res.status = FCT_OK;
        res.reason = FCT_TERM_NONE;
    } else {
        res.status = FCT_RECONSTRUCTION_FAILED;
    }
    return res;
}

void fct_factorization_free(fct_factorization *f)
{
    if (f == NULL)
        return;
    free(f->factors);
    free(f->unresolved_cofactors);
    f->factors = NULL;
    f->unresolved_cofactors = NULL;
    f->factor_count = 0;
    f->unresolved_count = 0;
}
