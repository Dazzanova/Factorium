/*
 * tests/test_pollard_rho.c — Pollard's Rho: small/large semiprimes,
 * determinism per seed, gcd edge cases, budgets, and counter invariants
 * (gcd_calls == iterations, modmul_ops == 3 * iterations).
 */
#include "factorium/pollard_rho.h"
#include "factorium/modular_arithmetic.h"
#include "test_harness.h"

/* Asserts the shared invariants of any successful rho result. */
static void check_valid_split(const fct_factor_result *r, uint64_t n)
{
    CHECK(r->status == FCT_OK);
    CHECK(r->reason == FCT_TERM_FOUND);
    if (r->status != FCT_OK)
        return;
    CHECK(r->factor > 1);
    CHECK(r->factor < n);
    CHECK(n % r->factor == 0);
    CHECK(r->cofactor == n / r->factor);
    CHECK(r->cofactor > 1);
    CHECK_EQ_U64(r->gcd_calls, r->iterations);
    CHECK_EQ_U64(r->modmul_ops, 3 * r->iterations);
}

static void test_invalid_and_even(void)
{
    fct_factor_result r;
    for (uint64_t n = 0; n <= 3; n++) {
        r = fct_pollard_rho(n, NULL);
        CHECK(r.status == FCT_INVALID_INPUT);
        CHECK(r.reason == FCT_TERM_INVALID_INPUT);
    }
    r = fct_pollard_rho(4, NULL);
    CHECK(r.status == FCT_OK && r.factor == 2 && r.cofactor == 2);
    r = fct_pollard_rho(18446744073709551614ULL, NULL); /* 2^64 - 2 */
    CHECK(r.status == FCT_OK);
    CHECK_EQ_U64(r.factor, 2);
    CHECK_EQ_U64(r.cofactor, 9223372036854775807ULL);
    CHECK_EQ_U64(r.iterations, 0);
    CHECK_EQ_U64(r.gcd_calls, 0);
    CHECK_EQ_U64(r.modmul_ops, 0);
}

static void test_small_composites(void)
{
    fct_factor_result r;
    /* n = 15: attempt 0 uses x0 = 2, c = 0; found on the first step. */
    r = fct_pollard_rho(15, NULL);
    check_valid_split(&r, 15);
    CHECK(r.factor == 3 || r.factor == 5);
    CHECK_EQ_U64(r.seed, 2);
    CHECK_EQ_U64(r.c, 0);
    CHECK_EQ_U64(r.iterations, 1);
    CHECK_EQ_U64(r.gcd_calls, 1);
    CHECK_EQ_U64(r.modmul_ops, 3);

    /* n = 9: x0 = 2, c = 0 -> gcd(3, 9) = 3 (gcd == n never happens). */
    r = fct_pollard_rho(9, NULL);
    check_valid_split(&r, 9);
    CHECK_EQ_U64(r.factor, 3);
    CHECK_EQ_U64(r.cofactor, 3);
    CHECK_EQ_U64(r.iterations, 1);

    /* Odd sweep: for every result either a valid split is reported or an
     * honest NO_FACTOR (primes land in the latter bucket). */
    for (uint64_t n = 21; n <= 1001; n += 2) {
        r = fct_pollard_rho(n, NULL);
        if (r.status == FCT_OK) {
            check_valid_split(&r, n);
        } else {
            CHECK(r.status == FCT_NO_FACTOR);
        }
        if (fh_failures > 0)
            break;
    }
}

static void test_large_semiprimes(void)
{
    fct_factor_result r;
    /* Close factors: 1000000007 * 1000000009. */
    r = fct_pollard_rho(1000000016000000063ULL, NULL);
    check_valid_split(&r, 1000000016000000063ULL);
    CHECK(r.factor == 1000000007ULL || r.factor == 1000000009ULL);
    /* Distant factors: largest 32-bit prime * 4000000007. */
    uint64_t n64;
    CHECK(fct_mul_checked(4294967291ULL, 4000000007ULL, &n64));
    r = fct_pollard_rho(n64, NULL);
    check_valid_split(&r, n64);
    CHECK(r.factor == 4294967291ULL || r.factor == 4000000007ULL);
    /* Small prime times a 30-bit prime. */
    r = fct_pollard_rho(3ULL * 998244353ULL, NULL);
    check_valid_split(&r, 3ULL * 998244353ULL);
    CHECK(r.factor == 3 || r.factor == 998244353ULL);
}

static void test_determinism(void)
{
    fct_rho_config c = fct_rho_config_default();
    c.seed = 42;
    fct_factor_result a = fct_pollard_rho(1000000016000000063ULL, &c);
    fct_factor_result b = fct_pollard_rho(1000000016000000063ULL, &c);
    check_valid_split(&a, 1000000016000000063ULL);
    CHECK_EQ_U64(a.factor, b.factor);
    CHECK_EQ_U64(a.iterations, b.iterations);
    CHECK_EQ_U64(a.c, b.c);
    CHECK_EQ_U64(a.seed, b.seed);

    /* A different seed must still produce a valid split. */
    fct_rho_config c2 = fct_rho_config_default();
    c2.seed = 7;
    fct_factor_result d = fct_pollard_rho(1000000016000000063ULL, &c2);
    check_valid_split(&d, 1000000016000000063ULL);
    /* c values are bounded by the number of retries (c = attempt). */
    CHECK(a.c < c.max_retries);
}

static void test_budgets_and_prime(void)
{
    fct_rho_config c = fct_rho_config_default();
    c.max_iterations_per_attempt = 10;
    c.max_retries = 3;
    fct_factor_result r = fct_pollard_rho(9999999967ULL, &c); /* prime */
    CHECK(r.status == FCT_NO_FACTOR);
    CHECK(r.reason == FCT_TERM_BUDGET_EXHAUSTED);
    /* 10 iterations x 3 retries. */
    CHECK_EQ_U64(r.iterations, 30);
    CHECK_EQ_U64(r.gcd_calls, 30);
    CHECK_EQ_U64(r.modmul_ops, 90);

    c.max_iterations_per_attempt = 0; /* 0 => default */
    c.max_retries = 1;
    r = fct_pollard_rho(9999999967ULL, &c);
    CHECK(r.status == FCT_NO_FACTOR);
    CHECK(r.iterations <= FCT_RHO_DEFAULT_MAX_ITERATIONS);

    /* Huge prime with a single tiny attempt: still never "prime". */
    fct_rho_config tiny = fct_rho_config_default();
    tiny.max_iterations_per_attempt = 1;
    tiny.max_retries = 1;
    r = fct_pollard_rho(18446744073709551557ULL, &tiny);
    CHECK(r.status == FCT_NO_FACTOR);
    CHECK_EQ_U64(r.iterations, 1);
}

static void test_prime_square_invariant(void)
{
    /* n = p^2: rho can surface either p or (rarely) collapse to n and
     * exhaust retries. Both outcomes are legal; assert only invariants. */
    uint64_t n = 4000000007ULL * 4000000007ULL; /* fits in uint64 */
    CHECK(n / 4000000007ULL == 4000000007ULL);
    fct_factor_result r = fct_pollard_rho(n, NULL);
    if (r.status == FCT_OK) {
        check_valid_split(&r, n);
        CHECK_EQ_U64(r.factor, 4000000007ULL);
    } else {
        CHECK(r.status == FCT_NO_FACTOR);
        CHECK(r.reason == FCT_TERM_BUDGET_EXHAUSTED);
    }
}

int main(void)
{
    fct_rho_config def = fct_rho_config_default();
    CHECK_EQ_U64(def.max_iterations_per_attempt, FCT_RHO_DEFAULT_MAX_ITERATIONS);
    CHECK_EQ_U64(def.max_retries, FCT_RHO_DEFAULT_MAX_RETRIES);
    CHECK_EQ_U64(def.seed, FCT_RHO_DEFAULT_SEED);
    test_invalid_and_even();
    test_small_composites();
    test_large_semiprimes();
    test_determinism();
    test_budgets_and_prime();
    test_prime_square_invariant();
    TEST_MAIN_END();
}
