/*
 * tests/test_pollard_rho.c — Pollard's Rho: small/large semiprimes,
 * determinism per seed, gcd edge cases, budgets, and counter invariants
 * (gcd_calls == iterations, modmul_ops == 3 * iterations).
 */
#include "factorium/pollard_rho.h"
#include "factorium/modular_arithmetic.h"
#include "test_harness.h"

/* Asserts the shared invariants of any successful rho result. */
static void check_valid_split(const FactorResult *r, uint64_t n)
{
    CHECK(r->code == RESULT_FOUND);
    if (r->code != RESULT_FOUND)
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
    FactorResult r;
    for (uint64_t n = 0; n <= 3; n++) {
        r = rho(n, RHO_DEFAULT_ITERATIONS, RHO_DEFAULT_RETRIES,
                RHO_DEFAULT_SEED);
        CHECK(r.code == RESULT_INVALID);
    }
    r = rho(4, RHO_DEFAULT_ITERATIONS, RHO_DEFAULT_RETRIES, RHO_DEFAULT_SEED);
    CHECK(r.code == RESULT_FOUND && r.factor == 2 && r.cofactor == 2);
    r = rho(18446744073709551614ULL, RHO_DEFAULT_ITERATIONS,
            RHO_DEFAULT_RETRIES, RHO_DEFAULT_SEED); /* 2^64 - 2 */
    CHECK(r.code == RESULT_FOUND);
    CHECK_EQ_U64(r.factor, 2);
    CHECK_EQ_U64(r.cofactor, 9223372036854775807ULL);
    CHECK_EQ_U64(r.iterations, 0);
    CHECK_EQ_U64(r.gcd_calls, 0);
    CHECK_EQ_U64(r.modmul_ops, 0);
}

static void test_small_composites(void)
{
    FactorResult r;
    /* n = 15: attempt 0 uses x0 = 2, c = 0; found on the first step. */
    r = rho(15, RHO_DEFAULT_ITERATIONS, RHO_DEFAULT_RETRIES, RHO_DEFAULT_SEED);
    check_valid_split(&r, 15);
    CHECK(r.factor == 3 || r.factor == 5);
    CHECK_EQ_U64(r.seed, 2);
    CHECK_EQ_U64(r.c, 0);
    CHECK_EQ_U64(r.iterations, 1);
    CHECK_EQ_U64(r.gcd_calls, 1);
    CHECK_EQ_U64(r.modmul_ops, 3);

    /* n = 9: x0 = 2, c = 0 -> gcd(3, 9) = 3 (gcd == n never happens). */
    r = rho(9, RHO_DEFAULT_ITERATIONS, RHO_DEFAULT_RETRIES, RHO_DEFAULT_SEED);
    check_valid_split(&r, 9);
    CHECK_EQ_U64(r.factor, 3);
    CHECK_EQ_U64(r.cofactor, 3);
    CHECK_EQ_U64(r.iterations, 1);

    /* Odd sweep: for every result either a valid split is reported or an
     * honest NO_FACTOR (primes land in the latter bucket). */
    for (uint64_t n = 21; n <= 1001; n += 2) {
        r = rho(n, RHO_DEFAULT_ITERATIONS, RHO_DEFAULT_RETRIES,
                RHO_DEFAULT_SEED);
        if (r.code == RESULT_FOUND) {
            check_valid_split(&r, n);
        } else {
            CHECK(r.code == RESULT_NO_FACTOR);
        }
        if (fh_failures > 0)
            break;
    }
}

static void test_large_semiprimes(void)
{
    FactorResult r;
    /* Close factors: 1000000007 * 1000000009. */
    r = rho(1000000016000000063ULL, RHO_DEFAULT_ITERATIONS,
            RHO_DEFAULT_RETRIES, RHO_DEFAULT_SEED);
    check_valid_split(&r, 1000000016000000063ULL);
    CHECK(r.factor == 1000000007ULL || r.factor == 1000000009ULL);
    /* Distant factors: largest 32-bit prime * 4000000007. */
    uint64_t n64;
    CHECK(mul_checked(4294967291ULL, 4000000007ULL, &n64));
    r = rho(n64, RHO_DEFAULT_ITERATIONS, RHO_DEFAULT_RETRIES,
            RHO_DEFAULT_SEED);
    check_valid_split(&r, n64);
    CHECK(r.factor == 4294967291ULL || r.factor == 4000000007ULL);
    /* Small prime times a 30-bit prime. */
    r = rho(3ULL * 998244353ULL, RHO_DEFAULT_ITERATIONS, RHO_DEFAULT_RETRIES,
            RHO_DEFAULT_SEED);
    check_valid_split(&r, 3ULL * 998244353ULL);
    CHECK(r.factor == 3 || r.factor == 998244353ULL);
}

static void test_determinism(void)
{
    FactorResult a = rho(1000000016000000063ULL, RHO_DEFAULT_ITERATIONS,
                         RHO_DEFAULT_RETRIES, 42);
    FactorResult b = rho(1000000016000000063ULL, RHO_DEFAULT_ITERATIONS,
                         RHO_DEFAULT_RETRIES, 42);
    check_valid_split(&a, 1000000016000000063ULL);
    CHECK_EQ_U64(a.factor, b.factor);
    CHECK_EQ_U64(a.iterations, b.iterations);
    CHECK_EQ_U64(a.c, b.c);
    CHECK_EQ_U64(a.seed, b.seed);

    /* A different seed must still produce a valid split. */
    FactorResult d = rho(1000000016000000063ULL, RHO_DEFAULT_ITERATIONS,
                         RHO_DEFAULT_RETRIES, 7);
    check_valid_split(&d, 1000000016000000063ULL);
    /* c values are bounded by the number of retries (c = attempt). */
    CHECK(a.c < RHO_DEFAULT_RETRIES);
}

static void test_budgets_and_prime(void)
{
    FactorResult r = rho(9999999967ULL, 10, 3, RHO_DEFAULT_SEED); /* prime */
    CHECK(r.code == RESULT_NO_FACTOR);
    /* 10 iterations x 3 retries. */
    CHECK_EQ_U64(r.iterations, 30);
    CHECK_EQ_U64(r.gcd_calls, 30);
    CHECK_EQ_U64(r.modmul_ops, 90);

    r = rho(9999999967ULL, RHO_DEFAULT_ITERATIONS, 1, RHO_DEFAULT_SEED);
    CHECK(r.code == RESULT_NO_FACTOR);
    CHECK(r.iterations <= RHO_DEFAULT_ITERATIONS);

    /* Huge prime with a single tiny attempt: still never "prime". */
    r = rho(18446744073709551557ULL, 1, 1, RHO_DEFAULT_SEED);
    CHECK(r.code == RESULT_NO_FACTOR);
    CHECK_EQ_U64(r.iterations, 1);
}

static void test_prime_square_invariant(void)
{
    /* n = p^2: rho can surface either p or (rarely) collapse to n and
     * exhaust retries. Both outcomes are legal; assert only invariants. */
    uint64_t n = 4000000007ULL * 4000000007ULL; /* fits in uint64 */
    CHECK(n / 4000000007ULL == 4000000007ULL);
    FactorResult r = rho(n, RHO_DEFAULT_ITERATIONS, RHO_DEFAULT_RETRIES,
                         RHO_DEFAULT_SEED);
    if (r.code == RESULT_FOUND) {
        check_valid_split(&r, n);
        CHECK_EQ_U64(r.factor, 4000000007ULL);
    } else {
        CHECK(r.code == RESULT_NO_FACTOR);
    }
}

int main(void)
{
    test_invalid_and_even();
    test_small_composites();
    test_large_semiprimes();
    test_determinism();
    test_budgets_and_prime();
    test_prime_square_invariant();
    TEST_MAIN_END();
}
