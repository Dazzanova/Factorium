/*
 * tests/test_factorization.c — complete factorization layer: multiplicities,
 * sorting, reconstruction via checked arithmetic, known 64-bit
 * factorizations, and explicit incomplete/limit outcomes.
 */
#include "factorium/factorization.h"
#include "factorium/modular_arithmetic.h"
#include "factorium/primality.h"
#include "test_harness.h"

/* Shared invariant: sorted, all prime, product reconstructs n. */
static void check_complete(const Factorization *f, uint64_t n)
{
    uint64_t prod = 1;

    CHECK(f->code == RESULT_COMPLETE);
    CHECK(f->n == n);
    if (f->code != RESULT_COMPLETE)
        return;
    CHECK(f->product_verified);
    CHECK_EQ_U64(f->product, n);
    CHECK(f->unresolved_count == 0);
    CHECK(f->factor_count > 0);
    for (int i = 0; i < f->factor_count; i++) {
        CHECK(is_prime(f->factors[i]));
        if (i > 0)
            CHECK(f->factors[i - 1] <= f->factors[i]); /* sorted */
        if (!mul_checked(prod, f->factors[i], &prod)) {
            CHECK(false); /* product overflowed a uint64 */
            break;
        }
    }
    CHECK_EQ_U64(prod, n);
}

static void check_sequence(const Factorization *f,
                           const uint64_t *want, int count)
{
    CHECK_EQ_U64((uint64_t)f->factor_count, (uint64_t)count);
    if (f->factor_count != count)
        return;
    for (int i = 0; i < count; i++)
        CHECK_EQ_U64(f->factors[i], want[i]);
}

static void test_small_composites(void)
{
    Factorization f;
    f = factorize(360, NULL);
    {
        const uint64_t want[] = {2, 2, 2, 3, 3, 5};
        check_sequence(&f, want, 6);
    }
    check_complete(&f, 360);

    f = factorize(6300, NULL); /* 2^2 * 3^2 * 5^2 * 7 */
    {
        const uint64_t want[] = {2, 2, 3, 3, 5, 5, 7};
        check_sequence(&f, want, 7);
    }
    check_complete(&f, 6300);

    f = factorize(59049, NULL); /* 3^10 */
    CHECK(f.code == RESULT_COMPLETE);
    CHECK_EQ_U64(f.factor_count, 10);
    for (int i = 0; i < f.factor_count; i++)
        CHECK_EQ_U64(f.factors[i], 3);
    check_complete(&f, 59049);

    f = factorize(97, NULL); /* prime */
    CHECK(f.code == RESULT_COMPLETE);
    check_sequence(&f, &(const uint64_t){97}, 1);
}

static void test_powers_of_two(void)
{
    Factorization f = factorize(1ULL << 40, NULL);
    CHECK(f.code == RESULT_COMPLETE);
    CHECK_EQ_U64(f.factor_count, 40);
    for (int i = 0; i < f.factor_count; i++)
        CHECK_EQ_U64(f.factors[i], 2);
    check_complete(&f, 1ULL << 40);
}

static void test_primes(void)
{
    Factorization f;
    f = factorize(9999999967ULL, NULL);
    check_sequence(&f, &(const uint64_t){9999999967ULL}, 1);
    check_complete(&f, 9999999967ULL);

    f = factorize(18446744073709551557ULL, NULL); /* largest prime < 2^64 */
    check_sequence(&f, &(const uint64_t){18446744073709551557ULL}, 1);
    check_complete(&f, 18446744073709551557ULL);

    for (uint64_t n = 2; n <= 50; n++) {
        f = factorize(n, NULL);
        check_complete(&f, n); /* primes -> [n], composites -> split */
    }
}

static void test_known_64bit_factorizations(void)
{
    /* 2^64 - 1 = 3 * 5 * 17 * 257 * 641 * 65537 * 6700417 */
    Factorization f = factorize(18446744073709551615ULL, NULL);
    {
        const uint64_t want[] = {3, 5, 17, 257, 641, 65537, 6700417};
        check_sequence(&f, want, 7);
    }
    check_complete(&f, 18446744073709551615ULL);

    /* 2^63 - 1 = 7^2 * 73 * 127 * 337 * 92737 * 649657 */
    f = factorize(9223372036854775807ULL, NULL);
    {
        const uint64_t want[] = {7, 7, 73, 127, 337, 92737, 649657};
        check_sequence(&f, want, 7);
    }
    check_complete(&f, 9223372036854775807ULL);
}

static void test_semiprimes(void)
{
    Factorization f;

    /* Distant factors: rho must split it. */
    uint64_t n64;
    CHECK(mul_checked(4294967291ULL, 4000000007ULL, &n64));
    f = factorize(n64, NULL);
    {
        const uint64_t want[] = {4000000007ULL, 4294967291ULL};
        check_sequence(&f, want, 2);
    }
    check_complete(&f, n64);

    /* Close factors: Fermat splits it on the first iteration. */
    f = factorize(1000000016000000063ULL, NULL);
    {
        const uint64_t want[] = {1000000007ULL, 1000000009ULL};
        check_sequence(&f, want, 2);
    }
    check_complete(&f, 1000000016000000063ULL);

    /* Prime square. */
    f = factorize(4000000007ULL * 4000000007ULL, NULL);
    CHECK(f.code == RESULT_COMPLETE);
    CHECK_EQ_U64(f.factor_count, 2);
    check_complete(&f, 4000000007ULL * 4000000007ULL);
}

static void test_invalid_inputs(void)
{
    Factorization f = factorize(0, NULL);
    CHECK(f.code == RESULT_INVALID);
    CHECK_EQ_U64(f.factor_count, 0);

    f = factorize(1, NULL);
    CHECK(f.code == RESULT_INVALID);
    CHECK_EQ_U64(f.factor_count, 0);
}

static void test_incomplete_tiny_budget(void)
{
    /* Distant 64-bit semiprime with every budget shrunk to the minimum.
     * Verified behavior of each stage on this input:
     *  - trial division: no factor <= 97,
     *  - Fermat: 1 iteration misses (a0 far from the solution),
     *  - Rho: 1 step with x0 = 2, c = 0 gives x = 4, y = 16, gcd(12, n) = 1,
     *  - p-1 with B = 2, one base: M = 2, gcd(2^M - 1, n) = gcd(3, n) = 1. */
    FactorizeOptions cfg = default_options();
    Factorization f;
    uint64_t n;
    uint64_t n2;
    uint64_t nsmall = 1000003ULL * 999999937ULL;

    cfg.max_iterations = 1;
    cfg.max_retries = 1;
    cfg.smoothness_bound = 2;
    cfg.max_bases = 1;
    CHECK(mul_checked(4294967291ULL, 4000000007ULL, &n));

    f = factorize(n, &cfg);
    CHECK(f.code == RESULT_INCOMPLETE);
    CHECK_EQ_U64(f.factor_count, 0);
    CHECK_EQ_U64(f.unresolved_count, 1);
    if (f.unresolved_count == 1)
        CHECK_EQ_U64(f.unresolved[0], n);
    CHECK(!f.product_verified);

    /* With a small prime peeled off first: partial factors + remainder.
     * A smaller distant semiprime keeps 2 * n inside uint64_t; same budget
     * recipe, stage outcomes verified independently: no factor <= 97,
     * Fermat's first b2 = 60735165 is not a square (the solution sits at
     * a = 500499970), rho's first gcd(12, n) = 1, p-1's gcd(3, n) = 1. */
    CHECK(mul_checked(2, nsmall, &n2));
    f = factorize(n2, &cfg);
    CHECK(f.code == RESULT_INCOMPLETE);
    check_sequence(&f, &(const uint64_t){2}, 1);
    CHECK_EQ_U64(f.unresolved_count, 1);
    if (f.unresolved_count == 1)
        CHECK_EQ_U64(f.unresolved[0], nsmall);
    CHECK(!f.product_verified);
    CHECK_EQ_U64(f.product, 2);
}

static void test_limits(void)
{
    /* Depth cap: with max_depth = 1 the split of 90 stalls at depth 2. */
    FactorizeOptions cfg = default_options();
    Factorization f;

    cfg.max_depth = 1;
    f = factorize(360, &cfg);
    CHECK(f.code == RESULT_DEPTH_LIMIT);
    {
        const uint64_t want[] = {2, 2};
        check_sequence(&f, want, 2);
    }
    CHECK_EQ_U64(f.unresolved_count, 1);
    if (f.unresolved_count == 1)
        CHECK_EQ_U64(f.unresolved[0], 90);

    /* Factor-count cap: the third factor exceeds max_factors = 2. */
    cfg = default_options();
    cfg.max_factors = 2;
    f = factorize(360, &cfg);
    CHECK(f.code == RESULT_FACTOR_LIMIT);
    CHECK_EQ_U64(f.factor_count, 2);
}

static void test_reconstruction_overflow_guard(void)
{
    /* n = 2^64 - 1: factors must reconstruct without checked-mul overflow. */
    Factorization f = factorize(18446744073709551615ULL, NULL);
    uint64_t prod = 1;

    CHECK(f.product_verified);
    CHECK(f.code == RESULT_COMPLETE);
    for (int i = 0; i < f.factor_count; i++)
        CHECK(mul_checked(prod, f.factors[i], &prod) && prod >= f.factors[i]);
    CHECK_EQ_U64(prod, 18446744073709551615ULL);
}

int main(void)
{
    FactorizeOptions def = default_options();
    CHECK_EQ_U64(def.max_iterations, FACTORIZE_DEFAULT_ITERATIONS);
    CHECK_EQ_U64(def.max_retries, FACTORIZE_DEFAULT_RETRIES);
    CHECK_EQ_U64(def.smoothness_bound, FACTORIZE_DEFAULT_BOUND);
    CHECK_EQ_U64(def.max_bases, FACTORIZE_DEFAULT_BASES);
    CHECK_EQ_U64(def.seed, FACTORIZE_DEFAULT_SEED);
    CHECK_EQ_U64(def.max_depth, FACTORIZE_DEFAULT_DEPTH);
    CHECK_EQ_U64(def.max_factors, FACTORIZE_DEFAULT_MAX_FACTORS);
    test_small_composites();
    test_powers_of_two();
    test_primes();
    test_known_64bit_factorizations();
    test_semiprimes();
    test_invalid_inputs();
    test_incomplete_tiny_budget();
    test_limits();
    test_reconstruction_overflow_guard();
    TEST_MAIN_END();
}
