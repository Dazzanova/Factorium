/*
 * tests/test_pollard_p_minus_one.c — stage 1 behavior: favorable smooth
 * p-1, unfavorable (all bases exhausted), d == n degeneration, prime and
 * even inputs, work limits, invalid inputs.
 *
 * Every expectation below was cross-checked at development time against
 * an independent Python stage-1 computation of M = prod p^k <= B.
 */
#include "factorium/pollard_p_minus_one.h"
#include "test_harness.h"

static void test_invalid_and_even(void)
{
    FactorResult r;
    for (uint64_t n = 0; n <= 3; n++) {
        r = pminus1(n, PM1_DEFAULT_BOUND, PM1_DEFAULT_START_BASE,
                    PM1_DEFAULT_BASES, PM1_DEFAULT_MAX_OPS);
        CHECK(r.code == RESULT_INVALID);
    }
    r = pminus1(15, PM1_MAX_BOUND + 1, PM1_DEFAULT_START_BASE,
                PM1_DEFAULT_BASES, PM1_DEFAULT_MAX_OPS);
    CHECK(r.code == RESULT_INVALID);

    /* Even inputs split immediately, before any sieve or exponent work. */
    r = pminus1(4, PM1_DEFAULT_BOUND, PM1_DEFAULT_START_BASE,
                PM1_DEFAULT_BASES, PM1_DEFAULT_MAX_OPS);
    CHECK(r.code == RESULT_FOUND && r.factor == 2 && r.cofactor == 2);
    CHECK_EQ_U64(r.iterations, 0);
    CHECK_EQ_U64(r.modmul_ops, 0);
    CHECK_EQ_U64(r.gcd_calls, 0);
}

static void test_favorable_smooth_p_minus_one(void)
{
    /* 29 - 1 = 2^2 * 7 <= 100, so base 2 finds 29 at once even though
     * 1000000007 - 1 = 10^9 + 6 is nowhere near 100-smooth. */
    FactorResult r = pminus1(29ULL * 1000000007ULL, 100,
                             PM1_DEFAULT_START_BASE, PM1_DEFAULT_BASES,
                             PM1_DEFAULT_MAX_OPS);
    CHECK(r.code == RESULT_FOUND);
    CHECK_EQ_U64(r.factor, 29);
    CHECK_EQ_U64(r.cofactor, 1000000007ULL);
    CHECK_EQ_U64(r.iterations, 1);   /* one base consumed */
    CHECK(r.modmul_ops > 0);
    CHECK_EQ_U64(r.gcd_calls, 2);    /* gcd(a, n) pre-check + gcd(x-1, n) */

    /* 3 - 1 = 2: tiny bound suffices. */
    r = pminus1(3ULL * 998244353ULL, 100, PM1_DEFAULT_START_BASE,
                PM1_DEFAULT_BASES, PM1_DEFAULT_MAX_OPS);
    CHECK(r.code == RESULT_FOUND);
    CHECK_EQ_U64(r.factor, 3);
    CHECK_EQ_U64(r.cofactor, 998244353ULL);
}

static void test_even_and_powers_of_two(void)
{
    FactorResult r = pminus1(1ULL << 40, PM1_DEFAULT_BOUND,
                             PM1_DEFAULT_START_BASE, PM1_DEFAULT_BASES,
                             PM1_DEFAULT_MAX_OPS);
    CHECK(r.code == RESULT_FOUND);
    CHECK_EQ_U64(r.factor, 2);
    CHECK_EQ_U64(r.cofactor, 1ULL << 39);
    CHECK_EQ_U64(r.modmul_ops, 0);   /* even input: returned pre-search */
    CHECK_EQ_U64(r.gcd_calls, 0);
    /* 3^10 with a 100-bound: 2 | M surfaces gcd(2^M-1, n) = 3^5.
     * (With the default B the run differs: 3^10 | M makes base 2
     * degenerate, and base 3 is caught by the gcd(a, n) pre-check.) */
    r = pminus1(59049, 100, PM1_DEFAULT_START_BASE, PM1_DEFAULT_BASES,
                PM1_DEFAULT_MAX_OPS);
    CHECK(r.code == RESULT_FOUND);
    CHECK_EQ_U64(r.factor, 243);     /* gcd = 3^5, still a valid split */
    CHECK_EQ_U64(r.cofactor, 243);
}

static void test_unfavorable_smoothness(void)
{
    /* 1000003 - 1 = 2 * 3 * 166667 and 1000033 - 1 has a 359-part; with
     * B = 100 neither factor is smooth, and none of bases 2..6 collides:
     * verified as d = 1 for every base. */
    FactorResult r = pminus1(1000003ULL * 1000033ULL, 100,
                             PM1_DEFAULT_START_BASE, PM1_DEFAULT_BASES,
                             PM1_DEFAULT_MAX_OPS);
    CHECK(r.code == RESULT_NO_BASE);
    CHECK_EQ_U64(r.iterations, 5);   /* default 5 bases, all tried */
    CHECK(r.factor == 0 && r.cofactor == 0);

    /* Power of a large prime: (4000000007)^2 with B = 100 — p-1 is not
     * 100-smooth, so no base surfaces the factor either. */
    r = pminus1(4000000007ULL * 4000000007ULL, 100, PM1_DEFAULT_START_BASE,
                PM1_DEFAULT_BASES, PM1_DEFAULT_MAX_OPS);
    CHECK(r.code == RESULT_NO_BASE);
    CHECK_EQ_U64(r.iterations, 5);
}

static void test_degenerate_d_equals_n(void)
{
    /* n = 11 * 31 = 341, B = 10: 10 | M and 30 | M, so every coprime
     * base a has a^M = 1 mod 341 -> gcd(x-1, n) = n. Bases are
     * abandoned (not misreported as factors). */
    FactorResult r = pminus1(341, 10, PM1_DEFAULT_START_BASE,
                             PM1_DEFAULT_BASES, PM1_DEFAULT_MAX_OPS);
    CHECK(r.code == RESULT_NO_BASE);
    CHECK_EQ_U64(r.iterations, 5);
    CHECK(r.factor == 0);
}

static void test_prime_inputs(void)
{
    FactorResult r = pminus1(9999999967ULL, 1000, PM1_DEFAULT_START_BASE,
                             PM1_DEFAULT_BASES, PM1_DEFAULT_MAX_OPS);
    CHECK(r.code == RESULT_NO_BASE);   /* never a primality claim */
    CHECK_EQ_U64(r.iterations, 5);

    /* Largest prime below 2^64 with default bounds. */
    r = pminus1(18446744073709551557ULL, PM1_DEFAULT_BOUND,
                PM1_DEFAULT_START_BASE, PM1_DEFAULT_BASES,
                PM1_DEFAULT_MAX_OPS);
    CHECK(r.code == RESULT_NO_BASE);
    CHECK_EQ_U64(r.iterations, PM1_DEFAULT_BASES);
}

static void test_work_limit_and_config(void)
{
    /* Cap smaller than one prime-power pass. */
    FactorResult r = pminus1(29ULL * 1000000007ULL, 100000,
                             PM1_DEFAULT_START_BASE, PM1_DEFAULT_BASES, 10);
    CHECK(r.code == RESULT_WORK_LIMIT);
    CHECK_EQ_U64(r.iterations, 1);  /* base started, not completed */
    CHECK(r.modmul_ops >= 10);
    CHECK(r.factor == 0);

    /* Configurable start base: base 3 finds a factor base 2 also would. */
    r = pminus1(29ULL * 1000000007ULL, 100, 3, PM1_DEFAULT_BASES,
                PM1_DEFAULT_MAX_OPS);
    CHECK(r.code == RESULT_FOUND);
    CHECK_EQ_U64(r.factor, 29);

    /* max_bases = 1: only the first base is spent. */
    r = pminus1(1000003ULL * 1000033ULL, 100, PM1_DEFAULT_START_BASE, 1,
                PM1_DEFAULT_MAX_OPS);
    CHECK(r.code == RESULT_NO_BASE);
    CHECK_EQ_U64(r.iterations, 1);
}

int main(void)
{
    test_invalid_and_even();
    test_favorable_smooth_p_minus_one();
    test_even_and_powers_of_two();
    test_unfavorable_smoothness();
    test_degenerate_d_equals_n();
    test_prime_inputs();
    test_work_limit_and_config();
    TEST_MAIN_END();
}
