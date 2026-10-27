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

static fct_pminus1_config cfg_b(uint64_t B)
{
    fct_pminus1_config c = fct_pminus1_config_default();
    c.smoothness_bound = B;
    return c;
}

static void test_invalid_and_even(void)
{
    fct_factor_result r;
    for (uint64_t n = 0; n <= 3; n++) {
        r = fct_pollard_p_minus_one(n, NULL);
        CHECK(r.status == FCT_INVALID_INPUT);
        CHECK(r.reason == FCT_TERM_INVALID_INPUT);
    }
    fct_pminus1_config huge = cfg_b(FCT_PM1_MAX_SMOOTHNESS_BOUND + 1);
    r = fct_pollard_p_minus_one(15, &huge);
    CHECK(r.status == FCT_INVALID_INPUT);

    /* Even inputs split immediately, before any sieve or exponent work. */
    r = fct_pollard_p_minus_one(4, NULL);
    CHECK(r.status == FCT_OK && r.factor == 2 && r.cofactor == 2);
    CHECK_EQ_U64(r.iterations, 0);
    CHECK_EQ_U64(r.modmul_ops, 0);
    CHECK_EQ_U64(r.gcd_calls, 0);
}

static void test_favorable_smooth_p_minus_one(void)
{
    /* 29 - 1 = 2^2 * 7 <= 100, so base 2 finds 29 at once even though
     * 1000000007 - 1 = 10^9 + 6 is nowhere near 100-smooth. */
    fct_pminus1_config c = cfg_b(100);
    fct_factor_result r = fct_pollard_p_minus_one(29ULL * 1000000007ULL, &c);
    CHECK(r.status == FCT_OK);
    CHECK(r.reason == FCT_TERM_FOUND);
    CHECK_EQ_U64(r.factor, 29);
    CHECK_EQ_U64(r.cofactor, 1000000007ULL);
    CHECK_EQ_U64(r.iterations, 1);   /* one base consumed */
    CHECK(r.modmul_ops > 0);
    CHECK_EQ_U64(r.gcd_calls, 2);    /* gcd(a, n) pre-check + gcd(x-1, n) */

    /* 3 - 1 = 2: tiny bound suffices. */
    r = fct_pollard_p_minus_one(3ULL * 998244353ULL, &c);
    CHECK(r.status == FCT_OK);
    CHECK_EQ_U64(r.factor, 3);
    CHECK_EQ_U64(r.cofactor, 998244353ULL);
}

static void test_even_and_powers_of_two(void)
{
    fct_factor_result r = fct_pollard_p_minus_one(1ULL << 40, NULL);
    CHECK(r.status == FCT_OK);
    CHECK_EQ_U64(r.factor, 2);
    CHECK_EQ_U64(r.cofactor, 1ULL << 39);
    CHECK_EQ_U64(r.modmul_ops, 0);   /* even input: returned pre-search */
    CHECK_EQ_U64(r.gcd_calls, 0);
    /* 3^10 with a 100-bound: 2 | M surfaces gcd(2^M-1, n) = 3^5.
     * (With the default B the run differs: 3^10 | M makes base 2
     * degenerate, and base 3 is caught by the gcd(a, n) pre-check.) */
    fct_pminus1_config c100 = cfg_b(100);
    r = fct_pollard_p_minus_one(59049, &c100);
    CHECK(r.status == FCT_OK);
    CHECK_EQ_U64(r.factor, 243);     /* gcd = 3^5, still a valid split */
    CHECK_EQ_U64(r.cofactor, 243);
}

static void test_unfavorable_smoothness(void)
{
    /* 1000003 - 1 = 2 * 3 * 166667 and 1000033 - 1 has a 359-part; with
     * B = 100 neither factor is smooth, and none of bases 2..6 collides:
     * verified as d = 1 for every base. */
    fct_pminus1_config c = cfg_b(100);
    fct_factor_result r =
        fct_pollard_p_minus_one(1000003ULL * 1000033ULL, &c);
    CHECK(r.status == FCT_NO_FACTOR);
    CHECK(r.reason == FCT_TERM_BASES_EXHAUSTED);
    CHECK_EQ_U64(r.iterations, 5);   /* default 5 bases, all tried */
    CHECK(r.factor == 0 && r.cofactor == 0);

    /* Power of a large prime: (4000000007)^2 with B = 100 — p-1 is not
     * 100-smooth, so no base surfaces the factor either. */
    r = fct_pollard_p_minus_one(4000000007ULL * 4000000007ULL, &c);
    CHECK(r.status == FCT_NO_FACTOR);
    CHECK(r.reason == FCT_TERM_BASES_EXHAUSTED);
    CHECK_EQ_U64(r.iterations, 5);
}

static void test_degenerate_d_equals_n(void)
{
    /* n = 11 * 31 = 341, B = 10: 10 | M and 30 | M, so every coprime
     * base a has a^M = 1 mod 341 -> gcd(x-1, n) = n. Bases are
     * abandoned (not misreported as factors). */
    fct_pminus1_config c = cfg_b(10);
    fct_factor_result r = fct_pollard_p_minus_one(341, &c);
    CHECK(r.status == FCT_NO_FACTOR);
    CHECK(r.reason == FCT_TERM_BASES_EXHAUSTED);
    CHECK_EQ_U64(r.iterations, 5);
    CHECK(r.factor == 0);
}

static void test_prime_inputs(void)
{
    fct_pminus1_config c = fct_pminus1_config_default();
    c.smoothness_bound = 1000;
    fct_factor_result r = fct_pollard_p_minus_one(9999999967ULL, &c);
    CHECK(r.status == FCT_NO_FACTOR);   /* never a primality claim */
    CHECK(r.reason == FCT_TERM_BASES_EXHAUSTED);
    CHECK_EQ_U64(r.iterations, 5);

    /* Largest prime below 2^64 with default bounds. */
    r = fct_pollard_p_minus_one(18446744073709551557ULL, NULL);
    CHECK(r.status == FCT_NO_FACTOR);
    CHECK_EQ_U64(r.iterations, FCT_PM1_DEFAULT_MAX_BASES);
}

static void test_work_limit_and_config(void)
{
    fct_pminus1_config c = fct_pminus1_config_default();
    c.smoothness_bound = 100000;
    c.max_modmul_ops = 10;          /* smaller than one prime-power pass */
    fct_factor_result r = fct_pollard_p_minus_one(29ULL * 1000000007ULL, &c);
    CHECK(r.status == FCT_WORK_LIMIT);
    CHECK(r.reason == FCT_TERM_WORK_LIMIT);
    CHECK_EQ_U64(r.iterations, 1);  /* base started, not completed */
    CHECK(r.modmul_ops >= 10);
    CHECK(r.factor == 0);

    /* Configurable start base: base 3 finds a factor base 2 also would. */
    fct_pminus1_config c3 = cfg_b(100);
    c3.start_base = 3;
    r = fct_pollard_p_minus_one(29ULL * 1000000007ULL, &c3);
    CHECK(r.status == FCT_OK);
    CHECK_EQ_U64(r.factor, 29);

    /* max_bases = 1: only the first base is spent. */
    fct_pminus1_config c1 = cfg_b(100);
    c1.max_bases = 1;
    r = fct_pollard_p_minus_one(1000003ULL * 1000033ULL, &c1);
    CHECK(r.status == FCT_NO_FACTOR);
    CHECK_EQ_U64(r.iterations, 1);
}

int main(void)
{
    fct_pminus1_config def = fct_pminus1_config_default();
    CHECK_EQ_U64(def.smoothness_bound, FCT_PM1_DEFAULT_SMOOTHNESS_BOUND);
    CHECK_EQ_U64(def.start_base, FCT_PM1_DEFAULT_START_BASE);
    CHECK_EQ_U64(def.max_bases, FCT_PM1_DEFAULT_MAX_BASES);
    CHECK_EQ_U64(def.max_modmul_ops, FCT_PM1_DEFAULT_MAX_MODMUL_OPS);
    test_invalid_and_even();
    test_favorable_smooth_p_minus_one();
    test_even_and_powers_of_two();
    test_unfavorable_smoothness();
    test_degenerate_d_equals_n();
    test_prime_inputs();
    test_work_limit_and_config();
    TEST_MAIN_END();
}
