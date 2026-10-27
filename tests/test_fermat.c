/*
 * tests/test_fermat.c — Fermat's method: edge cases, squares, close and
 * distant factors, budget accounting, and the 2^64 overflow boundary.
 */
#include "factorium/fermat.h"
#include "test_harness.h"

static fct_fermat_config cfg_with(uint64_t iters)
{
    fct_fermat_config c = fct_fermat_config_default();
    c.max_iterations = iters;
    return c;
}

static void test_invalid_and_even(void)
{
    fct_factor_result r;
    for (uint64_t n = 0; n <= 3; n++) {
        r = fct_fermat(n, NULL);
        CHECK(r.status == FCT_INVALID_INPUT);
        CHECK(r.reason == FCT_TERM_INVALID_INPUT);
    }
    r = fct_fermat(4, NULL);
    CHECK(r.status == FCT_OK);
    CHECK_EQ_U64(r.factor, 2);
    CHECK_EQ_U64(r.cofactor, 2);
    r = fct_fermat(100, NULL);
    CHECK(r.status == FCT_OK && r.factor == 2 && r.cofactor == 50);
    r = fct_fermat(18446744073709551614ULL, NULL); /* 2^64 - 2 */
    CHECK(r.status == FCT_OK);
    CHECK_EQ_U64(r.factor, 2);
    CHECK_EQ_U64(r.cofactor, 9223372036854775807ULL);
    CHECK_EQ_U64(r.iterations, 0);
    CHECK_EQ_U64(r.gcd_calls, 0);
    CHECK_EQ_U64(r.modmul_ops, 0);
}

static void test_squares_and_trivial(void)
{
    fct_factor_result r;
    r = fct_fermat(9, NULL);
    CHECK(r.status == FCT_OK && r.factor == 3 && r.cofactor == 3);
    r = fct_fermat(25, NULL);
    CHECK(r.status == FCT_OK && r.factor == 5 && r.cofactor == 5);
    r = fct_fermat(15, NULL);
    CHECK(r.status == FCT_OK);
    CHECK_EQ_U64(r.factor, 3);
    CHECK_EQ_U64(r.cofactor, 5);
    CHECK_EQ_U64(r.iterations, 1);
    r = fct_fermat(91, NULL); /* 7 * 13 */
    CHECK(r.status == FCT_OK && r.factor == 7 && r.cofactor == 13);
    CHECK_EQ_U64(r.iterations, 1);
    r = fct_fermat(77, NULL); /* 7 * 11 */
    CHECK(r.status == FCT_OK && r.factor == 7 && r.cofactor == 11);
    CHECK_EQ_U64(r.iterations, 1);
    r = fct_fermat(2015, NULL); /* 5*13*31; a=48, b=17 -> 31*65 */
    CHECK(r.status == FCT_OK);
    CHECK_EQ_U64(r.factor, 31);
    CHECK_EQ_U64(r.cofactor, 65);
    CHECK_EQ_U64(r.iterations, 4);
}

static void test_close_factors(void)
{
    fct_factor_result r;
    /* 1000000007 * 1000000009 = 1000000008^2 - 1 */
    r = fct_fermat(1000000016000000063ULL, NULL);
    CHECK(r.status == FCT_OK);
    CHECK_EQ_U64(r.factor, 1000000007ULL);
    CHECK_EQ_U64(r.cofactor, 1000000009ULL);
    CHECK_EQ_U64(r.iterations, 1);
    /* 2^64 - 1 = (2^32 - 1)(2^32 + 1): exercises a^2 == 2^64 exactly,
     * which does not fit in uint64_t. */
    r = fct_fermat(18446744073709551615ULL, NULL);
    CHECK(r.status == FCT_OK);
    CHECK_EQ_U64(r.factor, 4294967295ULL);
    CHECK_EQ_U64(r.cofactor, 4294967297ULL);
    CHECK_EQ_U64(r.iterations, 1);
}

static void test_distant_factors(void)
{
    fct_factor_result r;
    fct_fermat_config c = cfg_with(1ULL << 20);
    /* 3 * 1000003: solution at a = 500003, a0 = 1733. */
    r = fct_fermat(3000009, &c);
    CHECK(r.status == FCT_OK);
    CHECK_EQ_U64(r.factor, 3);
    CHECK_EQ_U64(r.cofactor, 1000003);
    CHECK_EQ_U64(r.iterations, 500003ULL - 1733ULL + 1ULL);
}

static void test_budget(void)
{
    fct_factor_result r;
    fct_fermat_config c1 = cfg_with(100000);
    fct_fermat_config c0 = cfg_with(0);
    /* Same input, small budget: must stop exactly at the budget. */
    r = fct_fermat(3000009, &c1);
    CHECK(r.status == FCT_NO_FACTOR);
    CHECK(r.reason == FCT_TERM_BUDGET_EXHAUSTED);
    CHECK_EQ_U64(r.iterations, 100000);
    /* Zero budget: stop immediately. */
    r = fct_fermat(3000009, &c0);
    CHECK(r.status == FCT_NO_FACTOR);
    CHECK_EQ_U64(r.iterations, 0);
    /* Prime input with a budget: NO_FACTOR (never a primality claim). */
    r = fct_fermat(9999999967ULL, &c1);
    CHECK(r.status == FCT_NO_FACTOR);
    CHECK(r.reason == FCT_TERM_BUDGET_EXHAUSTED);
    CHECK_EQ_U64(r.iterations, 100000);
    /* gcd/modmul counters stay zero: Fermat does no modular arithmetic. */
    CHECK_EQ_U64(r.gcd_calls, 0);
    CHECK_EQ_U64(r.modmul_ops, 0);
}

int main(void)
{
    fct_fermat_config def = fct_fermat_config_default();
    CHECK_EQ_U64(def.max_iterations, FCT_FERMAT_DEFAULT_MAX_ITERATIONS);
    test_invalid_and_even();
    test_squares_and_trivial();
    test_close_factors();
    test_distant_factors();
    test_budget();
    TEST_MAIN_END();
}
