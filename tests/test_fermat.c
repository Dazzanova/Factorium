/*
 * tests/test_fermat.c — Fermat's method: edge cases, squares, close and
 * distant factors, budget accounting, and the 2^64 overflow boundary.
 */
#include "factorium/fermat.h"
#include "test_harness.h"

static void test_invalid_and_even(void)
{
    FactorResult r;
    for (uint64_t n = 0; n <= 3; n++) {
        r = fermat(n, FERMAT_DEFAULT_ITERATIONS);
        CHECK(r.code == RESULT_INVALID);
    }
    r = fermat(4, FERMAT_DEFAULT_ITERATIONS);
    CHECK(r.code == RESULT_FOUND);
    CHECK_EQ_U64(r.factor, 2);
    CHECK_EQ_U64(r.cofactor, 2);
    r = fermat(100, FERMAT_DEFAULT_ITERATIONS);
    CHECK(r.code == RESULT_FOUND && r.factor == 2 && r.cofactor == 50);
    r = fermat(18446744073709551614ULL, FERMAT_DEFAULT_ITERATIONS); /* 2^64 - 2 */
    CHECK(r.code == RESULT_FOUND);
    CHECK_EQ_U64(r.factor, 2);
    CHECK_EQ_U64(r.cofactor, 9223372036854775807ULL);
    CHECK_EQ_U64(r.iterations, 0);
    CHECK_EQ_U64(r.gcd_calls, 0);
    CHECK_EQ_U64(r.modmul_ops, 0);
}

static void test_squares_and_trivial(void)
{
    FactorResult r;
    r = fermat(9, FERMAT_DEFAULT_ITERATIONS);
    CHECK(r.code == RESULT_FOUND && r.factor == 3 && r.cofactor == 3);
    r = fermat(25, FERMAT_DEFAULT_ITERATIONS);
    CHECK(r.code == RESULT_FOUND && r.factor == 5 && r.cofactor == 5);
    r = fermat(15, FERMAT_DEFAULT_ITERATIONS);
    CHECK(r.code == RESULT_FOUND);
    CHECK_EQ_U64(r.factor, 3);
    CHECK_EQ_U64(r.cofactor, 5);
    CHECK_EQ_U64(r.iterations, 1);
    r = fermat(91, FERMAT_DEFAULT_ITERATIONS); /* 7 * 13 */
    CHECK(r.code == RESULT_FOUND && r.factor == 7 && r.cofactor == 13);
    CHECK_EQ_U64(r.iterations, 1);
    r = fermat(77, FERMAT_DEFAULT_ITERATIONS); /* 7 * 11 */
    CHECK(r.code == RESULT_FOUND && r.factor == 7 && r.cofactor == 11);
    CHECK_EQ_U64(r.iterations, 1);
    r = fermat(2015, FERMAT_DEFAULT_ITERATIONS); /* 5*13*31; a=48, b=17 -> 31*65 */
    CHECK(r.code == RESULT_FOUND);
    CHECK_EQ_U64(r.factor, 31);
    CHECK_EQ_U64(r.cofactor, 65);
    CHECK_EQ_U64(r.iterations, 4);
}

static void test_close_factors(void)
{
    FactorResult r;
    /* 1000000007 * 1000000009 = 1000000008^2 - 1 */
    r = fermat(1000000016000000063ULL, FERMAT_DEFAULT_ITERATIONS);
    CHECK(r.code == RESULT_FOUND);
    CHECK_EQ_U64(r.factor, 1000000007ULL);
    CHECK_EQ_U64(r.cofactor, 1000000009ULL);
    CHECK_EQ_U64(r.iterations, 1);
    /* 2^64 - 1 = (2^32 - 1)(2^32 + 1): exercises a^2 == 2^64 exactly,
     * which does not fit in uint64_t. */
    r = fermat(18446744073709551615ULL, FERMAT_DEFAULT_ITERATIONS);
    CHECK(r.code == RESULT_FOUND);
    CHECK_EQ_U64(r.factor, 4294967295ULL);
    CHECK_EQ_U64(r.cofactor, 4294967297ULL);
    CHECK_EQ_U64(r.iterations, 1);
}

static void test_distant_factors(void)
{
    FactorResult r;
    /* 3 * 1000003: solution at a = 500003, a0 = 1733. */
    r = fermat(3000009, 1ULL << 20);
    CHECK(r.code == RESULT_FOUND);
    CHECK_EQ_U64(r.factor, 3);
    CHECK_EQ_U64(r.cofactor, 1000003);
    CHECK_EQ_U64(r.iterations, 500003ULL - 1733ULL + 1ULL);
}

static void test_budget(void)
{
    FactorResult r;
    /* Same input, small budget: must stop exactly at the budget. */
    r = fermat(3000009, 100000);
    CHECK(r.code == RESULT_NO_FACTOR);
    CHECK_EQ_U64(r.iterations, 100000);
    /* Zero budget: stop immediately. */
    r = fermat(3000009, 0);
    CHECK(r.code == RESULT_NO_FACTOR);
    CHECK_EQ_U64(r.iterations, 0);
    /* Prime input with a budget: NO_FACTOR (never a primality claim). */
    r = fermat(9999999967ULL, 100000);
    CHECK(r.code == RESULT_NO_FACTOR);
    CHECK_EQ_U64(r.iterations, 100000);
    /* gcd/modmul counters stay zero: Fermat does no modular arithmetic. */
    CHECK_EQ_U64(r.gcd_calls, 0);
    CHECK_EQ_U64(r.modmul_ops, 0);
}

int main(void)
{
    test_invalid_and_even();
    test_squares_and_trivial();
    test_close_factors();
    test_distant_factors();
    test_budget();
    TEST_MAIN_END();
}
