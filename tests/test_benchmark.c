/*
 * tests/test_benchmark.c — Unit and integration tests for benchmarking module.
 */
#define _POSIX_C_SOURCE 200809L

#include "test_harness.h"
#include "factorium/benchmark.h"
#include "factorium/modular_arithmetic.h"
#include "factorium/primality.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void test_dataset_integrity(void)
{
    size_t count = 0;
    const BenchCase *cases = get_dataset(false, &count);
    CHECK(cases != NULL);
    CHECK(count >= 20);

    for (size_t i = 0; i < count; i++) {
        const BenchCase *c = &cases[i];
        CHECK(strlen(c->id) > 0);
        CHECK(strlen(c->category) > 0);
        CHECK(strlen(c->description) > 0);
        CHECK(strlen(c->hypothesis) > 0);
        CHECK(c->n >= 2);

        /* Verify bit length */
        unsigned calculated_bits = 0;
        uint64_t tmp = c->n;
        while (tmp > 0) {
            calculated_bits++;
            tmp >>= 1;
        }
        CHECK_EQ_U64((uint64_t)calculated_bits, (uint64_t)c->bit_length);

        /* If ground truth p, q are supplied, check validity */
        if (c->known_factor_p > 1 && c->known_factor_p < c->n) {
            CHECK_EQ_U64(c->n % c->known_factor_p, 0);
        }
        if (c->known_factor_q > 1 && c->known_factor_q < c->n) {
            CHECK_EQ_U64(c->n % c->known_factor_q, 0);
        }
    }
}

static void test_smoke_dataset(void)
{
    size_t smoke_count = 0;
    const BenchCase *smoke = get_dataset(true, &smoke_count);
    CHECK(smoke != NULL);
    CHECK(smoke_count > 0 && smoke_count <= 25);

    /* Verify each category is represented */
    bool has_close = false, has_distant = false, has_small = false, has_p1 = false;
    for (size_t i = 0; i < smoke_count; i++) {
        if (strcmp(smoke[i].category, "close_semiprimes") == 0) has_close = true;
        if (strcmp(smoke[i].category, "distant_semiprimes") == 0) has_distant = true;
        if (strcmp(smoke[i].category, "small_prime_factors") == 0) has_small = true;
        if (strstr(smoke[i].category, "pminus1") != NULL) has_p1 = true;
    }
    CHECK(has_close);
    CHECK(has_distant);
    CHECK(has_small);
    CHECK(has_p1);
}

static void test_single_trial_run(void)
{
    BenchCase tc;
    strcpy(tc.id, "test_360");
    strcpy(tc.category, "small_factors");
    tc.n = 360;
    tc.bit_length = 9;
    tc.known_factor_p = 2;
    tc.known_factor_q = 3;
    strcpy(tc.description, "test case 360");
    strcpy(tc.hypothesis, "trivial factor");

    BenchParams p = default_params(BENCH_ALGO_FERMAT);
    BenchRecord rec = run_trial(42, &tc, BENCH_ALGO_FERMAT, &p);

    CHECK_EQ_U64(rec.trial_id, 42);
    CHECK_EQ_U64(rec.n, 360);
    CHECK_EQ_U64((uint64_t)rec.code, (uint64_t)RESULT_FOUND);
    CHECK_EQ_U64(rec.factor_found, 2);
    CHECK_EQ_U64(rec.cofactor_found, 180);
    CHECK(rec.factor_is_valid);
    CHECK(rec.runtime_ns > 0);
}

static void test_summary_calculation(void)
{
    BenchRecord recs[5];
    memset(recs, 0, sizeof(recs));

    for (size_t i = 0; i < 5; i++) {
        strcpy(recs[i].dataset_id, "test_summary");
        strcpy(recs[i].category, "unit_test");
        recs[i].algo = BENCH_ALGO_RHO;
        recs[i].runtime_ns = (i + 1) * 1000; /* 1000, 2000, 3000, 4000, 5000 */
        recs[i].iterations = (i + 1) * 10;
        recs[i].gcd_calls = (i + 1) * 5;
        recs[i].modmul_ops = (i + 1) * 20;
        recs[i].code = (i < 4) ? RESULT_FOUND : RESULT_NO_FACTOR;
    }

    BenchSummary sum = calculate_summary(recs, 5);

    CHECK_EQ_U64(sum.num_trials, 5);
    CHECK_EQ_U64(sum.num_success, 4);
    CHECK_EQ_U64(sum.num_timeout, 1);
    CHECK_EQ_U64(sum.num_fail, 0);

    /* Median of 1000, 2000, 3000, 4000, 5000 is 3000 */
    CHECK_EQ_U64(sum.min_runtime_ns, 1000);
    CHECK_EQ_U64(sum.max_runtime_ns, 5000);
    CHECK_EQ_U64(sum.median_runtime_ns, 3000);
    CHECK_EQ_U64(sum.median_iterations, 30);
    CHECK_EQ_U64(sum.median_gcd_calls, 15);
    CHECK_EQ_U64(sum.median_modmul_ops, 60);
    CHECK(sum.success_rate >= 0.79 && sum.success_rate <= 0.81);
}

static void test_csv_export_and_escaping(void)
{
    BenchRecord rec;
    memset(&rec, 0, sizeof(rec));
    rec.trial_id = 1;
    strcpy(rec.dataset_id, "escape,\"test\"");
    strcpy(rec.category, "cat,with,commas");
    rec.n = 15;
    rec.bit_length = 4;
    rec.ground_truth_p = 3;
    rec.ground_truth_q = 5;
    rec.algo = BENCH_ALGO_FERMAT;
    rec.code = RESULT_FOUND;
    rec.factor_found = 3;
    rec.cofactor_found = 5;
    rec.factor_is_valid = true;
    rec.runtime_ns = 1234;
    rec.iterations = 1;

    char buf[1024];
    FILE *mem = fmemopen(buf, sizeof(buf), "w");
    CHECK(mem != NULL);
    bool ok = export_raw_csv(&rec, 1, mem);
    CHECK(ok);
    fclose(mem);

    /* Check that commas and quotes inside field were properly escaped */
    CHECK(strstr(buf, "\"escape,\"\"test\"\"\"") != NULL);
    CHECK(strstr(buf, "\"cat,with,commas\"") != NULL);
    CHECK(strstr(buf, "1234") != NULL);
}

static void test_fast_smoke_suite(void)
{
    BenchOptions opts = default_suite_options(true);
    opts.rho_trials = 1;
    opts.max_iterations = 5000;
    opts.smoothness_bound = 5000;

    BenchSuite res = run_suite(&opts);
    CHECK(res.records != NULL);
    CHECK(res.record_count > 0);
    CHECK(res.summaries != NULL);
    CHECK(res.summary_count > 0);

    /* Verify all records have valid status and non-negative counters */
    for (size_t i = 0; i < res.record_count; i++) {
        CHECK(res.records[i].n >= 2);
        CHECK(res.records[i].trial_id > 0);
    }

    free_suite(&res);
}

int main(void)
{
    test_dataset_integrity();
    test_smoke_dataset();
    test_single_trial_run();
    test_summary_calculation();
    test_csv_export_and_escaping();
    test_fast_smoke_suite();
    TEST_MAIN_END();
}
