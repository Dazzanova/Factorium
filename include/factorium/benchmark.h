/*
 * factorium/benchmark.h — benchmarking harness (Phase 2).
 *
 * Runs each algorithm on a curated dataset, times every trial with
 * CLOCK_MONOTONIC, computes summary statistics and writes RFC 4180 CSV.
 */
#ifndef FACTORIUM_BENCHMARK_H
#define FACTORIUM_BENCHMARK_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "factorium/common.h"
#include "factorium/fermat.h"
#include "factorium/pollard_p_minus_one.h"
#include "factorium/pollard_rho.h"

#define BENCH_MAX_NAME 64
#define BENCH_MAX_DESC 256

/* One test case from the dataset. */
typedef struct {
    char id[BENCH_MAX_NAME];
    char category[BENCH_MAX_NAME];
    uint64_t n;
    unsigned bit_length;
    uint64_t known_factor_p;
    uint64_t known_factor_q; /* 0 when n has more than two factors */
    char description[BENCH_MAX_DESC];
    char hypothesis[BENCH_MAX_DESC];
} BenchCase;

typedef enum {
    BENCH_ALGO_FERMAT = 0,
    BENCH_ALGO_RHO,
    BENCH_ALGO_PMINUS1,
    BENCH_ALGO_COUNT
} BenchAlgo;

const char *algo_name(BenchAlgo algo);

/* Parameters for one trial. */
typedef struct {
    uint64_t max_iterations;
    uint64_t max_retries;
    uint64_t seed;
    uint64_t smoothness_bound;
    uint64_t start_base;
    uint64_t max_bases;
    uint64_t max_modmul_ops;
} BenchParams;

BenchParams default_params(BenchAlgo algo);

/* Result of one trial. */
typedef struct {
    uint64_t trial_id;
    char dataset_id[BENCH_MAX_NAME];
    char category[BENCH_MAX_NAME];
    uint64_t n;
    unsigned bit_length;
    uint64_t ground_truth_p;
    uint64_t ground_truth_q;
    BenchAlgo algo;
    BenchParams params;
    ResultCode code;
    uint64_t factor_found;
    uint64_t cofactor_found;
    bool factor_is_valid;
    uint64_t runtime_ns;
    uint64_t iterations;
    uint64_t gcd_calls;
    uint64_t modmul_ops;
} BenchRecord;

/* Aggregated statistics for a group of trials. */
typedef struct {
    char dataset_id[BENCH_MAX_NAME];
    char category[BENCH_MAX_NAME];
    BenchAlgo algo;
    uint64_t num_trials;
    uint64_t num_success;
    uint64_t num_timeout;
    uint64_t num_fail;
    double success_rate;
    uint64_t median_runtime_ns;
    uint64_t min_runtime_ns;
    uint64_t max_runtime_ns;
    uint64_t q1_runtime_ns;
    uint64_t q3_runtime_ns;
    uint64_t median_iterations;
    uint64_t median_gcd_calls;
    uint64_t median_modmul_ops;
} BenchSummary;

/* Dataset access. */
size_t get_dataset_count(bool smoke_only);
const BenchCase *get_dataset(bool smoke_only, size_t *count_out);

/* Run one trial. */
BenchRecord run_trial(uint64_t trial_id, const BenchCase *test_case,
                      BenchAlgo algo, const BenchParams *params);

typedef struct {
    bool smoke_only;
    uint64_t rho_trials;
    uint64_t max_iterations;
    uint64_t smoothness_bound;
} BenchOptions;

BenchOptions default_suite_options(bool smoke_only);

typedef struct {
    BenchRecord *records;
    size_t record_count;
    size_t record_cap;
    BenchSummary *summaries;
    size_t summary_count;
    size_t summary_cap;
} BenchSuite;

BenchSuite run_suite(const BenchOptions *opts);
void free_suite(BenchSuite *suite);

BenchSummary calculate_summary(const BenchRecord *records, size_t count);

bool export_raw_csv(const BenchRecord *records, size_t count, FILE *out);
bool export_summary_csv(const BenchSummary *summaries, size_t count, FILE *out);

#endif /* FACTORIUM_BENCHMARK_H */
