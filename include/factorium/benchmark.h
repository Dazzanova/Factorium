/*
 * factorium/benchmark.h — Benchmarking harness for comparative analysis.
 *
 * Provides controlled datasets, timing with CLOCK_MONOTONIC, raw trial records,
 * aggregate statistics (median, IQR, min, max, success rate), and RFC 4180 CSV export.
 */
#ifndef FACTORIUM_BENCHMARK_H
#define FACTORIUM_BENCHMARK_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "factorium/common.h"
#include "factorium/fermat.h"
#include "factorium/pollard_rho.h"
#include "factorium/pollard_p_minus_one.h"

#define FCT_BENCH_MAX_NAME_LEN 64
#define FCT_BENCH_MAX_DESC_LEN 256

/* A single test case in a benchmark dataset. */
typedef struct {
    char id[FCT_BENCH_MAX_NAME_LEN];            /* Unique identifier (e.g., "close_32_01") */
    char category[FCT_BENCH_MAX_NAME_LEN];      /* Category name (e.g., "close_semiprimes") */
    uint64_t n;                                 /* Input integer to factor */
    unsigned bit_length;                        /* Number of bits in n */
    uint64_t known_factor_p;                    /* Ground-truth factor p */
    uint64_t known_factor_q;                    /* Ground-truth factor q (or 0 if >2 factors) */
    char description[FCT_BENCH_MAX_DESC_LEN];   /* Short description */
    char hypothesis[FCT_BENCH_MAX_DESC_LEN];    /* Hypothesis being tested */
} fct_bench_case;

/* Algorithm identifier for benchmarks. */
typedef enum {
    FCT_BENCH_ALGO_FERMAT = 0,
    FCT_BENCH_ALGO_RHO,
    FCT_BENCH_ALGO_PMINUS1,
    FCT_BENCH_ALGO_COUNT
} fct_bench_algo;

const char *fct_bench_algo_name(fct_bench_algo algo);

/* Execution parameters for a benchmark trial. */
typedef struct {
    uint64_t max_iterations;    /* Fermat and Rho iteration cap */
    uint64_t max_retries;       /* Rho retry cap */
    uint64_t seed;              /* Rho initial seed */
    uint64_t smoothness_bound;  /* p-1 bound B */
    uint64_t start_base;        /* p-1 start base */
    uint64_t max_bases;         /* p-1 number of bases */
    uint64_t max_modmul_ops;    /* p-1 max modular multiplications */
} fct_bench_params;

fct_bench_params fct_bench_params_default(fct_bench_algo algo);

/* Result of a single benchmark trial. */
typedef struct {
    uint64_t trial_id;
    char dataset_id[FCT_BENCH_MAX_NAME_LEN];
    char category[FCT_BENCH_MAX_NAME_LEN];
    uint64_t n;
    unsigned bit_length;
    uint64_t ground_truth_p;
    uint64_t ground_truth_q;
    fct_bench_algo algo;
    fct_bench_params params;
    fct_status status;
    fct_termination reason;
    uint64_t factor_found;
    uint64_t cofactor_found;
    bool factor_is_valid;       /* True if factor divides n and matches ground truth */
    uint64_t runtime_ns;        /* Measured via CLOCK_MONOTONIC */
    uint64_t iterations;
    uint64_t gcd_calls;
    uint64_t modmul_ops;
} fct_bench_record;

/* Aggregate summary statistics for a group of trials. */
typedef struct {
    char dataset_id[FCT_BENCH_MAX_NAME_LEN];
    char category[FCT_BENCH_MAX_NAME_LEN];
    fct_bench_algo algo;
    uint64_t num_trials;
    uint64_t num_success;
    uint64_t num_timeout;       /* Budget / work limit reached */
    uint64_t num_fail;          /* Invalid input, memory, or error */
    double success_rate;        /* num_success / num_trials */
    uint64_t median_runtime_ns;
    uint64_t min_runtime_ns;
    uint64_t max_runtime_ns;
    uint64_t q1_runtime_ns;     /* 25th percentile */
    uint64_t q3_runtime_ns;     /* 75th percentile */
    uint64_t median_iterations;
    uint64_t median_gcd_calls;
    uint64_t median_modmul_ops;
} fct_bench_summary;

/* Dataset retrieval */
size_t fct_bench_get_dataset_count(bool smoke_only);
const fct_bench_case *fct_bench_get_dataset(bool smoke_only, size_t *count_out);

/* Single trial runner */
fct_bench_record fct_bench_run_trial(uint64_t trial_id,
                                     const fct_bench_case *test_case,
                                     fct_bench_algo algo,
                                     const fct_bench_params *params);

/* Batch runner options */
typedef struct {
    bool smoke_only;            /* Use fast, smaller subset */
    uint64_t rho_trials;        /* Number of randomized Rho seeds (default: 5, smoke: 1) */
    uint64_t max_iterations;    /* Override max iterations if > 0 */
    uint64_t smoothness_bound;  /* Override p-1 bound B if > 0 */
} fct_bench_suite_options;

fct_bench_suite_options fct_bench_suite_options_default(bool smoke_only);

/* Suite runner */
typedef struct {
    fct_bench_record *records;
    size_t record_count;
    size_t record_cap;
    fct_bench_summary *summaries;
    size_t summary_count;
    size_t summary_cap;
} fct_bench_suite_result;

fct_bench_suite_result fct_bench_run_suite(const fct_bench_suite_options *opts);
void fct_bench_suite_result_free(fct_bench_suite_result *res);

/* Summary statistics calculation */
fct_bench_summary fct_bench_calculate_summary(const fct_bench_record *records,
                                              size_t count);

/* CSV Export */
bool fct_bench_export_raw_csv(const fct_bench_record *records,
                              size_t count,
                              FILE *out);

bool fct_bench_export_summary_csv(const fct_bench_summary *summaries,
                                  size_t count,
                                  FILE *out);

#endif /* FACTORIUM_BENCHMARK_H */
