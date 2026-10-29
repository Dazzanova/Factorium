/*
 * factorium/src/bench_main.c — Benchmark runner CLI.
 *
 * Runs controlled comparative factorization experiments, measures runtimes
 * via CLOCK_MONOTONIC, computes summary statistics, and outputs raw/summary CSVs.
 */
#define _POSIX_C_SOURCE 200809L

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "factorium/benchmark.h"

#define U64_FMT PRIu64

static void print_usage(FILE *out)
{
    fputs(
        "Factorium Benchmark Runner (Phase 2)\n"
        "\n"
        "Usage:\n"
        "  factorium_bench [options]\n"
        "\n"
        "Options:\n"
        "  --smoke              run fast smoke benchmark subset (<2s)\n"
        "  --full               run full dataset benchmark (default)\n"
        "  --raw-out <file>     destination path for raw CSV records (default: bench_raw.csv)\n"
        "  --summary-out <file> destination path for summary CSV (default: bench_summary.csv)\n"
        "  --rho-trials <n>     number of randomized Rho seeds (default: 5, smoke: 1)\n"
        "  --iters <n>          override max iterations budget\n"
        "  --bound <B>          override p-1 smoothness bound\n"
        "  --quiet              suppress console status output\n"
        "  --help               display this help message\n",
        out);
}

int main(int argc, char **argv)
{
    bool smoke = false;
    const char *raw_path = "bench_raw.csv";
    const char *summary_path = "bench_summary.csv";
    uint64_t rho_trials = 0;
    uint64_t max_iters = 0;
    uint64_t smoothness = 0;
    bool quiet = false;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--smoke") == 0) {
            smoke = true;
        } else if (strcmp(argv[i], "--full") == 0) {
            smoke = false;
        } else if (strcmp(argv[i], "--raw-out") == 0 && i + 1 < argc) {
            raw_path = argv[++i];
        } else if (strcmp(argv[i], "--summary-out") == 0 && i + 1 < argc) {
            summary_path = argv[++i];
        } else if (strcmp(argv[i], "--rho-trials") == 0 && i + 1 < argc) {
            rho_trials = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--iters") == 0 && i + 1 < argc) {
            max_iters = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--bound") == 0 && i + 1 < argc) {
            smoothness = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--quiet") == 0) {
            quiet = true;
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            print_usage(stdout);
            return 0;
        } else {
            fprintf(stderr, "error: unrecognized option '%s'\n", argv[i]);
            print_usage(stderr);
            return 2;
        }
    }

    fct_bench_suite_options opts = fct_bench_suite_options_default(smoke);
    if (rho_trials > 0) opts.rho_trials = rho_trials;
    if (max_iters > 0) opts.max_iterations = max_iters;
    if (smoothness > 0) opts.smoothness_bound = smoothness;

    if (!quiet) {
        printf("Running Factorium Benchmark Suite [%s mode]...\n", smoke ? "SMOKE" : "FULL");
        printf("Configuration: rho_trials=%" U64_FMT ", max_iters=%" U64_FMT ", bound=%" U64_FMT "\n\n",
               opts.rho_trials, opts.max_iterations, opts.smoothness_bound);
    }

    fct_bench_suite_result res = fct_bench_run_suite(&opts);
    if (!res.records || res.record_count == 0) {
        fprintf(stderr, "error: benchmark execution failed or yielded zero records\n");
        fct_bench_suite_result_free(&res);
        return 1;
    }

    /* Export Raw CSV */
    FILE *raw_file = fopen(raw_path, "w");
    if (!raw_file) {
        fprintf(stderr, "error: cannot open raw output file '%s' for writing\n", raw_path);
        fct_bench_suite_result_free(&res);
        return 1;
    }
    fct_bench_export_raw_csv(res.records, res.record_count, raw_file);
    fclose(raw_file);

    /* Export Summary CSV */
    FILE *sum_file = fopen(summary_path, "w");
    if (!sum_file) {
        fprintf(stderr, "error: cannot open summary output file '%s' for writing\n", summary_path);
        fct_bench_suite_result_free(&res);
        return 1;
    }
    fct_bench_export_summary_csv(res.summaries, res.summary_count, sum_file);
    fclose(sum_file);

    if (!quiet) {
        printf("Completed %zu trials across %zu dataset entries.\n", res.record_count, res.summary_count / 3);
        printf("Raw records saved to:     %s\n", raw_path);
        printf("Summary metrics saved to: %s\n\n", summary_path);

        /* Print quick preview table */
        printf("%-20s %-8s %-10s %-12s %-12s\n", "Dataset", "Algo", "Success", "Median (us)", "Median Iters");
        printf("----------------------------------------------------------------------\n");
        for (size_t i = 0; i < res.summary_count; i++) {
            const fct_bench_summary *s = &res.summaries[i];
            printf("%-20s %-8s %3" U64_FMT "/%-4" U64_FMT " %10.2f us %12" U64_FMT "\n",
                   s->dataset_id, fct_bench_algo_name(s->algo),
                   s->num_success, s->num_trials,
                   (double)s->median_runtime_ns / 1000.0,
                   s->median_iterations);
        }
    }

    fct_bench_suite_result_free(&res);
    return 0;
}
