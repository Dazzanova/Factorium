/*
 * factorium/src/benchmark.c — trial runner, statistics and CSV export.
 */
#define _POSIX_C_SOURCE 199309L

#include "factorium/benchmark.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define U64_FMT PRIu64

/* ------------------------------------------------------------------ */
/* Algorithm Names & Parameters                                       */
/* ------------------------------------------------------------------ */

const char *algo_name(BenchAlgo algo)
{
    switch (algo) {
    case BENCH_ALGO_FERMAT: return "fermat";
    case BENCH_ALGO_RHO: return "rho";
    case BENCH_ALGO_PMINUS1: return "pminus1";
    default: return "unknown";
    }
}

BenchParams default_params(BenchAlgo algo)
{
    BenchParams p;
    memset(&p, 0, sizeof(p));
    switch (algo) {
    case BENCH_ALGO_FERMAT:
        p.max_iterations = 1000000ULL; /* 10^6 iterations */
        break;
    case BENCH_ALGO_RHO:
        p.max_iterations = 1000000ULL;
        p.max_retries = 8;
        p.seed = 1;
        break;
    case BENCH_ALGO_PMINUS1:
        p.smoothness_bound = 100000ULL;
        p.start_base = 2;
        p.max_bases = 5;
        p.max_modmul_ops = 1ULL << 22; /* 4M modmuls */
        break;
    default:
        break;
    }
    return p;
}

BenchOptions default_suite_options(bool smoke_only)
{
    BenchOptions o;
    o.smoke_only = smoke_only;
    o.rho_trials = smoke_only ? 1 : 5;
    o.max_iterations = smoke_only ? 20000ULL : 500000ULL;
    o.smoothness_bound = smoke_only ? 10000ULL : 100000ULL;
    return o;
}

/* ------------------------------------------------------------------ */
/* Monotonic Clock Timing Helper                                      */
/* ------------------------------------------------------------------ */

static inline uint64_t get_time_ns(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

/* ------------------------------------------------------------------ */
/* Single Trial Runner                                                */
/* ------------------------------------------------------------------ */

BenchRecord run_trial(uint64_t trial_id,
                                     const BenchCase *tc,
                                     BenchAlgo algo,
                                     const BenchParams *params_in)
{
    BenchRecord rec;
    memset(&rec, 0, sizeof(rec));

    rec.trial_id = trial_id;
    snprintf(rec.dataset_id, sizeof(rec.dataset_id), "%s", tc->id);
    snprintf(rec.category, sizeof(rec.category), "%s", tc->category);
    rec.n = tc->n;
    rec.bit_length = tc->bit_length;
    rec.ground_truth_p = tc->known_factor_p;
    rec.ground_truth_q = tc->known_factor_q;
    rec.algo = algo;

    BenchParams p = params_in ? *params_in : default_params(algo);
    rec.params = p;

    FactorResult res;
    memset(&res, 0, sizeof(res));

    uint64_t t0 = 0, t1 = 0;

    switch (algo) {
    case BENCH_ALGO_FERMAT:
        t0 = get_time_ns();
        res = fermat(tc->n, p.max_iterations ? p.max_iterations
                                             : FERMAT_DEFAULT_ITERATIONS);
        t1 = get_time_ns();
        break;
    case BENCH_ALGO_RHO:
        t0 = get_time_ns();
        res = rho(tc->n,
                  p.max_iterations ? p.max_iterations : RHO_DEFAULT_ITERATIONS,
                  p.max_retries ? p.max_retries : RHO_DEFAULT_RETRIES,
                  p.seed ? p.seed : RHO_DEFAULT_SEED);
        t1 = get_time_ns();
        break;
    case BENCH_ALGO_PMINUS1:
        t0 = get_time_ns();
        res = pminus1(tc->n,
                      p.smoothness_bound ? p.smoothness_bound : PM1_DEFAULT_BOUND,
                      p.start_base ? p.start_base : PM1_DEFAULT_START_BASE,
                      p.max_bases ? p.max_bases : PM1_DEFAULT_BASES,
                      p.max_modmul_ops ? p.max_modmul_ops : PM1_DEFAULT_MAX_OPS);
        t1 = get_time_ns();
        break;
    default:
        res.code = RESULT_INVALID;
        break;
    }

    rec.runtime_ns = (t1 >= t0) ? (t1 - t0) : 0;
    rec.code = res.code;
    rec.factor_found = res.factor;
    rec.cofactor_found = res.cofactor;
    rec.iterations = res.iterations;
    rec.gcd_calls = res.gcd_calls;
    rec.modmul_ops = res.modmul_ops;

    /* Accept the factor only if it really divides n. */
    rec.factor_is_valid = false;
    if (res.code == RESULT_FOUND && res.factor > 1 && res.factor < tc->n) {
        if (tc->n % res.factor == 0) {
            rec.factor_is_valid = true;
        }
    }

    return rec;
}

/* ------------------------------------------------------------------ */
/* Summary Calculation                                                */
/* ------------------------------------------------------------------ */

static int cmp_u64_asc(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a;
    uint64_t y = *(const uint64_t *)b;
    return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

BenchSummary calculate_summary(const BenchRecord *records,
                                              size_t count)
{
    BenchSummary s;
    memset(&s, 0, sizeof(s));

    if (count == 0)
        return s;

    snprintf(s.dataset_id, sizeof(s.dataset_id), "%s", records[0].dataset_id);
    snprintf(s.category, sizeof(s.category), "%s", records[0].category);
    s.algo = records[0].algo;
    s.num_trials = (uint64_t)count;

    uint64_t *runtimes = (uint64_t *)malloc(count * sizeof(uint64_t));
    uint64_t *iters = (uint64_t *)malloc(count * sizeof(uint64_t));
    uint64_t *gcds = (uint64_t *)malloc(count * sizeof(uint64_t));
    uint64_t *modmuls = (uint64_t *)malloc(count * sizeof(uint64_t));

    if (!runtimes || !iters || !gcds || !modmuls) {
        free(runtimes);
        free(iters);
        free(gcds);
        free(modmuls);
        return s;
    }

    for (size_t i = 0; i < count; i++) {
        switch (records[i].code) {
        case RESULT_FOUND:
            s.num_success++;
            break;
        case RESULT_NO_FACTOR:
        case RESULT_NO_BASE:
        case RESULT_WORK_LIMIT:
        case RESULT_FACTOR_LIMIT:
            s.num_timeout++;
            break;
        default:
            s.num_fail++;
            break;
        }
        runtimes[i] = records[i].runtime_ns;
        iters[i] = records[i].iterations;
        gcds[i] = records[i].gcd_calls;
        modmuls[i] = records[i].modmul_ops;
    }

    s.success_rate = (double)s.num_success / (double)s.num_trials;

    /* Sort arrays to compute percentiles */
    qsort(runtimes, count, sizeof(uint64_t), cmp_u64_asc);
    qsort(iters, count, sizeof(uint64_t), cmp_u64_asc);
    qsort(gcds, count, sizeof(uint64_t), cmp_u64_asc);
    qsort(modmuls, count, sizeof(uint64_t), cmp_u64_asc);

    s.min_runtime_ns = runtimes[0];
    s.max_runtime_ns = runtimes[count - 1];
    s.median_runtime_ns = runtimes[count / 2];
    s.q1_runtime_ns = runtimes[count / 4];
    s.q3_runtime_ns = runtimes[(3 * count) / 4];

    s.median_iterations = iters[count / 2];
    s.median_gcd_calls = gcds[count / 2];
    s.median_modmul_ops = modmuls[count / 2];

    free(runtimes);
    free(iters);
    free(gcds);
    free(modmuls);

    return s;
}

/* ------------------------------------------------------------------ */
/* Suite Runner                                                       */
/* ------------------------------------------------------------------ */

BenchSuite run_suite(const BenchOptions *opts_in)
{
    BenchOptions opts = opts_in ? *opts_in : default_suite_options(false);

    BenchSuite res;
    memset(&res, 0, sizeof(res));

    size_t case_count = 0;
    const BenchCase *cases = get_dataset(opts.smoke_only, &case_count);
    if (!cases || case_count == 0)
        return res;

    /* Capacity: Fermat (1) + Rho (rho_trials) + p-1 (1) per test case */
    size_t trials_per_case = 1 + opts.rho_trials + 1;
    res.record_cap = case_count * trials_per_case;
    res.records = (BenchRecord *)malloc(res.record_cap * sizeof(BenchRecord));

    res.summary_cap = case_count * BENCH_ALGO_COUNT;
    res.summaries = (BenchSummary *)malloc(res.summary_cap * sizeof(BenchSummary));

    if (!res.records || !res.summaries) {
        free(res.records);
        free(res.summaries);
        memset(&res, 0, sizeof(res));
        return res;
    }

    uint64_t trial_id = 1;

    for (size_t c = 0; c < case_count; c++) {
        const BenchCase *tc = &cases[c];

        /* --- 1. Fermat --- */
        BenchParams fp = default_params(BENCH_ALGO_FERMAT);
        if (opts.max_iterations > 0) fp.max_iterations = opts.max_iterations;
        size_t fermat_start = res.record_count;
        res.records[res.record_count++] = run_trial(trial_id++, tc, BENCH_ALGO_FERMAT, &fp);
        res.summaries[res.summary_count++] = calculate_summary(&res.records[fermat_start], 1);

        /* --- 2. Pollard Rho (Multiple trials with controlled seeds) --- */
        BenchParams rp = default_params(BENCH_ALGO_RHO);
        if (opts.max_iterations > 0) rp.max_iterations = opts.max_iterations;
        size_t rho_start = res.record_count;
        for (uint64_t r = 0; r < opts.rho_trials; r++) {
            rp.seed = r + 1; /* Reproducible seeds 1, 2, 3, ... */
            res.records[res.record_count++] = run_trial(trial_id++, tc, BENCH_ALGO_RHO, &rp);
        }
        res.summaries[res.summary_count++] = calculate_summary(&res.records[rho_start], opts.rho_trials);

        /* --- 3. Pollard p-1 --- */
        BenchParams pp = default_params(BENCH_ALGO_PMINUS1);
        if (opts.smoothness_bound > 0) pp.smoothness_bound = opts.smoothness_bound;
        size_t pm1_start = res.record_count;
        res.records[res.record_count++] = run_trial(trial_id++, tc, BENCH_ALGO_PMINUS1, &pp);
        res.summaries[res.summary_count++] = calculate_summary(&res.records[pm1_start], 1);
    }

    return res;
}

void free_suite(BenchSuite *res)
{
    if (!res) return;
    free(res->records);
    free(res->summaries);
    res->records = NULL;
    res->summaries = NULL;
    res->record_count = res->record_cap = 0;
    res->summary_count = res->summary_cap = 0;
}

/* ------------------------------------------------------------------ */
/* CSV Serialization Helper with RFC 4180 Escaping                    */
/* ------------------------------------------------------------------ */

static void print_csv_field(FILE *out, const char *s)
{
    if (s == NULL) {
        fputc(',', out);
        return;
    }
    bool need_quotes = false;
    for (const char *p = s; *p; p++) {
        if (*p == ',' || *p == '"' || *p == '\n' || *p == '\r') {
            need_quotes = true;
            break;
        }
    }
    if (!need_quotes) {
        fputs(s, out);
    } else {
        fputc('"', out);
        for (const char *p = s; *p; p++) {
            if (*p == '"') {
                fputc('"', out);
            }
            fputc(*p, out);
        }
        fputc('"', out);
    }
}

bool export_raw_csv(const BenchRecord *records,
                              size_t count,
                              FILE *out)
{
    if (!records || !out)
        return false;

    /* Header */
    fprintf(out,
            "trial_id,dataset_id,category,n,bit_length,ground_truth_p,ground_truth_q,"
            "algorithm,param_iters,param_retries,param_seed,param_b,param_bases,param_max_modmul,"
            "status,reason,factor_found,cofactor_found,factor_is_valid,runtime_ns,iterations,gcd_calls,modmul_ops\n");

    for (size_t i = 0; i < count; i++) {
        const BenchRecord *r = &records[i];
        fprintf(out, "%" U64_FMT ",", r->trial_id);
        print_csv_field(out, r->dataset_id);
        fputc(',', out);
        print_csv_field(out, r->category);
        fputc(',', out);
        fprintf(out, "%" U64_FMT ",%u,%" U64_FMT ",%" U64_FMT ",",
                r->n, r->bit_length, r->ground_truth_p, r->ground_truth_q);
        print_csv_field(out, algo_name(r->algo));
        fputc(',', out);
        fprintf(out, "%" U64_FMT ",%" U64_FMT ",%" U64_FMT ",%" U64_FMT ",%" U64_FMT ",%" U64_FMT ",",
                r->params.max_iterations, r->params.max_retries, r->params.seed,
                r->params.smoothness_bound, r->params.max_bases, r->params.max_modmul_ops);
        print_csv_field(out, status_string(r->code));
        fputc(',', out);
        print_csv_field(out, reason_string(r->code));
        fputc(',', out);
        fprintf(out, "%" U64_FMT ",%" U64_FMT ",%s,%" U64_FMT ",%" U64_FMT ",%" U64_FMT ",%" U64_FMT "\n",
                r->factor_found, r->cofactor_found, r->factor_is_valid ? "true" : "false",
                r->runtime_ns, r->iterations, r->gcd_calls, r->modmul_ops);
    }

    return true;
}

bool export_summary_csv(const BenchSummary *summaries,
                                  size_t count,
                                  FILE *out)
{
    if (!summaries || !out)
        return false;

    /* Header */
    fprintf(out,
            "dataset_id,category,algorithm,num_trials,num_success,num_timeout,num_fail,success_rate,"
            "median_runtime_ns,min_runtime_ns,max_runtime_ns,q1_runtime_ns,q3_runtime_ns,"
            "median_iterations,median_gcd_calls,median_modmul_ops\n");

    for (size_t i = 0; i < count; i++) {
        const BenchSummary *s = &summaries[i];
        print_csv_field(out, s->dataset_id);
        fputc(',', out);
        print_csv_field(out, s->category);
        fputc(',', out);
        print_csv_field(out, algo_name(s->algo));
        fputc(',', out);
        fprintf(out, "%" U64_FMT ",%" U64_FMT ",%" U64_FMT ",%" U64_FMT ",%.4f,"
                     "%" U64_FMT ",%" U64_FMT ",%" U64_FMT ",%" U64_FMT ",%" U64_FMT ","
                     "%" U64_FMT ",%" U64_FMT ",%" U64_FMT "\n",
                s->num_trials, s->num_success, s->num_timeout, s->num_fail, s->success_rate,
                s->median_runtime_ns, s->min_runtime_ns, s->max_runtime_ns, s->q1_runtime_ns, s->q3_runtime_ns,
                s->median_iterations, s->median_gcd_calls, s->median_modmul_ops);
    }

    return true;
}
