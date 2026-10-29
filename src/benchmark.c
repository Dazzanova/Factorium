/*
 * factorium/src/benchmark.c — Implementation of benchmarking harness and datasets.
 */
#define _POSIX_C_SOURCE 199309L

#include "factorium/benchmark.h"
#include "factorium/modular_arithmetic.h"
#include "factorium/primality.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define U64_FMT PRIu64

/* ------------------------------------------------------------------ */
/* Datasets with Ground-Truth Factorizations and Hypotheses           */
/* ------------------------------------------------------------------ */

static const fct_bench_case ALL_DATASETS[] = {
    /* --- 1. Close Semiprimes (Fermat Favorable) --- */
    {
        "close_16_twin", "close_semiprimes",
        64507ULL, 16, 251ULL, 257ULL,
        "16-bit twin-like semiprime (delta=6)",
        "Fermat finds factor in 3 iterations; dominates Rho and p-1"
    },
    {
        "close_32_twin", "close_semiprimes",
        4295229443ULL, 33, 65537ULL, 65539ULL,
        "33-bit twin prime product (delta=2)",
        "Fermat finishes in 1 iteration; Rho and p-1 require significantly more work"
    },
    {
        "close_32_d16", "close_semiprimes",
        4294049777ULL, 32, 65521ULL, 65537ULL,
        "32-bit close semiprime (delta=16)",
        "Fermat requires 8 iterations; extremely fast"
    },
    {
        "close_40_d30", "close_semiprimes",
        1000036000099ULL, 40, 1000003ULL, 1000033ULL,
        "40-bit close semiprime (delta=30)",
        "Fermat completes in 15 iterations; highly favorable compared to Rho"
    },
    {
        "close_40_d10k", "close_semiprimes",
        1010006030009ULL, 40, 1000003ULL, 1010003ULL,
        "40-bit semiprime with delta=10000",
        "Fermat requires 5000 iterations; still competitive"
    },
    {
        "close_63_d12", "close_semiprimes",
        4611686039902224373ULL, 63, 2147483647ULL, 2147483659ULL,
        "63-bit semiprime near 2^62 with delta=12",
        "Fermat succeeds in 6 iterations on 63-bit input; u128 prevents overflow"
    },
    {
        "close_64_d14", "close_semiprimes",
        9223372037000249951ULL, 64, 3037000493ULL, 3037000507ULL,
        "64-bit semiprime with delta=14",
        "Fermat finishes in 7 iterations on large 64-bit composite"
    },

    /* --- 2. Distant Semiprimes (Pollard Rho Favorable) --- */
    {
        "distant_32_p257", "distant_semiprimes",
        4294905101ULL, 32, 257ULL, 16711693ULL,
        "32-bit semiprime with small factor p=257 (8-bit)",
        "Pollard Rho finds factor p in O(sqrt(p)) ~ 16 steps; Fermat struggles"
    },
    {
        "distant_33_p1009", "distant_semiprimes",
        4294981039ULL, 33, 1009ULL, 4256671ULL,
        "33-bit semiprime with p=1009 and q=4256671",
        "Rho discovers p in ~30-50 iterations; Fermat requires millions"
    },
    {
        "distant_42_p17", "distant_semiprimes",
        4285662610591ULL, 42, 17ULL, 252097800623ULL,
        "42-bit semiprime with tiny factor 17",
        "Pollard Rho resolves factor in <10 steps; Fermat requires ~10^11 iterations"
    },
    {
        "distant_49_p104k", "distant_semiprimes",
        418916000733103ULL, 49, 104729ULL, 4000000007ULL,
        "49-bit semiprime with p=104729 (17 bits)",
        "Pollard Rho finds p in O(sqrt(p)) ~ 320 steps; Fermat exhausts 1M budget"
    },
    {
        "distant_63_p65k", "distant_semiprimes",
        9223372036854349817ULL, 63, 65537ULL, 140735340904441ULL,
        "63-bit semiprime with 16-bit factor p=65537",
        "Rho finds factor in ~256 steps; Fermat cannot succeed within practical budget"
    },

    /* --- 3. Composites Containing Small Prime Factors --- */
    {
        "small_fact_360", "small_prime_factors",
        360ULL, 9, 2ULL, 3ULL,
        "Highly composite number 360 = 2^3 * 3^2 * 5",
        "Even number split instantly by Fermat; Rho discovers small factor immediately"
    },
    {
        "small_fact_primorial", "small_prime_factors",
        30030ULL, 15, 2ULL, 3ULL,
        "Primorial 2 * 3 * 5 * 7 * 11 * 13",
        "Fermat handles even immediately; Rho splits small prime in 1-2 steps"
    },
    {
        "small_fact_power2_40", "small_prime_factors",
        1099511627776ULL, 41, 2ULL, 2ULL,
        "Power of 2: 2^40",
        "All algorithms detect factor 2 immediately"
    },
    {
        "small_fact_odd_smooth", "small_prime_factors",
        21523360575ULL, 35, 3ULL, 5ULL,
        "Odd 35-bit composite 3 * 5^2 * 1171 * 245071",
        "Rho discovers small factor 3 or 5 in 1-3 iterations; Fermat requires odd search"
    },

    /* --- 4. Inputs with Favorable p-1 Smoothness --- */
    {
        "p1_favorable_7681", "pminus1_favorable",
        7681023043ULL, 33, 7681ULL, 1000003ULL,
        "33-bit semiprime; p-1 = 7680 = 2^9 * 3 * 5 (5-smooth), q-1 = 2 * 3 * 166667 (166667 > 100k)",
        "Pollard p-1 stage 1 extracts p with B=100000 (d = p exactly, since q-1 is not B-smooth)"
    },
    {
        "p1_favorable_65537", "pminus1_favorable",
        65537000458759ULL, 46, 65537ULL, 1000000007ULL,
        "46-bit semiprime where p-1 = 65536 = 2^16 (2-smooth)",
        "Pollard p-1 stage 1 identifies p with bound B=2 in 1 step; Fermat fails"
    },
    {
        "p1_favorable_12289", "pminus1_favorable",
        1228900086023ULL, 41, 12289ULL, 100000007ULL,
        "41-bit semiprime where p-1 = 12288 = 2^12 * 3 (3-smooth)",
        "Pollard p-1 discovers factor p at bound B >= 3 with minimal modmuls"
    },
    {
        "p1_favorable_61441", "pminus1_favorable",
        61441000430087ULL, 46, 61441ULL, 1000000007ULL,
        "46-bit semiprime where p-1 = 61440 = 2^12 * 3 * 5 (5-smooth)",
        "Pollard p-1 succeeds with B=5; demonstrates power of smooth factorization"
    },

    /* --- 5. Inputs with Unfavorable p-1 Structure (Safe Primes) --- */
    {
        "p1_unfavorable_safe1", "pminus1_unfavorable",
        2001363019343ULL, 41, 1000667ULL, 2000029ULL,
        "Safe prime p=1000667 (p-1 has prime factor 500333 > 100k)",
        "Pollard p-1 stage 1 fails with B=100000; Rho succeeds in O(sqrt(p)) steps (varies by seed)"
    },
    {
        "p1_unfavorable_safe2", "pminus1_unfavorable",
        2002141278457ULL, 41, 1000919ULL, 2000303ULL,
        "Both factors safe primes: p-1 has prime factor 500459, q-1 has 1000151 (both > 100k)",
        "Pollard p-1 exhausts all bases without finding factor; Rho succeeds"
    },
    {
        "p1_unfavorable_safe3", "pminus1_unfavorable",
        4004091079237ULL, 42, 1001003ULL, 4000079ULL,
        "Both factors safe primes: p-1 has prime factor 500501, q-1 has 2000039 (both > 100k)",
        "Pollard p-1 exhausts bases; illustrates dependence on smoothness assumption"
    },

    /* --- 6. Bit-Length Scaling (Balanced Semiprimes) --- */
    {
        "scale_16_bit", "bit_length_scaling",
        64507ULL, 16, 251ULL, 257ULL,
        "16-bit balanced semiprime",
        "Base baseline for runtime and operation count scaling"
    },
    {
        "scale_24_bit", "bit_length_scaling",
        16777207ULL, 24, 4093ULL, 4099ULL,
        "24-bit balanced semiprime",
        "Rho scales with O(n^(1/4)) ~ 64 operations"
    },
    {
        "scale_32_bit", "bit_length_scaling",
        4294049777ULL, 32, 65521ULL, 65537ULL,
        "32-bit balanced semiprime",
        "Rho requires ~256 operations; Fermat depends strictly on factor difference"
    },
    {
        "scale_41_bit", "bit_length_scaling",
        1099515822059ULL, 41, 1048573ULL, 1048583ULL,
        "41-bit balanced semiprime",
        "Demonstrates runtime growth across 40-bit boundary"
    },
    {
        "scale_49_bit", "bit_length_scaling",
        281475647799167ULL, 49, 16777213ULL, 16777259ULL,
        "49-bit balanced semiprime",
        "Rho requires ~4096 operations; Fermat completes due to small delta"
    },
    {
        "scale_56_bit", "bit_length_scaling",
        72057579542413141ULL, 56, 268435399ULL, 268435459ULL,
        "56-bit balanced semiprime",
        "Evaluation of Rho and Fermat near 56 bits"
    },
    {
        "scale_60_bit", "bit_length_scaling",
        1152921470247108503ULL, 60, 1073741789ULL, 1073741827ULL,
        "60-bit balanced semiprime",
        "Evaluation of 60-bit scalability without integer overflow"
    },
    {
        "scale_63_bit", "bit_length_scaling",
        4611686039902224373ULL, 63, 2147483647ULL, 2147483659ULL,
        "63-bit balanced semiprime near uint64 limit",
        "Verifies full 63-bit execution with no arithmetic overflow"
    },

    /* --- 7. Edge Cases and Difficult Inputs --- */
    {
        "edge_square_32", "edge_cases",
        4295098369ULL, 33, 65537ULL, 65537ULL,
        "Perfect square of prime 65537^2",
        "Fermat detects perfect square on step 0 (a = sqrt(n), b = 0)"
    },
    {
        "edge_prime_32", "edge_cases",
        1000000007ULL, 30, 1000000007ULL, 0ULL,
        "Prime input 1000000007",
        "Algorithms must report no_factor/exhausted without hanging or false positive"
    },
    {
        "edge_carmichael_561", "edge_cases",
        561ULL, 10, 3ULL, 11ULL,
        "Carmichael number 561 = 3 * 11 * 17",
        "Fermat and Rho resolve without being misled by pseudoprimality"
    },
    {
        "edge_carmichael_41041", "edge_cases",
        41041ULL, 16, 7ULL, 11ULL,
        "Carmichael number 41041 = 7 * 11 * 13 * 41",
        "Tests composite splitting on multi-prime Carmichael number"
    },
    {
        "edge_max_uint64_minus_1", "edge_cases",
        18446744073709551614ULL, 64, 2ULL, 7ULL,
        "2^64 - 2 (even 64-bit integer)",
        "Verifies no overflow near UINT64_MAX boundary"
    }
};

static const fct_bench_case SMOKE_DATASETS[] = {
    {
        "close_16_twin", "close_semiprimes",
        64507ULL, 16, 251ULL, 257ULL,
        "16-bit twin-like semiprime (delta=6)",
        "Fermat finds factor in 3 iterations; dominates Rho and p-1"
    },
    {
        "close_32_twin", "close_semiprimes",
        4295229443ULL, 33, 65537ULL, 65539ULL,
        "33-bit twin prime product (delta=2)",
        "Fermat finishes in 1 iteration; Rho and p-1 require significantly more work"
    },
    {
        "close_40_d10k", "close_semiprimes",
        1010006030009ULL, 40, 1000003ULL, 1010003ULL,
        "40-bit semiprime with delta=10000",
        "Fermat requires 5000 iterations; still competitive"
    },
    {
        "close_63_d12", "close_semiprimes",
        4611686039902224373ULL, 63, 2147483647ULL, 2147483659ULL,
        "63-bit semiprime near 2^62 with delta=12",
        "Fermat succeeds in 6 iterations on 63-bit input; u128 prevents overflow"
    },
    {
        "distant_32_p257", "distant_semiprimes",
        4294905101ULL, 32, 257ULL, 16711693ULL,
        "32-bit semiprime with small factor p=257 (8-bit)",
        "Pollard Rho finds factor p in O(sqrt(p)) ~ 16 steps; Fermat struggles"
    },
    {
        "distant_33_p1009", "distant_semiprimes",
        4294981039ULL, 33, 1009ULL, 4256671ULL,
        "33-bit semiprime with p=1009 and q=4256671",
        "Rho discovers p in ~30-50 iterations; Fermat requires millions"
    },
    {
        "distant_49_p104k", "distant_semiprimes",
        418916000733103ULL, 49, 104729ULL, 4000000007ULL,
        "49-bit semiprime with p=104729 (17 bits)",
        "Pollard Rho finds p in O(sqrt(p)) ~ 320 steps; Fermat exhausts 1M budget"
    },
    {
        "small_fact_360", "small_prime_factors",
        360ULL, 9, 2ULL, 3ULL,
        "Highly composite number 360 = 2^3 * 3^2 * 5",
        "Even number split instantly by Fermat; Rho discovers small factor immediately"
    },
    {
        "small_fact_power2_40", "small_prime_factors",
        1099511627776ULL, 41, 2ULL, 2ULL,
        "Power of 2: 2^40",
        "All algorithms detect factor 2 immediately"
    },
    {
        "p1_favorable_7681", "pminus1_favorable",
        7681023043ULL, 33, 7681ULL, 1000003ULL,
        "33-bit semiprime; p-1 = 7680 = 2^9 * 3 * 5 (5-smooth), q-1 = 2 * 3 * 166667 (166667 > 100k)",
        "Pollard p-1 stage 1 extracts p with B=100000 (d = p exactly, since q-1 is not B-smooth)"
    },
    {
        "p1_favorable_65537", "pminus1_favorable",
        65537000458759ULL, 46, 65537ULL, 1000000007ULL,
        "46-bit semiprime where p-1 = 65536 = 2^16 (2-smooth)",
        "Pollard p-1 stage 1 identifies p with bound B=2 in 1 step; Fermat fails"
    },
    {
        "p1_unfavorable_safe1", "pminus1_unfavorable",
        2001363019343ULL, 41, 1000667ULL, 2000029ULL,
        "Safe prime p=1000667 (p-1 has prime factor 500333 > 100k)",
        "Pollard p-1 stage 1 fails with B=100000; Rho succeeds in O(sqrt(p)) steps (varies by seed)"
    },
    {
        "p1_unfavorable_safe3", "pminus1_unfavorable",
        4004091079237ULL, 42, 1001003ULL, 4000079ULL,
        "Both factors safe primes: p-1 has prime factor 500501, q-1 has 2000039 (both > 100k)",
        "Pollard p-1 exhausts bases; illustrates dependence on smoothness assumption"
    },
    {
        "scale_16_bit", "bit_length_scaling",
        64507ULL, 16, 251ULL, 257ULL,
        "16-bit balanced semiprime",
        "Base baseline for runtime and operation count scaling"
    },
    {
        "scale_32_bit", "bit_length_scaling",
        4294049777ULL, 32, 65521ULL, 65537ULL,
        "32-bit balanced semiprime",
        "Rho requires ~256 operations; Fermat depends strictly on factor difference"
    },
    {
        "scale_49_bit", "bit_length_scaling",
        281475647799167ULL, 49, 16777213ULL, 16777259ULL,
        "49-bit balanced semiprime",
        "Rho requires ~4096 operations; Fermat completes due to small delta"
    },
    {
        "scale_60_bit", "bit_length_scaling",
        1152921470247108503ULL, 60, 1073741789ULL, 1073741827ULL,
        "60-bit balanced semiprime",
        "Evaluation of 60-bit scalability without integer overflow"
    },
    {
        "edge_square_32", "edge_cases",
        4295098369ULL, 33, 65537ULL, 65537ULL,
        "Perfect square of prime 65537^2",
        "Fermat detects perfect square on step 0 (a = sqrt(n), b = 0)"
    },
    {
        "edge_prime_32", "edge_cases",
        1000000007ULL, 30, 1000000007ULL, 0ULL,
        "Prime input 1000000007",
        "Algorithms must report no_factor/exhausted without hanging or false positive"
    },
    {
        "edge_carmichael_561", "edge_cases",
        561ULL, 10, 3ULL, 11ULL,
        "Carmichael number 561 = 3 * 11 * 17",
        "Fermat and Rho resolve without being misled by pseudoprimality"
    },
    {
        "edge_max_uint64_minus_1", "edge_cases",
        18446744073709551614ULL, 64, 2ULL, 7ULL,
        "2^64 - 2 (even 64-bit integer)",
        "Verifies no overflow near UINT64_MAX boundary"
    }
};

static const size_t TOTAL_DATASETS = sizeof(ALL_DATASETS) / sizeof(ALL_DATASETS[0]);
static const size_t TOTAL_SMOKE_DATASETS = sizeof(SMOKE_DATASETS) / sizeof(SMOKE_DATASETS[0]);

/* ------------------------------------------------------------------ */
/* Algorithm Names & Parameters                                       */
/* ------------------------------------------------------------------ */

const char *fct_bench_algo_name(fct_bench_algo algo)
{
    switch (algo) {
    case FCT_BENCH_ALGO_FERMAT: return "fermat";
    case FCT_BENCH_ALGO_RHO: return "rho";
    case FCT_BENCH_ALGO_PMINUS1: return "pminus1";
    default: return "unknown";
    }
}

fct_bench_params fct_bench_params_default(fct_bench_algo algo)
{
    fct_bench_params p;
    memset(&p, 0, sizeof(p));
    switch (algo) {
    case FCT_BENCH_ALGO_FERMAT:
        p.max_iterations = 1000000ULL; /* 10^6 iterations */
        break;
    case FCT_BENCH_ALGO_RHO:
        p.max_iterations = 1000000ULL;
        p.max_retries = 8;
        p.seed = 1;
        break;
    case FCT_BENCH_ALGO_PMINUS1:
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

fct_bench_suite_options fct_bench_suite_options_default(bool smoke_only)
{
    fct_bench_suite_options o;
    o.smoke_only = smoke_only;
    o.rho_trials = smoke_only ? 1 : 5;
    o.max_iterations = smoke_only ? 20000ULL : 500000ULL;
    o.smoothness_bound = smoke_only ? 10000ULL : 100000ULL;
    return o;
}

/* ------------------------------------------------------------------ */
/* Dataset Access                                                     */
/* ------------------------------------------------------------------ */

size_t fct_bench_get_dataset_count(bool smoke_only)
{
    return smoke_only ? TOTAL_SMOKE_DATASETS : TOTAL_DATASETS;
}

const fct_bench_case *fct_bench_get_dataset(bool smoke_only, size_t *count_out)
{
    if (smoke_only) {
        if (count_out) *count_out = TOTAL_SMOKE_DATASETS;
        return SMOKE_DATASETS;
    }
    if (count_out) *count_out = TOTAL_DATASETS;
    return ALL_DATASETS;
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

fct_bench_record fct_bench_run_trial(uint64_t trial_id,
                                     const fct_bench_case *tc,
                                     fct_bench_algo algo,
                                     const fct_bench_params *params_in)
{
    fct_bench_record rec;
    memset(&rec, 0, sizeof(rec));

    rec.trial_id = trial_id;
    snprintf(rec.dataset_id, sizeof(rec.dataset_id), "%s", tc->id);
    snprintf(rec.category, sizeof(rec.category), "%s", tc->category);
    rec.n = tc->n;
    rec.bit_length = tc->bit_length;
    rec.ground_truth_p = tc->known_factor_p;
    rec.ground_truth_q = tc->known_factor_q;
    rec.algo = algo;

    fct_bench_params p = params_in ? *params_in : fct_bench_params_default(algo);
    rec.params = p;

    fct_factor_result res;
    memset(&res, 0, sizeof(res));

    uint64_t t0 = 0, t1 = 0;

    switch (algo) {
    case FCT_BENCH_ALGO_FERMAT: {
        fct_fermat_config cfg;
        cfg.max_iterations = p.max_iterations ? p.max_iterations : FCT_FERMAT_DEFAULT_MAX_ITERATIONS;
        t0 = get_time_ns();
        res = fct_fermat(tc->n, &cfg);
        t1 = get_time_ns();
        break;
    }
    case FCT_BENCH_ALGO_RHO: {
        fct_rho_config cfg;
        cfg.max_iterations_per_attempt = p.max_iterations ? p.max_iterations : FCT_RHO_DEFAULT_MAX_ITERATIONS;
        cfg.max_retries = p.max_retries ? p.max_retries : FCT_RHO_DEFAULT_MAX_RETRIES;
        cfg.seed = p.seed ? p.seed : FCT_RHO_DEFAULT_SEED;
        t0 = get_time_ns();
        res = fct_pollard_rho(tc->n, &cfg);
        t1 = get_time_ns();
        break;
    }
    case FCT_BENCH_ALGO_PMINUS1: {
        fct_pminus1_config cfg;
        cfg.smoothness_bound = p.smoothness_bound ? p.smoothness_bound : FCT_PM1_DEFAULT_SMOOTHNESS_BOUND;
        cfg.start_base = p.start_base ? p.start_base : FCT_PM1_DEFAULT_START_BASE;
        cfg.max_bases = p.max_bases ? p.max_bases : FCT_PM1_DEFAULT_MAX_BASES;
        cfg.max_modmul_ops = p.max_modmul_ops ? p.max_modmul_ops : FCT_PM1_DEFAULT_MAX_MODMUL_OPS;
        t0 = get_time_ns();
        res = fct_pollard_p_minus_one(tc->n, &cfg);
        t1 = get_time_ns();
        break;
    }
    default:
        res.status = FCT_INVALID_INPUT;
        res.reason = FCT_TERM_INVALID_INPUT;
        break;
    }

    rec.runtime_ns = (t1 >= t0) ? (t1 - t0) : 0;
    rec.status = res.status;
    rec.reason = res.reason;
    rec.factor_found = res.factor;
    rec.cofactor_found = res.cofactor;
    rec.iterations = res.iterations;
    rec.gcd_calls = res.gcd_calls;
    rec.modmul_ops = res.modmul_ops;

    /* Validate found factor */
    rec.factor_is_valid = false;
    if (res.status == FCT_OK && res.factor > 1 && res.factor < tc->n) {
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

fct_bench_summary fct_bench_calculate_summary(const fct_bench_record *records,
                                              size_t count)
{
    fct_bench_summary s;
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
        if (records[i].status == FCT_OK) {
            s.num_success++;
        } else if (records[i].status == FCT_NO_FACTOR || records[i].status == FCT_WORK_LIMIT) {
            s.num_timeout++;
        } else {
            s.num_fail++;
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

fct_bench_suite_result fct_bench_run_suite(const fct_bench_suite_options *opts_in)
{
    fct_bench_suite_options opts = opts_in ? *opts_in : fct_bench_suite_options_default(false);

    fct_bench_suite_result res;
    memset(&res, 0, sizeof(res));

    size_t case_count = 0;
    const fct_bench_case *cases = fct_bench_get_dataset(opts.smoke_only, &case_count);
    if (!cases || case_count == 0)
        return res;

    /* Capacity: Fermat (1) + Rho (rho_trials) + p-1 (1) per test case */
    size_t trials_per_case = 1 + opts.rho_trials + 1;
    res.record_cap = case_count * trials_per_case;
    res.records = (fct_bench_record *)malloc(res.record_cap * sizeof(fct_bench_record));

    res.summary_cap = case_count * FCT_BENCH_ALGO_COUNT;
    res.summaries = (fct_bench_summary *)malloc(res.summary_cap * sizeof(fct_bench_summary));

    if (!res.records || !res.summaries) {
        free(res.records);
        free(res.summaries);
        memset(&res, 0, sizeof(res));
        return res;
    }

    uint64_t trial_id = 1;

    for (size_t c = 0; c < case_count; c++) {
        const fct_bench_case *tc = &cases[c];

        /* --- 1. Fermat --- */
        fct_bench_params fp = fct_bench_params_default(FCT_BENCH_ALGO_FERMAT);
        if (opts.max_iterations > 0) fp.max_iterations = opts.max_iterations;
        size_t fermat_start = res.record_count;
        res.records[res.record_count++] = fct_bench_run_trial(trial_id++, tc, FCT_BENCH_ALGO_FERMAT, &fp);
        res.summaries[res.summary_count++] = fct_bench_calculate_summary(&res.records[fermat_start], 1);

        /* --- 2. Pollard Rho (Multiple trials with controlled seeds) --- */
        fct_bench_params rp = fct_bench_params_default(FCT_BENCH_ALGO_RHO);
        if (opts.max_iterations > 0) rp.max_iterations = opts.max_iterations;
        size_t rho_start = res.record_count;
        for (uint64_t r = 0; r < opts.rho_trials; r++) {
            rp.seed = r + 1; /* Reproducible seeds 1, 2, 3, ... */
            res.records[res.record_count++] = fct_bench_run_trial(trial_id++, tc, FCT_BENCH_ALGO_RHO, &rp);
        }
        res.summaries[res.summary_count++] = fct_bench_calculate_summary(&res.records[rho_start], opts.rho_trials);

        /* --- 3. Pollard p-1 --- */
        fct_bench_params pp = fct_bench_params_default(FCT_BENCH_ALGO_PMINUS1);
        if (opts.smoothness_bound > 0) pp.smoothness_bound = opts.smoothness_bound;
        size_t pm1_start = res.record_count;
        res.records[res.record_count++] = fct_bench_run_trial(trial_id++, tc, FCT_BENCH_ALGO_PMINUS1, &pp);
        res.summaries[res.summary_count++] = fct_bench_calculate_summary(&res.records[pm1_start], 1);
    }

    return res;
}

void fct_bench_suite_result_free(fct_bench_suite_result *res)
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

bool fct_bench_export_raw_csv(const fct_bench_record *records,
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
        const fct_bench_record *r = &records[i];
        fprintf(out, "%" U64_FMT ",", r->trial_id);
        print_csv_field(out, r->dataset_id);
        fputc(',', out);
        print_csv_field(out, r->category);
        fputc(',', out);
        fprintf(out, "%" U64_FMT ",%u,%" U64_FMT ",%" U64_FMT ",",
                r->n, r->bit_length, r->ground_truth_p, r->ground_truth_q);
        print_csv_field(out, fct_bench_algo_name(r->algo));
        fputc(',', out);
        fprintf(out, "%" U64_FMT ",%" U64_FMT ",%" U64_FMT ",%" U64_FMT ",%" U64_FMT ",%" U64_FMT ",",
                r->params.max_iterations, r->params.max_retries, r->params.seed,
                r->params.smoothness_bound, r->params.max_bases, r->params.max_modmul_ops);
        print_csv_field(out, fct_status_str(r->status));
        fputc(',', out);
        print_csv_field(out, fct_termination_str(r->reason));
        fputc(',', out);
        fprintf(out, "%" U64_FMT ",%" U64_FMT ",%s,%" U64_FMT ",%" U64_FMT ",%" U64_FMT ",%" U64_FMT "\n",
                r->factor_found, r->cofactor_found, r->factor_is_valid ? "true" : "false",
                r->runtime_ns, r->iterations, r->gcd_calls, r->modmul_ops);
    }

    return true;
}

bool fct_bench_export_summary_csv(const fct_bench_summary *summaries,
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
        const fct_bench_summary *s = &summaries[i];
        print_csv_field(out, s->dataset_id);
        fputc(',', out);
        print_csv_field(out, s->category);
        fputc(',', out);
        print_csv_field(out, fct_bench_algo_name(s->algo));
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
