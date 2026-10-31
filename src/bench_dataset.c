/*
 * factorium/src/bench_dataset.c — the benchmark datasets.
 */
#include "factorium/benchmark.h"

#include <stddef.h>

static const BenchCase ALL_DATASETS[] = {
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

static const BenchCase SMOKE_DATASETS[] = {
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

size_t get_dataset_count(bool smoke_only)
{
    return smoke_only ? TOTAL_SMOKE_DATASETS : TOTAL_DATASETS;
}

const BenchCase *get_dataset(bool smoke_only, size_t *count_out)
{
    if (smoke_only) {
        if (count_out) *count_out = TOTAL_SMOKE_DATASETS;
        return SMOKE_DATASETS;
    }
    if (count_out) *count_out = TOTAL_DATASETS;
    return ALL_DATASETS;
}
