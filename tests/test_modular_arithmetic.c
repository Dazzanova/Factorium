/*
 * tests/test_modular_arithmetic.c — gcd, modular mul/add/pow, isqrt,
 * is_square, checked multiplication. Hand-computed literals plus, where
 * __uint128_t is available, a cross-check of mod_mul against the compiler
 * extension.
 */
#include "factorium/modular_arithmetic.h"
#include "test_harness.h"

/* The cross-checks below use the documented __int128 extension directly. */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"

static void test_gcd(void)
{
    CHECK_EQ_U64(gcd(0, 0), 0);
    CHECK_EQ_U64(gcd(5, 0), 5);
    CHECK_EQ_U64(gcd(0, 7), 7);
    CHECK_EQ_U64(gcd(1, 1), 1);
    CHECK_EQ_U64(gcd(270, 192), 6);
    CHECK_EQ_U64(gcd(462, 105), 21);
    CHECK_EQ_U64(gcd(9223372036854775807ULL, 9223372036854775806ULL), 1);
    /* Consecutive Fibonacci numbers (F(89), F(90)) exercise the long chain. */
    CHECK_EQ_U64(gcd(1779979416004714189ULL, 2880067194370816120ULL), 1);
    CHECK_EQ_U64(gcd(4294967291ULL, 4294967291ULL), 4294967291ULL);
}

static void test_mod_mul_hand_computed(void)
{
    CHECK_EQ_U64(mul_mod(0, 12345, 97), 0);
    CHECK_EQ_U64(mul_mod(3, 4, 10), 2);
    CHECK_EQ_U64(mul_mod(5, 5, 7), 4);            /* 25 mod 7            */
    CHECK_EQ_U64(mul_mod(UINT64_MAX, UINT64_MAX, UINT64_MAX), 0); /* (m-1)^2 mod m */
    /* (2^64-1) mod (2^63+1) = m-3; (m-3)^2 mod m = 9 */
    CHECK_EQ_U64(mul_mod(UINT64_MAX, UINT64_MAX, UINT64_C(9223372036854775808) + 1), 9);
    /* 4000000007 = -2 (mod 4000000009) -> (-2)^2 = 4 */
    CHECK_EQ_U64(mul_mod(4000000007ULL, 4000000007ULL, 4000000009ULL), 4);
    CHECK_EQ_U64(mul_mod(4000000007ULL, 2ULL, 4000000009ULL), 4000000005ULL); /* -4 */
    CHECK_EQ_U64(mul_mod(0, 0, 1), 0);
    CHECK_EQ_U64(mul_mod(0, 5, 1), 0);
}

static void test_mod_mul_u128_crosscheck(void)
{
#if defined(__SIZEOF_INT128__)
    /* Fixed-seed LCG; verify against the 128-bit reference. */
    uint64_t s = 0x2545F4914F6CDD1DULL;
    for (int i = 0; i < 200000; i++) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        uint64_t m = 1 + (s >> 8) % UINT64_C(1000000007); /* 1..1e9+6 */
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        uint64_t a = s % m;
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        uint64_t b = s % m;
        uint64_t want = (uint64_t)(((unsigned __int128)a * b) % m);
        if (mul_mod(a, b, m) != want) {
            CHECK(false); /* report once; keep going for more coverage */
            break;
        }
        if (i == 199999) CHECK(true);
    }
#endif
    /* Near-2^64 cross-check (any implementation). */
    {
        /* 18446744073709551557 is prime (largest prime < 2^64);
         * 18446744073709551533 is prime as well. Product mod m below
         * is computed against the 128-bit reference when available. */
        uint64_t p = 18446744073709551557ULL;
        uint64_t q = 18446744073709551533ULL;
        uint64_t got = mul_mod(p, q, UINT64_C(1000000007));
#if defined(__SIZEOF_INT128__)
        uint64_t want = (uint64_t)(((unsigned __int128)p * (unsigned __int128)q) % 1000000007ULL);
        CHECK_EQ_U64(got, want);
#else
        (void)got;
#endif
    }
}

static void test_mod_add(void)
{
    CHECK_EQ_U64(add_mod(3, 4, 10), 7);
    CHECK_EQ_U64(add_mod(8, 9, 10), 7);          /* 17 mod 10           */
    CHECK_EQ_U64(add_mod(9, 9, 10), 8);          /* 18 mod 10           */
    CHECK_EQ_U64(add_mod(0, 0, 1), 0);
    /* Overflowing sum: a = b = 2^64-2, m = 2^64-1 -> (2^65-4) mod m = m-2 */
    CHECK_EQ_U64(add_mod(UINT64_MAX - 1, UINT64_MAX - 1, UINT64_MAX),
                 UINT64_MAX - 2);
    CHECK_EQ_U64(add_mod(123, 0, 1000), 123);
    CHECK_EQ_U64(add_mod(999, 999, 1000), 998);
}

static void test_mod_pow(void)
{
    CHECK_EQ_U64(pow_mod(0, 0, 7, NULL), 1);           /* convention: x^0 = 1 */
    CHECK_EQ_U64(pow_mod(5, 0, 7, NULL), 1);
    CHECK_EQ_U64(pow_mod(5, 0, 1, NULL), 0);           /* 1 mod 1 = 0         */
    CHECK_EQ_U64(pow_mod(2, 0, 1000, NULL), 1);
    CHECK_EQ_U64(pow_mod(2, 10, 1000, NULL), 24);
    CHECK_EQ_U64(pow_mod(3, 3, 5, NULL), 2);
    CHECK_EQ_U64(pow_mod(2, 63, 1000, NULL), 808);     /* 2^63 mod 1000      */
    /* Fermat: a^(p-1) == 1 mod p for prime p (base must not be 0 mod p). */
    CHECK_EQ_U64(pow_mod(5, 4000000006ULL, 4000000007ULL, NULL), 1);
    CHECK_EQ_U64(pow_mod(12345, 4294967290ULL, 4294967291ULL, NULL), 1);
    /* 3^13 = 1594323 -> 23 mod 100 */
    CHECK_EQ_U64(pow_mod(3, 13, 100, NULL), 23);
    /* base >= m is reduced: 10^3 mod 7 == 3^3 mod 7 == 6 */
    CHECK_EQ_U64(pow_mod(10, 3, 7, NULL), 6);
    /* exponent UINT64_MAX, small modulus, cross-checked below */
#if defined(__SIZEOF_INT128__)
    {
        uint64_t m = 1000000009ULL; /* prime */
        uint64_t got = pow_mod(123456789ULL, UINT64_MAX, m, NULL);
        /* Reference: 123456789^(2^64-1) mod p. 2^64 mod (p-1): compute. */
        uint64_t p1 = m - 1;
        uint64_t e2 = (uint64_t)(((unsigned __int128)1 << 64) % p1); /* = 2^64 mod (p-1) */
        uint64_t want = pow_mod(123456789ULL, e2 - 1, m, NULL); /* 2^64-1 mod (p-1) == e2-1 */
        CHECK_EQ_U64(got, want);
    }
#endif
}

static void test_mod_pow_counted(void)
{
    uint64_t ops = 99;
    pow_mod(3, 13, 100, &ops);                       /* bitlen 4, popcnt 3 */
    CHECK_EQ_U64(ops, 5);
    CHECK_EQ_U64(pow_mod(3, 13, 100, NULL), 23);
    ops = 99;
    pow_mod(7, 0, 100, &ops);
    CHECK_EQ_U64(ops, 0);
    ops = 99;
    pow_mod(7, 1, 100, &ops);                        /* bitlen 1, popcnt 1 */
    CHECK_EQ_U64(ops, 0);
    ops = 99;
    pow_mod(2, UINT64_MAX, 1000, &ops);              /* bitlen 64, popcnt 64 */
    CHECK_EQ_U64(ops, 63 + 63);
}

static void test_isqrt(void)
{
    CHECK_EQ_U64(isqrt(0), 0);
    CHECK_EQ_U64(isqrt(1), 1);
    CHECK_EQ_U64(isqrt(2), 1);
    CHECK_EQ_U64(isqrt(3), 1);
    CHECK_EQ_U64(isqrt(4), 2);
    CHECK_EQ_U64(isqrt(8), 2);
    CHECK_EQ_U64(isqrt(9), 3);
    /* (2^32-1)^2 = 18446744065119617025 */
    CHECK_EQ_U64(isqrt(18446744065119617025ULL), 4294967295ULL);
    CHECK_EQ_U64(isqrt(18446744065119617024ULL), 4294967294ULL);
    CHECK_EQ_U64(isqrt(18446744065119617026ULL), 4294967295ULL);
    /* 2^62 = (2^31)^2 */
    CHECK_EQ_U64(isqrt(4611686018427387904ULL), 2147483648ULL);
    CHECK_EQ_U64(isqrt(4611686018427387903ULL), 2147483647ULL);
    CHECK_EQ_U64(isqrt(4611686018427387905ULL), 2147483648ULL);
    /* 65536^2 = 4294967296 */
    CHECK_EQ_U64(isqrt(4294967296ULL), 65536ULL);
    CHECK_EQ_U64(isqrt(4294967297ULL), 65536ULL);
    CHECK_EQ_U64(isqrt(4294967295ULL), 65535ULL);
    /* top of the range */
    CHECK_EQ_U64(isqrt(UINT64_MAX), 4294967295ULL);
    CHECK_EQ_U64(isqrt(UINT64_MAX - 1), 4294967295ULL);
}

static void test_is_square(void)
{
    uint64_t r;
    CHECK(is_square(0, &r) && r == 0);
    CHECK(is_square(1, &r) && r == 1);
    CHECK(is_square(4, &r) && r == 2);
    CHECK(!is_square(2, &r));
    CHECK(!is_square(15, &r));
    CHECK(is_square(18446744065119617025ULL, &r) && r == 4294967295ULL);
    CHECK(is_square(4611686018427387904ULL, &r) && r == 2147483648ULL);
    CHECK(!is_square(UINT64_MAX, &r));
    CHECK(!is_square(18446744065119617026ULL, &r));
}

static void test_mul_checked(void)
{
    uint64_t out = 0;
    CHECK(mul_checked(0, UINT64_MAX, &out) && out == 0);
    CHECK(mul_checked(2, 3, &out) && out == 6);
    CHECK(mul_checked(UINT64_MAX, 1, &out) && out == UINT64_MAX);
    CHECK(!mul_checked(UINT64_MAX, 2, &out));
    CHECK(!mul_checked(4294967296ULL, 4294967296ULL, &out)); /* 2^64   */
    CHECK(mul_checked(4294967295ULL, 4294967295ULL, &out)
          && out == 18446744065119617025ULL);
    CHECK(!mul_checked(1000000000000000000ULL, 1000000000000000000ULL, &out));
    CHECK(mul_checked(4000000007ULL, 4000000007ULL, &out)
          && out == 16000000056000000049ULL);
}

int main(void)
{
    test_gcd();
    test_mod_mul_hand_computed();
    test_mod_mul_u128_crosscheck();
    test_mod_add();
    test_mod_pow();
    test_mod_pow_counted();
    test_isqrt();
    test_is_square();
    test_mul_checked();
    TEST_MAIN_END();
}
#pragma GCC diagnostic pop
