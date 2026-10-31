/*
 * tests/test_primality.c — deterministic Miller-Rabin edge cases.
 *
 * Large literals below were cross-checked at development time against an
 * independent Python Miller-Rabin implementation; 2^63-1 uses its known
 * factorization 7^2 * 73 * 127 * 337 * 92737 * 649657.
 */
#include "factorium/primality.h"
#include "factorium/modular_arithmetic.h"
#include "test_harness.h"

static void test_tiny(void)
{
    CHECK(is_prime(2));
    CHECK(is_prime(3));
    CHECK(!is_prime(0));
    CHECK(!is_prime(1));
    CHECK(!is_prime(4));
    CHECK(!is_prime(9));
    CHECK(is_prime(5));
    CHECK(is_prime(7));
    CHECK(!is_prime(8));
}

static void test_small_primes_and_composites(void)
{
    CHECK(is_prime(97));
    CHECK(is_prime(7919));
    CHECK(is_prime(1000000007ULL));
    CHECK(is_prime(1000000009ULL));
    CHECK(is_prime(998244353ULL));
    CHECK(is_prime(2147483647ULL));          /* M31                  */
    CHECK(is_prime(4000000007ULL));
    CHECK(is_prime(4294967291ULL));           /* largest 32-bit prime */
    CHECK(is_prime(2305843009213693951ULL));  /* M61                  */

    CHECK(!is_prime(97 * 101));
    uint64_t pq;
    CHECK(mul_checked(1000000007ULL, 1000000009ULL, &pq));
    CHECK(!is_prime(pq));
    CHECK(!is_prime(4000000007ULL * 4000000007ULL)); /* = 16000000056000000049 */
}

static void test_pseudoprimes(void)
{
    /* Carmichael numbers fool Fermat tests. */
    CHECK(!is_prime(561));
    CHECK(!is_prime(1105));
    CHECK(!is_prime(41041));
    CHECK(!is_prime(825265));
    /* Strong pseudoprimes to base 2. */
    CHECK(!is_prime(2047));
    CHECK(!is_prime(1373653));
}

static void test_near_64bit(void)
{
    /* Largest prime below 2^64. */
    CHECK(is_prime(18446744073709551557ULL));
    CHECK(is_prime(18446744073709551533ULL));
    /* 2^63 - 1 = 7^2 * 73 * 127 * 337 * 92737 * 649657 */
    CHECK(!is_prime(9223372036854775807ULL));
    /* 2^64 - 1 = 3 * 5 * 17 * 257 * 641 * 65537 * 6700417 */
    CHECK(!is_prime(UINT64_MAX));
    /* 2^64 - 2, even and huge. */
    CHECK(!is_prime(UINT64_MAX - 1));
    /* Odd neighbor of a huge prime is composite. */
    CHECK(!is_prime(18446744073709551555ULL));
    /* 64-bit semiprime with both factors ~ 2^32. */
    uint64_t n;
    CHECK(mul_checked(4294967291ULL, 4000000007ULL, &n));
    CHECK(!is_prime(n));
}

int main(void)
{
    test_tiny();
    test_small_primes_and_composites();
    test_pseudoprimes();
    test_near_64bit();
    TEST_MAIN_END();
}
