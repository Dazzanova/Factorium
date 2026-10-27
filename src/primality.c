/*
 * factorium/src/primality.c — deterministic Miller-Rabin for uint64_t.
 *
 * Witness set: first twelve primes. Proven deterministic for all
 * n < 318,665,857,834,031,151,167,461 (> 2^64), so the result is exact.
 */
#include "factorium/primality.h"
#include "factorium/modular_arithmetic.h"

static const uint64_t WITNESSES[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
};

bool fct_is_prime(uint64_t n)
{
    if (n < 2)
        return false;

    /* Small-prime pass: exact for the witnesses themselves and a cheap
     * filter; also guarantees n is odd and > 37 before the loop below. */
    for (size_t i = 0; i < sizeof(WITNESSES) / sizeof(WITNESSES[0]); i++) {
        uint64_t p = WITNESSES[i];
        if (n == p)
            return true;
        if (n % p == 0)
            return false;
    }

    /* Write n - 1 = d * 2^r with d odd. */
    uint64_t d = n - 1;
    int r = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        r++;
    }

    for (size_t i = 0; i < sizeof(WITNESSES) / sizeof(WITNESSES[0]); i++) {
        uint64_t a = WITNESSES[i]; /* a < n here, since n > 37 */
        uint64_t x = fct_mod_pow(a, d, n);
        if (x == 1 || x == n - 1)
            continue;
        bool passed = false;
        for (int i2 = 1; i2 < r; i2++) {
            x = fct_mod_mul(x, x, n);
            if (x == n - 1) {
                passed = true;
                break;
            }
        }
        if (!passed)
            return false;
    }
    return true;
}
