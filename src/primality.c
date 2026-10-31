/*
 * factorium/src/primality.c — deterministic Miller-Rabin for uint64_t.
 */
#include "factorium/primality.h"
#include "factorium/modular_arithmetic.h"

#include <stddef.h>

static const uint64_t WITNESSES[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
};

bool is_prime(uint64_t n)
{
    uint64_t d;
    int r;
    size_t i;

    if (n < 2)
        return false;

    /* Quick exact check for small numbers and small prime divisors. */
    for (i = 0; i < sizeof(WITNESSES) / sizeof(WITNESSES[0]); i++) {
        if (n == WITNESSES[i])
            return true;
        if (n % WITNESSES[i] == 0)
            return false;
    }

    /* Write n - 1 = d * 2^r with d odd. */
    d = n - 1;
    r = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        r++;
    }

    for (i = 0; i < sizeof(WITNESSES) / sizeof(WITNESSES[0]); i++) {
        uint64_t a = WITNESSES[i];
        uint64_t x = pow_mod(a, d, n, NULL);
        bool passed = false;

        if (x == 1 || x == n - 1)
            continue;
        for (int k = 1; k < r; k++) {
            x = mul_mod(x, x, n);
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
