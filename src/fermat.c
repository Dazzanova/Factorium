/*
 * factorium/src/fermat.c — Fermat's difference-of-squares method.
 *
 * a^2 is tracked in unsigned __int128 because it can be larger than 64
 * bits (a starts near sqrt(n), but the search may go past 2^32). On compilers
 * without that extension the loop stays below 2^32 and reports a work limit.
 */
#include "factorium/fermat.h"
#include "factorium/modular_arithmetic.h"

#if defined(__SIZEOF_INT128__)

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"

typedef unsigned __int128 u128;

/* floor(sqrt(x)) for a 128-bit value. */
static u128 isqrt128(u128 x)
{
    u128 r;
    if (x < 2)
        return x;
    r = x / 2 + 1;
    for (;;) {
        u128 next = (r + x / r) / 2;
        if (next >= r)
            break;
        r = next;
    }
    while (r > x / r)
        r--;
    return r;
}

static FactorResult fermat_wide(uint64_t n, uint64_t budget, uint64_t a0)
{
    FactorResult res = {0};
    uint64_t a = a0;
    u128 b2 = (u128)a * a - n;
    uint64_t tried;

    res.code = RESULT_NO_FACTOR;

    for (tried = 0; tried < budget; tried++) {
        u128 b = isqrt128(b2);
        if (b * b == b2) {
            res.iterations = tried + 1;
            if (b < a) {
                uint64_t p = a - (uint64_t)b;
                if (p > 1) {
                    res.code = RESULT_FOUND;
                    res.factor = p;
                    res.cofactor = n / p;
                    return res;
                }
            }
            /* b == a means the factor is 1: n is prime, so stop. */
            break;
        }
        b2 += (u128)2 * a + 1; /* (a + 1)^2 - n */
        a++;
    }

    res.iterations = budget;
    return res;
}

#pragma GCC diagnostic pop

#endif /* __SIZEOF_INT128__ */

FactorResult fermat(uint64_t n, uint64_t max_iterations)
{
    FactorResult res = {0};
    uint64_t a;

    res.code = RESULT_NO_FACTOR;

    if (n < 4) {
        res.code = RESULT_INVALID;
        return res;
    }
    if ((n & 1) == 0) {
        res.code = RESULT_FOUND;
        res.factor = 2;
        res.cofactor = n / 2;
        return res;
    }
    if (max_iterations == 0)
        return res;

    a = isqrt(n);
    if (a * a < n) /* a <= 2^32 - 1, so a * a cannot overflow here */
        a++;

#if defined(__SIZEOF_INT128__)

    return fermat_wide(n, max_iterations, a);

#else

    for (uint64_t tried = 0; tried < max_iterations; tried++) {
        uint64_t b;
        res.iterations = tried + 1;
        if (a >= (1ULL << 32)) {
            res.code = RESULT_WORK_LIMIT;
            return res;
        }
        if (is_square(a * a - n, &b) && b < a) {
            uint64_t p = a - b;
            if (p > 1) {
                res.code = RESULT_FOUND;
                res.factor = p;
                res.cofactor = n / p;
                return res;
            }
            break;
        }
        a++;
    }
    return res;

#endif
}
