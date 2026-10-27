/*
 * factorium/src/modular_arithmetic.c — overflow-safe uint64 modular
 * arithmetic. See include/factorium/modular_arithmetic.h for the compiler
 * assumption regarding unsigned __int128.
 */
#include "factorium/modular_arithmetic.h"

uint64_t fct_gcd(uint64_t a, uint64_t b)
{
    while (b != 0) {
        uint64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

/* (a + b) mod m, 0 <= a, b < m, m >= 1.
 * Rewrites the comparison as a >= m - b so that a + b is never formed when
 * it would wrap: if a >= m - b the answer is a - (m - b) (= a + b - m). */
uint64_t fct_mod_add(uint64_t a, uint64_t b, uint64_t m)
{
    if (a >= m - b)
        return a - (m - b);
    return a + b;
}

#if defined(__SIZEOF_INT128__)

/* __int128 is a documented compiler extension (GCC/Clang); keep -Wpedantic
 * clean by silencing it for the lines that use it. */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"

uint64_t fct_mod_mul(uint64_t a, uint64_t b, uint64_t m)
{
    return (uint64_t)(((unsigned __int128)a * b) % m);
}

#pragma GCC diagnostic pop

#else /* portable double-and-add; requires 0 <= a, b < m */

uint64_t fct_mod_mul(uint64_t a, uint64_t b, uint64_t m)
{
    if (m == 1)
        return 0;
    uint64_t r = 0;
    while (b != 0) {
        if (b & 1)
            r = fct_mod_add(r, a, m);
        b >>= 1;
        a = fct_mod_add(a, a, m);
    }
    return r;
}

#endif

/* Bit-scan square-and-multiply, MSB to LSB.
 * For exp > 0 performs exactly (bit_length(exp) - 1) squarings plus
 * (popcount(exp) - 1) multiplications. */
static uint64_t mod_pow_core(uint64_t base, uint64_t exp, uint64_t m,
                             uint64_t *ops_out)
{
    if (m == 1) {
        if (ops_out)
            *ops_out = 0;
        return 0;
    }
    if (exp == 0) {
        if (ops_out)
            *ops_out = 0;
        return 1;
    }

    uint64_t b = base % m;
    int highest = 0;
    for (uint64_t e = exp >> 1; e != 0; e >>= 1)
        highest++;

    uint64_t result = b;      /* leading (MSB) bit */
    uint64_t ops = 0;
    for (int i = highest - 1; i >= 0; i--) {
        result = fct_mod_mul(result, result, m);
        ops++;
        if ((exp >> i) & 1) {
            result = fct_mod_mul(result, b, m);
            ops++;
        }
    }
    if (ops_out)
        *ops_out = ops;
    return result;
}

uint64_t fct_mod_pow(uint64_t base, uint64_t exp, uint64_t m)
{
    return mod_pow_core(base, exp, m, NULL);
}

void fct_mod_pow_counted(uint64_t base, uint64_t exp, uint64_t m,
                         uint64_t *modmul_ops_out)
{
    (void)mod_pow_core(base, exp, m, modmul_ops_out);
}

/* floor(sqrt(n)), pure integer Newton.
 * x0 = n/2 + 1 is a strict upper bound on sqrt(n) for n >= 2, and the
 * iterates stay in [floor(sqrt(n)), x0]; every x + n/x sum that is formed
 * is < 2^64 (for x >= 2^32 the sum is maximized at the start, x = 2^63,
 * where it equals 2^63 + 1). */
uint64_t fct_isqrt(uint64_t n)
{
    if (n < 2)
        return n;
    uint64_t x = n / 2 + 1;
    for (;;) {
        uint64_t y = (x + n / x) / 2;
        if (y >= x)
            break;
        x = y;
    }
    while (x > n / x) /* defensive; the loop above already yields floor */
        x--;
    return x;
}

bool fct_is_square(uint64_t n, uint64_t *root_out)
{
    uint64_t r = fct_isqrt(n);
    uint64_t sq;
    if (!fct_mul_checked(r, r, &sq) || sq != n)
        return false;
    if (root_out)
        *root_out = r;
    return true;
}

bool fct_mul_checked(uint64_t a, uint64_t b, uint64_t *out)
{
    if (a == 0 || b == 0) {
        *out = 0;
        return true;
    }
    if (a > UINT64_MAX / b)
        return false;
    *out = a * b;
    return true;
}
