/*
 * factorium/src/modular_arithmetic.c — overflow-safe modular arithmetic.
 */
#include "factorium/modular_arithmetic.h"

#include <stddef.h>

uint64_t gcd(uint64_t a, uint64_t b)
{
    while (b != 0) {
        uint64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

/* (a + b) mod m without ever forming a + b when it would wrap. */
uint64_t add_mod(uint64_t a, uint64_t b, uint64_t m)
{
    if (a >= m - b)
        return a - (m - b);
    return a + b;
}

#if defined(__SIZEOF_INT128__)

/* 128-bit product, then one modulo: the simplest correct way. */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"

uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t m)
{
    return (uint64_t)(((unsigned __int128)a * b) % m);
}

#pragma GCC diagnostic pop

#else

/* Portable fallback: add-and-double ("Russian peasant") multiplication. */
uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t m)
{
    uint64_t result = 0;
    if (m == 1)
        return 0;
    while (b != 0) {
        if (b & 1)
            result = add_mod(result, a, m);
        b >>= 1;
        a = add_mod(a, a, m);
    }
    return result;
}

#endif

/* Square-and-multiply, most significant bit first. */
uint64_t pow_mod(uint64_t base, uint64_t exp, uint64_t m, uint64_t *ops)
{
    uint64_t result;
    uint64_t count = 0;
    int highest = 0;

    if (m == 1) {
        if (ops != NULL)
            *ops = 0;
        return 0;
    }
    if (exp == 0) {
        if (ops != NULL)
            *ops = 0;
        return 1;
    }

    base %= m;
    for (uint64_t e = exp >> 1; e != 0; e >>= 1)
        highest++;

    result = base; /* the leading bit is always 1 */
    for (int i = highest - 1; i >= 0; i--) {
        result = mul_mod(result, result, m);
        count++;
        if ((exp >> i) & 1) {
            result = mul_mod(result, base, m);
            count++;
        }
    }
    if (ops != NULL)
        *ops = count;
    return result;
}

/* Integer Newton iteration: x becomes (x + n/x) / 2 until it stops moving. */
uint64_t isqrt(uint64_t n)
{
    uint64_t x;
    if (n < 2)
        return n;
    x = n / 2 + 1;
    for (;;) {
        uint64_t next = (x + n / x) / 2;
        if (next >= x)
            break;
        x = next;
    }
    while (x > n / x)
        x--;
    return x;
}

bool is_square(uint64_t n, uint64_t *root)
{
    uint64_t r = isqrt(n);
    uint64_t square;
    if (!mul_checked(r, r, &square) || square != n)
        return false;
    if (root != NULL)
        *root = r;
    return true;
}

bool mul_checked(uint64_t a, uint64_t b, uint64_t *out)
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
