/*
 * factorium/modular_arithmetic.h — overflow-safe unsigned 64-bit helpers.
 *
 * mul_mod() uses the unsigned __int128 compiler extension when it is
 * available; otherwise it falls back to a portable double-and-add loop.
 * Both give the same answers, only speed differs.
 */
#ifndef FACTORIUM_MODULAR_ARITHMETIC_H
#define FACTORIUM_MODULAR_ARITHMETIC_H

#include <stdbool.h>
#include <stdint.h>

/* Euclid's algorithm. gcd(0, 0) = 0. */
uint64_t gcd(uint64_t a, uint64_t b);

/* (a * b) mod m, for 0 <= a, b < m and m >= 1. Never overflows. */
uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t m);

/* (a + b) mod m, for 0 <= a, b < m and m >= 1. Never overflows. */
uint64_t add_mod(uint64_t a, uint64_t b, uint64_t m);

/* (base ^ exp) mod m. base is reduced mod m first. If ops is not NULL it
 * receives the number of modular multiplications that were performed. */
uint64_t pow_mod(uint64_t base, uint64_t exp, uint64_t m, uint64_t *ops);

/* floor(sqrt(n)) for any uint64_t, using integer Newton iterations. */
uint64_t isqrt(uint64_t n);

/* True if n is a perfect square; the root is stored in *root. */
bool is_square(uint64_t n, uint64_t *root);

/* *out = a * b, or false if the product does not fit in a uint64_t. */
bool mul_checked(uint64_t a, uint64_t b, uint64_t *out);

#endif /* FACTORIUM_MODULAR_ARITHMETIC_H */
