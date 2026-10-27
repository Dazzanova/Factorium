/*
 * factorium/modular_arithmetic.h — overflow-safe uint64 modular arithmetic.
 *
 * Every function is total on its documented domain and free of undefined
 * behavior (no signed overflow, no unbounded shifts).
 *
 * Compiler assumption: on platforms where __SIZEOF_INT128__ is defined
 * (GCC/Clang on x86-64, aarch64, riscv64, ...), modular multiplication is
 * implemented with the unsigned __int128 extension. Otherwise a portable
 * double-and-add (Russian peasant) loop is used. Both give identical
 * results; only performance differs.
 */
#ifndef FACTORIUM_MODULAR_ARITHMETRIC_H
#define FACTORIUM_MODULAR_ARITHMETRIC_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* Euclidean GCD. gcd(0, 0) = 0. */
uint64_t fct_gcd(uint64_t a, uint64_t b);

/* (a * b) mod m. Requires 0 <= a, b < m and m >= 1. Overflow-safe. */
uint64_t fct_mod_mul(uint64_t a, uint64_t b, uint64_t m);

/* (a + b) mod m. Requires 0 <= a, b < m and m >= 1. Overflow-safe. */
uint64_t fct_mod_add(uint64_t a, uint64_t b, uint64_t m);

/* (base ^ exp) mod m. Requires m >= 1; base is reduced mod m first.
 * By convention exp == 0 returns 1 % m. */
uint64_t fct_mod_pow(uint64_t base, uint64_t exp, uint64_t m);

/* fct_mod_pow that also reports the number of modular multiplications
 * performed (bit-scan square-and-multiply: for exp > 0 the count is
 * (bit_length(exp) - 1) + (popcount(exp) - 1); for exp == 0 it is 0). */
void fct_mod_pow_counted(uint64_t base, uint64_t exp, uint64_t m,
                         uint64_t *modmul_ops_out);

/* floor(sqrt(n)) for the full uint64_t range. Pure integer Newton; no libm. */
uint64_t fct_isqrt(uint64_t n);

/* Perfect-square detection. On success sets *root_out = sqrt(n). */
bool fct_is_square(uint64_t n, uint64_t *root_out);

/* Sets *out = a * b and returns true; returns false if a * b exceeds
 * UINT64_MAX (in which case *out is unchanged). */
bool fct_mul_checked(uint64_t a, uint64_t b, uint64_t *out);

#endif /* FACTORIUM_MODULAR_ARITHMETRIC_H */
