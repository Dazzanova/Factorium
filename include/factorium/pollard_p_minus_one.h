/*
 * factorium/pollard_p_minus_one.h — Pollard's p-1, stage 1.
 *
 * If p - 1 is B-smooth (all its prime factors are <= B), then
 * M = product of all prime powers <= B is a multiple of p - 1. For a base
 * a coprime to p, Fermat's little theorem gives a^M = 1 (mod p), so
 * gcd(a^M - 1, n) reveals p.
 *
 * If p - 1 has a prime factor larger than B, stage 1 usually finds
 * nothing. Running out of bases is not a proof that n is prime.
 */
#ifndef FACTORIUM_POLLARD_P_MINUS_ONE_H
#define FACTORIUM_POLLARD_P_MINUS_ONE_H

#include "factorium/common.h"

#define PM1_DEFAULT_BOUND 100000ULL
#define PM1_DEFAULT_START_BASE 2
#define PM1_DEFAULT_BASES 5
#define PM1_DEFAULT_MAX_OPS (1ULL << 24)
#define PM1_MAX_BOUND 1000000000ULL /* one byte of sieve per number */

/* Tries max_bases consecutive bases (start_base, start_base + 1, ...).
 * max_modmul_ops caps the total number of modular multiplications. */
FactorResult pminus1(uint64_t n, uint64_t B, uint64_t start_base,
                     uint64_t max_bases, uint64_t max_modmul_ops);

#endif /* FACTORIUM_POLLARD_P_MINUS_ONE_H */
