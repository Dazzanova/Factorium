/*
 * factorium/primality.h — deterministic primality for the full uint64_t range.
 */
#ifndef FACTORIUM_PRIMALITY_H
#define FACTORIUM_PRIMALITY_H

#include <stdbool.h>
#include <stdint.h>

/* Deterministic Miller-Rabin.
 *
 * Witnesses: the first twelve primes {2,3,5,7,11,13,17,19,23,29,31,37}.
 * This set is proven to be deterministic for every n <
 * 318,665,857,834,031,151,167,461, which is far above UINT64_MAX
 * (18,446,744,073,709,551,615), so fct_is_prime is exact for all uint64_t:
 * it never reports a composite as prime and never reports a prime as
 * composite.
 *
 * n < 2 returns false. */
bool fct_is_prime(uint64_t n);

#endif /* FACTORIUM_PRIMALITY_H */
