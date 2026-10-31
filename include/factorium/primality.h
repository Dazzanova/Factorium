/*
 * factorium/primality.h — exact primality test for every uint64_t.
 */
#ifndef FACTORIUM_PRIMALITY_H
#define FACTORIUM_PRIMALITY_H

#include <stdbool.h>
#include <stdint.h>

/* Deterministic Miller-Rabin.
 *
 * The witnesses 2, 3, 5, ..., 37 are enough to decide primality for every
 * n < 318,665,857,834,031,151,167,461, which is larger than 2^64. So this
 * function is exact over the whole uint64_t range: no composite is ever
 * called prime and no prime is ever called composite. n < 2 returns false.
 */
bool is_prime(uint64_t n);

#endif /* FACTORIUM_PRIMALITY_H */
