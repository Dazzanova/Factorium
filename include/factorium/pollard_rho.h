/*
 * factorium/pollard_rho.h — Pollard's Rho with Floyd's cycle detection.
 *
 * f(x) = (x^2 + c) mod n. We walk x once and y twice, and take
 * d = gcd(|x - y|, n). d = 1 means "keep going", d = n means the walk
 * collapsed (try another c), and 1 < d < n is a factor.
 *
 * Attempt i uses x0 = (seed + i) mod (n - 1) + 1 and c = i mod n, so the
 * same seed always gives the same run. Running out of budget is not a
 * proof that n is prime.
 */
#ifndef FACTORIUM_POLLARD_RHO_H
#define FACTORIUM_POLLARD_RHO_H

#include "factorium/common.h"

#define RHO_DEFAULT_ITERATIONS (1ULL << 20)
#define RHO_DEFAULT_RETRIES 8
#define RHO_DEFAULT_SEED 1

/* Tries max_retries different (x0, c) pairs, each for at most
 * max_iterations steps. On success, seed and c describe the winning run. */
FactorResult rho(uint64_t n, uint64_t max_iterations,
                 uint64_t max_retries, uint64_t seed);

#endif /* FACTORIUM_POLLARD_RHO_H */
