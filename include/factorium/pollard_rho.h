/*
 * factorium/pollard_rho.h — Pollard's Rho factorization.
 *
 * Iterates f(x) = (x^2 + c) mod n with Floyd's cycle detection and uses
 * d = gcd(|x - y|, n) to surface a nontrivial factor. Expected running
 * time is heuristic: about O(n^{1/4}) modular multiplications per factor
 * for balanced semiprimes.
 *
 * Determinism: attempt i (0-based) uses x0 = (seed + i) % (n - 1) + 1 and
 * c = i % n, so a fixed seed plus a fixed config reproduce bit-identical
 * runs. On success `seed` reports the effective initial iterate x0 and `c`
 * its polynomial constant, so the winning run can be replayed directly as
 * x = y = seed, f(t) = (t*t + c) mod n.
 *
 * Handling of gcd results: d == 1 means "no factor yet" (keep looping);
 * d == n means the cycle collapsed (x == y mod n) — that attempt is
 * abandoned and the next (x0, c) pair is tried.
 *
 * FCT_NO_FACTOR means the attempt budget ran out — never a primality
 * proof. Only fct_is_prime() establishes primality.
 */
#ifndef FACTORIUM_POLLARD_RHO_H
#define FACTORIUM_POLLARD_RHO_H

#include "factorium/common.h"

#define FCT_RHO_DEFAULT_MAX_ITERATIONS (1ULL << 20)
#define FCT_RHO_DEFAULT_MAX_RETRIES 8
#define FCT_RHO_DEFAULT_SEED 1

typedef struct {
    uint64_t max_iterations_per_attempt; /* 0 => default (2^20)  */
    uint64_t max_retries;                /* 0 => default (8)     */
    uint64_t seed;                       /* initial x0 for attempt 0 */
} fct_rho_config;

fct_rho_config fct_rho_config_default(void);

/* Finds a nontrivial factor of n >= 4.
 *  - even n returns factor 2 immediately;
 *  - n < 4 is FCT_INVALID_INPUT;
 *  - FCT_OK with factor/cofactor on success (seed/c report the successful
 *    attempt);
 *  - FCT_NO_FACTOR / FCT_TERM_BUDGET_EXHAUSTED when every attempt hit its
 *    iteration budget. */
fct_factor_result fct_pollard_rho(uint64_t n, const fct_rho_config *cfg);

#endif /* FACTORIUM_POLLARD_RHO_H */
