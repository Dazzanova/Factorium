/*
 * factorium/pollard_p_minus_one.h — Pollard's p-1, stage 1 only.
 *
 * Idea: if p - 1 is B-smooth (all prime factors <= B), then
 * M = product of prime powers p^k <= B is a multiple of p - 1, so for
 * any base a coprime to p, a^M = 1 (mod p) by Fermat. Computing
 * x = a^M mod n and then d = gcd(x - 1, n) surfaces p with good
 * probability whenever p - 1 is B-smooth, even if q - 1 is not.
 *
 * Dependence on smoothness (documented limitation): if the unknown
 * factor p has a prime factor of p - 1 larger than B, stage 1 will
 * typically find nothing for that p. Larger B raises the hit rate at
 * linear-ish cost in sieve and exponent work; stage 2 (not implemented
 * in Phase 1) extends the bound cheaply for the next prime range.
 *
 * gcd outcomes: d in (1, n) -> factor found. d == 1 -> this base saw
 * neither factor; another base is tried. d == n -> a^M = 1 mod n (both
 * factors smooth minus one, or n prime); that base is abandoned.
 * A first check gcd(a, n) catches the case where the base shares a
 * factor with n outright.
 *
 * FCT_NO_FACTOR means every configured base was exhausted — never a
 * primality proof.
 */
#ifndef FACTORIUM_POLLARD_P_MINUS_ONE_H
#define FACTORIUM_POLLARD_P_MINUS_ONE_H

#include "factorium/common.h"

#define FCT_PM1_DEFAULT_SMOOTHNESS_BOUND 100000ULL
#define FCT_PM1_DEFAULT_MAX_BASES 5
#define FCT_PM1_DEFAULT_START_BASE 2
#define FCT_PM1_DEFAULT_MAX_MODMUL_OPS (1ULL << 20)
/* Stage 1 needs one byte per candidate up to B. */
#define FCT_PM1_MAX_SMOOTHNESS_BOUND 1000000000ULL

typedef struct {
    uint64_t smoothness_bound; /* B, primes <= B enter M; 0 => default  */
    uint64_t start_base;       /* first base a; 0 => default 2          */
    uint64_t max_bases;        /* how many consecutive bases to try;
                                * 0 => default 5                        */
    uint64_t max_modmul_ops;   /* hard work cap across all bases;
                                * 0 => default 2^20                     */
} fct_pminus1_config;

fct_pminus1_config fct_pminus1_config_default(void);

/* Finds a nontrivial factor of n >= 4 via stage 1.
 *  - even n returns factor 2 immediately;
 *  - n < 4 or B > FCT_PM1_MAX_SMOOTHNESS_BOUND => FCT_INVALID_INPUT;
 *  - sieve allocation failure => FCT_MEMORY;
 *  - work cap hit => FCT_WORK_LIMIT with the progress made so far;
 *  - all bases exhausted => FCT_NO_FACTOR / FCT_TERM_BASES_EXHAUSTED.
 * modmul_ops counts modular multiplications; iterations counts bases. */
fct_factor_result fct_pollard_p_minus_one(uint64_t n,
                                          const fct_pminus1_config *cfg);

#endif /* FACTORIUM_POLLARD_P_MINUS_ONE_H */
