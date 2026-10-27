/*
 * factorium/common.h — shared result types for the Factorium core library.
 *
 * Design rules:
 *  - Math modules return values (or heap-allocated, caller-freed arrays);
 *    they never print, never measure time, and never interpret a failed
 *    factorization attempt as proof of primality.
 *  - "No factor found within budget" (FCT_NO_FACTOR) is a normal outcome of
 *    a probabilistic/limited search, NOT evidence that n is prime. Only
 *    fct_is_prime() (deterministic Miller-Rabin) establishes primality.
 */
#ifndef FACTORIUM_COMMON_H
#define FACTORIUM_COMMON_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* Outcome of a factorization operation. */
typedef enum {
    FCT_OK = 0,                  /* factor found, or complete factorization */
    FCT_NO_FACTOR,               /* search exhausted within budget; no factor */
    FCT_INVALID_INPUT,           /* n out of domain (n < 4 for factor finding,
                                    n < 2 for factorization) or bad config */
    FCT_WORK_LIMIT,              /* hard resource cap reached mid-search */
    FCT_INCOMPLETE,              /* factorization finished but >= 1 cofactor
                                    could not be resolved within limits */
    FCT_RECONSTRUCTION_FAILED,   /* product of factors does not equal n */
    FCT_MEMORY,                  /* allocation failure */
    FCT_INTERNAL                 /* unexpected internal error (defensive) */
} fct_status;

const char *fct_status_str(fct_status s);

/* Why a search stopped. */
typedef enum {
    FCT_TERM_NONE = 0,
    FCT_TERM_FOUND,              /* nontrivial factor found */
    FCT_TERM_BUDGET_EXHAUSTED,   /* iteration/retry budget used up */
    FCT_TERM_BASES_EXHAUSTED,    /* p-1: every base tried gave gcd 1 or n */
    FCT_TERM_INVALID_INPUT,
    FCT_TERM_WORK_LIMIT,
    FCT_TERM_RECURSION_LIMIT,    /* factorization: depth cap */
    FCT_TERM_FACTOR_LIMIT        /* factorization: factor-count cap */
} fct_termination;

const char *fct_termination_str(fct_termination t);

/* Result of a single factor-discovery attempt (any algorithm). */
typedef struct {
    fct_status status;
    uint64_t factor;      /* nontrivial factor found (status == FCT_OK)   */
    uint64_t cofactor;    /* n / factor                                   */
    uint64_t iterations;  /* main-loop iterations executed                */
    uint64_t gcd_calls;   /* number of gcd evaluations                    */
    uint64_t modmul_ops;  /* modular multiplications performed            */
    uint64_t seed;        /* seed used (Pollard Rho; 0 for others)        */
    uint64_t c;           /* polynomial constant used (Pollard Rho)       */
    fct_termination reason;
} fct_factor_result;

#endif /* FACTORIUM_COMMON_H */
