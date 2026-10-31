/*
 * factorium/common.h — one result code shared by every algorithm.
 *
 * The code says both what happened (status) and why the search stopped
 * (reason). status_string()/reason_string() give the text used by the CLI,
 * the JSON output and the CSV export.
 *
 * Important: "no factor" only means the search ran out of budget. It is
 * never a proof that n is prime -- is_prime() is the only primality test.
 */
#ifndef FACTORIUM_COMMON_H
#define FACTORIUM_COMMON_H

#include <stdbool.h>
#include <stdint.h>

typedef enum {
    RESULT_FOUND = 0,    /* a nontrivial factor was found               */
    RESULT_COMPLETE,     /* full factorization, product verified        */
    RESULT_NO_FACTOR,    /* budget used up, nothing found               */
    RESULT_NO_BASE,      /* p-1: every base missed                      */
    RESULT_INVALID,      /* bad input or bad configuration              */
    RESULT_INCOMPLETE,   /* factorization stopped with cofactors left   */
    RESULT_DEPTH_LIMIT,  /* stopped by the recursion depth cap          */
    RESULT_WORK_LIMIT,   /* hard work cap reached                       */
    RESULT_FACTOR_LIMIT, /* the factor array is full                    */
    RESULT_MEMORY,       /* allocation failure                          */
    RESULT_INTERNAL,     /* defensive: should never happen              */
    RESULT_BAD_PRODUCT   /* factors do not multiply back to n           */
} ResultCode;

const char *status_string(ResultCode code);
const char *reason_string(ResultCode code);
int exit_code(ResultCode code);

/* Result of one factor-search attempt (Fermat, Rho or p-1). */
typedef struct {
    ResultCode code;
    uint64_t factor;     /* the nontrivial factor found                */
    uint64_t cofactor;   /* n / factor                                 */
    uint64_t iterations; /* main loop iterations / bases tried         */
    uint64_t gcd_calls;  /* number of gcd evaluations                  */
    uint64_t modmul_ops; /* modular multiplications performed         */
    uint64_t seed;       /* Rho: initial iterate of the winning attempt */
    uint64_t c;          /* Rho: polynomial constant of that attempt   */
} FactorResult;

#endif /* FACTORIUM_COMMON_H */
