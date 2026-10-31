/*
 * factorium/fermat.h — Fermat's difference-of-squares method.
 *
 * If n = a^2 - b^2 = (a - b)(a + b), then (a - b) is a factor of n.
 * The search starts at a = ceil(sqrt(n)) and increases a until a^2 - n is
 * a perfect square. It is fast when the two factors are close together.
 *
 * Running out of budget is not a proof that n is prime.
 */
#ifndef FACTORIUM_FERMAT_H
#define FACTORIUM_FERMAT_H

#include "factorium/common.h"

#define FERMAT_DEFAULT_ITERATIONS (1ULL << 20)

/* Tries at most max_iterations values of a. 0 means "no attempts". */
FactorResult fermat(uint64_t n, uint64_t max_iterations);

#endif /* FACTORIUM_FERMAT_H */
