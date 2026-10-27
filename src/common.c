/*
 * factorium/src/common.c — human-readable names for result codes.
 */
#include "factorium/common.h"

const char *fct_status_str(fct_status s)
{
    switch (s) {
    case FCT_OK: return "ok";
    case FCT_NO_FACTOR: return "no_factor";
    case FCT_INVALID_INPUT: return "invalid_input";
    case FCT_WORK_LIMIT: return "work_limit";
    case FCT_INCOMPLETE: return "incomplete";
    case FCT_RECONSTRUCTION_FAILED: return "reconstruction_failed";
    case FCT_MEMORY: return "memory";
    case FCT_INTERNAL: return "internal";
    }
    return "unknown";
}

const char *fct_termination_str(fct_termination t)
{
    switch (t) {
    case FCT_TERM_NONE: return "none";
    case FCT_TERM_FOUND: return "found";
    case FCT_TERM_BUDGET_EXHAUSTED: return "budget_exhausted";
    case FCT_TERM_BASES_EXHAUSTED: return "bases_exhausted";
    case FCT_TERM_INVALID_INPUT: return "invalid_input";
    case FCT_TERM_WORK_LIMIT: return "work_limit";
    case FCT_TERM_RECURSION_LIMIT: return "recursion_limit";
    case FCT_TERM_FACTOR_LIMIT: return "factor_limit";
    }
    return "unknown";
}
