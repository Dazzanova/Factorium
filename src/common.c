/*
 * factorium/src/common.c — text and exit code for each result code.
 */
#include "factorium/common.h"

#include <stddef.h>

typedef struct {
    ResultCode code;
    const char *status;
    const char *reason;
    int exit;
} CodeInfo;

static const CodeInfo CODES[] = {
    { RESULT_FOUND,        "ok",                    "found",             0 },
    { RESULT_COMPLETE,     "ok",                    "none",              0 },
    { RESULT_NO_FACTOR,    "no_factor",             "budget_exhausted",  1 },
    { RESULT_NO_BASE,      "no_factor",             "bases_exhausted",   1 },
    { RESULT_INVALID,      "invalid_input",         "invalid_input",     2 },
    { RESULT_INCOMPLETE,   "incomplete",            "budget_exhausted",  3 },
    { RESULT_DEPTH_LIMIT,  "incomplete",            "recursion_limit",   3 },
    { RESULT_WORK_LIMIT,   "work_limit",            "work_limit",        4 },
    { RESULT_FACTOR_LIMIT, "work_limit",            "factor_limit",      4 },
    { RESULT_MEMORY,       "memory",                "none",              5 },
    { RESULT_INTERNAL,     "internal",              "none",              5 },
    { RESULT_BAD_PRODUCT,  "reconstruction_failed", "none",              5 },
};

/* Returns the table row for a code; unknown codes behave as INTERNAL. */
static const CodeInfo *find(ResultCode code)
{
    size_t i;
    for (i = 0; i < sizeof(CODES) / sizeof(CODES[0]); i++) {
        if (CODES[i].code == code)
            return &CODES[i];
    }
    return &CODES[10]; /* RESULT_INTERNAL */
}

const char *status_string(ResultCode code)
{
    return find(code)->status;
}

const char *reason_string(ResultCode code)
{
    return find(code)->reason;
}

int exit_code(ResultCode code)
{
    return find(code)->exit;
}
