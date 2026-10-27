/*
 * tests/test_harness.h — minimal assertion harness for Factorium tests.
 *
 * Each test file is a standalone executable: it includes this header, runs
 * its checks, and ends with TEST_MAIN_END(), which reports and returns 1 if
 * any check failed. Keep tests deterministic: fixed seeds only.
 */
#ifndef FACTORIUM_TEST_HARNESS_H
#define FACTORIUM_TEST_HARNESS_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int fh_checks = 0;
static int fh_failures = 0;

#define CHECK(cond)                                                       \
    do {                                                                  \
        fh_checks++;                                                      \
        if (!(cond)) {                                                    \
            fh_failures++;                                                \
            fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        }                                                                 \
    } while (0)

#define CHECK_EQ_U64(a, b)                                                \
    do {                                                                  \
        fh_checks++;                                                      \
        uint64_t _a = (uint64_t)(a), _b = (uint64_t)(b);                  \
        if (_a != _b) {                                                   \
            fh_failures++;                                                \
            fprintf(stderr, "FAIL %s:%d: %s == %s (got %llu, want %llu)\n", \
                    __FILE__, __LINE__, #a, #b,                           \
                    (unsigned long long)_a, (unsigned long long)_b);      \
        }                                                                 \
    } while (0)

#define CHECK_EQ_STR(a, b)                                                \
    do {                                                                  \
        fh_checks++;                                                      \
        const char *_a = (a), *_b = (b);                                  \
        if (_a == NULL || _b == NULL || strcmp(_a, _b) != 0) {            \
            fh_failures++;                                                \
            fprintf(stderr, "FAIL %s:%d: %s == %s (got \"%s\", want \"%s\")\n", \
                    __FILE__, __LINE__, #a, #b,                           \
                    _a ? _a : "(null)", _b ? _b : "(null)");              \
        }                                                                 \
    } while (0)

#define TEST_MAIN_END()                                                   \
    do {                                                                  \
        fprintf(stderr, "%s: %d checks, %d failures\n",                   \
                "test", fh_checks, fh_failures);                          \
        return fh_failures != 0 ? 1 : 0;                                  \
    } while (0)

#endif /* FACTORIUM_TEST_HARNESS_H */
