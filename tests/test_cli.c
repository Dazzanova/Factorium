/*
 * tests/test_cli.c — CLI behavior: exit codes, JSON content, human output.
 *
 * Run as: test_cli <path-to-factorium-binary>
 * Each case captures stdout to a temp file in the working directory
 * (the build dir when driven by CTest).
 */
#define _POSIX_C_SOURCE 200809L

#include "test_harness.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/wait.h>

static const char *g_bin = NULL;
static const char *g_capture = "test_cli_capture.tmp";

/* Runs `binary <args...>` with stdout captured; returns the exit code.
 * On a launch failure returns -1. */
static int run(const char *args, char *out, size_t out_size)
{
    char cmd[4096];
    int n = snprintf(cmd, sizeof(cmd), "%s %s > %s 2>/dev/null", g_bin, args,
                     g_capture);
    if (n < 0 || (size_t)n >= sizeof(cmd)) {
        out[0] = '\0';
        return -1;
    }
    int st = system(cmd);
    if (st == -1) {
        out[0] = '\0';
        return -1;
    }
    int code = WIFEXITED(st) ? WEXITSTATUS(st) : -1;

    out[0] = '\0';
    FILE *f = fopen(g_capture, "r");
    if (f != NULL) {
        size_t got = fread(out, 1, out_size - 1, f);
        out[got] = '\0';
        fclose(f);
        remove(g_capture);
    }
    return code;
}

static void test_help_and_usage(void)
{
    char out[16384];
    int rc = run("help", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 0);
    CHECK(strstr(out, "Usage:") != NULL);
    CHECK(strstr(out, "Exit codes:") != NULL);

    rc = run("", out, sizeof(out)); /* no command at all */
    CHECK_EQ_U64((uint64_t)rc, 2);

    rc = run("nonsense", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 2);
    CHECK(strstr(out, "unknown command") == NULL); /* stderr only */
}

static void test_factor_complete(void)
{
    char out[16384];
    int rc = run("factor 360", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 0);
    CHECK(strstr(out, "status   : ok") != NULL);
    CHECK(strstr(out, "factors  : 2 2 2 3 3 5") != NULL);
    CHECK(strstr(out, "(verified)") != NULL);

    rc = run("factor 18446744073709551615", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 0);
    CHECK(strstr(out, "6700417") != NULL);
}

static void test_factor_single_algorithms(void)
{
    char out[16384];
    int rc = run("factor 91 --algo fermat --json", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 0);
    CHECK(strstr(out, "\"algorithm\":\"fermat\"") != NULL);
    CHECK(strstr(out, "\"status\":\"ok\"") != NULL);
    CHECK(strstr(out, "\"factor\":7") != NULL);
    CHECK(strstr(out, "\"cofactor\":13") != NULL);

    rc = run("factor 91 --algo rho --json", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 0);
    CHECK(strstr(out, "\"algorithm\":\"rho\"") != NULL);
    CHECK(strstr(out, "\"status\":\"ok\"") != NULL);

    rc = run("factor 15 --algo pminus1 --json", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 0);
    CHECK(strstr(out, "\"algorithm\":\"pminus1\"") != NULL);
    CHECK(strstr(out, "\"status\":\"ok\"") != NULL);

    /* rho seed reporting and reproducibility: identical args, identical
     * output. */
    char out2[16384];
    int rc_a = run("factor 1000000016000000063 --algo rho --seed 42 --json",
                   out, sizeof(out));
    int rc_b = run("factor 1000000016000000063 --algo rho --seed 42 --json",
                   out2, sizeof(out2));
    CHECK_EQ_U64((uint64_t)rc_a, 0);
    CHECK_EQ_U64((uint64_t)rc_b, 0);
    CHECK(strcmp(out, out2) == 0);
}

static void test_exit_codes(void)
{
    char out[16384];

    /* Budget exhausted: NO_FACTOR -> 1, iterations exactly as asked. */
    int rc = run("factor 9999999967 --algo fermat --iters 1000 --json",
                 out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 1);
    CHECK(strstr(out, "\"status\":\"no_factor\"") != NULL);
    CHECK(strstr(out, "\"iterations\":1000") != NULL);

    /* Invalid inputs -> 2. */
    rc = run("factor 1 --json", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 2);
    CHECK(strstr(out, "\"status\":\"invalid_input\"") != NULL);
    rc = run("factor 0 --json", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 2);
    rc = run("factor abc --json", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 2);
    CHECK(strstr(out, "\"status\":\"invalid_input\"") != NULL);
    rc = run("factor 15 --bogus --json", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 2);
    rc = run("factor 15 --algo bogus --json", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 2);
    rc = run("factor 15 --json --iters", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 2);
    CHECK(strstr(out, "\"status\":\"invalid_input\"") != NULL);
    rc = run("factor --json", out, sizeof(out)); /* missing n */
    CHECK_EQ_U64((uint64_t)rc, 2);

    /* Incomplete factorization -> 3 (same recipe as the library test). */
    rc = run("factor 17179869194064771037 --iters 1 --retries 1 --b 2 "
             "--bases 1 --json",
             out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 3);
    CHECK(strstr(out, "\"status\":\"incomplete\"") != NULL);
    CHECK(strstr(out, "\"unresolved\":[17179869194064771037]") != NULL);
}

static void test_compare(void)
{
    char out[16384];
    int rc = run("compare 15 --json", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 0);
    CHECK(strstr(out, "\"command\":\"compare\"") != NULL);
    CHECK(strstr(out, "\"algorithm\":\"fermat\"") != NULL);
    CHECK(strstr(out, "\"algorithm\":\"rho\"") != NULL);
    CHECK(strstr(out, "\"algorithm\":\"pminus1\"") != NULL);

    /* Human table lists all three algorithms. */
    rc = run("compare 15", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 0);
    CHECK(strstr(out, "algorithm  status") != NULL);
    CHECK(strstr(out, "fermat") != NULL);
    CHECK(strstr(out, "pminus1") != NULL);

    /* compare on an invalid input: at least one INVALID_INPUT -> 2. */
    rc = run("compare 3 --json", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 2);
    CHECK(strstr(out, "\"status\":\"invalid_input\"") != NULL);
}

static void test_verify(void)
{
    char out[16384];
    int rc = run("verify 360 2 2 2 3 3 5 --json", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 0);
    CHECK(strstr(out, "\"valid\":true") != NULL);
    CHECK(strstr(out, "\"product_matches\":true") != NULL);
    CHECK(strstr(out, "\"all_prime\":true") != NULL);

    /* Composite factor is not a valid prime factorization. */
    rc = run("verify 360 2 2 2 3 3 6 --json", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 1);
    CHECK(strstr(out, "\"valid\":false") != NULL);
    CHECK(strstr(out, "\"all_prime\":false") != NULL);

    /* Wrong product. */
    rc = run("verify 360 2 2 2 3 3 4 --json", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 1);
    CHECK(strstr(out, "\"product_matches\":false") != NULL);

    /* 2^64 - 1 with its full prime factorization. */
    rc = run("verify 18446744073709551615 3 5 17 257 641 65537 6700417 --json",
             out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 0);
    CHECK(strstr(out, "\"valid\":true") != NULL);

    /* Human form. */
    rc = run("verify 10 2 5", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 0);
    CHECK(strstr(out, "valid    : yes") != NULL);

    /* Invalid inputs -> 2. */
    rc = run("verify 1 1 --json", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 2);
    rc = run("verify 10 --json", out, sizeof(out)); /* no factors */
    CHECK_EQ_U64((uint64_t)rc, 2);
}

static void test_time_flag(void)
{
    char out[16384];
    int rc = run("factor 15 --time", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 0);
    CHECK(strstr(out, "elapsed_ms:") != NULL);
    /* --time must not appear in JSON output. */
    rc = run("factor 15 --time --json", out, sizeof(out));
    CHECK_EQ_U64((uint64_t)rc, 0);
    CHECK(strstr(out, "elapsed_ms") == NULL);
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: %s <factorium-binary>\n", argv[0]);
        return 2;
    }
    g_bin = argv[1];
    test_help_and_usage();
    test_factor_complete();
    test_factor_single_algorithms();
    test_exit_codes();
    test_compare();
    test_verify();
    test_time_flag();
    TEST_MAIN_END();
}
