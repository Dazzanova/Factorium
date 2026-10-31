/*
 * factorium/src/main.c — command line interface.
 *
 * Commands: factor, compare, verify, help. The human-readable output and
 * the --json output are both produced here; the algorithm code never
 * prints anything.
 */
#define _POSIX_C_SOURCE 200809L

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "factorium/common.h"
#include "factorium/factorization.h"
#include "factorium/fermat.h"
#include "factorium/modular_arithmetic.h"
#include "factorium/pollard_p_minus_one.h"
#include "factorium/pollard_rho.h"
#include "factorium/primality.h"

#define U64_FMT PRIu64

/* ------------------------------------------------------------------ */
/* usage                                                              */
/* ------------------------------------------------------------------ */

static void print_usage(FILE *out)
{
    fputs(
        "Factorium - comparative integer factorization laboratory (Phase 1)\n"
        "\n"
        "Usage:\n"
        "  factorium factor <n> [options]      factor <n> (complete by default)\n"
        "  factorium compare <n> [options]     run Fermat, Rho and p-1 on <n>\n"
        "  factorium verify <n> <f1>...<fk>    validate a prime factorization\n"
        "  factorium help                      show this message\n"
        "\n"
        "Options:\n"
        "  --algo <a>    factor only: auto (default), fermat, rho, pminus1\n"
        "  --iters <N>   iteration budget for Fermat/Rho (0 = default 1048576)\n"
        "  --retries <N> Rho attempts (0 = default 8)\n"
        "  --seed <S>    Rho seed (0 = default 1); same seed => same run\n"
        "  --b <N>       p-1 smoothness bound B (0 = default 100000)\n"
        "  --bases <K>   p-1 consecutive bases (0 = default 5)\n"
        "  --json        machine-readable output (fixed schema)\n"
        "  --time        print elapsed wall-clock time (human output only)\n"
        "\n"
        "Exit codes:\n"
        "  0  success: factor found / complete factorization / valid check\n"
        "  1  no factor found within budget (normal outcome, not an error)\n"
        "  2  invalid input or usage error\n"
        "  3  incomplete factorization (unresolved cofactors remain)\n"
        "  4  work/resource limit reached\n"
        "  5  memory or internal failure\n"
        "\n"
        "Note: \"no factor found\" is never a primality proof. Primality is\n"
        "established only by the deterministic Miller-Rabin test (verify).\n",
        out);
}

/* Strict uint64 parsing: digits only, nothing else. */
static bool parse_u64(const char *s, uint64_t *out)
{
    uint64_t value = 0;

    if (s == NULL || *s == '\0')
        return false;
    for (; *s != '\0'; s++) {
        if (*s < '0' || *s > '9')
            return false;
        if (value > (UINT64_MAX - (uint64_t)(*s - '0')) / 10)
            return false; /* would overflow */
        value = value * 10 + (uint64_t)(*s - '0');
    }
    *out = value;
    return true;
}

/* ------------------------------------------------------------------ */
/* JSON output                                                        */
/* ------------------------------------------------------------------ */

/* Prints one JSON object in a single line. `need_comma` mirrors the old
 * buffer: it is true when the next element needs a leading comma. */
typedef struct {
    bool need_comma;
} Json;

static void json_sep(Json *j)
{
    if (j->need_comma)
        putchar(',');
}

static void json_object_begin(Json *j)
{
    json_sep(j);
    putchar('{');
    j->need_comma = false;
}

static void json_object_end(Json *j)
{
    putchar('}');
    j->need_comma = true;
}

static void json_array_begin(Json *j)
{
    json_sep(j);
    putchar('[');
    j->need_comma = false;
}

static void json_array_end(Json *j)
{
    putchar(']');
    j->need_comma = true;
}

static void json_key(Json *j, const char *key)
{
    json_sep(j);
    printf("\"%s\":", key);
    j->need_comma = false;
}

static void json_string(Json *j, const char *key, const char *value)
{
    json_key(j, key);
    printf("\"%s\"", value);
    j->need_comma = true;
}

static void json_u64(Json *j, const char *key, uint64_t value)
{
    json_key(j, key);
    printf("%" U64_FMT, value);
    j->need_comma = true;
}

static void json_bool(Json *j, const char *key, bool value)
{
    json_key(j, key);
    fputs(value ? "true" : "false", stdout);
    j->need_comma = true;
}

static void json_u64_array(Json *j, const char *key,
                           const uint64_t *values, int count)
{
    int i;
    json_key(j, key);
    json_array_begin(j);
    for (i = 0; i < count; i++) {
        json_sep(j);
        printf("%" U64_FMT, values[i]);
        j->need_comma = true;
    }
    json_array_end(j);
}

/* Report a usage/input error on the requested channel; returns 2. */
static int report_error(const char *msg, bool json, const char *cmd)
{
    if (json) {
        Json j = {false};
        json_object_begin(&j);
        json_string(&j, "command", cmd ? cmd : "");
        json_string(&j, "status", "invalid_input");
        json_string(&j, "error", msg);
        json_object_end(&j);
        putchar('\n');
        return 2;
    }
    fprintf(stderr, "error: %s\n\n", msg);
    print_usage(stderr);
    return 2;
}

/* ------------------------------------------------------------------ */
/* options                                                            */
/* ------------------------------------------------------------------ */

typedef enum {
    ALGO_AUTO = 0,
    ALGO_FERMAT,
    ALGO_RHO,
    ALGO_PMINUS1
} Algo;

typedef struct {
    Algo algo;
    uint64_t iters;
    uint64_t retries;
    uint64_t seed;
    uint64_t b;
    uint64_t bases;
    bool json;
    bool show_time;
} Options;

static void options_init(Options *o)
{
    memset(o, 0, sizeof(*o));
    o->algo = ALGO_AUTO;
}

static bool json_flag_present(int argc, char **argv)
{
    int i;
    for (i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--json") == 0)
            return true;
    }
    return false;
}

static const char *const OPTION_HELP_HINT =
    "expected a value for an option (--iters, --retries, --seed, --b, --bases)";

/* Parses the arguments of factor/compare. Returns -1 when parsing is done,
 * or the exit code of an error that was already reported. */
static int parse_factor_args(int argc, char **argv, Options *o, uint64_t *n,
                             bool *have_n, const char *cmd, bool json)
{
    int i;
    for (i = 2; i < argc; i++) {
        const char *arg = argv[i];

        if (strncmp(arg, "--", 2) == 0) {
            if (strcmp(arg, "--json") == 0) {
                o->json = true;
            } else if (strcmp(arg, "--time") == 0) {
                o->show_time = true;
            } else if (strcmp(arg, "--algo") == 0) {
                const char *value;
                if (i + 1 >= argc)
                    return report_error("expected a value for --algo", json, cmd);
                value = argv[++i];
                if (strcmp(value, "auto") == 0)
                    o->algo = ALGO_AUTO;
                else if (strcmp(value, "fermat") == 0)
                    o->algo = ALGO_FERMAT;
                else if (strcmp(value, "rho") == 0)
                    o->algo = ALGO_RHO;
                else if (strcmp(value, "pminus1") == 0)
                    o->algo = ALGO_PMINUS1;
                else
                    return report_error("--algo must be auto, fermat, rho or pminus1",
                                        json, cmd);
            } else if (strcmp(arg, "--iters") == 0 ||
                       strcmp(arg, "--retries") == 0 ||
                       strcmp(arg, "--seed") == 0 ||
                       strcmp(arg, "--b") == 0 ||
                       strcmp(arg, "--bases") == 0) {
                uint64_t value;
                if (i + 1 >= argc)
                    return report_error(OPTION_HELP_HINT, json, cmd);
                if (!parse_u64(argv[++i], &value))
                    return report_error("option value is not a non-negative integer",
                                        json, cmd);
                if (strcmp(arg, "--iters") == 0)
                    o->iters = value;
                else if (strcmp(arg, "--retries") == 0)
                    o->retries = value;
                else if (strcmp(arg, "--seed") == 0)
                    o->seed = value;
                else if (strcmp(arg, "--b") == 0)
                    o->b = value;
                else
                    o->bases = value;
            } else {
                return report_error("unknown option", json, cmd);
            }
        } else {
            if (*have_n)
                return report_error("unexpected extra argument", json, cmd);
            if (!parse_u64(arg, n))
                return report_error("n is not a non-negative integer", json, cmd);
            *have_n = true;
        }
    }
    if (!*have_n)
        return report_error("missing required argument <n>", json, cmd);
    return -1;
}

/* ------------------------------------------------------------------ */
/* factor                                                             */
/* ------------------------------------------------------------------ */

static const char *algo_name(Algo algo)
{
    switch (algo) {
    case ALGO_AUTO: return "auto";
    case ALGO_FERMAT: return "fermat";
    case ALGO_RHO: return "rho";
    case ALGO_PMINUS1: return "pminus1";
    }
    return "unknown";
}

static FactorResult run_single(uint64_t n, const Options *o)
{
    switch (o->algo) {
    case ALGO_FERMAT:
        return fermat(n, o->iters ? o->iters : FERMAT_DEFAULT_ITERATIONS);
    case ALGO_RHO:
        return rho(n, o->iters ? o->iters : RHO_DEFAULT_ITERATIONS,
                   o->retries ? o->retries : RHO_DEFAULT_RETRIES,
                   o->seed ? o->seed : RHO_DEFAULT_SEED);
    case ALGO_PMINUS1:
        return pminus1(n, o->b ? o->b : PM1_DEFAULT_BOUND,
                       PM1_DEFAULT_START_BASE,
                       o->bases ? o->bases : PM1_DEFAULT_BASES,
                       PM1_DEFAULT_MAX_OPS);
    case ALGO_AUTO:
    default: {
        FactorResult bad = {0};
        bad.code = RESULT_INVALID;
        return bad;
    }
    }
}

static void print_single_human(uint64_t n, const Options *o,
                               const FactorResult *r)
{
    printf("command  : factor\n");
    printf("n        : %" U64_FMT "\n", n);
    printf("algorithm: %s\n", algo_name(o->algo));
    printf("status   : %s\n", status_string(r->code));
    printf("reason   : %s\n", reason_string(r->code));
    if (r->code == RESULT_FOUND)
        printf("factor   : %" U64_FMT "\n", r->factor);
    else
        printf("factor   : -\n");
    if (r->code == RESULT_FOUND)
        printf("cofactor : %" U64_FMT "\n", r->cofactor);
    else
        printf("cofactor : -\n");
    printf("iterations: %" U64_FMT "\n", r->iterations);
    printf("gcd_calls: %" U64_FMT "\n", r->gcd_calls);
    printf("modmul_ops: %" U64_FMT "\n", r->modmul_ops);
    if (o->algo == ALGO_RHO && r->code == RESULT_FOUND) {
        printf("seed     : %" U64_FMT "\n", r->seed);
        printf("c        : %" U64_FMT "\n", r->c);
    }
}

static int cmd_factor_single(uint64_t n, const Options *o)
{
    FactorResult r = run_single(n, o);

    if (o->json) {
        Json j = {false};
        json_object_begin(&j);
        json_string(&j, "command", "factor");
        json_u64(&j, "n", n);
        json_string(&j, "algorithm", algo_name(o->algo));
        json_string(&j, "status", status_string(r.code));
        json_string(&j, "reason", reason_string(r.code));
        json_u64(&j, "factor", r.factor);
        json_u64(&j, "cofactor", r.cofactor);
        json_u64(&j, "iterations", r.iterations);
        json_u64(&j, "gcd_calls", r.gcd_calls);
        json_u64(&j, "modmul_ops", r.modmul_ops);
        json_u64(&j, "seed", r.seed);
        json_u64(&j, "c", r.c);
        json_object_end(&j);
        putchar('\n');
    } else {
        print_single_human(n, o, &r);
    }
    return exit_code(r.code);
}

static int cmd_factor_auto(uint64_t n, const Options *o)
{
    FactorizeOptions cfg = default_options();
    Factorization f;

    if (o->iters)
        cfg.max_iterations = o->iters;
    if (o->retries)
        cfg.max_retries = o->retries;
    if (o->seed)
        cfg.seed = o->seed;
    if (o->b)
        cfg.smoothness_bound = o->b;
    if (o->bases)
        cfg.max_bases = o->bases;

    f = factorize(n, &cfg);

    if (o->json) {
        Json j = {false};
        json_object_begin(&j);
        json_string(&j, "command", "factor");
        json_u64(&j, "n", n);
        json_string(&j, "algorithm", "auto");
        json_string(&j, "status", status_string(f.code));
        json_string(&j, "reason", reason_string(f.code));
        json_u64_array(&j, "factors", f.factors, f.factor_count);
        json_u64_array(&j, "unresolved", f.unresolved, f.unresolved_count);
        json_u64(&j, "product", f.product);
        json_bool(&j, "product_verified", f.product_verified);
        json_u64(&j, "iterations", f.iterations);
        json_u64(&j, "gcd_calls", f.gcd_calls);
        json_u64(&j, "modmul_ops", f.modmul_ops);
        json_object_end(&j);
        putchar('\n');
    } else {
        int i;
        printf("command  : factor\n");
        printf("n        : %" U64_FMT "\n", n);
        printf("algorithm: auto\n");
        printf("status   : %s\n", status_string(f.code));
        printf("reason   : %s\n", reason_string(f.code));
        printf("factors  :");
        if (f.factor_count == 0) {
            printf(" (none)\n");
        } else {
            for (i = 0; i < f.factor_count; i++)
                printf(" %" U64_FMT, f.factors[i]);
            printf("\n");
        }
        if (f.unresolved_count > 0) {
            printf("unresolved:");
            for (i = 0; i < f.unresolved_count; i++)
                printf(" %" U64_FMT, f.unresolved[i]);
            printf("\n");
        }
        printf("product  : %" U64_FMT "%s\n", f.product,
               f.product_verified ? " (verified)" : " (not verified)");
        printf("iterations: %" U64_FMT "\n", f.iterations);
        printf("gcd_calls: %" U64_FMT "\n", f.gcd_calls);
        printf("modmul_ops: %" U64_FMT "\n", f.modmul_ops);
    }
    return exit_code(f.code);
}

/* ------------------------------------------------------------------ */
/* compare                                                            */
/* ------------------------------------------------------------------ */

static void print_result_json(Json *j, const char *algo, const FactorResult *r)
{
    json_object_begin(j);
    json_string(j, "algorithm", algo);
    json_string(j, "status", status_string(r->code));
    json_string(j, "reason", reason_string(r->code));
    json_u64(j, "factor", r->factor);
    json_u64(j, "cofactor", r->cofactor);
    json_u64(j, "iterations", r->iterations);
    json_u64(j, "gcd_calls", r->gcd_calls);
    json_u64(j, "modmul_ops", r->modmul_ops);
    json_u64(j, "seed", r->seed);
    json_u64(j, "c", r->c);
    json_object_end(j);
}

static int cmd_compare(uint64_t n, const Options *o)
{
    static const char *const names[] = {"fermat", "rho", "pminus1"};
    FactorResult rs[3];
    ResultCode overall = RESULT_COMPLETE;
    int i;

    rs[0] = fermat(n, o->iters ? o->iters : FERMAT_DEFAULT_ITERATIONS);
    rs[1] = rho(n, o->iters ? o->iters : RHO_DEFAULT_ITERATIONS,
                o->retries ? o->retries : RHO_DEFAULT_RETRIES,
                o->seed ? o->seed : RHO_DEFAULT_SEED);
    rs[2] = pminus1(n, o->b ? o->b : PM1_DEFAULT_BOUND,
                    PM1_DEFAULT_START_BASE,
                    o->bases ? o->bases : PM1_DEFAULT_BASES,
                    PM1_DEFAULT_MAX_OPS);

    /* Any invalid input wins; otherwise ok only if all three succeeded. */
    for (i = 0; i < 3; i++) {
        if (rs[i].code == RESULT_INVALID)
            overall = RESULT_INVALID;
        else if (rs[i].code != RESULT_FOUND && overall != RESULT_INVALID)
            overall = RESULT_NO_FACTOR;
    }

    if (o->json) {
        Json j = {false};
        json_object_begin(&j);
        json_string(&j, "command", "compare");
        json_u64(&j, "n", n);
        json_string(&j, "status", status_string(overall));
        json_key(&j, "results");
        json_array_begin(&j);
        for (i = 0; i < 3; i++)
            print_result_json(&j, names[i], &rs[i]);
        json_array_end(&j);
        json_object_end(&j);
        putchar('\n');
    } else {
        printf("n = %" U64_FMT "\n\n", n);
        printf("%-10s %-12s %-17s %-10s %-10s %11s %10s %11s\n",
               "algorithm", "status", "reason", "factor", "cofactor",
               "iterations", "gcd_calls", "modmul_ops");
        for (i = 0; i < 3; i++) {
            printf("%-10s %-12s %-17s ", names[i],
                   status_string(rs[i].code), reason_string(rs[i].code));
            if (rs[i].code == RESULT_FOUND)
                printf("%-10" U64_FMT " %-10" U64_FMT " ",
                       rs[i].factor, rs[i].cofactor);
            else
                printf("%-10s %-10s ", "-", "-");
            printf("%11" U64_FMT " %10" U64_FMT " %11" U64_FMT "\n",
                   rs[i].iterations, rs[i].gcd_calls, rs[i].modmul_ops);
        }
    }
    return exit_code(overall);
}

/* ------------------------------------------------------------------ */
/* verify                                                             */
/* ------------------------------------------------------------------ */

static int cmd_verify(uint64_t n, const uint64_t *factors, size_t count,
                      bool json)
{
    uint64_t product = 1;
    bool overflow = false;
    bool matches;
    bool all_prime = true;
    bool valid;
    size_t i;

    if (n < 2)
        return report_error("verify requires n >= 2", json, "verify");
    if (count == 0)
        return report_error("verify requires at least one factor", json,
                            "verify");

    for (i = 0; i < count; i++) {
        if (!mul_checked(product, factors[i], &product)) {
            overflow = true;
            break;
        }
    }
    matches = !overflow && product == n;

    for (i = 0; i < count; i++) {
        if (!is_prime(factors[i]))
            all_prime = false;
    }
    valid = matches && all_prime;

    if (json) {
        Json j = {false};
        json_object_begin(&j);
        json_string(&j, "command", "verify");
        json_u64(&j, "n", n);
        json_string(&j, "status", valid ? "ok" : "reconstruction_failed");
        json_u64_array(&j, "factors", factors, (int)count);
        json_u64(&j, "product", overflow ? 0 : product);
        json_bool(&j, "product_matches", matches);
        json_bool(&j, "all_prime", all_prime);
        json_bool(&j, "valid", valid);
        json_object_end(&j);
        putchar('\n');
    } else {
        printf("command  : verify\n");
        printf("n        : %" U64_FMT "\n", n);
        printf("factors  :");
        for (i = 0; i < count; i++)
            printf(" %" U64_FMT, factors[i]);
        printf("\n");
        printf("product  : %" U64_FMT "%s\n", overflow ? 0 : product,
               overflow ? " (overflow)" : "");
        printf("matches  : %s\n", matches ? "yes" : "no");
        printf("all prime: %s\n", all_prime ? "yes" : "no");
        printf("valid    : %s\n", valid ? "yes" : "no");
    }
    return valid ? 0 : 1;
}

/* ------------------------------------------------------------------ */
/* main                                                               */
/* ------------------------------------------------------------------ */

static double elapsed_ms(const struct timespec *a, const struct timespec *b)
{
    return (double)(b->tv_sec - a->tv_sec) * 1000.0 +
           (double)(b->tv_nsec - a->tv_nsec) / 1000000.0;
}

int main(int argc, char **argv)
{
    struct timespec start;
    const char *cmd;

    clock_gettime(CLOCK_MONOTONIC, &start);

    if (argc < 2) {
        fprintf(stderr, "error: missing command\n\n");
        print_usage(stderr);
        return 2;
    }
    cmd = argv[1];

    if (strcmp(cmd, "help") == 0 || strcmp(cmd, "--help") == 0 ||
        strcmp(cmd, "-h") == 0) {
        print_usage(stdout);
        return 0;
    }

    if (strcmp(cmd, "factor") == 0 || strcmp(cmd, "compare") == 0) {
        bool json = json_flag_present(argc, argv);
        Options o;
        uint64_t n = 0;
        bool have_n = false;
        int rc;

        options_init(&o);
        rc = parse_factor_args(argc, argv, &o, &n, &have_n, cmd, json);
        if (rc >= 0)
            return rc;

        if (strcmp(cmd, "factor") == 0)
            rc = (o.algo == ALGO_AUTO) ? cmd_factor_auto(n, &o)
                                       : cmd_factor_single(n, &o);
        else
            rc = cmd_compare(n, &o);

        if (o.show_time && !o.json) {
            struct timespec end;
            clock_gettime(CLOCK_MONOTONIC, &end);
            printf("elapsed_ms: %.3f\n", elapsed_ms(&start, &end));
        }
        return rc;
    }

    if (strcmp(cmd, "verify") == 0) {
        bool json = json_flag_present(argc, argv);
        uint64_t values[1024];
        size_t count = 0;
        bool have_n = false;
        uint64_t n = 0;
        int i;

        for (i = 2; i < argc; i++) {
            uint64_t value;
            if (strcmp(argv[i], "--json") == 0 ||
                strcmp(argv[i], "--time") == 0)
                continue;
            if (!parse_u64(argv[i], &value))
                return report_error("argument is not a non-negative integer",
                                    json, "verify");
            if (!have_n) {
                n = value;
                have_n = true;
                continue;
            }
            if (count >= sizeof(values) / sizeof(values[0]))
                return report_error("too many factors (max 1024)", json, "verify");
            values[count++] = value;
        }
        if (!have_n)
            return report_error("missing required argument <n>", json, "verify");
        return cmd_verify(n, values, count, json);
    }

    fprintf(stderr, "error: unknown command '%s'\n\n", cmd);
    print_usage(stderr);
    return 2;
}
