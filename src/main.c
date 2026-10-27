/*
 * factorium/src/main.c — command-line interface.
 *
 * The CLI owns parsing, presentation (human or JSON), timing, and exit
 * codes; it contains no factorization logic of its own.
 *
 * Exit codes:
 *   0  success: factor found / complete factorization / valid verification
 *   1  no factor found within the configured budget (a normal outcome of
 *      a limited search, never a primality claim)
 *   2  invalid input or usage error
 *   3  incomplete factorization (unresolved cofactors remain)
 *   4  work/resource limit reached
 *   5  memory or internal failure
 */
#define _POSIX_C_SOURCE 200809L

#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
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
/* status -> exit code                                                */
/* ------------------------------------------------------------------ */

static int exit_code_for_status(fct_status s)
{
    switch (s) {
    case FCT_OK: return 0;
    case FCT_NO_FACTOR: return 1;
    case FCT_INVALID_INPUT: return 2;
    case FCT_INCOMPLETE: return 3;
    case FCT_WORK_LIMIT: return 4;
    case FCT_MEMORY:
    case FCT_INTERNAL:
    case FCT_RECONSTRUCTION_FAILED: return 5;
    }
    return 5;
}

/* ------------------------------------------------------------------ */
/* strict uint64 parsing (digits only, no signs, no locale surprises) */
/* ------------------------------------------------------------------ */

static bool parse_u64(const char *s, uint64_t *out)
{
    if (s == NULL || *s == '\0')
        return false;
    uint64_t v = 0;
    for (const char *p = s; *p != '\0'; p++) {
        if (*p < '0' || *p > '9')
            return false;
        uint64_t digit = (uint64_t)(*p - '0');
        if (v > (UINT64_MAX - digit) / 10)
            return false; /* would overflow */
        v = v * 10 + digit;
    }
    *out = v;
    return true;
}

/* ------------------------------------------------------------------ */
/* minimal growable JSON writer with string escaping                  */
/* ------------------------------------------------------------------ */

typedef struct {
    char *data;
    size_t len;
    size_t cap;
    bool failed;     /* allocation failure: output is silently dropped  */
    bool need_comma; /* separator state for object/array elements       */
} jbuf;

static void jb_init(jbuf *b)
{
    b->data = NULL;
    b->len = 0;
    b->cap = 0;
    b->failed = false;
    b->need_comma = false;
}

static void jb_free(jbuf *b)
{
    free(b->data);
    b->data = NULL;
    b->len = b->cap = 0;
}

static bool jb_reserve(jbuf *b, size_t extra)
{
    if (b->failed)
        return false;
    if (b->len + extra + 1 <= b->cap)
        return true;
    size_t ncap = b->cap ? b->cap : 256;
    while (ncap < b->len + extra + 1) {
        if (ncap > SIZE_MAX / 2) {
            b->failed = true;
            return false;
        }
        ncap *= 2;
    }
    char *nd = (char *)realloc(b->data, ncap);
    if (nd == NULL) {
        b->failed = true;
        return false;
    }
    b->data = nd;
    b->cap = ncap;
    return true;
}

static void jb_putc(jbuf *b, char c)
{
    if (!jb_reserve(b, 1))
        return;
    b->data[b->len++] = c;
    b->data[b->len] = '\0';
}

static void jb_puts(jbuf *b, const char *s)
{
    size_t n = strlen(s);
    if (!jb_reserve(b, n))
        return;
    memcpy(b->data + b->len, s, n);
    b->len += n;
    b->data[b->len] = '\0';
}

static void jb_put_escaped(jbuf *b, const char *s)
{
    jb_putc(b, '"');
    for (const char *p = s; *p != '\0'; p++) {
        unsigned char ch = (unsigned char)*p;
        switch (ch) {
        case '"': jb_puts(b, "\\\""); break;
        case '\\': jb_puts(b, "\\\\"); break;
        case '\b': jb_puts(b, "\\b"); break;
        case '\f': jb_puts(b, "\\f"); break;
        case '\n': jb_puts(b, "\\n"); break;
        case '\r': jb_puts(b, "\\r"); break;
        case '\t': jb_puts(b, "\\t"); break;
        default:
            if (ch < 0x20) {
                char tmp[7];
                snprintf(tmp, sizeof(tmp), "\\u%04x", (unsigned)ch);
                jb_puts(b, tmp);
            } else {
                jb_putc(b, (char)ch);
            }
        }
    }
    jb_putc(b, '"');
}

static void jb_sep(jbuf *b)
{
    if (b->need_comma)
        jb_putc(b, ',');
}

static void jb_obj_begin(jbuf *b)
{
    jb_sep(b);
    jb_putc(b, '{');
    b->need_comma = false;
}

static void jb_obj_end(jbuf *b)
{
    jb_putc(b, '}');
    b->need_comma = true;
}

static void jb_arr_begin(jbuf *b)
{
    jb_sep(b);
    jb_putc(b, '[');
    b->need_comma = false;
}

static void jb_arr_end(jbuf *b)
{
    jb_putc(b, ']');
    b->need_comma = true;
}

static void jb_key(jbuf *b, const char *key)
{
    jb_sep(b);
    jb_put_escaped(b, key);
    jb_putc(b, ':');
    b->need_comma = false; /* value follows without a separator */
}

static void jb_kv_str(jbuf *b, const char *key, const char *val)
{
    jb_key(b, key);
    jb_put_escaped(b, val);
    b->need_comma = true;
}

static void jb_kv_u64(jbuf *b, const char *key, uint64_t val)
{
    char tmp[32];
    snprintf(tmp, sizeof(tmp), "%" U64_FMT, val);
    jb_key(b, key);
    jb_puts(b, tmp);
    b->need_comma = true;
}

static void jb_kv_bool(jbuf *b, const char *key, bool val)
{
    jb_key(b, key);
    jb_puts(b, val ? "true" : "false");
    b->need_comma = true;
}

static void jb_arr_u64(jbuf *b, const uint64_t *vals, size_t count)
{
    jb_arr_begin(b);
    for (size_t i = 0; i < count; i++) {
        char tmp[32];
        snprintf(tmp, sizeof(tmp), "%" U64_FMT, vals[i]);
        jb_sep(b);
        jb_puts(b, tmp);
        b->need_comma = true;
    }
    jb_arr_end(b);
}

static void jb_arr_result(jbuf *b, const char *const *algos,
                          const fct_factor_result *rs, size_t count)
{
    jb_arr_begin(b);
    for (size_t i = 0; i < count; i++) {
        jb_obj_begin(b);
        jb_kv_str(b, "algorithm", algos[i]);
        jb_kv_str(b, "status", fct_status_str(rs[i].status));
        jb_kv_str(b, "reason", fct_termination_str(rs[i].reason));
        jb_kv_u64(b, "factor", rs[i].factor);
        jb_kv_u64(b, "cofactor", rs[i].cofactor);
        jb_kv_u64(b, "iterations", rs[i].iterations);
        jb_kv_u64(b, "gcd_calls", rs[i].gcd_calls);
        jb_kv_u64(b, "modmul_ops", rs[i].modmul_ops);
        jb_kv_u64(b, "seed", rs[i].seed);
        jb_kv_u64(b, "c", rs[i].c);
        jb_obj_end(b);
    }
    jb_arr_end(b);
}

/* Emit a finished buffer plus newline, or fail loudly. */
static int jb_flush(const jbuf *b)
{
    if (b->failed) {
        fprintf(stderr, "error: out of memory while building JSON output\n");
        return 5;
    }
    fputs(b->data ? b->data : "{}", stdout);
    fputc('\n', stdout);
    return 0;
}

/* ------------------------------------------------------------------ */
/* usage / diagnostics                                                */
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

/* Report a usage/input error in the requested channel; returns exit code. */
static int report_error(const char *msg, bool json, const char *cmd)
{
    if (json) {
        jbuf b;
        jb_init(&b);
        jb_obj_begin(&b);
        jb_kv_str(&b, "command", cmd ? cmd : "");
        jb_kv_str(&b, "status", "invalid_input");
        jb_kv_str(&b, "error", msg);
        jb_obj_end(&b);
        int rc = jb_flush(&b);
        jb_free(&b);
        return rc == 0 ? 2 : rc;
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
} algo_t;

typedef struct {
    algo_t algo;
    uint64_t iters;
    uint64_t retries;
    uint64_t seed;
    uint64_t b;
    uint64_t bases;
    bool json;
    bool show_time;
} options;

static void options_init(options *o)
{
    o->algo = ALGO_AUTO;
    o->iters = 0;
    o->retries = 0;
    o->seed = 0;
    o->b = 0;
    o->bases = 0;
    o->json = false;
    o->show_time = false;
}

static bool json_flag_present(int argc, char **argv)
{
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--json") == 0)
            return true;
    }
    return false;
}

static const char *const OPTION_HELP_HINT =
    "expected a value for an option (--iters, --retries, --seed, --b, --bases)";

/* Parses flags/options for factor/compare; the first positional is *n.
 * Returns 0 on success, or an error code already reported. */
static int parse_factor_args(int argc, char **argv, options *o, uint64_t *n,
                             bool *have_n, const char *cmd, bool json)
{
    for (int i = 2; i < argc; i++) {
        const char *a = argv[i];
        if (strncmp(a, "--", 2) == 0) {
            if (strcmp(a, "--json") == 0) {
                o->json = true;
            } else if (strcmp(a, "--time") == 0) {
                o->show_time = true;
            } else if (strcmp(a, "--algo") == 0) {
                if (i + 1 >= argc)
                    return report_error("expected a value for --algo", json, cmd);
                const char *val = argv[++i];
                if (strcmp(val, "auto") == 0) {
                    o->algo = ALGO_AUTO;
                } else if (strcmp(val, "fermat") == 0) {
                    o->algo = ALGO_FERMAT;
                } else if (strcmp(val, "rho") == 0) {
                    o->algo = ALGO_RHO;
                } else if (strcmp(val, "pminus1") == 0) {
                    o->algo = ALGO_PMINUS1;
                } else {
                    return report_error(
                        "--algo must be auto, fermat, rho or pminus1",
                        json, cmd);
                }
            } else if (strcmp(a, "--iters") == 0 || strcmp(a, "--retries") == 0 ||
                       strcmp(a, "--seed") == 0 || strcmp(a, "--b") == 0 ||
                       strcmp(a, "--bases") == 0) {
                if (i + 1 >= argc)
                    return report_error(OPTION_HELP_HINT, json, cmd);
                const char *val = argv[++i];
                uint64_t nv;
                if (!parse_u64(val, &nv))
                    return report_error("option value is not a non-negative integer",
                                        json, cmd);
                if (strcmp(a, "--iters") == 0) {
                    o->iters = nv;
                } else if (strcmp(a, "--retries") == 0) {
                    o->retries = nv;
                } else if (strcmp(a, "--seed") == 0) {
                    o->seed = nv;
                } else if (strcmp(a, "--b") == 0) {
                    o->b = nv;
                } else {
                    o->bases = nv;
                }
            } else {
                return report_error("unknown option", json, cmd);
            }
        } else {
            if (*have_n)
                return report_error("unexpected extra argument", json, cmd);
            if (!parse_u64(a, n))
                return report_error("n is not a non-negative integer", json, cmd);
            *have_n = true;
        }
    }
    if (!*have_n)
        return report_error("missing required argument <n>", json, cmd);
    return -1; /* caller continues */
}

/* ------------------------------------------------------------------ */
/* factor                                                             */
/* ------------------------------------------------------------------ */

static const char *algo_name(algo_t a)
{
    switch (a) {
    case ALGO_AUTO: return "auto";
    case ALGO_FERMAT: return "fermat";
    case ALGO_RHO: return "rho";
    case ALGO_PMINUS1: return "pminus1";
    }
    return "unknown";
}

static fct_factor_result run_single(uint64_t n, const options *o)
{
    switch (o->algo) {
    case ALGO_FERMAT: {
        fct_fermat_config c;
        c.max_iterations = o->iters ? o->iters
                                    : FCT_FERMAT_DEFAULT_MAX_ITERATIONS;
        return fct_fermat(n, &c);
    }
    case ALGO_RHO: {
        fct_rho_config c;
        c.max_iterations_per_attempt =
            o->iters ? o->iters : FCT_RHO_DEFAULT_MAX_ITERATIONS;
        c.max_retries = o->retries ? o->retries : FCT_RHO_DEFAULT_MAX_RETRIES;
        c.seed = o->seed ? o->seed : FCT_RHO_DEFAULT_SEED;
        return fct_pollard_rho(n, &c);
    }
    case ALGO_PMINUS1: {
        fct_pminus1_config c;
        c.smoothness_bound = o->b ? o->b : FCT_PM1_DEFAULT_SMOOTHNESS_BOUND;
        c.start_base = 0; /* 0 => default 2 */
        c.max_bases = o->bases ? o->bases : FCT_PM1_DEFAULT_MAX_BASES;
        c.max_modmul_ops = 0; /* 0 => default cap */
        return fct_pollard_p_minus_one(n, &c);
    }
    case ALGO_AUTO:
    default:
        break;
    }
    fct_factor_result invalid;
    memset(&invalid, 0, sizeof(invalid));
    invalid.status = FCT_INVALID_INPUT;
    invalid.reason = FCT_TERM_INVALID_INPUT;
    return invalid;
}

static void print_single_human(uint64_t n, const options *o,
                               const fct_factor_result *r)
{
    printf("command  : factor\n");
    printf("n        : %" U64_FMT "\n", n);
    printf("algorithm: %s\n", algo_name(o->algo));
    printf("status   : %s\n", fct_status_str(r->status));
    printf("reason   : %s\n", fct_termination_str(r->reason));
    printf("factor   : ");
    if (r->status == FCT_OK)
        printf("%" U64_FMT "\n", r->factor);
    else
        printf("-\n");
    printf("cofactor : ");
    if (r->status == FCT_OK)
        printf("%" U64_FMT "\n", r->cofactor);
    else
        printf("-\n");
    printf("iterations: %" U64_FMT "\n", r->iterations);
    printf("gcd_calls: %" U64_FMT "\n", r->gcd_calls);
    printf("modmul_ops: %" U64_FMT "\n", r->modmul_ops);
    if (o->algo == ALGO_RHO && r->status == FCT_OK) {
        printf("seed     : %" U64_FMT "\n", r->seed);
        printf("c        : %" U64_FMT "\n", r->c);
    }
}

static int cmd_factor_single(uint64_t n, const options *o)
{
    fct_factor_result r = run_single(n, o);
    if (o->json) {
        jbuf b;
        jb_init(&b);
        jb_obj_begin(&b);
        jb_kv_str(&b, "command", "factor");
        jb_kv_u64(&b, "n", n);
        jb_kv_str(&b, "algorithm", algo_name(o->algo));
        jb_kv_str(&b, "status", fct_status_str(r.status));
        jb_kv_str(&b, "reason", fct_termination_str(r.reason));
        jb_kv_u64(&b, "factor", r.factor);
        jb_kv_u64(&b, "cofactor", r.cofactor);
        jb_kv_u64(&b, "iterations", r.iterations);
        jb_kv_u64(&b, "gcd_calls", r.gcd_calls);
        jb_kv_u64(&b, "modmul_ops", r.modmul_ops);
        jb_kv_u64(&b, "seed", r.seed);
        jb_kv_u64(&b, "c", r.c);
        jb_obj_end(&b);
        int rc = jb_flush(&b);
        jb_free(&b);
        if (rc != 0)
            return rc;
    } else {
        print_single_human(n, o, &r);
    }
    return exit_code_for_status(r.status);
}

static int cmd_factor_auto(uint64_t n, const options *o)
{
    fct_factorization_config cfg = fct_factorization_config_default();
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

    fct_factorization f = fct_factorize(n, &cfg);

    if (o->json) {
        jbuf b;
        jb_init(&b);
        jb_obj_begin(&b);
        jb_kv_str(&b, "command", "factor");
        jb_kv_u64(&b, "n", n);
        jb_kv_str(&b, "algorithm", "auto");
        jb_kv_str(&b, "status", fct_status_str(f.status));
        jb_kv_str(&b, "reason", fct_termination_str(f.reason));
        jb_key(&b, "factors");
        jb_arr_u64(&b, f.factors, f.factor_count);
        jb_key(&b, "unresolved");
        jb_arr_u64(&b, f.unresolved_cofactors, f.unresolved_count);
        jb_kv_u64(&b, "product", f.product);
        jb_kv_bool(&b, "product_verified", f.product_verified);
        jb_kv_u64(&b, "iterations", f.iterations);
        jb_kv_u64(&b, "gcd_calls", f.gcd_calls);
        jb_kv_u64(&b, "modmul_ops", f.modmul_ops);
        jb_obj_end(&b);
        int rc = jb_flush(&b);
        jb_free(&b);
        if (rc != 0) {
            fct_factorization_free(&f);
            return rc;
        }
    } else {
        printf("command  : factor\n");
        printf("n        : %" U64_FMT "\n", n);
        printf("algorithm: auto\n");
        printf("status   : %s\n", fct_status_str(f.status));
        printf("reason   : %s\n", fct_termination_str(f.reason));
        printf("factors  :");
        if (f.factor_count == 0) {
            printf(" (none)\n");
        } else {
            for (size_t i = 0; i < f.factor_count; i++)
                printf(" %" U64_FMT, f.factors[i]);
            printf("\n");
        }
        if (f.unresolved_count > 0) {
            printf("unresolved:");
            for (size_t i = 0; i < f.unresolved_count; i++)
                printf(" %" U64_FMT, f.unresolved_cofactors[i]);
            printf("\n");
        }
        printf("product  : %" U64_FMT "%s\n", f.product,
               f.product_verified ? " (verified)" : " (not verified)");
        printf("iterations: %" U64_FMT "\n", f.iterations);
        printf("gcd_calls: %" U64_FMT "\n", f.gcd_calls);
        printf("modmul_ops: %" U64_FMT "\n", f.modmul_ops);
    }

    fct_status st = f.status;
    fct_factorization_free(&f);
    return exit_code_for_status(st);
}

/* ------------------------------------------------------------------ */
/* compare                                                            */
/* ------------------------------------------------------------------ */

static int cmd_compare(uint64_t n, const options *o)
{
    static const char *const names[] = {"fermat", "rho", "pminus1"};
    fct_factor_result rs[3];
    rs[0] = fct_fermat(n, &(fct_fermat_config){
                               .max_iterations = o->iters
                                                     ? o->iters
                                                     : FCT_FERMAT_DEFAULT_MAX_ITERATIONS});
    rs[1] = fct_pollard_rho(n,
                            &(fct_rho_config){
                                .max_iterations_per_attempt =
                                    o->iters ? o->iters
                                             : FCT_RHO_DEFAULT_MAX_ITERATIONS,
                                .max_retries = o->retries ? o->retries
                                                          : FCT_RHO_DEFAULT_MAX_RETRIES,
                                .seed = o->seed ? o->seed : FCT_RHO_DEFAULT_SEED});
    rs[2] = fct_pollard_p_minus_one(n,
                                    &(fct_pminus1_config){
                                        .smoothness_bound =
                                            o->b ? o->b
                                                 : FCT_PM1_DEFAULT_SMOOTHNESS_BOUND,
                                        .start_base = 0,
                                        .max_bases = o->bases ? o->bases
                                                              : FCT_PM1_DEFAULT_MAX_BASES,
                                        .max_modmul_ops = 0});

    /* Overall status: any invalid => invalid; all found => ok; else none. */
    fct_status overall = FCT_OK;
    for (int i = 0; i < 3; i++) {
        if (rs[i].status == FCT_INVALID_INPUT)
            overall = FCT_INVALID_INPUT;
        else if (rs[i].status != FCT_OK && overall != FCT_INVALID_INPUT)
            overall = FCT_NO_FACTOR;
    }

    if (o->json) {
        jbuf b;
        jb_init(&b);
        jb_obj_begin(&b);
        jb_kv_str(&b, "command", "compare");
        jb_kv_u64(&b, "n", n);
        jb_kv_str(&b, "status", fct_status_str(overall));
        jb_key(&b, "results");
        jb_arr_result(&b, names, rs, 3);
        jb_obj_end(&b);
        int rc = jb_flush(&b);
        jb_free(&b);
        if (rc != 0)
            return rc;
    } else {
        printf("n = %" U64_FMT "\n\n", n);
        printf("%-10s %-12s %-17s %-10s %-10s %11s %10s %11s\n",
               "algorithm", "status", "reason", "factor", "cofactor",
               "iterations", "gcd_calls", "modmul_ops");
        for (int i = 0; i < 3; i++) {
            printf("%-10s %-12s %-17s ", names[i],
                   fct_status_str(rs[i].status),
                   fct_termination_str(rs[i].reason));
            if (rs[i].status == FCT_OK) {
                printf("%-10" U64_FMT " %-10" U64_FMT " ", rs[i].factor,
                       rs[i].cofactor);
            } else {
                printf("%-10s %-10s ", "-", "-");
            }
            printf("%11" U64_FMT " %10" U64_FMT " %11" U64_FMT "\n",
                   rs[i].iterations, rs[i].gcd_calls, rs[i].modmul_ops);
        }
    }
    return exit_code_for_status(overall);
}

/* ------------------------------------------------------------------ */
/* verify                                                             */
/* ------------------------------------------------------------------ */

static int cmd_verify(uint64_t n, const uint64_t *fs, size_t k, bool json)
{
    if (n < 2)
        return report_error("verify requires n >= 2", json, "verify");
    if (k == 0)
        return report_error("verify requires at least one factor", json, "verify");

    uint64_t prod = 1;
    bool overflow = false;
    for (size_t i = 0; i < k; i++) {
        if (!fct_mul_checked(prod, fs[i], &prod)) {
            overflow = true;
            break;
        }
    }
    bool matches = !overflow && prod == n;
    bool all_prime = true;
    for (size_t i = 0; i < k; i++) {
        if (!fct_is_prime(fs[i]))
            all_prime = false;
    }
    bool valid = matches && all_prime;

    if (json) {
        jbuf b;
        jb_init(&b);
        jb_obj_begin(&b);
        jb_kv_str(&b, "command", "verify");
        jb_kv_u64(&b, "n", n);
        jb_kv_str(&b, "status", valid ? "ok" : "reconstruction_failed");
        jb_key(&b, "factors");
        jb_arr_u64(&b, fs, k);
        jb_kv_u64(&b, "product", overflow ? 0 : prod);
        jb_kv_bool(&b, "product_matches", matches);
        jb_kv_bool(&b, "all_prime", all_prime);
        jb_kv_bool(&b, "valid", valid);
        jb_obj_end(&b);
        int rc = jb_flush(&b);
        jb_free(&b);
        if (rc != 0)
            return rc;
    } else {
        printf("command  : verify\n");
        printf("n        : %" U64_FMT "\n", n);
        printf("factors  :");
        for (size_t i = 0; i < k; i++)
            printf(" %" U64_FMT, fs[i]);
        printf("\n");
        printf("product  : %" U64_FMT "%s\n", overflow ? 0 : prod,
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
    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    if (argc < 2) {
        fprintf(stderr, "error: missing command\n\n");
        print_usage(stderr);
        return 2;
    }
    const char *cmd = argv[1];

    if (strcmp(cmd, "help") == 0 || strcmp(cmd, "--help") == 0 ||
        strcmp(cmd, "-h") == 0) {
        print_usage(stdout);
        return 0;
    }

    if (strcmp(cmd, "factor") == 0 || strcmp(cmd, "compare") == 0) {
        bool json = json_flag_present(argc, argv);
        options o;
        options_init(&o);
        uint64_t n = 0;
        bool have_n = false;
        int rc = parse_factor_args(argc, argv, &o, &n, &have_n, cmd, json);
        if (rc >= 0)
            return rc;

        rc = (strcmp(cmd, "factor") == 0)
                 ? (o.algo == ALGO_AUTO ? cmd_factor_auto(n, &o)
                                        : cmd_factor_single(n, &o))
                 : cmd_compare(n, &o);

        if (o.show_time && !o.json) {
            struct timespec t1;
            clock_gettime(CLOCK_MONOTONIC, &t1);
            printf("elapsed_ms: %.3f\n", elapsed_ms(&t0, &t1));
        }
        return rc;
    }

    if (strcmp(cmd, "verify") == 0) {
        bool json = json_flag_present(argc, argv);
        uint64_t values[1024];
        size_t count = 0;
        bool have_n = false;
        uint64_t n = 0;
        for (int i = 2; i < argc; i++) {
            const char *a = argv[i];
            if (strcmp(a, "--json") == 0 || strcmp(a, "--time") == 0)
                continue;
            uint64_t v;
            if (!parse_u64(a, &v))
                return report_error("argument is not a non-negative integer",
                                    json, "verify");
            if (!have_n) {
                n = v;
                have_n = true;
                continue;
            }
            if (count >= sizeof(values) / sizeof(values[0]))
                return report_error("too many factors (max 1024)", json, "verify");
            values[count++] = v;
        }
        if (!have_n)
            return report_error("missing required argument <n>", json, "verify");
        return cmd_verify(n, values, count, json);
    }

    fprintf(stderr, "error: unknown command '%s'\n\n", cmd);
    print_usage(stderr);
    return 2;
}
