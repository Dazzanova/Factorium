#!/usr/bin/env python3
"""Validate Factorium CLI JSON output with the Python 3 standard library.

Usage: check_cli_json.py <path-to-factorium-binary>

Checks:
  - `factor`/`compare`/`verify` JSON is syntactically valid and carries the
    fixed schema keys (unknown keys are rejected, so schema drift is caught).
  - Exit codes: 0 ok, 1 no factor, 2 invalid input, 3 incomplete.
  - Numeric fields are JSON integers (not floats/strings) and non-negative.
  - verify 17179869194064771037 = 4294967291 * 4000000007 validates exactly.
"""
import json
import subprocess
import sys

FAILURES = 0


def fail(msg: str) -> None:
    global FAILURES
    FAILURES += 1
    print(f"FAIL: {msg}", file=sys.stderr)


def run(binary, args):
    proc = subprocess.run([binary, *args], capture_output=True, text=True)
    return proc.returncode, proc.stdout, proc.stderr


def parse_json(binary, args, expect_rc):
    rc, out, err = run(binary, args)
    cmd = " ".join(args)
    if rc != expect_rc:
        fail(f"`{cmd}` exit code {rc}, expected {expect_rc}; stderr={err!r}")
        return None
    try:
        doc = json.loads(out)
    except json.JSONDecodeError as e:
        fail(f"`{cmd}` produced invalid JSON: {e}; stdout={out!r}")
        return None
    if not isinstance(doc, dict):
        fail(f"`{cmd}` JSON root is not an object: {type(doc).__name__}")
        return None
    return doc


def check_keys(doc, required, allowed, label):
    for k in required:
        if k not in doc:
            fail(f"{label}: missing required key {k!r}")
    for k in doc:
        if k not in allowed:
            fail(f"{label}: unexpected key {k!r} (schema drift?)")


def check_u64(doc, key, label):
    v = doc.get(key)
    if not isinstance(v, int) or isinstance(v, bool) or v < 0:
        fail(f"{label}: field {key!r} must be a non-negative integer, got {v!r}")
        return None
    if v >= 2**64:
        fail(f"{label}: field {key!r} exceeds uint64 range: {v}")
    return v


RESULT_KEYS = {
    "algorithm", "status", "reason", "factor", "cofactor",
    "iterations", "gcd_calls", "modmul_ops", "seed", "c",
}
STATUSES = {"ok", "no_factor", "invalid_input", "work_limit", "incomplete",
            "reconstruction_failed", "memory", "internal"}


def main(binary):
    # --- help ------------------------------------------------------------
    rc, out, _ = run(binary, ["help"])
    if rc != 0 or "Usage:" not in out:
        fail(f"help: exit={rc}, Usage banner present={'Usage:' in out}")

    # --- factor: complete factorization ----------------------------------
    doc = parse_json(binary, ["factor", "360", "--json"], 0)
    if doc:
        check_keys(
            doc,
            {"command", "n", "algorithm", "status", "reason", "factors",
             "unresolved", "product", "product_verified", "iterations",
             "gcd_calls", "modmul_ops"},
            {"command", "n", "algorithm", "status", "reason", "factors",
             "unresolved", "product", "product_verified", "iterations",
             "gcd_calls", "modmul_ops"},
            "factor/360",
        )
        if doc.get("factors") != [2, 2, 2, 3, 3, 5]:
            fail(f"factor/360 factors wrong: {doc.get('factors')!r}")
        if doc.get("product_verified") is not True or doc.get("product") != 360:
            fail(f"factor/360 product verification wrong: {doc!r}")
        check_u64(doc, "product", "factor/360")
        if doc.get("status") != "ok" or doc.get("algorithm") != "auto":
            fail(f"factor/360 status/algorithm wrong: {doc!r}")

    # --- factor: single algorithm ----------------------------------------
    doc = parse_json(binary, ["factor", "91", "--algo", "fermat", "--json"], 0)
    if doc:
        check_keys(
            doc,
            {"command", "n", "algorithm", "status", "reason", "factor",
             "cofactor", "iterations", "gcd_calls", "modmul_ops", "seed", "c"},
            {"command", "n", "algorithm", "status", "reason", "factor",
             "cofactor", "iterations", "gcd_calls", "modmul_ops", "seed", "c"},
            "factor/91-fermat",
        )
        if doc.get("factor") != 7 or doc.get("cofactor") != 13:
            fail(f"factor/91-fermat split wrong: {doc!r}")
        if doc.get("algorithm") != "fermat" or doc.get("status") != "ok":
            fail(f"factor/91-fermat status/algorithm wrong: {doc!r}")
        for k in ("n", "factor", "cofactor", "iterations", "gcd_calls",
                  "modmul_ops", "seed", "c"):
            check_u64(doc, k, "factor/91-fermat")

    # --- factor: no factor within budget -> exit 1 ------------------------
    doc = parse_json(
        binary,
        ["factor", "9999999967", "--algo", "fermat", "--iters", "1000", "--json"],
        1,
    )
    if doc:
        if doc.get("status") != "no_factor":
            fail(f"budget case status wrong: {doc!r}")
        if doc.get("iterations") != 1000:
            fail(f"budget case iterations wrong: {doc!r}")
        if doc.get("factor") != 0 or doc.get("cofactor") != 0:
            fail(f"no-factor result must report factor=0: {doc!r}")

    # --- factor: invalid input -> exit 2 ----------------------------------
    doc = parse_json(binary, ["factor", "abc", "--json"], 2)
    if doc and doc.get("status") != "invalid_input":
        fail(f"invalid input status wrong: {doc!r}")
    parse_json(binary, ["factor", "1", "--json"], 2)
    parse_json(binary, ["factor", "15", "--algo", "bogus", "--json"], 2)
    parse_json(binary, ["factor", "--json"], 2)

    # --- factor: incomplete -> exit 3 -------------------------------------
    doc = parse_json(
        binary,
        ["factor", "17179869194064771037", "--iters", "1", "--retries", "1",
         "--b", "2", "--bases", "1", "--json"],
        3,
    )
    if doc:
        if doc.get("status") != "incomplete":
            fail(f"incomplete case status wrong: {doc!r}")
        if doc.get("unresolved") != [17179869194064771037]:
            fail(f"incomplete case unresolved wrong: {doc.get('unresolved')!r}")
        if doc.get("product_verified") is not False:
            fail(f"incomplete case must not claim verification: {doc!r}")

    # --- compare ----------------------------------------------------------
    doc = parse_json(binary, ["compare", "15", "--json"], 0)
    if doc:
        check_keys(doc, {"command", "n", "status", "results"},
                   {"command", "n", "status", "results"}, "compare/15")
        results = doc.get("results")
        if not isinstance(results, list) or len(results) != 3:
            fail(f"compare/15 results must be a list of 3, got {results!r}")
        else:
            algos = [r.get("algorithm") for r in results]
            if algos != ["fermat", "rho", "pminus1"]:
                fail(f"compare/15 algorithm order wrong: {algos!r}")
            for r in results:
                if not isinstance(r, dict):
                    fail(f"compare/15 result not an object: {r!r}")
                    continue
                for k in RESULT_KEYS:
                    if k not in r:
                        fail(f"compare/15 result missing key {k!r}: {r!r}")
                if r.get("status") not in STATUSES:
                    fail(f"compare/15 bad status {r.get('status')!r}")
                for k in ("factor", "cofactor", "iterations", "gcd_calls",
                          "modmul_ops", "seed", "c"):
                    if k in r:
                        v = r[k]
                        if not isinstance(v, int) or isinstance(v, bool) or v < 0:
                            fail(f"compare/15 {k} must be uint, got {v!r}")

    parse_json(binary, ["compare", "3", "--json"], 2)  # prime input

    # --- verify -----------------------------------------------------------
    doc = parse_json(
        binary,
        ["verify", "17179869194064771037", "4294967291", "4000000007", "--json"],
        0,
    )
    if doc:
        check_keys(doc,
                   {"command", "n", "status", "factors", "product",
                    "product_matches", "all_prime", "valid"},
                   {"command", "n", "status", "factors", "product",
                    "product_matches", "all_prime", "valid"},
                   "verify/semiprime")
        if doc.get("valid") is not True:
            fail(f"verify/semiprime should be valid: {doc!r}")
        if doc.get("product") != 17179869194064771037:
            fail(f"verify/semiprime product wrong: {doc!r}")

    doc = parse_json(binary, ["verify", "360", "2", "2", "2", "3", "3", "6", "--json"], 1)
    if doc and (doc.get("valid") is not False or doc.get("all_prime") is not False):
        fail(f"verify with composite factor wrong: {doc!r}")

    parse_json(binary, ["verify", "10", "--json"], 2)  # no factors given

    if FAILURES:
        print(f"check_cli_json: {FAILURES} failure(s)", file=sys.stderr)
        return 1
    print("check_cli_json: all checks passed")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <factorium-binary>", file=sys.stderr)
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
