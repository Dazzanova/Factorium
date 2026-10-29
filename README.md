# Factorium — Comparative Integer Factorization Laboratory (Phase 1 & Phase 2)

**Factorium** is an integer factorization engine and experimental benchmarking laboratory implemented in C17, designed for a 5th-semester *Design and Analysis of Algorithms* course. It compares three classical factorization algorithms across varying composite structures:

1. **Fermat's Factorization Method** (Difference of Squares)
2. **Pollard's $\rho$ (Rho) Algorithm** (Floyd's Cycle Finding)
3. **Pollard's $p-1$ Algorithm** (Stage 1 Smoothness Method)

The project emphasizes mathematical correctness, strict overflow safety across the full 64-bit domain ($0 \le n < 2^{64}$), zero undefined behavior, determinism, modular architecture, automated test coverage, and reproducible experimental analysis.

---

## Architecture and Repository Structure

The core library (`factorium_core`) contains no I/O, no timing logic, and no memory leaks. The CLI (`factorium`) handles interactive commands and JSON output. The benchmark runner (`factorium_bench`) coordinates controlled dataset trials, monotonic timing (`CLOCK_MONOTONIC`), statistics aggregation, and RFC 4180 CSV export. The visualization pipeline (`scripts/plot_benchmarks.py`) plots empirical results.

```
Factorium/
├── CMakeLists.txt                      # Build configuration (C17, strict warnings, CTest)
├── README.md                           # Comprehensive documentation and technical manual
├── .gitignore                          # Clean build, temporary, and virtualenv exclusions
├── include/
│   └── factorium/
│       ├── common.h                    # Result structures, status codes, and termination reasons
│       ├── modular_arithmetic.h        # Overflow-safe modular operations and square roots
│       ├── primality.h                 # Deterministic 64-bit Miller-Rabin test
│       ├── fermat.h                    # Fermat difference-of-squares interface
│       ├── pollard_rho.h               # Pollard Rho interface
│       ├── pollard_p_minus_one.h       # Pollard p-1 stage 1 interface
│       ├── factorization.h             # Recursive complete-factorization engine
│       └── benchmark.h                 # Benchmark suite, dataset definitions, and CSV export
├── src/
│   ├── common.c                        # String conversions for statuses and reasons
│   ├── modular_arithmetic.c            # GCD, modmul, modpow, isqrt, checked multiplication
│   ├── primality.c                     # Miller-Rabin implementation (12 proven witnesses)
│   ├── fermat.c                        # Fermat's algorithm implementation
│   ├── pollard_rho.c                   # Pollard Rho with Floyd cycle detection
│   ├── pollard_p_minus_one.c           # Pollard p-1 stage 1 with Eratosthenes sieve
│   ├── factorization.c                 # Recursive composition pipeline and verification
│   ├── benchmark.c                     # Controlled datasets, trial runner, stats, and CSV writer
│   ├── bench_main.c                    # Dedicated benchmark runner binary (factorium_bench)
│   └── main.c                          # Core CLI binary: factor, compare, verify, JSON mode
├── tests/
│   ├── test_harness.h                  # Lightweight C test harness macro suite
│   ├── test_modular_arithmetic.c       # Tests for modular arithmetic boundary cases
│   ├── test_primality.c                # Miller-Rabin verification against known primes/composites
│   ├── test_fermat.c                   # Fermat tests (close factors, squares, budget caps)
│   ├── test_pollard_rho.c              # Pollard Rho tests (cycles, retries, reproducibility)
│   ├── test_pollard_p_minus_one.c      # Pollard p-1 tests (smooth/unfavorable p-1, bases)
│   ├── test_factorization.c            # Complete recursive factorization test suite
│   ├── test_cli.c                      # Integration tests for CLI arguments and exit codes
│   └── test_benchmark.c                # Unit and integration tests for benchmark harness & CSV
├── scripts/
│   ├── check_cli_json.py               # Python 3 stdlib JSON schema validator test
│   └── plot_benchmarks.py              # Publication-grade plotting script (pandas + matplotlib)
└── benchmarks/
    └── plots/                          # Generated empirical benchmark figures (PNG)
        ├── runtime_vs_bit_length.png
        ├── fermat_vs_factor_distance.png
        ├── rho_iteration_distribution.png
        ├── pminus1_vs_smoothness_bound.png
        └── success_timeout_rates.png
```

---

## Mathematical Foundations & Algorithms

### 1. Fermat's Factorization
- **Theory**: Based on the identity:
  $$n = a^2 - b^2 = (a - b)(a + b)$$
  Starting from $a = \lceil \sqrt{n} \rceil$, compute $b^2 = a^2 - n$. If $b^2$ is a perfect square, then $p = a - b$ and $q = a + b$ are nontrivial factors of $n$. At each step, $a \leftarrow a + 1$ and $b^2$ is updated incrementally as $b^2_{new} = b^2 + 2a + 1$ without division.
- **Factor Closeness**: For $n = p \cdot q$ with $p \le q$, the search succeeds at $a = \frac{p + q}{2}$ after $\frac{q - p}{2} + 1$ iterations. When factors are close ($p \approx q$, e.g. difference of a few thousands), Fermat finds the split almost instantly. When factors are distant (e.g. $p = 3$ and $q$ large), Fermat degenerates into trial division taking $O(n)$ steps.

### 2. Pollard's Rho Algorithm
- **Theory**: Pollard's $\rho$ algorithm models pseudorandom walks generated by a polynomial $f(x) = (x^2 + c) \pmod n$. By the Birthday Paradox, the sequence modulo a divisor $p$ of $n$ is expected to enter a cycle after $O(\sqrt{p}) \le O(n^{1/4})$ iterations.
- **Cycle Detection**: Floyd's cycle-finding algorithm maintains a tortoise $x \leftarrow f(x)$ and a hare $y \leftarrow f(f(y))$. At each step, $d = \gcd(|x - y|, n)$ is evaluated.
  - $d = 1$: No factor found yet; continue iteration.
  - $1 < d < n$: Nontrivial factor found!
  - $d = n$: The cycle collapsed ($x \equiv y \pmod n$); the current polynomial failed. Retries are triggered deterministically with shifted seed and updated $c$.

### 3. Pollard's $p-1$ Algorithm (Stage 1)
- **Theory**: Derived from Fermat's Little Theorem: if $p \mid n$ and $a$ is coprime to $p$, then $a^{p-1} \equiv 1 \pmod p$. If all prime factors of $p - 1$ are $\le B$ ($B$-smooth), then $p - 1 \mid M$, where:
  $$M = \prod_{q \le B, q \text{ prime}} q^{\lfloor \log_q B \rfloor}$$
  Consequently, $a^M \equiv 1 \pmod p$, implying $p \mid \gcd(a^M - 1, n)$.
- **Implementation**: Primes up to smoothness bound $B$ are generated using a byte Sieve of Eratosthenes. The accumulator $x \leftarrow x^{q^k} \pmod n$ is computed incrementally per prime power.
- **Dependence on Smoothness**: If $p - 1$ has even a single large prime factor exceeding $B$, Stage 1 will not identify $p$. If both $p-1$ and $q-1$ are smooth, $a^M \equiv 1 \pmod n$, leading to $\gcd(a^M - 1, n) = n$; in this case, another base $a$ is attempted.

### 4. Deterministic Miller-Rabin Primality Test
- Rather than probabilistic testing, Factorium uses deterministic Miller-Rabin with the 12 prime witnesses:
  $$\{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37\}$$
- **Correctness**: Proven (Pomerance, Selfridge, Wagstaff, Sorenson, Webster) to be deterministic for all $n < 318{,}665{,}857{,}834{,}031{,}151{,}167{,}461 > 2^{64}$. Thus, primality testing is exact across the complete 64-bit integer domain.
- **Separation of Concerns**: A failed search in Fermat, Rho, or $p-1$ is reported as `FCT_NO_FACTOR`, never as proof that $n$ is prime.

### 5. Overflow Safety and 64-bit Arithmetic
- **Compiler Extension**: On platforms where `__SIZEOF_INT128__` is supported (GCC/Clang on x86-64, aarch64), modular multiplication $(a \cdot b) \pmod m$ and Fermat's quadratic accumulator use `unsigned __int128` to prevent overflow when $a \cdot b \ge 2^{64}$.
- **Portable Fallback**: If compiled on an architecture without `__int128`, a portable double-and-add (Russian peasant) modular multiplication loop is used, and Fermat checks bounds before overflow.
- **Factor Reconstruction**: Factors are multiplied back using `fct_mul_checked` (`uint64_t` checked multiplication) to rigorously verify $n = \prod p_i^{k_i}$.

---

## Controlled Benchmark Datasets & Hypotheses

Factorium defines 36 curated, ground-truth-verified dataset items across 7 algorithmic categories:

| Category | Count | Ground-Truth Characteristics | Tested Hypothesis |
|:---|:---:|:---|:---|
| **`close_semiprimes`** | 7 | $|q - p| \in [2, 10000]$, $16 \le \text{bits} \le 64$ | Fermat resolves in $O((q-p)/2)$ steps; vastly outperforms Rho and $p-1$ when factors are close. |
| **`distant_semiprimes`** | 5 | $p \ll q$ (e.g. $p$ is 8-16 bit, $q$ is large) | Pollard Rho discovers $p$ in $O(\sqrt{p})$ steps; Fermat degrades toward $O(n)$ trial division. |
| **`small_prime_factors`** | 4 | Multiples of $2, 3, 5$, Primorial, $2^{40}$, smooth | Trivial factors are resolved immediately; Fermat handles evens at step 0, Rho finds small primes in $\le 3$ steps. |
| **`pminus1_favorable`** | 4 | $p-1$ is $B$-smooth ($B \le 5$, e.g. $2^{16}+1$) | Pollard $p-1$ extracts factor in a single pass with minimal operations; Fermat fails on distant factors. |
| **`pminus1_unfavorable`** | 3 | Every factor's $r - 1$ has a prime factor $> 100{,}000$; $p$ is a safe prime $2p' + 1$ with $p' > 100{,}000$ in all three, and $q$ is safe in two | Pollard $p-1$ stage 1 exhausts all bases (0% success); demonstrates failure under non-smooth factors. Rho succeeds (100%). |
| **`bit_length_scaling`** | 8 | Balanced semiprimes from 16 to 63 bits | Pollard Rho execution time scales with $O(n^{1/4})$; Fermat requires exponential iterations unless factors are adjacent. |
| **`edge_cases`** | 5 | Perfect squares ($p^2$), Carmichael numbers, prime, $2^{64}-2$ | Algorithms handle boundary conditions without integer overflow, false positive primality, or infinite loops. |

---

## Building and Testing

### Prerequisites
- C17 compatible compiler (GCC $\ge$ 9 or Clang $\ge$ 10)
- CMake $\ge$ 3.16
- Python $\ge$ 3.8 (for validation and plotting)

### Standard Build & CTest
```bash
# Configure and build
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build

# Execute automated test suite (10 test targets including smoke benchmark)
ctest --test-dir build --output-on-failure
```

### Sanitizers Build (ASan + UBSan)
```bash
cmake -B build-san -S . -DFACTORIUM_SANITIZERS=ON
cmake --build build-san
ctest --test-dir build-san --output-on-failure
```

---

## CLI Usage (`factorium`)

The `factorium` binary provides four commands: `factor`, `compare`, `verify`, and `help`.

```
Usage:
  factorium factor <n> [options]      factor <n> (complete recursive by default)
  factorium compare <n> [options]     run Fermat, Rho, and p-1 on <n>
  factorium verify <n> <f1>...<fk>    validate a prime factorization
  factorium help                      show help banner

Options:
  --algo <a>    factor only: auto (default), fermat, rho, pminus1
  --iters <k>   iteration cap for Fermat/Rho (0 = default 1048576)
  --retries <r> Rho retry cap (0 = default 8)
  --b <B>       p-1 smoothness bound B (0 = default 100000)
  --bases <m>   p-1 max consecutive bases (0 = default 5)
  --seed <s>    Rho initial seed (0 = default 1); same seed => same run
  --json        machine-readable JSON output (fixed schema)
  --time        display elapsed wall time (human output only)
```

### Examples

#### 1. Complete Recursive Factorization
```bash
$ ./build/factorium factor 360
command  : factor
n        : 360
algorithm: auto
status   : ok
reason   : none
factors  : 2 2 2 3 3 5
product  : 360 (verified)
iterations: 0
gcd_calls: 0
modmul_ops: 0
```

#### 2. Single Algorithm Mode
```bash
$ ./build/factorium factor 9999999967 --algo fermat --iters 1000
command  : factor
n        : 9999999967
algorithm: fermat
status   : no_factor
reason   : budget_exhausted
factor   : -
cofactor : -
iterations: 1000
gcd_calls: 0
modmul_ops: 0
```
Exit code `1`: no factor within the budget — a normal outcome of a bounded search, **not** a primality claim.

#### 3. Comparing Algorithms on a Distant-Factor Semiprime
```bash
$ ./build/factorium compare 17179869194064771037 --json
{"command":"compare","n":17179869194064771037,"status":"no_factor","results":[{"algorithm":"fermat","status":"no_factor","reason":"budget_exhausted","factor":0,"cofactor":0,"iterations":1048576,"gcd_calls":0,"modmul_ops":0,"seed":0,"c":0},{"algorithm":"rho","status":"ok","reason":"found","factor":4000000007,"cofactor":4294967291,"iterations":1119523,"gcd_calls":1119523,"modmul_ops":3358569,"seed":3,"c":1},{"algorithm":"pminus1","status":"ok","reason":"found","factor":4000000007,"cofactor":4294967291,"iterations":1,"gcd_calls":2,"modmul_ops":212451,"seed":0,"c":0}]}
```
Here Fermat exhausts its budget (the factors differ by ~2.9×10⁸), Pollard Rho succeeds, and Pollard p−1 succeeds in a single base pass (4000000006 = 4000000007 − 1 is 100000-smooth). Overall status is `no_factor` (exit code `1`) because not every algorithm succeeded.

#### 4. Validating a Factorization
```bash
$ ./build/factorium verify 17179869194064771037 4294967291 4000000007
command  : verify
n        : 17179869194064771037
factors  : 4294967291 4000000007
product  : 17179869194064771037
matches  : yes
all prime: yes
valid    : yes
```

---

## Exit Codes

| Code | Status | Meaning |
|:---:|:---|:---|
| `0` | Success | Factor discovered, complete factorization verified, or valid verification. |
| `1` | No Factor | Normal negative outcome: search exhausted budget without finding a factor (or not all algorithms succeeded in `compare`). |
| `2` | Invalid Input | Input outside domain ($n < 2$ for factor/verify, $n < 4$ for single-algorithm split), invalid syntax, or argument error. |
| `3` | Incomplete | Complete factorization ended with unresolved composite cofactors. |
| `4` | Work Limit | Hard resource cap reached mid-search (e.g. max modmuls in $p-1$). |
| `5` | Internal / Error | Memory allocation failure, reconstruction mismatch, or internal error. |

---

## JSON Schemas

### `factor` (Auto / Complete Mode)
```json
{
  "command": "factor",
  "n": 360,
  "algorithm": "auto",
  "status": "ok",
  "reason": "none",
  "factors": [2, 2, 2, 3, 3, 5],
  "unresolved": [],
  "product": 360,
  "product_verified": true,
  "iterations": 0,
  "gcd_calls": 0,
  "modmul_ops": 0
}
```

### `factor` (Single Algorithm) / `compare` result entries
```json
{
  "command": "factor",
  "n": 91,
  "algorithm": "fermat",
  "status": "ok",
  "reason": "found",
  "factor": 7,
  "cofactor": 13,
  "iterations": 1,
  "gcd_calls": 0,
  "modmul_ops": 0,
  "seed": 0,
  "c": 0
}
```
`compare` wraps the three single-algorithm result objects in a `results` array (order: `fermat`, `rho`, `pminus1`).

### `verify`
```json
{
  "command": "verify",
  "n": 17179869194064771037,
  "status": "ok",
  "factors": [4294967291, 4000000007],
  "product": 17179869194064771037,
  "product_matches": true,
  "all_prime": true,
  "valid": true
}
```

---

## Benchmark Runner (`factorium_bench`)

The `factorium_bench` binary runs comparative factorization experiments, measures runtimes with nanosecond resolution via `CLOCK_MONOTONIC`, repeats randomized Pollard Rho trials with controlled seeds, computes summary statistics (min, max, median, Q1, Q3, success rate), and exports RFC 4180 compliant CSV files.

### CLI Syntax
```bash
factorium_bench [options]

Options:
  --smoke              run fast smoke benchmark subset (<0.2s)
  --full               run full dataset benchmark (default)
  --raw-out <file>     destination path for raw CSV records (default: bench_raw.csv)
  --summary-out <file> destination path for summary CSV (default: bench_summary.csv)
  --rho-trials <n>     number of randomized Rho seeds (default: 5, smoke: 1)
  --iters <n>          override max iterations budget
  --bound <B>          override p-1 smoothness bound
  --quiet              suppress console status output
  --help               display help message
```

### Running Full Benchmarks
```bash
# Run full benchmark suite with 5 Rho seeds per input
./build/factorium_bench --full --raw-out bench_raw.csv --summary-out bench_summary.csv
```

---

## CSV Export Schemas

### 1. Raw Trials CSV (`bench_raw.csv`)
Records every individual execution trial:
- `trial_id`: Sequential trial counter
- `dataset_id`: Unique identifier (e.g. `close_32_twin`)
- `category`: Dataset category
- `n`: Input integer
- `bit_length`: Bit length of $n$
- `ground_truth_p`, `ground_truth_q`: True factors
- `algorithm`: `fermat`, `rho`, or `pminus1`
- `param_iters`, `param_retries`, `param_seed`, `param_b`, `param_bases`, `param_max_modmul`: Algorithm parameters
- `status`: `ok`, `no_factor`, `invalid_input`, `work_limit`
- `reason`: `found`, `budget_exhausted`, `bases_exhausted`, etc.
- `factor_found`, `cofactor_found`: Factors identified by algorithm
- `factor_is_valid`: Boolean verification against $n$ and ground truth
- `runtime_ns`: Elapsed monotonic clock time in nanoseconds
- `iterations`, `gcd_calls`, `modmul_ops`: Operation counters

### 2. Summary Statistics CSV (`bench_summary.csv`)
Aggregates trials grouped by `(dataset_id, algorithm)`:
- `dataset_id`, `category`, `algorithm`: Identifying keys
- `num_trials`: Number of executed trials
- `num_success`, `num_timeout`, `num_fail`: Outcome counters
- `success_rate`: Ratio of successes ($\in [0.0, 1.0]$)
- `median_runtime_ns`, `min_runtime_ns`, `max_runtime_ns`: Runtime percentiles
- `q1_runtime_ns`, `q3_runtime_ns`: Interquartile range (25th and 75th percentiles)
- `median_iterations`, `median_gcd_calls`, `median_modmul_ops`: Median workload metrics

---

## Visualization & Experimental Plots

The Python visualization script (`scripts/plot_benchmarks.py`) reads the generated CSV datasets and produces 5 figures saved to `benchmarks/plots/`:

```bash
# Set up Python environment
python3 -m venv .venv
.venv/bin/pip install pandas matplotlib

# Generate all 5 analysis plots
.venv/bin/python3 scripts/plot_benchmarks.py
```

### Generated Figures
1. **`runtime_vs_bit_length.png`**:
   - Plots median runtime ($\mu$s, logarithmic scale) against input bit length ($16 \le \text{bits} \le 63$).
   - Shows asymptotic growth rates across balanced semiprimes.
2. **`fermat_vs_factor_distance.png`**:
   - Compares Fermat iterations and runtime against factor distance $|q - p|$.
   - Verifies the empirical $O((q-p)/2)$ theoretical line against experimental observations.
3. **`rho_iteration_distribution.png`**:
   - Boxplots displaying Pollard Rho iteration variability across seeds $s \in \{1, 2, 3, 4, 5\}$.
   - Illustrates variance inherent in pseudorandom cycle detection.
4. **`pminus1_vs_smoothness_bound.png`**:
   - Plots runtime and modular multiplication work as smoothness bound $B$ sweeps from $10$ to $600{,}000$.
   - Contrasts $B$-smooth inputs (instant resolution) with safe primes (linear work growth without factor discovery).
5. **`success_timeout_rates.png`**:
   - Categorical bar chart illustrating success vs. timeout/exhaustion percentages for each algorithm across the 7 dataset categories.

---

## Experimental Observations & Summary

1. **Fermat Dominance on Close Semiprimes**:
   When factors are adjacent ($|q - p| \le 30$), Fermat finds the nontrivial factor in $1$ to $15$ iterations ($< 2\,\mu\text{s}$), outperforming both Pollard Rho and $p-1$ by multiple orders of magnitude regardless of bit length.
2. **Pollard Rho Superiority on Unbalanced Semiprimes**:
   When $p \ll q$, Fermat exhausts even a $500{,}000$ iteration budget, whereas Pollard Rho resolves the small factor in $O(\sqrt{p})$ steps ($< 10\,\mu\text{s}$).
3. **Pollard $p-1$ Precision**:
   On numbers where $p-1$ is $B$-smooth (e.g. $p = 65537 = 2^{16}+1$), Pollard $p-1$ succeeds with $B=2$ in $1$ step. On safe primes where $p-1$ contains a large prime cofactor ($> 100{,}000$), Stage 1 systematically exhausts all bases without false factorization.
4. **Primality and Edge Resiliency**:
   Prime inputs ($1000000007$) correctly report `no_factor` within budget limits without hanging or generating false claims, and perfect squares are resolved immediately.

---

## Known Scope & Limitations

1. **Domain**: Inputs are unsigned 64-bit integers ($0 \le n < 2^{64}$).
2. **Algorithm Boundaries**:
   - Fermat requires $|p - q| \ll \sqrt{n}$ to finish within standard budgets.
   - Pollard $p-1$ implements Stage 1 only; it requires $p-1$ to be $B$-smooth for some factor $p \mid n$.
   - Hard semiprimes with two distant random 32-bit primes whose $p-1$ values have large prime factors may require $O(2^{16})$ to $O(2^{32})$ work; they need realistic search budgets for Pollard Rho or advanced algorithms (ECM, Quadratic Sieve) not covered here.
3. **No Primality by Omission**: Failing to find a factor in a bounded search is never reported as primality. Deterministic Miller-Rabin is the sole arbiter of primality.
4. **Benchmark Statistics**: Median/Q1/Q3 use nearest-rank selection on sorted trial vectors; per-dataset $n$ is fixed, so runtime differences reflect algorithmic work, not input variation. Timing uses `CLOCK_MONOTONIC` on the host and is not portable across machines; the CSVs in the repository are reference runs from a single environment.
