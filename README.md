# Sums of two dependent uniform random variables: a characterization by joint tail bounds

This is the computational artifact for the paper *Sums of two dependent uniform random variables: a characterization
by joint tail bounds* by P. Chaudhary and K. Khandelwal, arXiv:2610.24197 [math.PR], 2026
(https://arxiv.org/abs/2610.24197; citation: `CITATION.cff`). The paper determines which laws the sum of two
standard uniform random variables can have under unknown dependence, by joint tail-sum (staircase) inequalities. This
archive holds code and data for every computational statement of the paper:

- the fibered reduction and two uniforms at every spacing; the corrections to the two- and three-point theorems of
  Mao, Wang and Wang (2019);
- the grid (1/k)Z, lattice laws as doubly stochastic X-rays, and the three readings of potential duality (convex
  order, folded order, single step), with the two-point example;
- windows and the three- and four-atom laws;
- the small-atoms corollary: the statement for k <= 7 and every step of its self-contained proof (grid duality, the
  convex-minorant lemma and the chain of inequalities) for k <= 12;
- tail-sum inequalities: repeated thresholds, the counting and reflection lemmas, and the step counts for two- and
  three-point laws;
- laws with a density: densities at most one (the construction of the proof), sharpness, the tightness of the lattice
  bound, exact certificates at density two (single steps, nonincreasing gaps, equally spaced thresholds, two steps),
  an explicit boundary coupling, the failure of summed single-step certificates, deposited certificate catalogues,
  and a worked example with density four.

Every computational claim of the paper is listed in `CLAIM_EVIDENCE_MAP.md` with the script that checks it, its data
file, the tier that reaches the stated range and its status (full stream verified / small range / index only / not
checked). Results are cited by their printed numbers; numbers with the prefix S refer to the paper's supplementary material.

## Identifiers

- Cite this version (1.0.0): https://doi.org/10.5281/zenodo.23231387
- All versions: https://doi.org/10.5281/zenodo.23231386
- Preprint: arXiv:2610.24197, https://doi.org/10.48550/arXiv.2610.24197
- Published version: not yet accepted
- Repository: https://github.com/methodtrace/2026-two-uniform-sums, tag `v1.0.0`

## Requirements

- CPython 3.12 or newer, with `numpy`, `scipy` and `pytest`. `pytest` is required: `run.py` runs the test suite
  first, and the run fails without it.
- A C compiler is not required. The only C program, `tools/c/grid_steps.c`, is used by `check_repeated_thresholds.py`
  (Example 4.2(b)), which looks for `cc`, `gcc` or `clang`; if none is found it skips the C helper and runs the
  exact numpy implementation alone (with a compiler, both run and must agree). The check passes either way.
- Optional: `sympy` for the symbolic identities in `check_density_sharpness.py` and `check_two_step_density.py`
  (skipped and reported as `[SKIP]` without it), and `pycddlib` (package `pycddlib-standalone`, version 3 API) for the
  floating-point vertex enumeration in `check_small_atoms.py`. Without pycddlib the quick tier of
  `check_small_atoms.py` still runs in full (its vertices come from an exhaustive exact enumeration); the full tier
  skips k = 7, which needs pycddlib, and reports it as `[SKIP]`.

```sh
python3 -m pip install -r requirements.txt    # the versions the results were produced with
```

The scripts and tests run from a checkout without installation (they add `src/` to the import path).

## Verifying the results byte for byte

```sh
python3 -m pip install -r requirements.txt
python3 run_all.py
```

`run_all.py` runs the test suite and every script at its quick range (through `run.py`, in the scratch folder
`runs/run_all/`), writes each script's results to `results/<step>.json` in a canonical form, and compares their
SHA-256 digests with `MANIFEST.sha256`. A correct run ends with

```
all results match MANIFEST.sha256 byte for byte
```

The canonical form depends only on the code: keys are sorted, run times, output paths and worker counts are removed,
the record of which optional engines ran is removed (the C helper for Example 4.2(b) and pycddlib for the small-atoms
vertices; each is a cross-check that must agree with the exact computation when it runs), and floating-point numbers are written as 0.0 below 1e-9 in absolute value and to six
significant digits otherwise (floating point enters only in cross-checks; every statement of the paper is checked
exactly). `run_all.py` needs every package in `requirements.txt` (without `sympy` the symbolic checks are skipped,
which changes the results); `pycddlib` is not needed for it.

## Running

```sh
python3 run.py                    # quick tier: tests + every script at its default range
python3 run.py --tier medium      # larger ranges for some steps
python3 run.py --tier full        # every script with --full (the ranges stated in the paper)
python3 run.py --only small_atoms,worked_example --tier full     # selected steps
```

Measured on a 2-core Linux machine (`nproc` = 2) with CPython 3.13: the quick tier takes about 2.5 minutes of wall
time (145 to 147 s in two runs; the test suite about 40 s of it), the full tier about 10.5 minutes without the test suite (630 s; `--skip-tests`), about 11 minutes with it.

Every step writes `runs/<step>/<name>.json` with its counts and a pass/fail entry per check, and `run.py` writes
`runs/run_summary.json` (`runs/` is not part of the archive). Each script can also be run alone; `python3 scripts/<script>.py --help` lists its
options; every script takes `--out DIR`, small default ranges, and `--full` for the paper's range.

| Script | Reproduces | Quick tier | Full tier |
|---|---|---|---|
| `check_two_point_example.py` | Example 3.1 | seconds | same |
| `check_mww_steps.py` | Proposition 6.1 and the law after it | seconds | same |
| `check_window_corollary.py` | Corollary S3.1 against a generic LP (FLOAT), Lemma S3.2 | seconds | same |
| `count_window_vertices.py` | five-point window vertex counts | seconds | same |
| `check_small_atoms.py` | Corollary 3.2 (statement), Lemma 3.2, Lemma 3.3, the chain of the proof | statement k <= 5, chain k <= 8 (20 s) | statement k <= 7, chain k <= 12 (3 min) |
| `check_tail_bounds.py` | Lemma 4.3, Theorem 4.2, Theorem 4.3 | 30 s | 1 min |
| `check_staircase_lp.py` | Proposition S4.1, Lemma 4.1, Lemma 4.2, Lemma 4.3, Remark 4.2, Theorem 4.2, Theorem 4.3 by the exact staircase program | 10 s | 30 s |
| `check_repeated_thresholds.py` | Example 4.2 | 7 s | same |
| `check_density_construction.py` | Theorem 5.1 (the construction of the proof) | seconds | same |
| `check_density_sharpness.py` | Proposition 5.1; a FLOAT experiment with at most three thresholds | 15 s | 2 min |
| `check_lattice_tightness.py` | Remark S4.1 | 10 s | same |
| `check_single_step_density.py` | Proposition 5.2; Lemma 5.1 certificates for K = 8, 12, 16; FLOAT searches after Problem 5.1 | K = 8 (3 s) | K = 8, 12, 16, float search r <= 8 (3 min on 2 cores) |
| `check_two_step_density.py` | Theorem 5.4; symbolic case identities (sympy) | grid 1/40 (1 s) | grid 1/120 |
| `check_nonincreasing_gaps_density.py` | Theorem 5.2 | r <= 5 (1 s) | r <= 8 |
| `check_equally_spaced_density.py` | Theorem 5.3 | r <= 8 (1 s) | r <= 24 |
| `check_density_coupling.py` | Proposition 5.3 | seconds | same |
| `check_hinge_sum_failure.py` | after Proposition 5.4: the sum of single-step certificates fails | seconds | same |
| `check_density_catalogues.py` | after Problem 5.1: deposited certificates at density 2 | grids 8, 12, two-step 48 (5 s) | adds grid 16 |
| `check_worked_example.py` | Example 4.1 | seconds | same |

Times are indicative.

## Exactness

Every count and every verdict is computed in exact integer or rational arithmetic, with certificates:

- the facet test (used for single steps on hat-P(10), the facet of Example 4.2(a) and the facets of
  Remark S4.1) reports for each potential one of: no tight permutation; two anti-diagonals without allowed
  edge; an equation; a rank modulo the prime 2^31 - 1 reaching the facet dimension (a modular rank never exceeds the
  rational rank, and the incidence part is totally unimodular); or an exact rational recomputation;
- in `check_small_atoms.py` every vertex of Q(k) and every doubly stochastic matrix is exact and checked exactly; a
  floating-point linear program may propose the support of a matrix, which is then solved and verified in rational
  arithmetic (with a fallback to an exact simplex method); separating potentials come from an exact simplex method;
- the density certificates are checked in exact rational arithmetic, and the deposited catalogues record by record.

Floating point enters only in cross-checks that are labeled as such: the coupling and fibered linear programs
(`couplings.py`, `fibered.py`, `check_window_corollary.py`), the search for lattice witnesses (each witness is
verified exactly), the density experiments after Problem 5.1 (numerical evidence only), the linear programs
that propose the convex functions of the Lemma 5.1 certificates (each certificate is verified exactly), and
the cddlib vertex enumeration of Q(k) in `check_small_atoms.py` (every vertex is recovered and checked exactly; for
k <= 6 the vertex lists are also obtained by an exhaustive exact enumeration and compared, while for k = 7 the
completeness of the list rests on cddlib).

## Contents

- `src/dsxray/`: the package (one module per topic; see the module docstrings).
- `scripts/`: one script per group of claims, named by what it reproduces; `scripts/common.py` holds shared helpers.
  The docstring of each script lists the claims it checks.
- `tests/`: unit and regression tests (`python3 -m pytest -q`, about 40 seconds).
- `tools/c/grid_steps.c`: exhaustive check of all step sets on a grid (Example 4.2(b)), compiled on demand.
- `data/facets/facets_kNN.jsonl` (k = 9, 10, 11): the certified classification of every ray of the normal fan of
  hat-P(k) (one JSON object per ray: ray, dim Z, facet, certificate, total support, edges outside every tight
  permutation, empty anti-diagonals). Used only by the floating-point witness search of `check_lattice_tightness.py`.
- `data/density/density_two_certificates_K08.jsonl.gz`, `..._K12.jsonl.gz`, `..._K16.jsonl.gz` and
  `..._two_step_K48.jsonl.gz`: one exact certificate per staircase (fields `J`, `K`, `m`, node values `phi` of a convex
  function on the grid 1/(Km) as rational strings, `margin`), for every staircase with at least two distinct
  thresholds on the grids 1/8, 1/12, 1/16 and every two-step staircase on 1/48 (density cap 2).
- `results/`: the canonical results of `run_all.py`, one file per step.
- `MANIFEST.sha256`: SHA-256 digests of the files in `results/`, checked by `run_all.py`.
- `requirements.txt`: the package versions the results were produced with.

## License

Code (`src/`, `scripts/`, `tests/`, `tools/`, `run.py`, `run_all.py`): Apache License 2.0 (`LICENSE`).
Data and results (`data/`, `results/`): Creative Commons Attribution 4.0 International (`LICENSE-DATA`). See `NOTICE`.
