# Claim-to-evidence map

Every computational statement of the paper *Sums of two dependent uniform random variables: a characterization by
joint tail bounds*, keyed by the result it belongs to, with the script or test that checks it, the data it reads,
the tier at which the stated range is reached, and its status. Results are cited by their printed numbers (prefix S: the supplementary material);
paths are relative to this archive; all scripts are in `scripts/`, all tests in `tests/`.

Tiers: **quick** = `python3 run.py` (the test suite and every step at its default range); **medium** =
`python3 run.py --tier medium`; **full** = `python3 run.py --tier full` (every script with `--full`). Measured times
are in `README.md`.

Status:
- **full stream verified**: the quick tier checks the whole stated range, or every record of the deposited data;
- **small range**: the quick tier checks a smaller range and the full tier the stated range (see the Tier column);
- **index only**: only an index or summary of a deposit is checked;
- **not checked**: the statement is not reproduced by the code.

Finite computations do not prove statements for all parameters; where the paper proves a statement, the computation
checks the proof or its finite input. All checks are exact (integer or rational arithmetic, with certificates) unless
marked FLOAT; FLOAT items are cross-checks or numerical evidence.

## The fibered reduction and two uniforms at every spacing

| Claim | Script / test | Data | Tier | Status |
|---|---|---|---|---|
| Theorem 2.2 and Remark 2.1: two-point laws {a, a+b} with mean one are attainable iff 1/b is an integer; the support condition of (i) (fibered LP; FLOAT) | `tests/test_fibered_reduction.py` (test_biatomic, test_support_condition) | | quick | full stream verified |
| Theorem 2.2: a non-constant support pattern at b = 2/5 is feasible; the three-point minimum for non-integer 1/(2b) (FLOAT) | `tests/test_fibered_reduction.py` (test_nonconstant_example, test_three_point_minimum) | | quick | full stream verified |
| Example S3.1: b = 2/5 gives the law (3/5, 1/5, 1/5) on {3/5, 7/5, 9/5} by the stated fiber permutations (exact construction) | `tests/test_two_uniform_reduction.py::test_exact_construction_spacing_two_fifths` | | quick | full stream verified |

## The grid and potential duality

| Claim | Script / test | Data | Tier | Status |
|---|---|---|---|---|
| Corollary 3.1: the tail form of a step potential equals its grid inequality; the single step decides symmetric three-point laws | `tests/test_potential_duality_readings.py` | | quick | full stream verified |
| Example 3.1: F = (delta_0.6 + delta_1.4)/2 has mean one, satisfies the convex (gap x^2/4 on [0, 0.6], at least 0.05 on [0.6, 1.4]) and folded conditions; the step at 0.3 gives 3/2 > 1.4; the single step 1{i >= 3} is a facet of hat-P(10) and is violated | `check_two_point_example.py`; `check_mww_steps.py` (example) | | quick | full stream verified |

## Windows and the three- and four-atom laws

| Claim | Script / test | Data | Tier | Status |
|---|---|---|---|---|
| Lemma S3.2: the vertices of the width-four window polytope lie among the five candidate points (m <= 60) | `tests/test_windows.py::test_vertex_rule`; `check_window_corollary.py` (vertex_rule) | | quick | full stream verified |
| Corollary 6.1: three-point intervals against the generic fibered LP, including a + b = 1 (FLOAT) | `tests/test_windows.py::test_three_point_interval`, `tests/test_fibered_reduction.py::test_three_point_minimum` | | quick | full stream verified |
| Corollary S3.1 against a generic fibered LP: 215 windows, 77 without a mean-one law (prediction agrees), 1,104 support-function comparisons, largest discrepancy below 10^-14 (FLOAT; Section S5) | `check_window_corollary.py`; `tests/test_windows.py::test_four_atom_corollary`, `::test_window_starting_at_one` | | quick | full stream verified |
| After Corollary S3.1: the centered five-point window gives polytopes with 1, 2, 5, 5, 9, 6, 9, 9, 9, 10 vertices for m = 1..10 | `count_window_vertices.py` (m = 1..10, exact integer hulls); `tests/test_windows.py::test_five_point_window_counts` (m = 1..8) | | quick | full stream verified |
| Proposition 6.1(a): two-point laws with 1/2 < b < 1 violate the single step at min(a/2, a+b-1) (3,626 laws on the grid 1/100) | `check_mww_steps.py` (two_point) | | quick | full stream verified |
| Proposition 6.1(b): symmetric three-point laws, 1/4 <= b < 1: the single step decides membership (7,575 laws) | `check_mww_steps.py` (three_point); `tests/test_potential_duality_readings.py::test_symmetric_three_point_single_step` | | quick | full stream verified |
| After Proposition 6.1: the law (1/160, 1/40, 31/32) on {123/280, 29/40, 283/280} has mean one and satisfies every single-step inequality | `check_mww_steps.py` (asymmetric) | | quick | full stream verified |

## Small atoms

| Claim | Script / test | Data | Tier | Status |
|---|---|---|---|---|
| Proposition 3.3(i): on the affine hull the k - 2 hinge inequalities give T <=cx 2U_k at every integer hinge (checked on every vertex of Q(k)) | `check_small_atoms.py` (statement) | | quick (k <= 5), full (k <= 7) | small range |
| Corollary 3.2, statement, both directions: permutation X-rays satisfy T <=cx 2U_k; every vertex of Q(k) (2, 7, 44, 258, 1,493, 8,300 vertices for k = 2..7) has an exact doubly stochastic matrix with its X-ray. Vertices by exhaustive exact enumeration (k <= 5 quick, k <= 6 full) and by pycddlib with exact recovery (FLOAT enumeration, exact vertices; k = 7 relies on its completeness) | `check_small_atoms.py` (statement); `tests/test_small_atoms.py` | | quick (k <= 5), full (k <= 7, about 3 minutes; needs pycddlib for k = 7) | small range |
| Corollary 3.2, negative control: without x <= 1 there are vertices outside hat-P(k) (1 of 5, 4 of 13, 15 of 33 for k = 3, 4, 5), e.g. (0, 3/2, 0, 3/2, 0) | `check_small_atoms.py` (negative_control); `tests/test_small_atoms.py::test_negative_control` | | quick | full stream verified |
| Lemma 3.2: necessity on random doubly stochastic matrices; sufficiency: each vertex outside hat-P(k) in the negative control gets an exact separating potential, the others none | `check_small_atoms.py` (negative_control, chain) | | quick | full stream verified |
| Lemma 3.3 with its pointwise bounds on each piece of the convex minorant: random u (k = 1..40) and every signed staircase | `check_small_atoms.py` (minorant_random, minorant_staircases); `tests/test_small_atoms.py` | | quick (4,000 u; staircases k <= 10), full (16,000 u; all 32,766 staircases k <= 14) | small range |
| Proof of Corollary 3.2: the three-line chain for random points of Q(k) (12 vertices and 6 mixtures per k) and 200 random u each, with an exact doubly stochastic matrix for every point | `check_small_atoms.py` (chain) | | quick (k <= 8), full (k <= 12) | small range |
| Remark S4.1: exact witnesses with atoms <= 101/100 for k = 9, 10, 11 (bounds, affine hull, all hinge inequalities, convex and folded order directly, a facet I(u) violated by 6/25, 9/25, 3/10); found by a FLOAT LP and verified exactly | `check_lattice_tightness.py` | `data/facets/facets_k09.jsonl`, `facets_k10.jsonl`, `facets_k11.jsonl` (search) | quick | full stream verified |

## Tail sums and step counts

| Claim | Script / test | Data | Tier | Status |
|---|---|---|---|---|
| Lemma 4.1: (u box u)(s) = #{l : xi_l <= s} with repeated thresholds and thresholds at 1 (against direct minimization) | `check_staircase_lp.py` (lem_tailsgen) | | quick (3,000 trials), full (12,000) | small range |
| Theorem 4.1: the implication (ii) => (iii) is trivial and the others are proved analytically; the computations check the ingredients (Lemma 4.1, Example 3.1 for strictness of convex order) | see the rows for those results | | | not checked |
| Proposition 4.1 (the tail-sum characterization for a general marginal F) | proved analytically | | | not checked |
| Remark 4.2: the tail-sum inequality with multiplicities equals I(u)/k on the grid | `check_staircase_lp.py` (rem_multiplicity) | | quick | full stream verified |
| Example 4.2(a): xi values, 53/12 > 52/12; u = (0,1,2,2,4,4) is a ray of height two giving a facet of hat-P(6); 25 tight permutations spanning dimension 8 | `check_repeated_thresholds.py` (part_a) | | quick | full stream verified |
| Example 4.2(b): F_W satisfies I(u_J) for all 2^(K-1) step sets, K = 6, 12, 18, 24 (8,388,608 at K = 24) | `check_repeated_thresholds.py` (part_b; C helper `tools/c/grid_steps.c` and numpy, which must agree; numpy alone without a compiler) | | quick | full stream verified |
| Example 4.2(c): J = {6, 12, 23, 24} on K = 36 excludes F_W by 1/36 | `check_repeated_thresholds.py` (part_c) | | quick | full stream verified |
| Remark 4.1: the thresholds (1/5, 2/5, 1/2, 4/5) and (1/5, 2/5, 3/5, 4/5) give the same xi = (2/5, ..., 9/5), with right sides 21/5 and 4 | `check_tail_bounds.py` (remark_notsharp) | | quick | full stream verified |
| Lemma 4.2 (reflection), random exact tests | `check_staircase_lp.py` (lem_reflectsteps) | | quick | full stream verified |
| Lemma 4.3 (counting), random exact tests and all pair sums | `check_tail_bounds.py` (counting_lemma); `check_staircase_lp.py` (lem_counting) | | quick (3,000 trials), full (222,702 pairs) | small range |
| Theorem 4.2: the explicit thresholds and the violation bound; fewer steps never exclude (random and adversarial staircases); the least number of steps by the exact staircase program | `check_tail_bounds.py` (two_point_*); `check_staircase_lp.py` (thm_twopoint_lp) | | quick (16,000 lower-bound staircases), full (40,000) | small range |
| Theorem 4.3: the explicit r-step staircase (values r, r, r+1; inequality t >= dist(L, 2Z)/L), the exact smallest middle mass, the least number of steps by the exact program | `check_tail_bounds.py` (three_point_family, t_r); `check_staircase_lp.py` (thm_threepoint_lp) | | quick | full stream verified |
| Corollary 4.1: the laws F_{w,beta} with r+1 < 1/beta < r+2 satisfy all potentials with at most r steps | `check_tail_bounds.py` (two_point_lower, two_point_exact_lp); `check_staircase_lp.py` | | quick | full stream verified |
| Proposition S4.1: the thresholds (3/280, 29/80, 363/560), xi, 317/80 - 277/70 = 3/560; the exact maximum is 0 for r = 1, 2 and 3/560 for r = 3, for F_* and its reflection (two independent programs; Section S5) | `check_staircase_lp.py` (prop_asym*); `check_mww_steps.py` (asymmetric) | | quick | full stream verified |

## Laws with a density

| Claim | Script / test | Data | Tier | Status |
|---|---|---|---|---|
| Theorem 5.1: the construction of the proof (margin, rounding, lattice convex order) on four laws and k = 24, 30, 40; the example with density 1 on [1/3, 5/6] u [7/6, 5/3] satisfies S <=cx 2U | `check_density_construction.py` | | quick | full stream verified |
| Proposition 5.1: mean, convex order with D(1) = (2c-1)/(24c), folded order with H(1/3) = H(2/3) = (4c+1)/(144c), violation (c-2)/(6c); violation 0 at c = 2 (exact on twelve values of c; symbolic in c with sympy, skipped without it); a FLOAT experiment with at most three thresholds | `check_density_sharpness.py` (prop_densitysharp) | | quick, full (larger FLOAT experiment) | full stream verified |
| Lemma 5.1 and Proposition 5.2: with t = (8 eta - 2)/3 and phi(s) = (s - t)^+/(1 - eta), (1/2) int phi + 2 int (N - phi)^+ = 2(1 - eta) - (t^-)^2/(4(1 - eta)), and every intermediate value of the proof (eta = j/240) | `check_single_step_density.py` (prop_singlestep); `tests/test_density_certificates.py` | | quick | full stream verified |
| Theorem 5.2: thresholds with nonincreasing gaps: the hinge certificate and its closed-form margin; r = 1 is Proposition 5.2 (494 gap patterns, 7,410 cases with r <= 8) | `check_nonincreasing_gaps_density.py`; `tests/test_density_two.py` | | quick (r <= 5), full (r <= 8) | small range |
| Theorem 5.3 (certificate): equally spaced thresholds, the certificate works for c <= r + 1 (960 cases with r <= 24) | `check_equally_spaced_density.py` (certificate) | | quick (r <= 8), full (r <= 24) | small range |
| Theorem 5.3 (sharpness): for c > r + 1 the explicit laws have mass and mean one, S <=cx 2U, and violate the staircase by (c - r - 1)/(2c(2r + 1)) (72 laws with r <= 24) | `check_equally_spaced_density.py` (sharpness) | | quick (r <= 8), full (r <= 24) | small range |
| Sentence after Theorem 5.3: at c = 2 and r >= 2 the certificate keeps margin >= m - 10 r eps >= 0 for thresholds within eps <= m/(10r) | `check_equally_spaced_density.py` (perturbation, 2 <= r <= 6 at every tier) | | quick | small range |
| Theorem 5.4: two-step staircases at density 2, slack at least min(b - a, 1 - b)/21, constant attained at rho = 3/4; 7,021 pairs on the grid 1/120, smallest margin 1/1080; case identities symbolically (sympy, skipped without it) | `check_two_step_density.py`; `tests/test_density_two.py` | | quick (grid 1/40, 741 pairs), full (grid 1/120) | small range |
| Proposition 5.3: an explicit piecewise affine coupling with sum density 2 on [1/2, 3/4] u [5/4, 3/2]; S <=cx 2U, equality in the single step at eta = 1/4 | `check_density_coupling.py`; `tests/test_density_two.py::test_coupling` | | quick | full stream verified |
| After Proposition 5.4: for thresholds 1/8, 1/4, 1/2, 5/8, 3/4 the sum of the Proposition 5.2 certificates exceeds the budget 11/2 by 45889/39312 (the staircase itself is certified by the catalogue) | `check_hinge_sum_failure.py`; `tests/test_density_two.py::test_hinge_sum_excess` | `data/density/density_two_certificates_K08.jsonl.gz` | quick | full stream verified |
| After Problem 5.1: for every staircase with at least two distinct thresholds on the grids 1/8, 1/12, 1/16 (120, 2,036, 32,752), a convex piecewise linear phi satisfies the hypothesis of Lemma 5.1 with c = 2 (found by a FLOAT LP, verified exactly) | `check_single_step_density.py` (certificates_K8, K12, K16) | | quick (K = 8), medium (K = 12), full (K = 16) | small range |
| After Problem 5.1: deposited exact certificates (c = 2) for the grids 1/8, 1/12, 1/16 and every two-step staircase on 1/48 (1,081); minimum margins 7/512, 7/768, 1/144, 7/3072; each record verified exactly | `check_density_catalogues.py`; `tests/test_density_two.py` | `data/density/density_two_certificates_K08.jsonl.gz`, `..._K12.jsonl.gz`, `..._K16.jsonl.gz`, `..._two_step_K48.jsonl.gz` | quick (grids 8, 12 and two-step 48), full (also grid 16) | small range |
| After Problem 5.1: floating-point searches over staircases with up to eight steps on finer grids find largest violation 0 at c = 2, attained only by single steps (FLOAT, numerical evidence) | `check_single_step_density.py` (float_search) | | quick (r <= 4, grid 1/20), full (r <= 8, grid 1/40) | small range |
| Proposition 5.4: duality and the reduction to bang-bang laws | proved analytically | | | not checked |
| Example 4.1: density 4 on [7/20, 33/80] u [61/80, 33/40] u [47/40, 99/80] u [127/80, 33/20]: mean one; S <=cx 2U (exact piecewise minimum); every single-step inequality holds, largest left minus right over 0 < eta <= 33/40 is -3/40 (at eta = 47/80); the staircase (7/40, 47/80) has xi = (7/20, 61/80, 47/40, 127/80), 5/2 > 99/40, violation 1/40; the the folded condition holds (exact, two implementations agree) | `check_worked_example.py` | | quick | full stream verified |

## Certified computations

| Claim | Script / test | Data | Tier | Status |
|---|---|---|---|---|
| Section S5: the certified facet test used for single steps on hat-P(10), the facet of Example 4.2(a) and the facets of Remark S4.1: every outcome carries a certificate (no tight permutation, two anti-diagonals without allowed edge, rank modulo 2^31 - 1, exact recomputation) | `src/dsxray/facet_test.py`; `check_two_point_example.py`, `check_repeated_thresholds.py`, `check_lattice_tightness.py` | | quick | full stream verified |
| Section S5: floating point enters only in labeled cross-checks (fibered LPs, the window corollary, the LPs that propose certificates and lattice witnesses, the cddlib vertex enumeration of Corollary 3.2); every certificate, witness and vertex is verified exactly | the scripts above | | quick | full stream verified |
| Cost figures (seconds, hours, machines) | not reproduced: machine-dependent | | | not checked |

## Corrections to Mao, Wang and Wang (Appendix A)

| Claim | Script / test | Data | Tier | Status |
|---|---|---|---|---|
| Item 1: the law (1/18, 1/18, 8/9) on {7/18, 13/18, 19/18} belongs to the set (exact construction) but violates the printed bound; the lattice reduction reproduces their three-point theorem with the bound a + b - 1 | `tests/test_two_uniform_reduction.py` (test_tri_atomic, test_printed_statement_fails_on_an_attainable_law, test_exact_construction_of_the_law) | | quick | full stream verified |
| Item 2: the printed parameter is 2.2 at b = 1/4, a = 0.8 (outside [0, 1]); the corrected (a + b - 1)/b lies in (0, 1) | `tests/test_two_uniform_reduction.py::test_printed_parameter_leaves_the_simplex` | | quick | full stream verified |
| Item 3: the printed vertex of case (iii) sums to 1 + bc; the corrected vertex is (1/18, 1/18, 8/9) at a = 13/18, b = 1/3 | `tests/test_two_uniform_reduction.py::test_vertex_normalization` | | quick | full stream verified |
| Supporting cross-checks: the coupling LPs reproduce published bounds (Lu, Ding and Dasgupta; lattice Makarov bounds; Mao and Wang) (FLOAT) | `tests/test_couplings_vs_published.py` | | quick | full stream verified |

## Additional checks not cited in the paper

`check_mww_steps.py --search` runs an arrangement search (bisection on the middle mass with two steps; not a
certificate). `check_density_sharpness.py` runs a floating-point search over piecewise-constant densities and
staircases with at most three thresholds.
