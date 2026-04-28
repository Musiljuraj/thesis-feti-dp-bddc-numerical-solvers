# bc_ddm - FEM, FETI-DP and BDDC experiments for a bachelor thesis

This repository contains the computational implementation accompanying the bachelor thesis

**Comparison of the solvers of systems of linear equations: FETI-DP vs BDDC**

The code is written mainly in **GNU Octave / MATLAB-style `.m` files** and implements a controlled sequential framework for studying two dual-primal domain decomposition methods:

- **FETI-DP** - Finite Element Tearing and Interconnecting - Dual-Primal
- **BDDC** - Balancing Domain Decomposition by Constraints

The repository is not meant to be a general-purpose PDE solver package. Its purpose is to make the numerical experiments in the thesis reproducible and to expose the algebraic structure of the methods: finite element assembly, substructuring, Schur complement reduction, primal constraints, FETI-DP and BDDC operator application, PCG convergence, and spectral diagnostics.

---

## Thesis and numerical scope

The thesis studies linear systems arising from finite element discretizations of elliptic boundary value problems. The custom Octave implementation uses the model Poisson problem

```text
-Delta u = f  in Omega = [0,1] x [0,1]
u = 0         on the left and right boundary
du/dn = 0     on the bottom and top boundary
```

The problem is discretized by continuous piecewise linear finite elements on a structured triangular mesh. The resulting stiffness system is then reformulated by non-overlapping domain decomposition. Interior subdomain degrees of freedom are eliminated locally, producing interface Schur complement problems. FETI-DP and BDDC are then implemented on top of the same substructuring data.

The repository covers three computational layers:

1. **Sequential Octave implementation**
   - P1 finite element mesh generation, element matrices and global assembly.
   - Homogeneous Dirichlet elimination.
   - Structured rectangular non-overlapping subdomain partitioning.
   - Local Schur complements and product/assembled interface spaces.
   - Corner-type primal degrees of freedom.
   - Matrix-free FETI-DP and BDDC operators and preconditioners.
   - PCG convergence diagnostics.

2. **Sequential spectral experiments**
   - Explicit assembly of small matrix-free operators by probing basis vectors.
   - Spectra of non-preconditioned reduced operators.
   - Spectra of preconditioned PCG operators.
   - Residual history plots.
   - LaTeX tables used in the thesis.

3. **PETSc benchmark post-processing**
   - Stored logs and processed outputs from PETSc benchmark runs.
   - Python scripts for extracting PETSc iteration counts, spectral indicators, timing summaries and residual plots.
   - These PETSc results are complementary to the Octave implementation; they are not a direct timing comparison with the Octave code.

---

## What is implemented

### FEM layer

Located mainly in `src/fem/`.

The finite element part implements:

- structured unit-square mesh generation,
- P1 triangular elements,
- local stiffness matrix computation,
- local load vector computation,
- global sparse stiffness matrix assembly,
- global load vector assembly,
- homogeneous Dirichlet elimination on the left and right boundary.

Main files:

```text
src/fem/mesh/mesh_unit_square_P1.m
src/fem/elements/triP1_stiffness.m
src/fem/elements/triP1_load.m
src/fem/assembly/assemble_stiffness_P1.m
src/fem/assembly/assemble_load_P1.m
src/fem/bc/apply_dirichlet_elimination.m
```

### Domain decomposition layer

Located mainly in `src/ddm/`.

This layer builds the non-overlapping substructuring framework used by both solvers:

- structured rectangular subdomain partitioning,
- local subdomain matrices,
- interior/interface splitting,
- local Schur complements,
- product interface space,
- assembled interface space,
- jump operator,
- assembly operator,
- primal corner selection,
- maps between primal and non-primal interface unknowns.

Main files:

```text
src/ddm/partition/build_subdomains_structured.m
src/ddm/partition/identify_interface_dofs.m
src/ddm/local/assemble_subdomain_matrices_P1.m
src/ddm/local/extract_subdomain_blocks.m
src/ddm/local/setup_local_schur.m
src/ddm/local/apply_local_schur.m
src/ddm/interface/build_product_interface.m
src/ddm/interface/build_jump_operator_B.m
src/ddm/interface/build_assembly_operator_R.m
src/ddm/interface/apply_blockdiag_S.m
src/ddm/coarse/select_primal_dofs.m
src/ddm/coarse/build_primal_maps.m
```

### Common infrastructure

Located mainly in `src/common/`.

This contains shared routines used by both FETI-DP and BDDC:

- construction of the full problem-data structure,
- multiplicity scaling,
- deterministic random seed helper,
- matrix-free PCG wrapper,
- spectral assembly and eigenvalue diagnostics,
- plotting routines,
- LaTeX table exporters.

Main files:

```text
src/common/diagnostics/build_problem_data.m
src/common/krylov/pcg_wrap.m
src/common/scaling/multiplicity_scaling.m
src/common/spectra/assemble_from_apply.m
src/common/spectra/full_spectrum_operator.m
src/common/spectra/full_spectrum_precond.m
src/common/spectra/export_ch6_table_tex.m
src/common/spectra/export_ch6_raw_table_tex.m
src/common/spectra/plot_spectrum_sorted_overlay.m
src/common/spectra/plot_spectrum_histogram_overlay.m
src/common/spectra/plot_residual_history_overlay.m
src/common/utils/rng_deterministic.m
```

### FETI-DP solver

Located in `src/feti_dp/`.

The FETI-DP implementation uses a multiplier-space formulation. Continuity of non-primal interface unknowns is enforced by Lagrange multipliers, while primal corner values are treated through a coarse problem.

Main files:

```text
src/feti_dp/setup/setup_fetidp.m
src/feti_dp/operators/applyA_lambda.m
src/feti_dp/operators/applyM_lambda.m
src/feti_dp/operators/solve_tildeS.m
src/feti_dp/solve/solve_fetidp.m
src/feti_dp/solve/reconstruct_fetidp_solution.m
```

Conceptually:

- `setup_fetidp` prepares the multiplier-space operator, primal/delta splitting, local factorizations, jump constraints and coarse data.
- `applyA_lambda` applies the FETI-DP operator in matrix-free form.
- `applyM_lambda` applies the FETI-DP preconditioner.
- `solve_fetidp` runs PCG in multiplier space.
- `reconstruct_fetidp_solution` reconstructs the global free-DOF solution after the multiplier solve.

### BDDC solver

Located in `src/bddc/`.

The BDDC implementation uses an assembled-interface formulation. Interface corrections are combined through scaling/averaging and a primal coarse correction.

Main files:

```text
src/bddc/setup/setup_bddc.m
src/bddc/operators/applyA_hat.m
src/bddc/operators/applyM_bddc.m
src/bddc/solve/solve_bddc.m
src/bddc/solve/reconstruct_bddc_solution.m
```

Conceptually:

- `setup_bddc` prepares the assembled interface operator, local Schur data, scaling weights and coarse basis.
- `applyA_hat` applies the assembled BDDC operator.
- `applyM_bddc` applies the BDDC preconditioner.
- `solve_bddc` runs PCG in assembled interface space.
- `reconstruct_bddc_solution` reconstructs the global free-DOF solution after the interface solve.

---

## Repository layout

```text
bc_ddm/
├── README.md
├── main/
│   ├── setup_paths.m
│   ├── ch6_define_cases.m
│   ├── ch6_run_case.m
│   ├── run_ch6_spectral_experiments.m
│   ├── figures/
│   │   ├── run_all_figures.m
│   │   └── fig_*.m
│   └── archive/
│       └── older demo and integration runners
├── src/
│   ├── fem/
│   ├── ddm/
│   ├── feti_dp/
│   ├── bddc/
│   └── common/
├── tests/
│   ├── ch2_fem/
│   ├── ch3_ddm/
│   ├── ch4_common/
│   ├── ch4_fetidp/
│   ├── ch4_bddc/
│   └── ch6_spectra/
├── output/
│   ├── figures/
│   ├── mats/
│   └── tables/
└── petsc-benchmarks-rtol1e10/
    ├── logs/
    ├── figures/
    ├── tables/
    ├── final/
    └── scripts/
```

Notes:

- `src/` is the main library code.
- `main/` contains current entry points for the Chapter 6 sequential experiments.
- `main/archive/` contains older development scripts. They are useful as references, but the current primary workflow is in `main/run_ch6_spectral_experiments.m` and `main/ch6_run_case.m`.
- `tests/` contains topic-based tests. Some folder names preserve earlier development numbering, so they do not perfectly match the final thesis chapter numbering. For example, FEM tests are under `tests/ch2_fem/`, while the final thesis presents FEM discretization in Chapter 3.
- `output/` contains generated Octave experiment artifacts. It can be deleted and regenerated.
- `petsc-benchmarks-rtol1e10/` contains stored PETSc benchmark logs and post-processed results.

---

## Requirements

### Required for the Octave implementation

- GNU Octave, recommended.
- MATLAB may work for most code paths, but the project was developed in an Octave-style workflow.
- Standard sparse linear algebra support.
- No additional Octave packages are expected for the core solver workflow.

### Optional for PETSc post-processing

Only needed if you want to rerun the Python scripts in `petsc-benchmarks-rtol1e10/scripts/`:

- Python 3
- `matplotlib`

The PETSc benchmark itself was produced externally with PETSc and MPI. This repository stores logs and processed outputs; it does not contain a full PETSc build or a complete PETSc execution wrapper.

---

## Quick start

Run all commands from the repository root.

### 1. Start Octave and load paths

```matlab
addpath('main');
setup_paths();
```

or from the shell:

```bash
octave --eval "addpath('main'); setup_paths();"
```

### 2. Run a small sanity check

```matlab
addpath('main');
setup_paths();

data = build_problem_data(8, 2, 2, @(x,y) 1.0);

disp(data.ddm.Nsub);
disp(size(data.Kff));
```

Expected meaning:

- `data.ddm.Nsub` is the number of subdomains.
- `data.Kff` is the reduced finite element stiffness matrix after Dirichlet elimination.

---

## Main workflow: Chapter 6 sequential spectral experiments

The main current runner is:

```text
main/run_ch6_spectral_experiments.m
```

Run it with:

```matlab
addpath('main');
run_ch6_spectral_experiments();
```

or from the shell:

```bash
octave --eval "addpath('main'); run_ch6_spectral_experiments();"
```

This runs the case set defined in:

```text
main/ch6_define_cases.m
```

The current case set is:

```text
n = 16, subdomains = 2x2
n = 24, subdomains = 2x2
n = 32, subdomains = 2x2
n = 32, subdomains = 4x4
n = 48, subdomains = 2x2
n = 48, subdomains = 4x4
n = 64, subdomains = 8x8
```

The batch runner performs, for each case:

1. build the common Poisson/FEM/DDM data,
2. set up FETI-DP,
3. solve the FETI-DP multiplier problem by PCG,
4. compute FETI-DP spectra if the dimension is below the configured limit,
5. set up BDDC,
6. solve the BDDC assembled interface problem by PCG,
7. compute BDDC spectra if the dimension is below the configured limit,
8. export `.mat` results,
9. export spectral plots,
10. export residual-history plots,
11. export LaTeX summary tables.

Generated files are written to:

```text
output/mats/ch6/
output/figures/ch6/
output/tables/ch6/
```

Important outputs include:

```text
output/tables/ch6/table_ch6_summary.tex
output/tables/ch6/table_ch6_unpreconditioned_summary.tex
output/figures/ch6/fig_spec_sorted_*.pdf
output/figures/ch6/fig_spec_hist_*.pdf
output/figures/ch6/fig_spec_raw_sorted_*.pdf
output/figures/ch6/fig_spec_raw_hist_*.pdf
output/figures/ch6/fig_residual_*.pdf
```

The full spectral computation can be expensive because the matrix-free operators are assembled explicitly by applying them to basis vectors. This is intentional: the thesis uses this only for small and moderate sequential diagnostic cases.

---

## Running a single FETI-DP vs BDDC case

Use `main/ch6_run_case.m` for one controlled case.

```matlab
addpath('main');
setup_paths();

cfg = struct();
cfg.n = 32;
cfg.nSubX = 4;
cfg.nSubY = 4;
cfg.seed = 1;
cfg.tol = 1e-10;
cfg.maxit = 300;
cfg.do_spectra = true;
cfg.nmax = 1200;
cfg.verbose = true;

out = ch6_run_case(cfg);

disp(out.fetidp.stats);
disp(out.bddc.stats);
```

The returned structure contains:

```text
out.cfg
out.case_id
out.base_meta
out.fetidp.n
out.fetidp.stats
out.fetidp.spec
out.fetidp.spec_A
out.fetidp.spec_skipped
out.fetidp.spec_skip_reason
out.bddc.n
out.bddc.stats
out.bddc.spec
out.bddc.spec_A
out.bddc.spec_skipped
out.bddc.spec_skip_reason
```

Here:

- `spec` is the spectrum of the preconditioned PCG-relevant operator.
- `spec_A` is the spectrum of the corresponding non-preconditioned reduced operator.
- `stats.resvec` contains the PCG residual history.
- `stats.iter` and `stats.relres` summarize PCG convergence.

---

## Running the solvers directly

The direct solver-level workflow is useful when debugging or extending the implementation.

```matlab
addpath('main');
setup_paths();

% Build common FEM + DDM data.
data = build_problem_data(32, 4, 4, @(x,y) 1.0);

% FETI-DP
data_f = setup_fetidp(data);
[lambda, stats_f] = solve_fetidp(data_f, 1e-10, 300);
[w_c, w_d, u_feti, diag_f] = reconstruct_fetidp_solution(lambda, data_f);

% BDDC
data_b = setup_bddc(data);
sol_b = solve_bddc(data_b, struct('tol', 1e-10, 'maxit', 300));
rec_b = reconstruct_bddc_solution(sol_b.u_hat, data_b);
u_bddc = rec_b.u_free;

% Compare with the monolithic FEM reference solve.
u_ref = data.Kff \ data.Ff;

fprintf('FETI-DP relative error: %.3e\n', norm(u_feti - u_ref) / norm(u_ref));
fprintf('BDDC    relative error: %.3e\n', norm(u_bddc - u_ref) / norm(u_ref));
```

---

## Configuration parameters

The most important configuration fields for `ch6_run_case` are:

```text
cfg.n              mesh parameter; unit square has (n+1)^2 nodes
cfg.nSubX          number of subdomains in x-direction
cfg.nSubY          number of subdomains in y-direction
cfg.seed           deterministic RNG seed
cfg.f_handle       right-hand side f(x,y), default @(x,y) 1.0
cfg.tol            PCG tolerance
cfg.maxit          maximum PCG iterations
cfg.do_spectra     true/false; whether to compute full spectra
cfg.nmax           maximum operator dimension allowed for full spectra
cfg.verbose        true/false; print progress and summaries
cfg.store_base     true/false; store common data in the output struct
cfg.store_data     true/false; store full solver data in the output struct
```

Constraints:

- `n` must be divisible by `nSubX`.
- `n` must be divisible by `nSubY`.
- `nSubX * nSubY` must be at least 2 for the domain decomposition workflow.
- The implemented decomposition is structured and aligned with the finite element grid.

Example with a different right-hand side:

```matlab
cfg = struct();
cfg.n = 32;
cfg.nSubX = 2;
cfg.nSubY = 2;
cfg.f_handle = @(x,y) 1.0 + x - 0.5*y;
out = ch6_run_case(cfg);
```

---

## Spectral analysis controls

The spectral analysis is intentionally explicit and diagnostic.

The code uses:

```text
assemble_from_apply.m
full_spectrum_operator.m
full_spectrum_precond.m
```

The reduced operators are applied in matrix-free form, but for spectral analysis they are explicitly assembled column by column by applying the operator to standard basis vectors.

Full spectra are computed only if the relevant PCG-space dimension is not larger than `cfg.nmax`:

```text
FETI-DP dimension: nLambda
BDDC dimension:    nHat
```

If the dimension is too large, the result is still valid as a PCG run, but the full spectrum is skipped. The output then contains:

```text
out.fetidp.spec_skipped
out.fetidp.spec_skip_reason
out.bddc.spec_skipped
out.bddc.spec_skip_reason
```

To disable spectra completely:

```matlab
cfg.do_spectra = false;
```

To allow larger explicit spectra:

```matlab
cfg.nmax = 2000;
```

Be careful with large values of `cfg.nmax`; explicit dense spectral analysis becomes expensive quickly.

---

## PETSc benchmark data

The folder

```text
petsc-benchmarks-rtol1e10/
```

contains benchmark logs and post-processed outputs for the PETSc part of the thesis. These results complement the Octave implementation. They compare PETSc's official BDDC and FETI-DP related implementations on a larger benchmark set.

Important files and folders:

```text
petsc-benchmarks-rtol1e10/logs/
petsc-benchmarks-rtol1e10/tables/raw_results.csv
petsc-benchmarks-rtol1e10/tables/summary_results.csv
petsc-benchmarks-rtol1e10/final/tables/petsc_main_results.csv
petsc-benchmarks-rtol1e10/final/tables/table_petsc_main_results.tex
petsc-benchmarks-rtol1e10/final/figures/fig_petsc_residual_CM.pdf
petsc-benchmarks-rtol1e10/final/figures/fig_petsc_residual_CL.pdf
petsc-benchmarks-rtol1e10/scripts/process_petsc_outputs.py
petsc-benchmarks-rtol1e10/scripts/export_petsc_thesis_outputs.py
```

The PETSc benchmark should be interpreted carefully:

- It uses PETSc's official benchmark example, not the custom Octave FEM code.
- It belongs to the same general class of structured scalar elliptic problems.
- It is useful for comparing PETSc iteration counts, spectral indicators, residual histories and representative timings.
- It is not a direct timing comparison between Octave and PETSc.
- It is not a complete strong-scaling or weak-scaling study.

The Python post-processing scripts currently assume the path

```text
~/bc/petsc-benchmarks-rtol1e10
```

through a hard-coded `OUT` variable. If the repository is located elsewhere, edit the `OUT` path in the scripts before rerunning them.

---

## Tests

The repository contains many unit and integration-style tests under `tests/`.

Representative tests include:

```text
tests/ch2_fem/test_mesh_unit_square_P1.m
tests/ch2_fem/test_triP1_stiffness.m
tests/ch2_fem/test_assemble_stiffness_P1.m
tests/ch3_ddm/test_build_subdomains_structured.m
tests/ch3_ddm/test_setup_local_schur.m
tests/ch4_common/test_pcg_wrap.m
tests/ch4_fetidp/test_solve_fetidp.m
tests/ch4_bddc/test_solve_bddc.m
tests/ch6_spectra/test_ch6_run_case.m
tests/ch6_spectra/test_full_spectrum_precond.m
```

Run individual tests after setting the paths:

```matlab
addpath('main');
setup_paths();

test_mesh_unit_square_P1();
test_triP1_stiffness();
test_build_problem_data();
test_solve_fetidp();
test_solve_bddc();
test_ch6_run_case();
```

A number of older runners exist under `main/archive/`. Treat these as development history and reference scripts unless you have checked that the particular archived runner still matches the current folder layout.

---

## Generated thesis figures

Figure-generation scripts are located in:

```text
main/figures/
```

The helper

```text
main/figures/run_all_figures.m
```

generates selected explanatory figures and writes them to `output/figures/...`.

Run from the repository root with:

```matlab
run('main/figures/run_all_figures.m');
```

The Chapter 6 spectral and residual figures are generated by the main Chapter 6 batch runner, not by `run_all_figures.m`.

---

## Development guide

Recommended files to read first:

```text
src/common/diagnostics/build_problem_data.m
main/ch6_run_case.m
src/feti_dp/setup/setup_fetidp.m
src/feti_dp/operators/applyA_lambda.m
src/feti_dp/operators/applyM_lambda.m
src/bddc/setup/setup_bddc.m
src/bddc/operators/applyA_hat.m
src/bddc/operators/applyM_bddc.m
```

Typical dependency chain:

```text
mesh_unit_square_P1
  -> assemble_stiffness_P1 / assemble_load_P1
  -> apply_dirichlet_elimination
  -> build_subdomains_structured
  -> identify_interface_dofs
  -> assemble_subdomain_matrices_P1
  -> extract_subdomain_blocks
  -> setup_local_schur
  -> build_product_interface
  -> select_primal_dofs
  -> build_primal_maps
  -> setup_fetidp / setup_bddc
  -> solve_fetidp / solve_bddc
  -> reconstruct_fetidp_solution / reconstruct_bddc_solution
```

When modifying code:

- FEM changes should be checked against the tests in `tests/ch2_fem/`.
- DDM interface and Schur complement changes should be checked against `tests/ch3_ddm/`.
- Common solver infrastructure changes should be checked against `tests/ch4_common/`.
- FETI-DP changes should be checked against `tests/ch4_fetidp/`.
- BDDC changes should be checked against `tests/ch4_bddc/`.
- Spectral-analysis changes should be checked against `tests/ch6_spectra/`.

---

## Known limitations

The implementation intentionally keeps the mathematical setting narrow and transparent.

Current limitations:

- The custom Octave implementation is sequential.
- The domain is the unit square.
- The mesh is a structured triangular P1 mesh.
- The subdomain decomposition is structured, rectangular and aligned with the grid.
- Homogeneous Dirichlet conditions are implemented by elimination.
- Homogeneous Neumann conditions are natural and do not add a boundary load term.
- General non-homogeneous boundary data are not the main supported workflow.
- General unstructured meshes are not supported.
- Heterogeneous coefficients and anisotropic diffusion are not part of the main implementation.
- Primal constraints are selected as structured subdomain corner/cross-point DOFs.
- Full spectra are computed only for small or moderate dimensions.
- PETSc data are stored and post-processed, but PETSc itself is not built or executed by the Octave code.

These restrictions are appropriate for the thesis goal: comparing the algebraic and spectral behavior of FETI-DP and BDDC in a controlled and reproducible setting.

---

## Cleaning generated outputs

The `output/` directory contains generated artifacts and can be removed if you want to rerun experiments from scratch:

```bash
rm -rf output/mats/ch6 output/figures/ch6 output/tables/ch6
```

Then regenerate with:

```bash
octave --eval "addpath('main'); run_ch6_spectral_experiments();"
```

The `.gitignore` is configured to ignore typical generated outputs such as figures, logs and `.mat` files.

---

## Main conclusion represented by the code

The repository is designed to support the thesis conclusion that, on the tested model problems, FETI-DP and BDDC behave as closely related dual-primal domain decomposition methods.

In the sequential Octave experiments, the non-preconditioned reduced operators may have substantially different and less favorable spectra, but after applying the corresponding dual-primal preconditioners, both methods produce compact spectra and fast PCG convergence on the tested cases.

The PETSc benchmark data provide a complementary view using official PETSc implementations on larger problems. Those results also show small iteration counts and favorable spectral indicators for both methods, while remaining a compact benchmark rather than a full parallel scaling study.
