# MATLAB to Python Migration Plan

## Goal

Create a Python implementation of the active, non-GUI MATLAB workflows while
keeping the MATLAB implementation unchanged and available as the reference
implementation.

The migration targets numerical and experiment workflows only. The MATLAB GUI,
App Designer files, and archived MATLAB code are out of scope for the initial
port.

## Reference Workflows

The active MATLAB workflows to reproduce are:

1. **Shared frame and numerical utilities**
   - Code construction through `FramesTOOLBOX/getCode.mlx`.
   - Frame parameter handling through `FrameParameters.m`.
   - Matrix encoding and dimension padding.
   - Gram-matrix eigenvalues and condition numbers.
   - Noise-amplification and error measurements.
   - MANOVA, Marchenko-Pastur, and empirical eigenvalue distributions.

2. **Coding-scheme comparison**
   - MATLAB entry point: `Matlab/Coding Scheme/compareCodes.m`.
   - Encoding: `EncodingScheme.m`.
   - Decoding under stragglers and noise: `DecodingScheme.m`.
   - Condition-number statistics.
   - MSE, Frobenius-norm error, and LaTeX-oriented result tables.

3. **Article and distribution figures**
   - MATLAB entry point: `Matlab/Results/CreateFiguresForArticle.m`.
   - Theoretical MANOVA and Marchenko-Pastur distributions.
   - Empirical eigenvalue simulations for coded and Wishart matrices.

## Compatibility Rules

- Do not delete, move, or modify the MATLAB implementation as part of the
  Python port.
- Treat MATLAB outputs and behavior as the compatibility reference.
- Preserve matrix dimensions, row/column conventions, code transposition, and
  normalization behavior exactly before improving the design.
- Preserve the current default parameters from
  `Matlab/Coding Scheme/compareCodesConfig.m`.
- Do not introduce a random seed into the default workflow unless it is an
  explicit opt-in. Random-number generation is currently part of the existing
  behavior.
- Compare floating-point results with documented tolerances rather than exact
  equality.
- Keep generated Python outputs separate from MATLAB outputs.

## Proposed Python Layout

```text
Python/
  plan.md
  pyproject.toml
  src/
    frame_codes/
      __init__.py
      config.py
      codes.py
      frames.py
      encoding.py
      decoding.py
      linear_algebra.py
      distributions.py
      metrics.py
      results.py
  scripts/
    compare_codes.py
    create_article_figures.py
  tests/
    test_codes.py
    test_encoding.py
    test_decoding.py
    test_distributions.py
    test_metrics.py
  outputs/
```

The exact module boundaries may be adjusted during implementation, but shared
numerical utilities should remain separate from experiment drivers.

## Implementation Stages

### 1. Establish the MATLAB reference contract

- Record supported code types and their accepted dimensions.
- Extract the exact input/output behavior of `getCode.mlx`.
- Document whether returned matrices are generators, transposes, or frame
  matrices at each call site.
- Record default parameters, normalization rules, padding rules, and output
  units.
- Create small reference cases with saved matrices and numerical outputs.

### 2. Set up the Python project

- Add a minimal `pyproject.toml`.
- Use NumPy for arrays and linear algebra.
- Use SciPy where MATLAB behavior depends on specialized numerical routines.
- Use Matplotlib for figures.
- Use pandas only where table export is useful.
- Keep runtime dependencies minimal and document supported Python versions.

### 3. Port shared utilities

Port and test the utilities in this order:

1. Matrix dimensions and partitioning.
2. Code construction and normalization.
3. Encoding.
4. Gram-matrix eigenvalues and condition numbers.
5. Error and norm metrics.
6. PDF and empirical-distribution helpers.
7. Frame parameter/statistics behavior.

Each port should have a focused test against a MATLAB reference case before
the next layer is implemented.

### 4. Port coding-scheme experiments

- Port the defaults from `compareCodesConfig.m` into Python configuration.
- Implement encoding and decoding without changing the current supported
  code paths.
- Preserve straggler selection, noise calculation, matrix padding, and code
  orientation.
- Produce Python-specific output files under `Python/outputs/`.
- Keep result columns and metric definitions compatible with the MATLAB tables.

### 5. Port article/distribution figures

- Reproduce theoretical MANOVA and Marchenko-Pastur curves.
- Reproduce coded-frame eigenvalue experiments.
- Reproduce Wishart eigenvalue experiments.
- Compare figures visually and compare sampled distributions numerically.
- Do not overwrite figures under `Matlab/Results/`.

### 6. Add parity and regression tests

Tests should cover:

- Code dimensions and orientation.
- Deterministic code construction.
- Matrix partitioning and zero-padding.
- Encoding and decoding without noise or stragglers.
- Decoding with representative straggler counts and SNR values.
- Gram eigenvalues and condition numbers.
- Noise-amplification metrics.
- Distribution normalization and support.
- Result-table fields and output paths.

Use small matrices for fast tests and retain a separate, slower integration
test for the full comparison workflow.

## Known Risks

- `getCode.mlx` is the central dependency but is not currently a plain MATLAB
  function, so its complete behavior must be characterized before porting.
- MATLAB and NumPy differ in default random-number generators and random
  sampling behavior; numerical parity for stochastic experiments will require
  controlled reference inputs or statistical tolerances.
- MATLAB uses 1-based indexing while Python uses 0-based indexing.
- MATLAB's `transpose` and conjugate-transpose operations must remain distinct
  for complex-valued code matrices.
- The MATLAB implementation uses explicit matrix inverses in several places;
  the first Python port should reproduce results before considering more stable
  linear solves.
- Existing MATLAB scripts contain workspace and current-directory assumptions;
  the Python implementation should replace these with explicit function
  inputs and output paths without changing the MATLAB behavior.

## Completion Criteria

The migration is complete when:

- The MATLAB workflows remain intact and documented.
- Python reproduces representative MATLAB numerical results within declared
  tolerances.
- Python reproduces the active comparison metrics and distribution figures.
- Python outputs are separated from MATLAB outputs.
- Tests cover the shared utilities and both non-GUI workflows.
- The root README links to the Python plan and usage documentation.
