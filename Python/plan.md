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

## Porting Principles

- Do not delete, move, or modify the MATLAB implementation as part of the
  Python port.
- Use MATLAB as an initial behavioral reference for understanding workflows,
  not as a requirement for ongoing backwards compatibility.
- Preserve mathematical intent and workflow requirements. Correct or improve
  legacy implementation details when appropriate, documenting meaningful
  behavioral differences.
- Preserve the current default parameters from
  `Matlab/Coding Scheme/compareCodesConfig.m`.
- Do not introduce a random seed into the default workflow unless it is an
  explicit opt-in. Random-number generation is currently part of the existing
  behavior.
- Keep generated Python outputs separate from MATLAB outputs.
- MATLAB-generated numerical fixtures are useful when available, but are not a
  prerequisite for migration. Validate components through logical properties,
  dimensional invariants, analytical expectations, and established results in
  the literature. Any MATLAB comparison is a one-time migration check, not an
  ongoing compatibility test.

### MATLAB implementation observations

- `getCode(n, m, codeType, CodeSubType, ...)` returns `transpose(genMatrix)`;
  its documented output shape is `m x n`, while the generator used to encode a
  length-`m` vector is `n x m`. The optional `permIdx` selects a construction
  where supported; the implementation defaults it to zero (random selection).
- The implemented code types are `Simple`, `BCH`, `LPF`, `BPF`, `Wishart`,
  `Reed Solomon`, `Omitted Vandermonde`, `Non consecutive powers`,
  `OrthoMatDot`, and `Circulant Permutation`. Reed-Solomon subtypes include
  `None`, `Uniform sample`, `Random Set`, `Non-uniform sample`, and
  `Multiple non-uniform sample`. Dimensions and toolbox availability are
  type-specific; for example, Circulant Permutation requires even `n` and
  dimensions paired in twos. The `.mlx` source contains no uniform validation
  enforcing `m <= n` beyond its documentation.
- `EncodingScheme` obtains `Code` in the orientation above, then sets
  `EncodingMatA = transpose(Code)`, restoring generator orientation
  (`nNodes x mForCodeA`). `encodeMatrix` partitions matrix rows or columns
  evenly and forms each encoded part as a weighted sum using one encoding
  matrix row. `FrameParameters` instead transposes `Code` and stores an
  explicit inverse-based decoder derived from that frame matrix.
- `FrameParameters` can normalize rows, columns, or neither. Its constructor
  accepts scalar positive `M`, `N`, and `Gamma` (`Gamma <= 1`) and string
  `Type`, `SubType`, and `normDim` values (`None`, `Row`, `Column`).
- `FixMatricesDimensions` pads with zeros to achieve divisibility. OrthoMatDot
  pads columns of A and rows of B based on their respective partition counts;
  other active encoders pad rows of A and columns of B. Circulant Permutation
  follows the latter rule with paired row/column encoding.
- Gram eigenvalues are computed from `H' * H` when H is tall, otherwise
  `H * H'` (`'` is MATLAB's conjugate transpose); the helper also returns
  `cond(A)`. Coding-scheme condition-number statistics use an inverse for
  square selected matrices; otherwise they calculate an explicit
  pseudoinverse expression and take the square root of its Gram condition
  number.
- `compareCodesConfig.m` defaults are `numBits=8`, `mDataSets=[29]`,
  `nNodes=[31]`, `SNR=[80]`, `numTrials=10000`, `N=2800`, `n2=1`, and code
  types `Non consecutive powers` and `Circulant Permutation`. `compareCodes.m`
  generates integer A and B entries in `[-51, 49]`, with A having
  `22*mDataSets*30` rows, then measures MSE, relative Frobenius error, and
  mean/min/max decoder condition numbers.
- Decoding adds independent real Gaussian noise with variance
  `10^(-SNR/10)` to each surviving coded result. Reported MSE is the mean of
  `sum(abs(error).^2) / sum(abs(reconstruction).^2)`; relative Frobenius error
  is `||error||_F / ||C||_F`.
- `manovaPDF(x, beta, gamma)` and `marchenkoPasturPDF(x, beta)` normalize
  sampled density values using `normPDF`; MANOVA may additionally transfer
  mass to the grid point nearest `1/gamma`. Article-figure defaults include
  `gamma=0.5`, `beta=0.8`, and 400-dimensional coded/Wishart examples.

## Proposed Python Layout

```text
Python/
  plan.md
  pyproject.toml
  src/
    frame_codes/
      __init__.py
      numerics/
        __init__.py
        metrics.py
        linear_algebra.py
        distributions.py
      frames/
        __init__.py
        codes.py
        parameters.py
      coding/
        __init__.py
        encoding.py
        decoding.py
      experiments/
        __init__.py
        compare_codes.py
        article_figures.py
      results.py
  tests/
    unit_tests/
      test_metrics.py
      test_codes.py
      test_encoding.py
      test_decoding.py
      test_distributions.py
  outputs/
```

Modules are grouped by responsibility; shared numerical utilities remain
separate from coding workflows and experiment drivers.

## Implementation Stages

### 1. Port shared utilities

Port utilities in dependency order, starting with the modules that can be
checked independently:

1. **Standalone linear algebra:** Gram-matrix eigenvalues and condition
   numbers.
2. **Standalone distributions:** PDF normalization, MANOVA and
   Marchenko-Pastur densities, empirical eigenvalue histograms, and sampling
   from a supplied distribution.
3. **Matrix shape utilities:** zero-padding/dimension fixing and equal row or
   column partitioning.
4. **Code construction:** supported code matrices and normalization.
5. **Encoding:** coded matrix partitions using the code matrices and shape
   utilities.
6. **Frame parameters/statistics:** this layer depends on code construction,
   Gram-matrix calculations, and distribution/statistics helpers.

Each port should have focused checks for its logical properties and
dimensional invariants. Validate theoretical numerical routines against
established literature results. If a utility has no dependency on other
project modules, it may be moved earlier after confirming that its validation
does not rely on unported behavior. MATLAB comparisons may be performed once
during migration, but should not be repeated as backwards-compatibility tests.

### 2. Port coding-scheme experiments

- Port the defaults from `compareCodesConfig.m` into Python configuration.
- Implement the intended encoding and decoding workflows for the currently
  supported code types.
- Preserve straggler selection, noise calculation, matrix padding, and code
  orientation.
- Produce Python-specific output files under `Python/outputs/`.
- Keep result fields and metric definitions consistent across Python workflows.

### 3. Port article/distribution figures

- Reproduce theoretical MANOVA and Marchenko-Pastur curves.
- Reproduce coded-frame eigenvalue experiments.
- Reproduce Wishart eigenvalue experiments.
- Compare figures visually and compare sampled distributions numerically.
- Do not overwrite figures under `Matlab/Results/`.

### 4. Add mathematical and workflow tests

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

Use small matrices for fast unit tests and retain a separate, slower
integration test for the full comparison workflow. Tests should assert
mathematical invariants and intended workflow behavior. Do not encode legacy
MATLAB quirks as ongoing compatibility requirements.

## Known Risks

- `getCode.mlx` is the central dependency but is not currently a plain MATLAB
  function, so its complete behavior must be characterized before porting.
- MATLAB and NumPy differ in default random-number generators and random
  sampling behavior; stochastic outputs should be validated statistically,
  not expected to match sample-by-sample.
- MATLAB uses 1-based indexing while Python uses 0-based indexing.
- MATLAB's `transpose` and conjugate-transpose operations must remain distinct
  for complex-valued code matrices.
- The MATLAB implementation uses explicit matrix inverses in several places;
  prefer numerically stable linear solves where they preserve the intended
  mathematical behavior.
- Existing MATLAB scripts contain workspace and current-directory assumptions;
  the Python implementation should replace these with explicit function
  inputs and output paths.

## Completion Criteria

The migration is complete when:

- The MATLAB workflows remain intact and documented.
- Python's numerical components satisfy documented mathematical properties
  and established analytical/theoretical expectations.
- Python reproduces the active comparison metrics and distribution figures.
- Python outputs are separated from MATLAB outputs.
- Tests cover the shared utilities and both non-GUI workflows.
- The root README links to the Python plan and usage documentation.
