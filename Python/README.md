# Python implementation

The Python package implements non-GUI numerical and coded-computation workflows
for Frame-codes, including code construction, matrix-multiplication experiments,
and article distribution figures.

## Introduction

Frame-codes provides concrete code families and experiment-specific operations.
The encoding, local worker simulation, pseudoinverse decoder, and Hydra runner
are provided by
[`codedcompute`](https://github.com/Royee-Yosibash/codedcompute), installed from
[PyPI](https://pypi.org/project/codedcompute/).

Workers execute locally and sequentially. Stragglers are simulated as missing
worker responses.

## Use cases

- Compare code families and worker precisions in coded matrix multiplication.
- Generate theoretical and empirical eigenvalue-distribution figures.
- Explore and Evaluate the Frame codes suggested in the paper
"IRREGULAR POLYNOMIAL CODES FOR CODED COMPUTATION WITH NUMERICAL STABILITY AND GRACEFUL DEGRADATION"

## Features

- Code families implementation include LPF, BPF, Wishart, Vandermonde, OrthoMatDot,
  Circulant Permutation, and SimpleBandCirculantCode.
- Composable Hydra configurations define experiment inputs, worker dtypes,
  straggler trials, and result collection.
- Numerical modules provide matrix partitioning, Gram-matrix calculations,
  distributions, and error metrics.

## Usage examples

Use the existing workflows and tests below as examples. Run them from the
`Python/` directory after installing the package as described below.

### Create a code matrix

See [`test_codes.py`](tests/unit_tests/test_codes.py) for factory selection,
code generation, and seeded random construction.
[`test_simple_band_circulant.py`](tests/unit_tests/test_simple_band_circulant.py)
demonstrates supported band configurations and recovery properties.

### Run a Hydra experiment

The configurations in [`experiments/`](experiments/) initialize the code family,
encoder, dispatcher, and decoder, then execute the configured pipeline.

| Configuration | Purpose | Results in the returned workspace |
| --- | --- | --- |
| [`simple_band_circulant_matrix_multiplication`](experiments/simple_band_circulant_matrix_multiplication.yaml) | Run matrix multiplication over every selected-size straggler combination. | `sweep_results` |
| [`matrix_multiplication_precision`](experiments/matrix_multiplication_precision.yaml) | Compare float64 and float16 workers on the same matrix product. | `sweep_results` |
| [`ncp_vs_circulant_matrix_multiplication`](experiments/ncp_vs_circulant_matrix_multiplication.yaml) | Compare non-consecutive-powers and Circulant Permutation codes over seeded trials. | `comparison_results` |

Shared initialization and pipeline steps are in
[`experiments/setup/`](experiments/setup/). Results are returned in memory;
the runner does not automatically write them to disk.

For configuration composition, execution, and result checks, see
[`test_hydra_distributed_matrix_multiplication.py`](tests/pipeline_tests/test_hydra_distributed_matrix_multiplication.py).
[`test_hydra_matrix_multiplication_precision.py`](tests/pipeline_tests/test_hydra_matrix_multiplication_precision.py)
demonstrates comparing the configured worker precisions.

For direct Python execution without Hydra, see
[`test_coding_pipeline.py`](tests/pipeline_tests/test_coding_pipeline.py).

### Generate article distribution figures

See [`figures/article_figures.py`](figures/article_figures.py) for
`create_article_figures()` and its `ArticleFigureConfig` parameters. The function
writes four PNG files to `Python/outputs/article_figures/` by default.

## Versioning

Latest version: **0.2.0**, released **November 15, 2026**.

Version 0.2.0 extracts the shared coded-computation framework into `codedcompute`.
Frame-codes now depends on `codedcompute>=0.5.0` from PyPI. Shared components are
imported from `codedcompute`; concrete code families remain under
`frame_codes.code_families`.

## Requirements and installation

Python 3.12 or newer is required. Runtime dependencies are `codedcompute`, Hydra
Core, NumPy, SciPy, and Matplotlib. OmegaConf is installed through `codedcompute`.

From the `Python/` directory, create an environment and install the package with
the test and lint tools:

```bash
python -m venv .venv
.venv/bin/python -m pip install -e ".[test,lint]"
```

The installation downloads `codedcompute` from PyPI. A separate local checkout
of that repository is not required.

### Tests and lint

```bash
.venv/bin/python -m unittest discover -s tests/unit_tests -p 'test_*.py' -v
.venv/bin/python -m unittest discover -s tests/pipeline_tests -p 'test_*.py' -v
.venv/bin/ruff check src tests
```

GitHub Actions runs Ruff and both test groups on Python 3.12. Unit tests cover
Frame-codes implementations; pipeline tests cover the code families, operations,
and experiment configurations working with `codedcompute`.
