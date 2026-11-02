"""Mathematical error and norm metrics for numerical workflows."""

import numpy as np
from numpy.typing import ArrayLike

from frame_codes.utils.numpy_types import as_inexact_array


def frobenius_norm(matrix: ArrayLike) -> float:
    """Calculate the Frobenius norm of a matrix.

    Args:
        matrix: A two-dimensional array-like matrix.

    Returns:
        The square root of the sum of the squared element magnitudes.
    """
    values = as_inexact_array(matrix)
    return float(np.sqrt(np.sum(np.abs(values) ** 2)))


def relative_frobenius_error(actual: ArrayLike, expected: ArrayLike) -> float:
    """Calculate the Frobenius error relative to an expected result.

    Args:
        actual: Result produced by a computation.
        expected: Reference result.

    Returns:
        `||actual - expected||_F / ||expected||_F`.

    Raises:
        ValueError: If the arrays have different shapes or the expected result
            has zero Frobenius norm.
        TypeError: If either input has a non-numeric dtype.
    """
    actual_values = as_inexact_array(actual)
    expected_values = as_inexact_array(expected)
    if actual_values.shape != expected_values.shape:
        raise ValueError("actual and expected must have matching shapes")
    expected_norm = frobenius_norm(expected_values)
    if expected_norm == 0:
        raise ValueError("expected result must have nonzero Frobenius norm")
    return frobenius_norm(actual_values - expected_values) / expected_norm


def error_measurement(differences: ArrayLike, metric: str) -> float:
    """Calculate an error summary using a supported metric.

    Args:
        differences: Array-like error values.
        metric: Metric name, case-insensitively. Supported values are `L1`,
            `L2`, and `MMSE`.

    Returns:
        The sum of absolute values for `L1`, or the sum of squared magnitudes
        for `L2` and `MMSE`.

    Raises:
        ValueError: If the metric name is not supported.
    """
    values = as_inexact_array(differences)
    normalized_metric = metric.upper()

    if normalized_metric == "L1":
        return float(np.sum(np.abs(values)))
    if normalized_metric in {"L2", "MMSE"}:
        return float(np.sum(np.abs(values) ** 2))
    raise ValueError(f"Unsupported error metric: {metric}")
