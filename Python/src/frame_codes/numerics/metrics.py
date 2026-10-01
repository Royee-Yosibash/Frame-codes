"""Mathematical error and norm metrics for numerical workflows."""

import numpy as np
from numpy.typing import ArrayLike


def frobenius_norm(matrix: ArrayLike) -> float:
    """Calculate the Frobenius norm of a matrix.

    Args:
        matrix: A two-dimensional array-like matrix.

    Returns:
        The square root of the sum of the squared element magnitudes.
    """
    values = np.asarray(matrix)
    return float(np.sqrt(np.sum(np.abs(values) ** 2)))


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
    values = np.asarray(differences)
    normalized_metric = metric.upper()

    if normalized_metric == "L1":
        return float(np.sum(np.abs(values)))
    if normalized_metric in {"L2", "MMSE"}:
        return float(np.sum(np.abs(values) ** 2))
    raise ValueError(f"Unsupported error metric: {metric}")
