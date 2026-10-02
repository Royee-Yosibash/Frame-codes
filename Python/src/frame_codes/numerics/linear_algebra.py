"""Linear algebra utilities for Gram matrices."""

import numpy as np
from numpy.typing import ArrayLike, NDArray

from frame_codes.utils.numpy_types import as_inexact_array


def gram_matrix_eigenvalues(matrix: ArrayLike) -> NDArray[np.float64]:
    """Calculate the eigenvalues of the smaller Gram matrix.

    Integer and lower-precision floating-point inputs are promoted to double
    precision before matrix products.

    Args:
        matrix: A non-empty two-dimensional real or complex numeric matrix.

    Returns:
        Ascending eigenvalues of the smaller Gram matrix.

    Raises:
        ValueError: If the input is not a non-empty two-dimensional matrix.
        TypeError: If the input dtype is not an integer, floating-point, or
            complex numeric type.
    """
    values = _prepare_matrix(matrix)
    gram_matrix = _calculate_gram_matrix(values)
    return np.linalg.eigvalsh(gram_matrix)


def gram_matrix_condition_number(matrix: ArrayLike) -> float:
    """Calculate the condition number of the smaller Gram matrix.

    Integer and lower-precision floating-point inputs are promoted to double
    precision before calculations.

    Args:
        matrix: A non-empty two-dimensional real or complex numeric matrix.

    Returns:
        The Gram-matrix condition number, or infinity when the input is
        numerically rank deficient.

    Raises:
        ValueError: If the input is not a non-empty two-dimensional matrix.
        TypeError: If the input dtype is not an integer, floating-point, or
            complex numeric type.
    """
    values = _prepare_matrix(matrix)
    singular_values = np.linalg.svd(values, compute_uv=False)
    if np.linalg.matrix_rank(values) < singular_values.size:
        return float("inf")
    return float((singular_values[0] / singular_values[-1]) ** 2)


def _prepare_matrix(matrix: ArrayLike) -> NDArray[np.float64] | NDArray[np.complex128]:
    """Validate a matrix and promote its dtype for numerical calculations.

    Args:
        matrix: A non-empty two-dimensional real or complex numeric matrix.

    Returns:
        The matrix promoted to `float64` or `complex128`.

    Raises:
        ValueError: If the input is not a non-empty two-dimensional matrix.
        TypeError: If the input dtype is not an integer, floating-point, or
            complex numeric type.
    """
    values = np.asarray(matrix)
    if values.ndim != 2 or 0 in values.shape:
        raise ValueError("matrix must be a non-empty two-dimensional array")
    return as_inexact_array(values)


def _calculate_gram_matrix(
    matrix: NDArray[np.float64] | NDArray[np.complex128],
) -> NDArray[np.float64] | NDArray[np.complex128]:
    """Calculate the smaller Gram matrix using a conjugate transpose.

    Args:
        matrix: A promoted two-dimensional real or complex matrix.

    Returns:
        The smaller of `matrix.conj().T @ matrix` and
        `matrix @ matrix.conj().T`.
    """
    if matrix.shape[0] > matrix.shape[1]:
        return matrix.conj().T @ matrix
    return matrix @ matrix.conj().T
