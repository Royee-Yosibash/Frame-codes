"""Unit tests for Gram-matrix linear algebra utilities."""

import unittest

import numpy as np

from frame_codes.numerics.linear_algebra import (
    gram_matrix_condition_number,
    gram_matrix_eigenvalues,
    normalize_code,
)


class TestGramMatrixEigenvalues(unittest.TestCase):
    """Test Gram matrix eigenvalue calculations across shapes and dtypes."""

    def test_handles_tall_integer_matrix_without_integer_overflow(self) -> None:
        """Integer matrix products are promoted before forming the Gram matrix."""
        matrix = np.array([[50_000, 0], [0, 60_000], [0, 0]], dtype=np.int32)

        eigenvalues = gram_matrix_eigenvalues(matrix)

        np.testing.assert_allclose(eigenvalues, [2.5e9, 3.6e9])

    def test_handles_wide_float32_matrix(self) -> None:
        """Wide matrices use the row Gram matrix and return its spectrum."""
        matrix = np.array([[1, 0, 0], [0, 2, 0]], dtype=np.float32)

        eigenvalues = gram_matrix_eigenvalues(matrix)

        np.testing.assert_allclose(eigenvalues, [1, 4])

    def test_handles_complex64_matrix_with_conjugate_transpose(self) -> None:
        """Complex inputs use conjugate transpose in the Gram matrix."""
        matrix = np.array([[1 + 1j, 0], [0, 2j]], dtype=np.complex64)

        eigenvalues = gram_matrix_eigenvalues(matrix)

        np.testing.assert_allclose(eigenvalues, [2, 4])

    def test_rejects_non_matrix_inputs(self) -> None:
        """Eigenvalue calculation requires a non-empty matrix."""
        with self.assertRaisesRegex(ValueError, "non-empty two-dimensional"):
            gram_matrix_eigenvalues([1, 2, 3])


class TestGramMatrixConditionNumber(unittest.TestCase):
    """Test Gram matrix condition-number calculations."""

    def test_handles_integer_matrix(self) -> None:
        """Integer input is promoted before calculating the condition number."""
        matrix = np.array([[50_000, 0], [0, 60_000], [0, 0]], dtype=np.int32)

        self.assertAlmostEqual(gram_matrix_condition_number(matrix), 1.44)

    def test_rank_deficient_matrix_has_infinite_condition_number(self) -> None:
        """Rank-deficient inputs have an infinite Gram condition number."""
        matrix = np.array([[1.0, 2.0], [2.0, 4.0]])

        self.assertTrue(np.isinf(gram_matrix_condition_number(matrix)))

    def test_rejects_non_matrix_inputs(self) -> None:
        """Condition number calculation requires a non-empty matrix."""
        with self.assertRaisesRegex(ValueError, "non-empty two-dimensional"):
            gram_matrix_condition_number([1, 2, 3])


class TestNormalizeCode(unittest.TestCase):
    """Test row and column normalization across numeric dtypes."""

    def test_normalizes_integer_rows(self) -> None:
        """Integer matrices are promoted and each row has unit norm."""
        matrix = np.array([[3, 4], [5, 12]], dtype=np.int32)

        normalized = normalize_code(matrix, "row")

        self.assertEqual(normalized.dtype, np.float64)
        np.testing.assert_allclose(np.linalg.norm(normalized, axis=1), 1)

    def test_normalizes_complex_columns(self) -> None:
        """Complex columns have unit norm after normalization."""
        matrix = np.array([[3 + 4j, 0], [0, 5j]], dtype=np.complex64)

        normalized = normalize_code(matrix, "column")

        self.assertEqual(normalized.dtype, np.complex128)
        np.testing.assert_allclose(np.linalg.norm(normalized, axis=0), 1)

    def test_rejects_zero_norm_vectors(self) -> None:
        """A zero vector cannot be normalized."""
        with self.assertRaisesRegex(ValueError, "zero-norm"):
            normalize_code(np.zeros((2, 2)), "row")


if __name__ == "__main__":
    unittest.main()
