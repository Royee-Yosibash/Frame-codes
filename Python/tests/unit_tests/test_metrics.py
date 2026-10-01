"""Unit tests for the standalone metric utilities."""

import unittest

import numpy as np

from frame_codes.numerics.metrics import error_measurement, frobenius_norm


class TestFrobeniusNorm(unittest.TestCase):
    """Test Frobenius norm behavior."""

    def test_calculates_magnitude_for_real_matrix(self) -> None:
        """The norm equals the square root of the sum of squares."""
        matrix = np.array([[3, 4], [0, 12]])

        self.assertEqual(frobenius_norm(matrix), 13.0)

    def test_calculates_magnitude_for_complex_matrix(self) -> None:
        """Complex entries contribute by magnitude."""
        matrix = np.array([[3 + 4j, 0]])

        self.assertEqual(frobenius_norm(matrix), 5.0)


class TestErrorMeasurement(unittest.TestCase):
    """Test mathematical error metric behavior."""

    def test_l1_sums_absolute_values_case_insensitively(self) -> None:
        """L1 sums magnitudes and accepts case-insensitive names."""
        self.assertEqual(error_measurement([-2, 3, -4], "l1"), 9.0)

    def test_l2_and_mmse_sum_squares(self) -> None:
        """L2 and MMSE both sum squared error magnitudes."""
        differences = np.array([-2, 3, -4])

        self.assertEqual(error_measurement(differences, "L2"), 29)
        self.assertEqual(error_measurement(differences, "MMSE"), 29)

    def test_l2_sums_complex_error_magnitudes(self) -> None:
        """Complex errors contribute their squared magnitudes."""
        self.assertAlmostEqual(error_measurement([1 + 1j], "L2"), 2.0)

    def test_rejects_unknown_metric(self) -> None:
        """Unknown metrics raise a clear value error."""
        with self.assertRaisesRegex(ValueError, "Unsupported error metric"):
            error_measurement([1, 2], "L-infinity")


if __name__ == "__main__":
    unittest.main()
