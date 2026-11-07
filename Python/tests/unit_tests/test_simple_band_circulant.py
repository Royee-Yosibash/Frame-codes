"""Tests for the SimpleBandCirculant code family."""

import unittest
from itertools import combinations

import numpy as np

from frame_codes.code_families.factory import create_code_family
from frame_codes.code_families.simple_band_circulant import (
    SimpleBandCirculantCode,
)
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters


class TestSimpleBandCirculantCode(unittest.TestCase):
    """Verify band construction, supported worker counts, and erasure rank."""

    def test_generates_each_supported_three_element_band_prefix(self) -> None:
        """The factory builds the identity, first band, and terminal band cases."""
        expected_codes = {
            3: np.array(
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            4: np.array(
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [1.0, 1.0, 1.0],
                ]
            ),
            6: np.array(
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [1.0, 1.0, 0.0],
                    [0.0, 1.0, 1.0],
                    [1.0, 0.0, 1.0],
                ]
            ),
            7: np.array(
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [1.0, 1.0, 0.0],
                    [0.0, 1.0, 1.0],
                    [1.0, 0.0, 1.0],
                    [1.0, 1.0, 1.0],
                ]
            ),
        }

        for worker_count, expected in expected_codes.items():
            with self.subTest(worker_count=worker_count):
                family = create_code_family("SimpleBandCirculantCode")
                self.assertIsInstance(family, SimpleBandCirculantCode)

                family.generate_code(CodeParameters(m=3, n=worker_count))

                np.testing.assert_array_equal(family.get_code(), expected)

    def test_rejects_unsupported_worker_counts_and_single_element_codes(self) -> None:
        """Only the prescribed band-prefix row counts are supported."""
        family = SimpleBandCirculantCode()

        for worker_count in (5, 8):
            with (
                self.subTest(worker_count=worker_count),
                self.assertRaisesRegex(ValueError, "supported"),
            ):
                family.generate_code(CodeParameters(m=3, n=worker_count))

        with self.assertRaisesRegex(ValueError, "m >= 2"):
            family.generate_code(CodeParameters(m=1, n=1))

    def test_generates_multiple_regular_bands_before_the_terminal_band(self) -> None:
        """Band indices advance the run length while each regular band has m rows."""
        family = SimpleBandCirculantCode()
        family.generate_code(CodeParameters(m=4, n=13))
        expected = np.array(
            [
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
                [1, 1, 0, 0],
                [0, 1, 1, 0],
                [0, 0, 1, 1],
                [1, 0, 0, 1],
                [1, 1, 1, 0],
                [0, 1, 1, 1],
                [1, 0, 1, 1],
                [1, 1, 0, 1],
                [1, 1, 1, 1],
            ],
            dtype=np.float64,
        )

        np.testing.assert_array_equal(family.get_code(), expected)

    def test_final_band_recovers_from_any_four_worker_rows(self) -> None:
        """Every four-row survivor matrix has full column rank."""
        family = SimpleBandCirculantCode()
        family.generate_code(CodeParameters(m=3, n=7))
        code = family.get_code()

        for row_indices in combinations(range(7), 4):
            with self.subTest(row_indices=row_indices):
                self.assertEqual(np.linalg.matrix_rank(code[list(row_indices)]), 3)

        self.assertEqual(np.linalg.matrix_rank(code[[0, 1, 3]]), 2)


if __name__ == "__main__":
    unittest.main()
